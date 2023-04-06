"""
node.py

This model computes equilibrium temperatures and thermal balances for the solar orbiter spacecraft. The code is
just a loop that goes through each shielding condition, each thermal condition, and each node of the spacecraft
to print out equilibrium temperatures and thermal balances.

Some key assumptions:
    - This is technically 2 separate single node analyses for the spacecraft, rather than a 2 node analysis. This is
    because the view factor between the solar sails and the spacecraft bus is assumed to be negligible, as is conductive
    contact due to the lack of mechanical booms.
        - The two separate nodes are thus the spacecraft bus (with heat shield), and the solar sails.
    - Solar cells were not considered at this phase of design in this model.
    - Conductive transfer in the interface between the shield and spacecraft was also not considered, and was assumed small
    through the use of potential thermal isolation techniques.

TODO:

"""

import pandas as pd
import numpy as np
from scipy import optimize
from scipy import integrate
import time
from tqdm import tqdm
import matplotlib.pyplot as plt

SB = 5.67e-8
AU = 1.496e11
solar_luminosity = 3.83e26


class Node:
    def __init__(self, name, properties, n_shield=0):
        """
        This sets up each Node class with the surface properties and the area.

        Inputs:

        name [str]: The name of the node
        skin_properties [arr]: Input an array consisting of [emissivity, absorptivity, reflectivity, density, area, sun view factor, space view factor, temperature range, and internal heat].
        area [arr]: Input an array consisting of [area that the sun sees, area that deep space sees, area covered by radiators.]
        n_shield [int]: The number of layers for the heat shield.
        """
        self.name = name
        self.emissivity = properties[0]
        self.reflectivity = properties[1]
        self.absorptivity = properties[2]
        self.density = properties[3]
        self.area = properties[4]
        self.sun_vf = properties[5]
        self.space_vf = properties[6]
        self.temp_range = properties[7]
        self.internal_heat = properties[8]
        self.n_shield = n_shield

        # Initialization temperature
        self.temp = 263

    def solar_heat_in(self, thermal_case, sail_deployed, panel_duty_cycle, payload_duty_cycle):
        """
        This computes the solar heat that the outer surface sees.

        The formula is Q = A * (L * a * cos(theta) * F / (4 * PI * R^2))
            - A = area of surface accepting Sun's heat [m^2]
            - L = luminosity [W]
            - a = absorptivity
            - F = view factor based on if the part of the spacecraft faces the sun or not.
            - theta = angle to the Sun (0 is such that the spacecraft is facing the sun directly)
            - R = distance from Sun [m]

        Outputs:

        self.solar_q_in [float]: Incoming heat from the sun in W.
        """
        if self.name == "Solar Sail" or self.name == "Boom":
            self.solar_flux = (sail_deployed * np.cos(thermal_case[1]) * solar_luminosity * self.absorptivity * self.sun_vf
            ) / (4 * np.pi * (thermal_case[0] * AU) ** 2)
        elif self.name == "Solar Panel":
            self.solar_flux = (panel_duty_cycle * np.cos(thermal_case[1]) * solar_luminosity * self.absorptivity * self.sun_vf
            ) / (4 * np.pi * (thermal_case[0] * AU) ** 2)
        elif self.name == "Doppler Magnetograph" or self.name == "Coronagraph":
            self.solar_flux = (payload_duty_cycle * np.cos(thermal_case[1]) * solar_luminosity * self.absorptivity * self.sun_vf) / (4 * np.pi * (thermal_case[0] * AU) ** 2)
        else:
            self.solar_flux = (np.cos(thermal_case[1]) * solar_luminosity * self.absorptivity * self.sun_vf) / (4 * np.pi * (thermal_case[0] * AU) ** 2)
        self.solar_q_in = self.solar_flux * self.area
        return self.solar_q_in

    def radiated_heat_coeff(self, sail_deployed):
        """
        Computes the coefficient for heat radiated by a certain surface.

        """
        if self.name == "Solar Sail" or self.name == "Boom":
            q_rad_coeff = -sail_deployed * self.area * self.emissivity * self.space_vf * SB
        else:
            q_rad_coeff = -self.area * self.emissivity * self.space_vf * SB
        return q_rad_coeff

    def heat_radiated(self, sail_deployed):
        """
        Computes the heat radiated by a certain surface.

        The formula is Q = A * e * sigma * T^4
            - A = area of surface
            - e = emissivity of the surface
            - sigma = stefan-boltzmann constant
            - T = surface temperature

        Inputs:
            - eq_temp [K]: The equilibrium temperature of the surface (calculated using temperatures())

        Outputs:

            self.q_radiated [float]: The total heat radiated by the node.
        """
        if self.name == "Solar Sail" or self.name == "Boom":
            self.q_radiated = (
                sail_deployed * self.emissivity * self.space_vf * SB * self.area * self.temp**4
            )
        else:
            self.q_radiated = (self.emissivity * self.space_vf * SB * self.area * self.temp**4
            )
        return self.q_radiated


def equations(node_temps, rad_matrix, cond_matrix, q_extra):
    '''
    System of equations used for calculating temperatures in steady-state analysis.
    
    --- DEPRECATED CODE ---
    '''
    system_of_equations = (
        np.dot(rad_matrix, node_temps**4)
        + np.dot(cond_matrix, node_temps)
        + q_extra
    )
    return system_of_equations

def heat(nodes_examined, R_ij, G_ij):
    '''
    Adds radiative and conductive heat between two nodes.

    --- DEPRECATED CODE ---
    '''
    q_rad = R_ij*(nodes_examined[1].temp**4 - nodes_examined[0].temp**4)
    q_cond = G_ij*(nodes_examined[1].temp - nodes_examined[0].temp)
    q_nodes = q_rad + q_cond
    return q_nodes


def temp_update(node_temperatures, rad_matrix, cond_matrix, capacities, q_extra):
    '''
    Recalculates the radiative, conductive, solar, and internal heat contributions. 

    Inputs:
        node_temperatures (np.array): temperatures at each node (K)
        rad_matrix (np.array): N x N matrix of radiative coefficients for each node where N = number of spacecraft nodes.
        cond_matrix (np.array): N x N matrix of conductive coefficients for each node where N = number of spacecraft nodes.
        capacities (np.array): N x 1 array of heat capacities for each node where N = number of spacecraft nodes.
        q_extra (np.array): N x 1 array of solar heat input + internal heat.
    '''
    q_rad = np.dot(rad_matrix, node_temperatures**4)
    q_cond = np.dot(cond_matrix, node_temperatures)
    q_total = (q_rad + q_cond + q_extra)
    return np.divide(q_total, capacities)

def runge_kutta(node_temperatures, time_step, capacities, rad_matrix, cond_matrix, q_extra):
    '''
    Computes one RK4 integration step for the temperatures.

    Inputs:
        node_temperatures (np.array): temperatures at each node (K)
        rad_matrix (np.array): N x N matrix of radiative coefficients for each node where N = number of spacecraft nodes.
        cond_matrix (np.array): N x N matrix of conductive coefficients for each node where N = number of spacecraft nodes.
        capacities (np.array): N x 1 array of heat capacities for each node where N = number of spacecraft nodes.
        q_extra (np.array): N x 1 array of solar heat input + internal heat.
    '''
    k1 = time_step * temp_update(node_temperatures, rad_matrix, cond_matrix, capacities, q_extra)
    k2 = time_step * temp_update(node_temperatures + (k1/2), rad_matrix, cond_matrix, capacities, q_extra)
    k3 = time_step * temp_update(node_temperatures + (k2/2), rad_matrix, cond_matrix, capacities, q_extra)
    k4 = time_step * temp_update(node_temperatures + k3, rad_matrix, cond_matrix, capacities, q_extra)
    k = (k1 + (2*k2) + (2*k3) + k4) / 6
    node_temperatures = node_temperatures + k
    return node_temperatures



def steady_state(nodes, relationships, sail_deployed, duty_cycle, thermal_case):
    '''
    Simultaneously solves temperatures at each time step using a steady-state thermal
    analysis.

    --- DEPRECATED CODE ---
    '''
    node_initial = np.asarray([nodes[i].temp for i in range(0, len(nodes))])
    radiative = np.triu(relationships, 1)
    radiative = np.where(radiative, radiative, radiative.T)
    conductive = np.tril(relationships, -1)
    conductive = np.where(conductive, conductive, conductive.T)
    capacities = np.diagonal(relationships)
    q_extra = []
    q_rad_coeff = []
    for i in range(0, len(nodes)):
        q_extra.append(
            nodes[i].solar_heat_in(thermal_case, sail_deployed, duty_cycle)
            + nodes[i].internal_heat
        )
        q_rad_coeff.append(nodes[i].radiated_heat_coeff(sail_deployed))

    cond_matrix = conductive + np.diag(-np.sum(conductive, axis=1))


    rad_matrix = np.diag(q_rad_coeff - np.sum(radiative, axis=1)) + radiative
    node_temperatures = optimize.newton(
        equations, node_initial, args=(rad_matrix, cond_matrix, q_extra), full_output=True, maxiter=1000
    )
    for idx, node in enumerate(nodes):
        node.temp = node_temperatures[idx]
    return node_temperatures


def time_variant(
    nodes, relationships, dt, sail_deployed, thermal_case, temp_ranges, verbose, plot):
    '''
    Calculates the temperature at each node for a single orbital time step by dividing 
    the orbit model time step into time steps required to keep the thermal model stable.

    Inputs:
        - nodes (array): input array containing each node object.
        - relationships (matrix): matrix from config.py containing conductive, radiative, and capacity information.
        - dt (int): time step in seconds from the orbital model.
        - sail_deployed (int): 1 or 0 flag to determine if the solar sail is deployed.
        - thermal_case (arr): array of [distance from sun (AU), angle from the sun (rad)]
        - temp_ranges (arr): array of temperature ranges for each component.
        - verbose (bool): True to print out temperature information, False to suppress printouts.
        - plot (bool): True to print active thermal control plots.

    Returns:
        - node_temperatures (arr): Array containing temperatures at each node.
    '''
    t_start = time.perf_counter()
    dt_original = dt
    node_temperatures = np.asarray([nodes[i].temp for i in range(0, len(nodes))])
    radiative = np.triu(relationships, 1)
    radiative = np.where(radiative, radiative, radiative.T)
    conductive = np.tril(relationships, -1)
    conductive = np.where(conductive, conductive, conductive.T)
    capacities = np.diagonal(relationships)
    panel_duty_cycle = 0.01
    payload_duty_cycle = 0.3

    q_rad_coeff = []
    q_extra = []
    for i in range(0, len(nodes)):
        q_extra.append(
            nodes[i].solar_heat_in(thermal_case, sail_deployed, panel_duty_cycle, payload_duty_cycle)
            + nodes[i].internal_heat
        )
        q_rad_coeff.append(nodes[i].radiated_heat_coeff(sail_deployed))

    rad_matrix = np.diag(q_rad_coeff - np.sum(radiative, axis=1)) + radiative
    cond_matrix = conductive + np.diag(-np.sum(conductive, axis=1))
    
    # Interpolate time steps to ensure thermal model is stable
    dt_interp = np.min(0.001*capacities)
    if dt_interp > 0:
        dt_arr = [dt_interp]*int(dt_original/dt_interp)
    else:
        dt_arr = [0]

    heater_cooler_power = []
    temp_steps = []
    for time_idx, time_step in enumerate(dt_arr):
        # Update node temperatures
        node_temperatures = runge_kutta(node_temperatures, time_step, capacities, rad_matrix, cond_matrix, q_extra)
        # Active thermal control adjustments
        for idx, temp in enumerate(node_temperatures):
            if idx == 12:
                tempC = temp - 273.15
                if (tempC < temp_ranges[idx][0]+10):
                    nodes[idx].internal_heat += 0.02
                elif (tempC > temp_ranges[idx][1]-10):
                    nodes[idx].internal_heat += 0.02
                elif abs(tempC - (temp_ranges[idx][1] + temp_ranges[idx][0])/2) < 1:
                    nodes[idx].internal_heat = 0
            if idx == 13:
                tempC = temp - 273.15
                if (tempC < temp_ranges[idx][0]+10):
                    nodes[idx].internal_heat += 0.05
                elif (tempC > temp_ranges[idx][1]-10):
                    nodes[idx].internal_heat -= 0.05
                elif abs(tempC - (temp_ranges[idx][1] + temp_ranges[idx][0])/2) < 1:
                    nodes[idx].internal_heat = 0
            if idx == 14:
                tempC = temp - 273.15
                if (tempC < temp_ranges[idx][0]+10):
                    nodes[idx].internal_heat += 0.5
                elif (tempC > temp_ranges[idx][1]-10):
                    nodes[idx].internal_heat -= 0.5
                elif abs(tempC - ((temp_ranges[idx][1] + temp_ranges[idx][0])/2)) < 1:
                    nodes[idx].internal_heat = 0
            if idx == 15:
                tempC = temp - 273.15
                if (tempC < temp_ranges[idx][0]+10):
                    nodes[idx].internal_heat += 0.01
                elif (tempC > temp_ranges[idx][1]-10):
                    nodes[idx].internal_heat -= 0.01
                elif abs(tempC - (temp_ranges[idx][1] + temp_ranges[idx][0])/2) < 1:
                    nodes[idx].internal_heat = 0
            q_extra = []
            for i in range(0, len(nodes)):
                q_extra.append(
                    nodes[i].solar_heat_in(thermal_case, sail_deployed, panel_duty_cycle, payload_duty_cycle)
                    + nodes[i].internal_heat
                )

        temp_steps.append([node_temps - 273.15 for node_temps in node_temperatures])
        heater_cooler_power.append(abs(nodes[12].internal_heat) + abs(nodes[13].internal_heat) + abs(nodes[14].internal_heat) + abs(nodes[15].internal_heat))

        
    for i in range(0, len(nodes)):
        nodes[i].temp = node_temperatures[i]

    disp_temps = np.round(node_temperatures - 273.15, 2)
    dist = np.round(thermal_case[0], 2)

    t_stop = time.perf_counter()

    dt_arr_plot = np.multiply(np.asarray(dt_arr), np.asarray(range(0, len(dt_arr))))

    if plot:
        plt.rcParams.update({'font.size': 32})
        plt.figure()
        plt.title(f'Active Thermal Control at 0.38 AU for METIS')
        plt.plot(dt_arr_plot, np.asarray(temp_steps)[:, 14], linewidth=4, label='METIS')
        plt.plot(dt_arr_plot, [nodes[14].temp_range[0]]*len(dt_arr_plot), linestyle='dashed', color='red', linewidth=4, label='METIS Limits')
        plt.plot(dt_arr_plot, [nodes[14].temp_range[1]]*len(dt_arr_plot), linestyle='dashed', color='red', linewidth=4)
        plt.legend(loc='center right', bbox_to_anchor=(0.95,0.3))
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (°C)')


    if verbose:
        print("=======================================================")      
        print(f"distance: {dist} AU")
        print(f"Spacecraft +Z (sun-facing): {disp_temps[0]} C, Spacecraft -Z (space-facing): {disp_temps[1]} C")
        print(f"Spacecraft +X: {disp_temps[2]} C, Spacecraft -X: {disp_temps[3]} C")
        print(f"Spacecraft +Y: {disp_temps[4]} C, Spacecraft -Y: {disp_temps[5]} C")
        print(f"Sails: {disp_temps[6]} C, Booms: {disp_temps[7]} C")
        print(f"Panels: {disp_temps[8]} C, Antenna: {disp_temps[9]} C")
        print(f"Outer Shield: {disp_temps[10]} C, Inner Shield: {disp_temps[11]} C")
        print(f"Batteries: {disp_temps[12]} C, Hydrazine: {disp_temps[13]} C")
        print(f"METIS: {disp_temps[14]} C, CDM: {disp_temps[15]} C")
        print(f"Average Power:  {np.round(np.mean(np.asarray(heater_cooler_power)), 2)} W, Maximum Power:  {np.round(np.max(np.asarray(heater_cooler_power)), 2)} W")
        print(f"Runtime: {round(t_stop - t_start, 2)} s")
    
    return node_temperatures
