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

SB = 5.67e-8
AU = 1.496e11
solar_luminosity = 3.83e26


class Node:
    def __init__(self, name, properties, n_shield=0):
        """
        This sets up each Node class with the surface properties and the area.

        Inputs:

        name [str]: The name of the node
        skin_properties [arr]: Input an array consisting of [emissivity, absorptivity, reflectivity, radiator_emissivity], all scaled from 0 to 1.
        area [arr]: Input an array consisting of [area that the sun sees, area that deep space sees, area covered by radiators.]
        case [arr]: Input an array consisting of [distance in AU, angle from sun in rad]
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
        self.temp = 273.0

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
    system_of_equations = (
        np.dot(rad_matrix, node_temps**4)
        + np.dot(cond_matrix, node_temps)
        + q_extra
    )
    # for i in range(0, len(node_temperatures)):
    # system_of_equations.append(np.dot(rad_matrix[i], node_temperatures**4) + np.dot(cond_matrix[i], node_temperatures) + q_extra[i])
    return system_of_equations

def heat(nodes_examined, R_ij, G_ij):
    q_rad = R_ij*(nodes_examined[1].temp**4 - nodes_examined[0].temp**4)
    q_cond = G_ij*(nodes_examined[1].temp - nodes_examined[0].temp)
    q_nodes = q_rad + q_cond
    return q_nodes


def temp_update(node_temperatures, rad_matrix, cond_matrix, capacities, q_extra):
    q_rad = np.dot(rad_matrix, node_temperatures**4)
    q_cond = np.dot(cond_matrix, node_temperatures)
    q_total = (q_rad + q_cond + q_extra)
    return np.divide(q_total, capacities)

def runge_kutta(node_temperatures, time_step, capacities, rad_matrix, cond_matrix, q_extra):
    k1 = time_step * temp_update(node_temperatures, rad_matrix, cond_matrix, capacities, q_extra)
    k2 = time_step * temp_update(node_temperatures + (k1/2), rad_matrix, cond_matrix, capacities, q_extra)
    k3 = time_step * temp_update(node_temperatures + (k2/2), rad_matrix, cond_matrix, capacities, q_extra)
    k4 = time_step * temp_update(node_temperatures + k3, rad_matrix, cond_matrix, capacities, q_extra)
    k = (k1 + (2*k2) + (2*k3) + k4) / 6
    node_temperatures = node_temperatures + k
    return node_temperatures



def steady_state(nodes, relationships, sail_deployed, duty_cycle, thermal_case):
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
    nodes, relationships, dt, sail_deployed, thermal_case, heater_power=400, electrical_power=800
):
    """
    Computes the equilibrium temperatures for each node.

    For the shield, it does this using the concept of thermal resistance. In the absence
    of shielding, the heat added to the spacecraft wall is equal to the inputted solar heat
    plus the internal heat of the spacecraft, while the heat lost is from radiation to deep
    space.

    q_internal + solar_q_in = A * e * sigma * T_wall^4.

    T_wall can be solved in this equation to get the wall equilibrium temperature.

    When shielding is added, there is also an interaction between the shield wall and the
    spacecraft wall. This shielding reduces the solar heat felt by the spacecraft wall by
    a factor of:

    1 / ((n_shields / e_shield) + (1 / e_wall) - n_shields)

    Information for this part of the model was modified from information about the radiation
    felt between a series of parallel plates. It had to be modified as the sun is not a parallel
    plate. However, the heat shield and spacecraft were assumed to be parallel plates for this
    analysis (which may not be true in practice).

    Sources:
        https://thermopedia.com/content/69/
        https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node136.html
        https://www.engineeringenotes.com/thermal-engineering/heat-transfer/use-of-radiation-shields-to-reduce-heat-transfer-thermal-engineering/30727

    Inputs:
    - bus [class]: The Node object corresponding to the spacecraft bus.
    - sail [class]: The Node object corresponding to the solar sail.
    - shield [class]: The Node object corresponding to the heat shield

    Outputs:

        self.q_absorbed - self.q_radiated [float]: The thermal balance of the node.
    """
    t_start = time.perf_counter()
    dt_original = dt
    node_temperatures = np.asarray([nodes[i].temp for i in range(0, len(nodes))])
    radiative = np.triu(relationships, 1)
    radiative = np.where(radiative, radiative, radiative.T)
    conductive = np.tril(relationships, -1)
    conductive = np.where(conductive, conductive, conductive.T)
    capacities = np.diagonal(relationships)
    panel_duty_cycle = 0.01
    payload_duty_cycle = 0

    q_extra = []
    q_rad_coeff = []
    for i in range(0, len(nodes)):
        q_extra.append(
            nodes[i].solar_heat_in(thermal_case, sail_deployed, panel_duty_cycle, payload_duty_cycle)
            + nodes[i].internal_heat
        )
        q_rad_coeff.append(nodes[i].radiated_heat_coeff(sail_deployed))

    rad_matrix = np.diag(q_rad_coeff - np.sum(radiative, axis=1)) + radiative
    cond_matrix = conductive + np.diag(-np.sum(conductive, axis=1))
    
    dt_interp = np.min(0.001*capacities)
    if dt_interp > 0:
        dt_arr = [dt_interp]*int(dt_original/dt_interp)
    else:
        dt_arr = [0]

    # integrator = integrate.RK45(fun = lambda self, node_temps: temp_update(node_temps, capacities, rad_matrix, cond_matrix, q_extra), t0 = dt_interp, y0 = node_initial, t_bound = dt_original, vectorized=True)
    for idx, time_step in enumerate(dt_arr):
        # q_rad = np.dot(rad_matrix, node_temperatures**4)
        # q_cond = np.dot(cond_matrix, node_temperatures)
        # q_total = q_rad + q_cond + q_extra
        node_temperatures = runge_kutta(node_temperatures, time_step, capacities, rad_matrix, cond_matrix, q_extra)
        # node_temperatures = node_temperatures + np.divide(time_step, capacities) * q_total
        # node_temperatures = node_temperatures + np.divide(time_step, capacities) * q_extra
        # q_cond = np.dot(cond_matrix, node_temperatures)
        # node_temperatures = node_temperatures + np.divide(time_step, capacities) * q_cond
        # q_rad = np.dot(rad_matrix, node_temperatures**4)
        # node_temperatures = node_temperatures + np.divide(time_step, capacities) * q_rad
        
    for i in range(0, len(nodes)):
        nodes[i].temp = node_temperatures[i]

    disp_temps = np.round(node_temperatures - 273.15, 2)
    dist = np.round(thermal_case[0], 2)

    t_stop = time.perf_counter()

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
    print(f"Runtime: {round(t_stop - t_start, 2)} s")
    return node_temperatures
