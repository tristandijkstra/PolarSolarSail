'''
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

'''

import pandas as pd
import numpy as np

SB = 5.67e-8
AU = 1.496e11
solar_luminosity = 3.83e26


class Node():
    def __init__(self, name, properties, n_shield=0):
        '''
        This sets up each Node class with the surface properties and the area.
        
        Inputs:
        
        name [str]: The name of the node
        skin_properties [arr]: Input an array consisting of [emissivity, absorptivity, reflectivity, radiator_emissivity], all scaled from 0 to 1.
        area [arr]: Input an array consisting of [area that the sun sees, area that deep space sees, area covered by radiators.]
        case [arr]: Input an array consisting of [distance in AU, angle from sun in rad]
        n_shield [int]: The number of layers for the heat shield.
        '''
        self.name = name
        self.emissivity = properties[0]
        self.absorptivity = properties[1]
        self.reflectivity = properties[2]
        self.density = properties[3]
        self.area = properties[4]
        self.sun_vf = properties[5]
        self.space_vf = properties[6]
        self.temp_range = properties[7]
        self.n_shield = n_shield
        self.temp = 300
    
    def solar_heat_in(self, thermal_case):
        '''
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
        '''
        self.solar_flux = (np.cos(thermal_case[1]) * solar_luminosity * self.absorptivity * self.sun_vf) / (4 * np.pi * (thermal_case[0] * AU)**2)
        self.solar_q_in = self.solar_flux * self.area
        return self.solar_q_in
    
    def internal_heat(self, heater_power, electrical_heat):
        '''
        This computes the heat generated inside the spacecraft.
        Inputs:
            - heater_power [W]
            - electrical_heat [W]
        
        Outputs:
        
            self.q_internal [float]: Internal heat generated by the spacecraft in W
        '''
        self.q_internal = heater_power + electrical_heat
        return self.q_internal
    
    def heat_radiated(self):
        '''
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
        '''
        self.q_radiated = self.emissivity * self.space_vf * SB * self.area * self.temp**4
        return self.q_radiated

def heat(nodes, R_ij, G_ij, solar_q, internal_q, radiated_q):
    q_rad = R_ij * SB * (nodes[1].temp**4 - nodes[0].temp**4)
    q_cond = G_ij * (nodes[1].temp - nodes[0].temp)
    q_total = q_rad + q_cond + solar_q + internal_q - radiated_q
    return q_total

def temperatures(nodes, relationships, dt, thermal_case, heater_power=400, electrical_power=800):
    '''
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
        '''
    node_temperatures = [nodes[i].temp for i in range(0, len(nodes))]
    radiative = np.triu(relationships, 1)
    radiative = np.where(radiative, radiative, radiative.T)
    conductive = np.tril(relationships, 1)
    conductive = np.where(conductive, conductive, conductive.T)
    capacities = np.diagonal(relationships)
    for i in range(0, len(nodes)):
        q_total = 0
        for j in range(0, len(nodes)):
            if i != j:
                capacity = capacities[i]
                G_ij = conductive[i, j]
                R_ij = radiative[i, j]
                solar_q = nodes[i].solar_heat_in(thermal_case)
                internal_q = nodes[i].internal_heat(heater_power, electrical_power)
                radiated_q = nodes[i].heat_radiated()
                nodes_examined = [nodes[i], nodes[j]]
                q_total += heat(nodes_examined, R_ij, G_ij, solar_q, internal_q, radiated_q)
        
        node_temperatures[i] = node_temperatures[i] + (dt / capacity) * q_total
        nodes[i].temp = node_temperatures[i]
        
    return node_temperatures

