'''
thermal_model.py

This code presents a callable version of thermal_script.py that the orbital simulation can
use to properly optimize the orbit.

-TODO:
    - Containerize outputs as dictionary containing each node and whether it has passed or 
      failed the thermal constraints.
    - Also consider sending temperature outputs to the orbital simulation as well.
'''

import node as nd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm
from budget_properties import Properties
import materials
import config

class Thermal():
    def __init__(self): 
        '''
        Necessary inputs:

            - Thermal case
            - Time step
        '''
        self.AU = 1.496e11
        self.SB = 5.67e-8
        self.yearInSeconds = 365 * 24 * 3600
        self.solar_luminosity = 3.83e26
        initial_node_temp = 300

        bus_materials = []
        self.bus_properties = []
        sail_materials = []
        self.sail_properties = []
        self.total_properties = []

        self.shield_layers = config.shield['layers']

        for i, node in enumerate(config.bus_nodes):
            bus_materials.append(materials.bus_material(node['external']))
            self.bus_properties.append([bus_materials[i]['emissivity'], bus_materials[i]['reflectivity'], 
                          bus_materials[i]['absorptivity'], bus_materials[i]['density'], 
                          node['area'], node['sun_vf'], node['space_vf'], node['temp_range']])
            
        self.total_properties.append(self.bus_properties)
        for i, node in enumerate(config.sail_nodes):
            sail_materials.append(materials.sail_material(node['external']))
            self.sail_properties.append([sail_materials[i]['emissivity'], sail_materials[i]['reflectivity'], 
                            sail_materials[i]['absorptivity'], sail_materials[i]['density'], 
                            node['area'], node['sun_vf'], node['space_vf'], node['temp_range']])

        self.total_properties.append(self.sail_properties)
        panel = materials.panel_material(config.solar_panels['external'])
        self.panel_properties = [panel['emissivity'], panel['reflectivity'], panel['absorptivity'], 
                    panel['density'], config.solar_panels['area'], config.solar_panels['sun_vf'],
                    config.solar_panels['space_vf'], config.solar_panels['temp_range']]
        
        self.total_properties.append(self.panel_properties)
        boom = materials.boom_material(config.booms['external'])
        self.boom_properties = [boom['emissivity'], boom['reflectivity'], boom['absorptivity'], boom['density'], 
                   config.booms['area'], config.booms['sun_vf'], config.booms['space_vf'], config.booms['temp_range']]
        self.total_properties.append(self.boom_properties)
        
        self.total_properties.append(self.boom_properties)
        shield = materials.shield_material(config.shield['external'])
        self.shield_properties = [shield['emissivity'], shield['reflectivity'], shield['absorptivity'],
                     shield['density'], config.shield['area'], config.shield['sun_vf'], 
                     config.shield['space_vf'], config.shield['temp_range']]
        self.total_properties.append(self.shield_properties)

        self.relationships = np.asarray(config.node_relationship)
        self.total_nodes = len(self.relationships)

        self.spacecraft = []
        self.spacecraft_bus = []
        self.solar_sail = [] 

        for node, node_properties in enumerate(self.bus_properties):
            self.spacecraft_bus.append(nd.Node("Spacecraft Bus", node_properties))
        self.spacecraft.append(self.spacecraft_bus)
        for node, node_properties in enumerate(self.sail_properties):
            self.solar_sail.append(nd.Node("Solar Sail", node_properties))
        self.spacecraft.append(self.solar_sail)
        self.spacecraft = [item for sublist in self.spacecraft for item in sublist]
        solar_panel = nd.Node("Solar Panel", self.panel_properties)
        self.spacecraft.append(solar_panel)
        boom = nd.Node("Boom", self.boom_properties)
        self.spacecraft.append(boom)
        if self.shield_layers > 0:
            heat_shield = nd.Node("Heat Shield", self.shield_properties, self.shield_layers)
        else:
            heat_shield = None
        self.spacecraft.append(heat_shield)


        self.node_keys = [config.nodes[i]['name'] for i in range(0, len(config.nodes))]
        self.node_temp_ranges = [config.nodes[i]['temp_range'] for i in range(0, len(config.nodes))]
        self.node_fail_step = [False]*len(self.node_keys)

        self.node_failure = []
        self.node_temperatures = []


        
    def step(self):
        '''
        Steps forward in time and runs the thermal nodal model again.

        Necessary inputs:
            - time_step: time step in seconds.
            - thermal_case [arr]: includes distance from the sun and angle to the sun as [dist, angle].
        '''
        for node in self.spacecraft:
            if node != None:
                flux_in = node.solar_heat_in(self.thermal_case)
        node_temp_step = nd.temperatures(self.spacecraft, self.relationships, self.time_step, self.thermal_case)
        node_fail_step = [False]*self.total_nodes

        for idx, temp in enumerate(node_temp_step):
            if temp not in range(self.node_temp_ranges[idx][0], self.node_temp_ranges[idx][1]):
                node_fail_step[idx] = True
        
        self.node_failure.append(self.node_fail_step)
        self.node_temperatures.append(node_temp_step)
        

    def stop_propagation(self, time_step, thermal_case):
        self.time_step = time_step
        self.thermal_case = thermal_case
        self.step()
        if all(self.node_fail_step) == False:
            return False
        else:
            return True
        
    def optimize_output(self):
        self.temps_dict = {}
        self.fail_dict = {}
        self.temps_dict = dict(zip(self.node_keys, self.node_temperatures))
        self.fail_dict = dict(zip(self.node_keys, self.node_failure))
        return self.fail_dict, self.temps_dict