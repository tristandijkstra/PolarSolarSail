"""
thermal_model.py

This code presents a callable version of thermal_script.py that the orbital simulation can
use to properly optimize the orbit.

-TODO:
    - Containerize outputs as dictionary containing each node and whether it has passed or 
      failed the thermal constraints.
    - Also consider sending temperature outputs to the orbital simulation as well.
"""

# import node as nd
import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import logging
# from tqdm import tqdm
from .budget_properties import Properties
from . import materials
from . import config
from . import node as nd

class Thermal:
    def __init__(self):
        """
        Necessary inputs:

        """
        self.AU = 1.496e11
        self.SB = 5.67e-8
        self.yearInSeconds = 365 * 24 * 3600
        self.solar_luminosity = 3.83e26

        bus_materials = []
        self.bus_properties = []
        sail_materials = []
        self.sail_properties = []
        self.total_properties = []
        shield_materials = []
        self.shield_properties = []

        self.shield_layers = config.shield["layers"]

        for i, node in enumerate(config.bus_nodes):
            bus_materials.append(materials.bus_material(node["external"]))
            self.bus_properties.append(
                [
                    bus_materials[i]["emissivity"],
                    bus_materials[i]["reflectivity"],
                    bus_materials[i]["absorptivity"],
                    bus_materials[i]["density"],
                    node["area"],
                    node["sun_vf"],
                    node["space_vf"],
                    node["temp_range"],
                    node["internal_heat"],
                ]
            )

        self.total_properties.append(self.bus_properties)

        solar_sail = materials.sail_material(config.sail["external"])
        self.sail_properties = [
            solar_sail["emissivity"],
            solar_sail["reflectivity"],
            solar_sail["absorptivity"],
            solar_sail["density"],
            config.sail["area"],
            config.sail["sun_vf"],
            config.sail["space_vf"],
            config.sail["temp_range"],
            config.sail["internal_heat"],
        ]

        self.total_properties.append(self.sail_properties)
        boom = materials.boom_material(config.booms["external"])
        self.boom_properties = [
            boom["emissivity"],
            boom["reflectivity"],
            boom["absorptivity"],
            boom["density"],
            config.booms["area"],
            config.booms["sun_vf"],
            config.booms["space_vf"],
            config.booms["temp_range"],
            config.booms["internal_heat"],
        ]
        self.total_properties.append(self.boom_properties)

        panel = materials.panel_material(config.solar_panels["external"])
        self.panel_properties = [
            panel["emissivity"],
            panel["reflectivity"],
            panel["absorptivity"],
            panel["density"],
            config.solar_panels["area"],
            config.solar_panels["sun_vf"],
            config.solar_panels["space_vf"],
            config.solar_panels["temp_range"],
            config.solar_panels["internal_heat"],
        ]

        self.total_properties.append(self.panel_properties)
        heat_shield = materials.shield_material(config.shield["external"])
        self.shield_properties = [
            heat_shield["emissivity"],
            heat_shield["reflectivity"],
            heat_shield["absorptivity"],
            heat_shield["density"],
            config.shield["area"],
            config.shield["sun_vf"],
            config.shield["space_vf"],
            config.shield["temp_range"],
            config.shield["internal_heat"],
        ]
        self.total_properties.append(self.shield_properties)


        self.spacecraft = []
        self.spacecraft_bus = []

        for node, node_properties in enumerate(self.bus_properties):
            self.spacecraft_bus.append(nd.Node("Spacecraft Bus", node_properties))
            self.spacecraft.append(self.spacecraft_bus[-1])
        self.solar_sail = nd.Node("Solar Sail", self.sail_properties)
        self.spacecraft.append(self.solar_sail)
        boom = nd.Node("Boom", self.boom_properties)
        self.spacecraft.append(boom)
        solar_panel = nd.Node("Solar Panel", self.panel_properties)
        self.spacecraft.append(solar_panel)
        if self.shield_layers > 0:
            self.heat_shield = nd.Node("Heat Shield", self.shield_properties, self.shield_layers)
            self.spacecraft.append(self.heat_shield)
        else:
            self.heat_shield = None

        self.node_keys = [config.nodes[i]["name"] for i in range(0, len(config.nodes))]
        self.node_temp_ranges = [
            config.nodes[i]["temp_range"] for i in range(0, len(config.nodes))
        ]

        self.node_failure = []
        self.node_temperatures = []
        self.start_time = False


    def __repr__(self) -> str:
        return "Thermal V1"
    
    def __str__(self) -> str:
        return "Thermal V1"


    def step(self, current_time, alt, coneAngle):
        """
        Steps forward in time and runs the thermal nodal model again.

        Necessary inputs:
            - time_step: time step in seconds.
            - thermal_case [arr]: includes distance from the sun and angle to the sun as [dist, angle].
        """
        if self.start_time == False:
            self.start_time = current_time
            dt = 0
            self.node_fail_step = [False] * len(self.spacecraft)
            self.node_initial = np.asarray([self.spacecraft[i].temp for i in range(0, len(self.spacecraft))])
            self.node_failure.append(self.node_fail_step)
            self.node_temperatures.append(self.node_initial)
        else:
            dt = current_time - self.start_time
            self.start_time = current_time

            self.total_nodes = len(self.spacecraft)

            sail_deployed = 0
            self.relationships = np.asarray(config.node_relationship_deployed)

            # node_temp_step = nd.steady_state(self.spacecraft, self.relationships, dt,
            #                                  sail_deployed, [alt, coneAngle])
            node_temp_step = nd.time_variant(
                self.spacecraft, self.relationships, dt, sail_deployed, [alt, coneAngle]
            )
            self.node_fail_step = [False] * self.total_nodes

            for idx, temp in enumerate(node_temp_step):
                if (temp < self.node_temp_ranges[idx][0]) and (temp > self.node_temp_ranges[idx][1]):
                    self.node_fail_step[idx] = True

            self.node_failure.append(self.node_fail_step)
            self.node_temperatures.append(node_temp_step)


    def stopPropagation(self, time_step):
        if any(self.node_fail_step) == True:
            print(f"Stopping Propagation => Heat")
            print(self.node_failure[-1])
            return True
        else:
            return False

    def optimize_output(self):
        self.temps_dict = dict(zip(self.node_keys, self.node_temperatures))
        self.fail_dict = dict(zip(self.node_keys, self.node_failure))
        return self.fail_dict, self.temps_dict
