"""
thermal_model.py

This code presents a callable version of thermal_script.py that the orbital simulation can
use to properly optimize the orbit.

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
    def __init__(self, dt=None, verbose=False):
        """
        Initializes the thermal model with the spacecraft nodes and their properties.

        Necessary inputs:

        dt (float): The time step as dictated by the orbital model.

        """
        self.AU = 1.496e11
        self.SB = 5.67e-8
        self.yearInSeconds = 365 * 24 * 3600
        self.solar_luminosity = 3.83e26

        bus_materials = []
        self.bus_properties = []
        self.sail_properties = []
        self.total_properties = []
        self.shield_properties = []

        self.shield_layers = config.shield_inner["layers"]

        # Set up the spacecraft nodes and their properties

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

        hga = materials.antenna_material(config.antenna["external"])
        self.antenna_properties = [
            hga["emissivity"],
            hga["reflectivity"],
            hga["absorptivity"],
            hga["density"],
            config.antenna["area"],
            config.antenna["sun_vf"],
            config.antenna["space_vf"],
            config.antenna["temp_range"],
            config.antenna["internal_heat"],
        ]

        self.total_properties.append(self.antenna_properties)
        
        front_shield = materials.shield_material(config.shield_front["external"])
        self.shield_properties = [
            front_shield["emissivity"],
            front_shield["reflectivity"],
            front_shield["absorptivity"],
            front_shield["density"],
            config.shield_front["area"],
            config.shield_front["sun_vf"],
            config.shield_front["space_vf"],
            config.shield_front["temp_range"],
            config.shield_front["internal_heat"],
        ]
        self.total_properties.append(self.shield_properties)

        inner_shield = materials.shield_material(config.shield_inner["external"])
        self.innershield_properties = [
            inner_shield["emissivity"],
            inner_shield["reflectivity"],
            inner_shield["absorptivity"],
            inner_shield["density"],
            config.shield_inner["area"],
            config.shield_inner["sun_vf"],
            config.shield_inner["space_vf"],
            config.shield_inner["temp_range"],
            config.shield_inner["internal_heat"],
        ]
        self.total_properties.append(self.innershield_properties)

        battery = materials.internal_material(config.batteries["external"])
        self.battery_properties = [
            battery["emissivity"],
            battery["reflectivity"],
            battery["absorptivity"],
            battery["density"],
            config.batteries["area"],
            config.batteries["sun_vf"],
            config.batteries["space_vf"],
            config.batteries["temp_range"],
            config.batteries["internal_heat"],
        ]
        self.total_properties.append(self.battery_properties)

        propellant = materials.internal_material(config.hydrazine["external"])
        self.propellant_properties = [
            propellant["emissivity"],
            propellant["reflectivity"],
            propellant["absorptivity"],
            propellant["density"],
            config.hydrazine["area"],
            config.hydrazine["sun_vf"],
            config.hydrazine["space_vf"],
            config.hydrazine["temp_range"],
            config.hydrazine["internal_heat"],
        ]
        self.total_properties.append(self.propellant_properties)

        coronagraph = materials.internal_material(config.metis["external"])
        self.coronagraph_properties = [
            coronagraph["emissivity"],
            coronagraph["reflectivity"],
            coronagraph["absorptivity"],
            coronagraph["density"],
            config.metis["area"],
            config.metis["sun_vf"],
            config.metis["space_vf"],
            config.metis["temp_range"],
            config.metis["internal_heat"],
        ]
        self.total_properties.append(self.coronagraph_properties)

        doppler = materials.internal_material(config.cdm["external"])
        self.doppler_properties = [
            doppler["emissivity"],
            doppler["reflectivity"],
            doppler["absorptivity"],
            doppler["density"],
            config.cdm["area"],
            config.cdm["sun_vf"],
            config.cdm["space_vf"],
            config.cdm["temp_range"],
            config.cdm["internal_heat"],
        ]
        self.total_properties.append(self.doppler_properties)

        self.spacecraft = []
        self.spacecraft_bus = []

        for node, node_properties in enumerate(self.bus_properties):
            self.spacecraft_bus.append(nd.Node("Spacecraft Bus", node_properties))
            self.spacecraft.append(self.spacecraft_bus[-1])
        self.solar_sail = nd.Node("Solar Sail", self.sail_properties)
        self.spacecraft.append(self.solar_sail)
        boom = nd.Node("Boom", self.boom_properties)
        self.spacecraft.append(boom)
        self.solar_panel = nd.Node("Solar Panel", self.panel_properties)
        self.spacecraft.append(self.solar_panel)
        antenna = nd.Node("Antenna", self.antenna_properties)
        self.spacecraft.append(antenna)
        if self.shield_layers > 0:
            self.heat_shield = nd.Node("Heat Shield", self.shield_properties, self.shield_layers)
            self.shield_inner = nd.Node("Heat Shield", self.innershield_properties, self.shield_layers)
            self.spacecraft.append(self.heat_shield)
            self.spacecraft.append(self.shield_inner)
        else:
            self.heat_shield = None
        self.batteries = nd.Node("Batteries", self.battery_properties)
        self.hydrazine = nd.Node("Propellant", self.propellant_properties)
        self.metis = nd.Node("Coronagraph", self.coronagraph_properties)
        self.cdm = nd.Node("Doppler Magnetograph", self.doppler_properties)
        self.spacecraft.append(self.batteries)
        self.spacecraft.append(self.hydrazine)
        self.spacecraft.append(self.metis)
        self.spacecraft.append(self.cdm)

        self.node_keys = [config.nodes[i]["name"] for i in range(0, len(config.nodes))]
        self.node_temp_ranges = [
            config.nodes[i]["temp_range"] for i in range(0, len(config.nodes))
        ]

        self.node_failure = []
        self.node_temperatures = []
        self.node_temp_state = [0 for _ in self.node_keys]
        self.start_time = False

        if dt is not None:
            self.customDT = True
            self.dt = dt
        else:
            self.customDT = False

        self.thermalFailure = False

        self.verbose = verbose

    def __repr__(self) -> str:
        return "Thermal V1"
    
    def __str__(self) -> str:
        return "Thermal V1"


    def step(self, current_time, alt, coneAngle):
        """
        Steps forward in time and runs the thermal nodal model.

        Necessary inputs:
            - current_time (int): current time in seconds.
            - alt (float): distance from the sun, inputted in AU.
            - coneAngle (float): cone angle of the solar sail from the sun, in rad.
        """
        self.total_nodes = len(self.spacecraft)

        # Used to limit the amount of simulations that are completed.
        if alt < 1.0:
            if self.start_time == False:
                self.start_time = current_time
                dt = 0
                self.node_fail_step = [False] * len(self.spacecraft)
                self.node_initial = np.asarray([self.spacecraft[i].temp for i in range(0, len(self.spacecraft))])
                self.node_failure.append(self.node_fail_step)
                self.node_temperatures.append(self.node_initial)
            else:
                if self.customDT:
                    dt = self.dt
                else:
                    dt = current_time - self.start_time

                self.start_time = current_time

                sail_deployed = 1
                self.relationships = np.asarray(config.node_relationship_deployed)

                node_temp_step = nd.time_variant(
                    self.spacecraft, self.relationships, dt, sail_deployed, [alt, coneAngle],
                    self.node_temp_ranges, verbose=self.verbose, plot=False
                )
                self.node_fail_step = [False] * self.total_nodes

                # Check if temperatures fall within allowable ranges.
                for idx, temp in enumerate(node_temp_step):
                    tempC = temp - 273.15
                    if (tempC < self.node_temp_ranges[idx][0]) or (tempC > self.node_temp_ranges[idx][1]):
                        print(self.node_keys[idx], round(tempC, 2), self.node_temp_ranges[idx])
                        self.node_fail_step[idx] = True
                self.node_temp_state = list(node_temp_step)
        else:
            self.node_fail_step = [False] * self.total_nodes
            node_temp_step = [273.0]*self.total_nodes
            self.node_failure.append(self.node_fail_step)
            self.node_temperatures.append(node_temp_step)
            self.node_temp_state = list(node_temp_step)


    def stopPropagation(self, time_step):
        """
        Checks if the thermal boundaries are broken, and if the propagation needs to be stopped.

        Necessary inputs:
            - time_step (int): The orbital time step
        Returns:
            - True if thermal boundaries are broken and propagation needs to be stopped, False otherwise.
        """
        if sum(self.node_fail_step) > 0:
            print(f"Stopping Propagation => Heat")
            print(self.node_fail_step)
            self.thermalFailure = True
            return True
        else:
            return False

    def optimize_output(self):
        """
        
        Outputs of the thermal model for the optimizer code.

        Necessary inputs:
            - time_step (int): The orbital time step

        Returns:
            - self.temps_dict (dict): dictionary containing each node and their temperatures over the simulation.
            - self.fail_dict (dict): dictionary containing each node and if their limits have been violated.
        """
        self.temps_dict = dict(zip(self.node_keys, self.node_temperatures))
        self.fail_dict = dict(zip(self.node_keys, self.node_failure))
        return self.fail_dict, self.temps_dict
