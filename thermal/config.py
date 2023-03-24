import numpy as np
import pandas as pd
from . import materials

"""
config.py

This file contains the inputs for the thermal model. The inputs will be automatically logged.

TODO:
    - Input correct values for configuration matrix.
"""

### General Parameters
spacecraft_radius = 1.5
spacecraft_length = 5
shell_thickness = 0.002
sail_thickness = 0.000004
sail_density = 1540
boom_length = 70
boom_radius = 0.029
boom_thickness = 0.000011
sail_area = 10000
panel_area = 18
heaters_power = 0
coolers_power = 0
electrical_heat = 800
log_file = "data/thermal.log"

### Structure Nodes
node_1 = {
    "name": "Spacecraft +X",
    "area": np.pi * spacecraft_radius * spacecraft_length,
    "sun_vf": 0.5,
    "space_vf": 0.5,
    "external": "multi-layer insulation",
    "temp_range": [-100, 100],
}
node_2 = {
    "name": "Spacecraft -X",
    "area": np.pi * spacecraft_radius * spacecraft_length,
    "sun_vf": 0.5,
    "space_vf": 0.5,
    "external": "az93 white paint",
    "temp_range": [-100, 100],
}
node_3 = {
    "name": "Spacecraft +Z",
    "area": np.pi * spacecraft_radius**2,
    "sun_vf": 0,
    "space_vf": 0,
    "external": "multi-layer insulation",
    "temp_range": [-100, 100],
}
node_4 = {
    "name": "Spacecraft -Z",
    "area": np.pi * spacecraft_radius**2,
    "sun_vf": 0,
    "space_vf": 1,
    "external": "multi-layer insulation",
    "temp_range": [-100, 100],
}
bus_nodes = [node_1, node_2, node_3, node_4]

### Sail Nodes
node_front = {
    "name": "Sail Front",
    "area": sail_area,
    "sun_vf": 1,
    "space_vf": 0,
    "external": "aluminum",
    "temp_range": [-100, 100],
}
node_back = {
    "name": "Sail Back",
    "area": sail_area,
    "sun_vf": 0,
    "space_vf": 1,
    "external": "chromium",
    "temp_range": [-100, 100],
}
sail_nodes = [node_front, node_back]

### External Nodes
booms = {
    "name": "Booms",
    "area": 4 * np.pi * (boom_radius) * boom_length,
    "sun_vf": 0.5,
    "space_vf": 0.5,
    "external": "carbon fibre",
    "temp_range": [-100, 100],
}
solar_panels = {
    "name": "Solar Panels",
    "area": panel_area,
    "sun_vf": 1,
    "space_vf": 0,
    "external": "osr",
    "temp_range": [-100, 100],
}
shield = {
    "name": "Heat Shield",
    "area": node_1["area"],
    "sun_vf": 1,
    "space_vf": 0,
    "external": "ceramic cloth",
    "temp_range": [-100, 100],
    "layers": 7,
}


nodes = [
    node_1,
    node_2,
    node_3,
    node_4,
    node_front,
    node_back,
    booms,
    solar_panels,
    shield,
]


### Node Relationships Matrix.
"""
    The setup for this matrix is as follows.
    i/j     1       2       3       4   
    1     c1*m1    F_12    F_13    F_14
    2      G_21    c2*m2   F_23    F_24
    3      G_31    G_32    c3*m3   F_34
    4      G_41    G_42    G_43   c4*m4  

    Here, c_i*m_i are the thermal capacity of each node (the specific heat capacity c times the mass m.)
    Then, G_ij are the conductive coefficients, equal to k_ij * A_ij / L_ij
    Then, R_ij are the radiative coefficients, equal to e_ij * A_ij * F_ij.
"""
densities = []
capacities = []
conductivities = []
emissivities = []

for node in bus_nodes:
    densities.append(materials.bus_material(node["external"])["density"])
    capacities.append(materials.bus_material(node["external"])["specific_heat"])
    conductivities.append(materials.bus_material(node["external"])["conductivity"])
    emissivities.append(materials.bus_material(node["external"])["emissivity"])

for node in sail_nodes:
    densities.append(materials.sail_material(node["external"])["density"])
    capacities.append(materials.sail_material(node["external"])["specific_heat"])
    conductivities.append(materials.sail_material(node["external"])["conductivity"])
    emissivities.append(materials.sail_material(node["external"])["emissivity"])


densities.append(materials.boom_material(booms["external"])["density"])
capacities.append(materials.boom_material(booms["external"])["specific_heat"])
conductivities.append(materials.boom_material(booms["external"])["conductivity"])
emissivities.append(materials.boom_material(booms["external"])["emissivity"])

densities.append(materials.panel_material(solar_panels["external"])["density"])
capacities.append(materials.panel_material(solar_panels["external"])["specific_heat"])
conductivities.append(
    materials.panel_material(solar_panels["external"])["conductivity"]
)
emissivities.append(materials.panel_material(solar_panels["external"])["emissivity"])

densities.append(materials.shield_material(shield["external"])["density"])
capacities.append(materials.shield_material(shield["external"])["specific_heat"])
conductivities.append(materials.shield_material(shield["external"])["conductivity"])
emissivities.append(materials.shield_material(shield["external"])["emissivity"])

masses = [
    node_1["area"] * shell_thickness * densities[0],
    node_2["area"] * shell_thickness * densities[1],
    node_3["area"] * shell_thickness * densities[2],
    node_4["area"] * shell_thickness * densities[3],
    node_front["area"] * sail_thickness * sail_density * 0.5,
    node_back["area"] * sail_thickness * sail_density * 0.5,
    densities[6]
    * np.pi
    * ((boom_radius) ** 2 - (boom_radius - boom_thickness) ** 2)
    * boom_length,
    densities[7] * solar_panels["area"],
    densities[8] * shield["area"],
]

view_factors = [
    [0, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1],
    [0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1],
    [0, 0, 0, 0, 0, 0, 0.1, 0.1, 0.02],
    [0, 0, 0, 0, 0, 0, 0.25, 0, 0.02],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0.1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]

radiative_area = [
    [
        0,
        node_1["area"],
        node_1["area"],
        node_1["area"],
        node_1["area"],
        node_1["area"],
        node_1["area"],
        node_1["area"],
        node_1["area"],
    ],
    [
        0,
        0,
        node_2["area"],
        node_2["area"],
        node_2["area"],
        node_2["area"],
        node_2["area"],
        node_2["area"],
        node_2["area"],
    ],
    [
        0,
        0,
        0,
        node_3["area"],
        node_3["area"],
        node_3["area"],
        node_3["area"],
        node_3["area"],
        node_3["area"],
    ],
    [
        0,
        0,
        0,
        0,
        node_4["area"],
        node_4["area"],
        node_4["area"],
        node_4["area"],
        node_4["area"],
    ],
    [
        0,
        0,
        0,
        0,
        0,
        node_front["area"],
        node_front["area"],
        node_front["area"],
        node_front["area"],
    ],
    [0, 0, 0, 0, 0, 0, node_back["area"], node_back["area"], node_back["area"]],
    [0, 0, 0, 0, 0, 0, 0, booms["area"], booms["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, solar_panels["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]

effective_emissivities = [
    [
        0,
        1 / ((1 / emissivities[0]) + (1 / emissivities[1])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[2])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[3])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[4])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[5])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[6])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[0]) + (1 / emissivities[8])),
    ],
    [
        0,
        0,
        1 / ((1 / emissivities[1]) + (1 / emissivities[2])),
        1 / ((1 / emissivities[1]) + (1 / emissivities[3])),
        1 / ((1 / emissivities[1]) + (1 / emissivities[4])),
        1 / ((1 / emissivities[1]) + (1 / emissivities[5])),
        1 / ((1 / emissivities[1]) + (1 / emissivities[6])),
        1 / ((1 / emissivities[1]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[1]) + (1 / emissivities[8])),
    ],
    [
        0,
        0,
        0,
        1 / ((1 / emissivities[2]) + (1 / emissivities[3])),
        1 / ((1 / emissivities[2]) + (1 / emissivities[4])),
        1 / ((1 / emissivities[2]) + (1 / emissivities[5])),
        1 / ((1 / emissivities[2]) + (1 / emissivities[6])),
        1 / ((1 / emissivities[2]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[2]) + (1 / emissivities[8])),
    ],
    [
        0,
        0,
        0,
        0,
        1 / ((1 / emissivities[3]) + (1 / emissivities[4])),
        1 / ((1 / emissivities[3]) + (1 / emissivities[5])),
        1 / ((1 / emissivities[3]) + (1 / emissivities[6])),
        1 / ((1 / emissivities[3]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[3]) + (1 / emissivities[8])),
    ],
    [
        0,
        0,
        0,
        0,
        0,
        1 / ((1 / emissivities[4]) + (1 / emissivities[5])),
        1 / ((1 / emissivities[4]) + (1 / emissivities[6])),
        1 / ((1 / emissivities[4]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[4]) + (1 / emissivities[8])),
    ],
    [
        0,
        0,
        0,
        0,
        0,
        0,
        1 / ((1 / emissivities[5]) + (1 / emissivities[6])),
        1 / ((1 / emissivities[5]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[5]) + (1 / emissivities[8])),
    ],
    [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1 / ((1 / emissivities[6]) + (1 / emissivities[7])),
        1 / ((1 / emissivities[6]) + (1 / emissivities[8])),
    ],
    [0, 0, 0, 0, 0, 0, 0, 0, 1 / ((1 / emissivities[7]) + (1 / emissivities[8]))],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]

radiative_relationship = np.multiply(
    np.multiply(view_factors, radiative_area), effective_emissivities
)

conductive_coeff = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [400, 400, 0, 0, 0, 0, 0, 0, 0],
    [400, 400, 400, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]

conductive_area = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [
        0.5
        * np.pi
        * (spacecraft_radius**2 - (spacecraft_radius - shell_thickness) ** 2),
        0.5
        * np.pi
        * (spacecraft_radius**2 - (spacecraft_radius - shell_thickness) ** 2),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ],
    [
        0.5
        * np.pi
        * (spacecraft_radius**2 - (spacecraft_radius - shell_thickness) ** 2),
        0.5
        * np.pi
        * (spacecraft_radius**2 - (spacecraft_radius - shell_thickness) ** 2),
        spacecraft_length * shell_thickness,
        0,
        0,
        0,
        0,
        0,
        0,
    ],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]

conductive_length = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 0, 0, 0, 0, 0, 0, 0],
    [
        1 / (2 * shell_thickness),
        1 / (2 * shell_thickness),
        1 / (2 * shell_thickness),
        0,
        0,
        0,
        0,
        0,
        0,
    ],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
]

conductive_relationship = np.multiply(
    np.multiply(conductive_coeff, conductive_area), conductive_length
)
conductive_relationship[conductive_relationship == np.nan] = 0

capacities_matrix = np.diag(np.multiply(masses, capacities))

node_relationship = conductive_relationship + capacities_matrix + radiative_relationship
