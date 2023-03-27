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
SB = 5.67e-8
spacecraft_radius = 1.5
spacecraft_length = 5
shell_thickness = 0.002
sail_thickness = 0.000004
shield_thickness = 0.005
sail_density = 1540
panel_density = 1.76
boom_length = 70
boom_radius = 0.229
boom_thickness = 0.0003
sail_area = 10000
panel_area = 18
heaters_power = 0
coolers_power = 0
electrical_heat = 800
log_file = "data/thermal.log"

### Structure Nodes
node_1 = {
    "name": "Spacecraft +Z",
    "area": np.pi * spacecraft_radius**2,
    "sun_vf": 0.125,
    "space_vf": 0.0025,
    "external": "multi-layer insulation",
    "internal": "carbon fibre",
    "temp_range": [0, 10000],
    "internal_heat": 0,
}
node_2 = {
    "name": "Spacecraft -Z",
    "area": 2*np.pi * spacecraft_radius**2,
    "sun_vf": 0,
    "space_vf": 0.5,
    "external": "az93 white paint",
    "internal": "carbon fibre",
    "temp_range": [0, 10000],
    "internal_heat": 800,
}
node_3 = {
    "name": "Spacecraft +X",
    "area": np.pi * spacecraft_radius * spacecraft_length,
    "sun_vf": 0.0189,
    "space_vf": 0.3599,
    "external": "aluminum",
    "internal": "carbon fibre",
    "temp_range": [0, 10000],
    "internal_heat": 0,
}
node_4 = {
    "name": "Spacecraft -X",
    "area": np.pi * spacecraft_radius * spacecraft_length,
    "sun_vf": 0.0189,
    "space_vf": 0.3599,
    "external": "aluminum",
    "internal": "carbon fibre",
    "temp_range": [0, 10000],
    "internal_heat": 0,
}
bus_nodes = [node_1, node_2, node_3, node_4]

### Sail Nodes
sail = {
    "name": "Solar Sail",
    "area": 2*sail_area,
    "sun_vf": 0.4875,
    "space_vf": 0.4875,
    "external": "standard",
    "internal": "cp-1",
    "temp_range": [0, 10000],
    "internal_heat": 0
}

### External Nodes
booms = {
    "name": "Booms",
    "area": 4 * np.pi * (boom_radius) * boom_length,
    "sun_vf": 0.4874,
    "space_vf": 0.4874,
    "external": "carbon fibre",
    "temp_range": [0, 10000],
    "internal_heat": 0,
}
solar_panels = {
    "name": "Solar Panels",
    "area": panel_area,
    "sun_vf": 0.7182,
    "space_vf": 0,
    "external": "osr",
    "temp_range": [0, 10000],
    "internal_heat": 0,
}

shield = {
    "name": "Heat Shield",
    "area": 2*node_1["area"],
    "sun_vf": 0.5,
    "space_vf": 0.15,
    "external": "ceramic cloth",
    "temp_range": [0, 10000],
    "internal_heat": 0,
    "layers": 7,
}



nodes = [
    node_1,
    node_2,
    node_3,
    node_4,
    sail,
    booms,
    solar_panels,
    shield
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
    densities.append(materials.bus_material(node["internal"])["density"])
    capacities.append(materials.bus_material(node["internal"])["specific_heat"])
    conductivities.append(materials.bus_material(node["external"])["conductivity"])
    emissivities.append(materials.bus_material(node["external"])["emissivity"])

densities.append(materials.sail_material(sail["external"])["density"])
capacities.append(materials.sail_material(sail["external"])["specific_heat"])
conductivities.append(materials.sail_material(sail["external"])["conductivity"])
emissivities.append(materials.sail_material(sail["external"])["emissivity"])


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
    sail["area"] * sail_thickness * densities[4] * 0.5,
    densities[5] * 4
    * np.pi
    * ((boom_radius) ** 2 - (boom_radius - boom_thickness) ** 2)
    * boom_length,
    panel_density * solar_panels["area"],
    shield["layers"] * densities[7] * shield["area"] * shield_thickness,
]


view_factors = [
    [0, 0.33*0.5, 0.33*0.5, 0.33*0.5, 0, 0, 0, 0.7],
    [0, 0, 0.33*0.5, 0.33*0.5, 0, 0, 0, 0],
    [0, 0, 0, 0.33*0.5, 0, 0, 0.1, 0.02],
    [0, 0, 0, 0, 0, 0, 0.1, 0.02],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

view_factors_deployed = [
    [0, 0, 0, 0, 0, 0, 0, 0.7],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0.0424413, 0.1, 0.1, 0.02],
    [0, 0, 0, 0, 0, 0.0424413, 0.1, 0.02],
    [0, 0, 0, 0, 0, 0.002, 0, 0.0001],
    [0, 0, 0, 0, 0, 0, 0.00178714, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# view_factors_deployed = [
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.02],
#     [0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.02],
#     [0, 0, 0, 0, 0, 0, 0.002014, 0.00018, 0],
#     [0, 0, 0, 0, 0, 0, 0.002014, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0.00178714, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0],
# ]



radiative_area = [
    [
        0,
        node_1["area"]*node_2["area"],
        node_1["area"]*node_3["area"],
        node_1["area"]*node_4["area"],
        node_1["area"]*sail["area"],
        node_1["area"]*booms["area"],
        node_1["area"]*solar_panels["area"],
        node_1["area"]*shield["area"],
    ],
    [
        0,
        0,
        node_2["area"]*node_3["area"],
        node_2["area"]*node_4["area"],
        node_2["area"]*sail["area"],
        node_2["area"]*booms["area"],
        node_2["area"]*solar_panels["area"],
        node_2["area"]*shield["area"],
    ],
    [
        0,
        0,
        0,
        node_3["area"]*node_4["area"],
        node_3["area"]*sail["area"],
        node_3["area"]*booms["area"],
        node_3["area"]*solar_panels["area"],
        node_3["area"]*shield["area"],
    ],
    [
        0,
        0,
        0,
        0,
        node_4["area"]*sail["area"],
        node_4["area"]*booms["area"],
        node_4["area"]*solar_panels["area"],
        node_4["area"]*shield["area"]
    ],
    [
        0,
        0,
        0,
        0,
        0,
        sail["area"]*booms["area"],
        sail["area"]*solar_panels["area"],
        sail["area"]*shield["area"],
    ],
    [0, 0, 0, 0, 0, 0, booms["area"]*solar_panels["area"], booms["area"]*shield["area"]],
    [0, 0, 0, 0, 0, 0, 0, solar_panels["area"]*shield["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0]
]

# effective_emissivities = [
#     [
#         0,
#         1 / ((1 / emissivities[0]) + (1 / emissivities[1])-1),
#         1 / ((1 / emissivities[0]) + (1 / emissivities[2])-1),
#         1 / ((1 / emissivities[0]) + (1 / emissivities[3])-1),
#         1 / ((1 / emissivities[0]) + (1 / emissivities[4])-1),
#         1 / ((1 / emissivities[0]) + (1 / emissivities[5])-1),
#         1 / ((1 / emissivities[0]) + (1 / emissivities[6])-1),
#         1 / ((1 / emissivities[0]) + (1 / emissivities[7])-1),
#         1 / ((shield_front["layers"] + 1)*(1 / emissivities[0]) + (1 / emissivities[8])-1),
#         1 / ((shield_front["layers"] + 1)*(1 / emissivities[0]) + (1 / emissivities[9])-1)
#     ],
#     [
#         0,
#         0,
#         1 / ((1 / emissivities[1]) + (1 / emissivities[2])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[3])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[4])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[5])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[6])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[7])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[8])-1),
#         1 / ((1 / emissivities[1]) + (1 / emissivities[9])-1)
#     ],
#     [
#         0,
#         0,
#         0,
#         1 / ((1 / emissivities[2]) + (1 / emissivities[3])-1),
#         1 / ((1 / emissivities[2]) + (1 / emissivities[4])-1),
#         1 / ((1 / emissivities[2]) + (1 / emissivities[5])-1),
#         1 / ((1 / emissivities[2]) + (1 / emissivities[6])-1),
#         1 / ((1 / emissivities[2]) + (1 / emissivities[7])-1),
#         1 / ((1 / emissivities[2]) + (1 / emissivities[8])-1),
#         1 / ((1 / emissivities[2]) + (1 / emissivities[9])-1)

#     ],
#     [
#         0,
#         0,
#         0,
#         0,
#         1 / ((1 / emissivities[3]) + (1 / emissivities[4])-1),
#         1 / ((1 / emissivities[3]) + (1 / emissivities[5])-1),
#         1 / ((1 / emissivities[3]) + (1 / emissivities[6])-1),
#         1 / ((1 / emissivities[3]) + (1 / emissivities[7])-1),
#         1 / ((1 / emissivities[3]) + (1 / emissivities[8])-1),
#         1 / ((1 / emissivities[3]) + (1 / emissivities[9])-1)
#     ],
#     [
#         0,
#         0,
#         0,
#         0,
#         0,
#         1 / ((1 / emissivities[4]) + (1 / emissivities[5])-1),
#         1 / ((1 / emissivities[4]) + (1 / emissivities[6])-1),
#         1 / ((1 / emissivities[4]) + (1 / emissivities[7])-1),
#         1 / ((1 / emissivities[4]) + (1 / emissivities[8])-1),
#         1 / ((1 / emissivities[4]) + (1 / emissivities[9])-1)
#     ],
#     [
#         0,
#         0,
#         0,
#         0,
#         0,
#         0,
#         1 / ((1 / emissivities[5]) + (1 / emissivities[6])-1),
#         1 / ((1 / emissivities[5]) + (1 / emissivities[7])-1),
#         1 / ((1 / emissivities[5]) + (1 / emissivities[8])-1),
#         1 / ((1 / emissivities[5]) + (1 / emissivities[9])-1)
#     ],
#     [
#         0,
#         0,
#         0,
#         0,
#         0,
#         0,
#         0,
#         1 / ((1 / emissivities[6]) + (1 / emissivities[7])-1),
#         1 / ((1 / emissivities[6]) + (1 / emissivities[8])-1),
#         1 / ((1 / emissivities[6]) + (1 / emissivities[9])-1)
#     ],
#     [0, 0, 0, 0, 0, 0, 0, 0, 1 / ((1 / emissivities[7]) + (1 / emissivities[8])-1), 1 / ((1 / emissivities[7]) + (1 / emissivities[9])-1)],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1 / ((1 / emissivities[8]) + (1 / emissivities[9])-1)],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# ]

effective_emissivities = [[0,
                           emissivities[0]*emissivities[1], 
                           emissivities[0]*emissivities[2],
                           emissivities[0]*emissivities[3],
                           emissivities[0]*emissivities[4],
                           emissivities[0]*emissivities[5],
                           emissivities[0]*emissivities[6],
                           1/(shield["layers"]+1)*emissivities[0]*emissivities[7]], 
                           [0, 
                            0,
                           emissivities[1]*emissivities[2],
                           emissivities[1]*emissivities[3],
                           emissivities[1]*emissivities[4],
                           emissivities[1]*emissivities[5],
                           emissivities[1]*emissivities[6],
                           emissivities[1]*emissivities[7]],
                           [0, 
                            0,
                            0,
                           emissivities[2]*emissivities[3],
                           emissivities[2]*emissivities[4],
                           emissivities[2]*emissivities[5],
                           emissivities[2]*emissivities[6],
                           emissivities[2]*emissivities[7]],
                           [0, 
                            0,
                            0,
                            0,
                           emissivities[3]*emissivities[4],
                           emissivities[3]*emissivities[5],
                           emissivities[3]*emissivities[6],
                           emissivities[3]*emissivities[7]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                           emissivities[4]*emissivities[5],
                           emissivities[4]*emissivities[6],
                           emissivities[4]*emissivities[7]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                            0,
                           emissivities[5]*emissivities[6],
                           emissivities[5]*emissivities[7]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                           emissivities[6]*emissivities[7]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0]]

radiative_relationship = SB * np.multiply(
    np.multiply(view_factors, radiative_area), effective_emissivities
)
conductive_coeff = [[0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [400, 400, 0, 0, 0, 0, 0, 0],
                          [400, 400, 400, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0],
                          [400, 0, 0, 0, 0, 0, 0, 0]]

# conductive_coeff = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                           [13.7, 13.7, 0, 0, 0, 0, 0, 0, 0, 0],
#                           [13.7, 13.7, 13.7, 0, 0, 0, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0.155, 0, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0.155, 0.155, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                           [0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0]]


conductive_area = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
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
        0
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
        0
    ],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [4*np.pi*0.01**2, 0, 0, 0, 0, 0, 0, 0]
]

conductive_length = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 0, 0, 0, 0, 0, 0],
    [
        1 / (2 * shell_thickness),
        1 / (2 * shell_thickness),
        1 / (2 * shell_thickness),
        0,
        0,
        0,
        0,
        0
    ],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0.5, 0, 0, 0, 0, 0, 0, 0]
]

conductive_relationship = np.multiply(
    np.multiply(conductive_coeff, conductive_area), conductive_length
)
conductive_relationship[conductive_relationship == np.nan] = 0

capacities_matrix = np.diag(np.multiply(masses, capacities))

node_relationship = conductive_relationship + capacities_matrix + radiative_relationship

radiative_relationship_deployed = SB * np.multiply(
    np.multiply(view_factors_deployed, radiative_area), effective_emissivities
)
node_relationship_deployed = conductive_relationship + capacities_matrix + radiative_relationship_deployed
