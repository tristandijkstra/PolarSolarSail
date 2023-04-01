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
spacecraft_width = 3
spacecraft_length = 3
shell_thickness = 0.002
sail_thickness = 0.000004
shield_thickness = 0.005
sail_density = 1540
panel_density = 1.76
boom_length = 70
boom_radius = 0.229
boom_thickness = 0.0003
antenna_radius = 0.55
sail_area = 10000
antenna_mass = 28.0
panel_area = 18
heaters_power = 0
coolers_power = 0
electrical_heat = 800
battery_mass = 7.5
hydrazine_mass = 36
metis_mass = 19
cdm_mass = 19
log_file = "data/thermal.log"
internal = True

### Structure Nodes
node_1 = {
    "name": "Spacecraft +Z",
    "area": spacecraft_width**2,
    "sun_vf": 0.0,
    "space_vf": 0.200,
    "external": "multi-layer insulation",
    "internal": "carbon fibre",
    "temp_range": [-200, 600],
    "internal_heat": 400,
}
node_2 = {
    "name": "Spacecraft -Z",
    "area": 1.75*spacecraft_width**2,
    "sun_vf": 0,
    "space_vf": 1,
    "external": "az93 white paint",
    "internal": "carbon fibre",
    "temp_range": [-200, 600],
    "internal_heat": 1200,
}
node_3 = {
    "name": "Spacecraft +X",
    "area": spacecraft_width*spacecraft_length,
    "sun_vf": 0.1,
    "space_vf": 0.7098,
    "external": "multi-layer insulation",
    "internal": "carbon fibre",
    "temp_range": [-200, 600],
    "internal_heat": 400,
}
node_4 = {
    "name": "Spacecraft -X",
    "area": spacecraft_width*spacecraft_length,
    "sun_vf": 0.1,
    "space_vf": 0.7098,
    "external": "multi-layer insulation",
    "internal": "carbon fibre",
    "temp_range": [-200, 600],
    "internal_heat": 400,
}
node_5 = {
    "name": "Spacecraft +Y",
    "area": spacecraft_width*spacecraft_length,
    "sun_vf": 0.1,
    "space_vf": 0.7098,
    "external": "az93 white paint",
    "internal": "carbon fibre",
    "temp_range": [-200, 600],
    "internal_heat": 400,
}
node_6 = {
    "name": "Spacecraft -Y",
    "area": spacecraft_width*spacecraft_length,
    "sun_vf": 0.1,
    "space_vf": 0.7098,
    "external": "az93 white paint",
    "internal": "carbon fibre",
    "temp_range": [-200, 600],
    "internal_heat": 400,
}


bus_nodes = [node_1, node_2, node_3, node_4, node_5, node_6]

### Sail Nodes
sail = {
    "name": "Solar Sail",
    "area": 2*sail_area,
    "sun_vf": 0.4875,
    "space_vf": 0.4875,
    "external": "standard",
    "internal": "cp-1",
    "temp_range": [-100, 260],
    "internal_heat": 0
}

### External Nodes
booms = {
    "name": "Booms",
    "area": 2 * 4 * np.pi * (boom_radius) * boom_length,
    "sun_vf": 0.4874,
    "space_vf": 0.4874,
    "external": "aluminum",
    "internal": "carbon fibre",
    "temp_range": [-100, 260],
    "internal_heat": 0,
}
solar_panels = {
    "name": "Solar Panels",
    "area": 2*panel_area,
    "sun_vf": 0.5,
    "space_vf": 0.2,
    "external": "osr",
    "temp_range": [-65, 200],
    "internal_heat": 0,
}

antenna = {
    "name": "Antenna",
    "area": np.pi*antenna_radius**2,
    "sun_vf": 0.5,
    "space_vf": 0.4,
    "external": "solar black",
    "internal": "titanium",
    "temp_range": [-200, 580],
    "internal_heat": 0,
}

shield_front = {
    "name": "Heat Shield",
    "area": 2*node_1["area"],
    "sun_vf": 0.50,
    "space_vf": 0.10,
    "external": "ceramic cloth",
    "internal": "multi-layer insulation",
    "temp_range": [-200, 630],
    "internal_heat": 0,
}

shield_inner = {
    "name": "Heat Shield Inner",
    "area": 2*node_1["area"],
    "sun_vf": 0,
    "space_vf": 0.33,
    "external": "multi-layer insulation",
    "internal": "multi-layer insulation",
    "temp_range": [-200, 630],
    "internal_heat": 0,
    "layers": 10,
}

# Internal Nodes

batteries = {
    "name": "Batteries",
    "area": 0.240,
    "sun_vf": 0,
    "space_vf": 0,
    "external": "battery",
    "internal": "battery",
    "temp_range": [0, 45],
    "internal_heat": 0.02,
}

hydrazine = {
    "name": "Hydrazine",
    "area": 0.348,
    "sun_vf": 0,
    "space_vf": 0,
    "external": "hydrazine",
    "internal": "hydrazine",
    "temp_range": [15, 40],
    "internal_heat": 1.87,
}

metis = {
    "name": "METIS",
    "area": 2.38,
    "sun_vf": 0.056,
    "space_vf": 0,
    "external": "metis",
    "internal": "metis",
    "temp_range": [-45, -5],
    "internal_heat": 0.02,
}

cdm = {
    "name": "CDM",
    "area": 0.395,
    "sun_vf": 0.0286,
    "space_vf": 0,
    "external": "cdm",
    "internal": "cdm",
    "temp_range": [-10, 30],
    "internal_heat": 0.01,
}


nodes = [
    node_1,
    node_2,
    node_3,
    node_4,
    node_5,
    node_6,
    sail,
    booms,
    solar_panels,
    antenna,
    shield_front,
    shield_inner
]

if internal:
    nodes = [
        node_1,
        node_2,
        node_3,
        node_4,
        node_5,
        node_6,
        sail,
        booms,
        solar_panels,
        antenna,
        shield_front,
        shield_inner,
        batteries,
        hydrazine,
        metis,
        cdm
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
    conductivities.append(materials.bus_material(node["internal"])["conductivity"])
    emissivities.append(materials.bus_material(node["external"])["emissivity"])

densities.append(materials.sail_material(sail["external"])["density"])
capacities.append(materials.sail_material(sail["external"])["specific_heat"])
conductivities.append(materials.sail_material(sail["external"])["conductivity"])
emissivities.append(materials.sail_material(sail["external"])["emissivity"])


densities.append(materials.boom_material(booms["internal"])["density"])
capacities.append(materials.boom_material(booms["internal"])["specific_heat"])
conductivities.append(materials.boom_material(booms["internal"])["conductivity"])
emissivities.append(materials.boom_material(booms["external"])["emissivity"])

densities.append(materials.panel_material(solar_panels["external"])["density"])
capacities.append(materials.panel_material(solar_panels["external"])["specific_heat"])
conductivities.append(
    materials.panel_material(solar_panels["external"])["conductivity"]
)
emissivities.append(materials.panel_material(solar_panels["external"])["emissivity"])

densities.append(materials.antenna_material(antenna["internal"])["density"])
capacities.append(materials.antenna_material(antenna["internal"])["specific_heat"])
conductivities.append(materials.antenna_material(antenna["internal"])["conductivity"])
emissivities.append(materials.antenna_material(antenna["external"])["emissivity"])

densities.append(materials.shield_material(shield_front["internal"])["density"])
capacities.append(materials.shield_material(shield_front["internal"])["specific_heat"])
conductivities.append(materials.shield_material(shield_front["internal"])["conductivity"])
emissivities.append(materials.shield_material(shield_front["external"])["emissivity"])

densities.append(materials.shield_material(shield_inner["internal"])["density"])
capacities.append(materials.shield_material(shield_inner["internal"])["specific_heat"])
conductivities.append(materials.shield_material(shield_inner["internal"])["conductivity"])
emissivities.append(materials.shield_material(shield_inner["external"])["emissivity"])

if internal:
    densities.append(materials.internal_material(batteries["internal"])["density"])
    capacities.append(materials.internal_material(batteries["internal"])["specific_heat"])
    conductivities.append(materials.internal_material(batteries["internal"])["conductivity"])
    emissivities.append(materials.internal_material(batteries["external"])["emissivity"])

    densities.append(materials.internal_material(hydrazine["internal"])["density"])
    capacities.append(materials.internal_material(hydrazine["internal"])["specific_heat"])
    conductivities.append(materials.internal_material(hydrazine["internal"])["conductivity"])
    emissivities.append(materials.internal_material(hydrazine["external"])["emissivity"])

    densities.append(materials.internal_material(metis["internal"])["density"])
    capacities.append(materials.internal_material(metis["internal"])["specific_heat"])
    conductivities.append(materials.internal_material(metis["internal"])["conductivity"])
    emissivities.append(materials.internal_material(metis["external"])["emissivity"])

    densities.append(materials.internal_material(cdm["internal"])["density"])
    capacities.append(materials.internal_material(cdm["internal"])["specific_heat"])
    conductivities.append(materials.internal_material(cdm["internal"])["conductivity"])
    emissivities.append(materials.internal_material(cdm["external"])["emissivity"])


masses = [
    node_1["area"] * shell_thickness * densities[0],
    node_2["area"] * shell_thickness * densities[1],
    node_3["area"] * shell_thickness * densities[2],
    node_4["area"] * shell_thickness * densities[3],
    node_5["area"] * shell_thickness * densities[4],
    node_6["area"] * shell_thickness * densities[5],
    sail["area"] * sail_thickness * densities[6] * 0.5,
    densities[7] * 4
    * np.pi
    * ((boom_radius) ** 2 - (boom_radius - boom_thickness) ** 2)
    * boom_length,
    panel_density * solar_panels["area"],
    antenna_mass,
    densities[9] * shield_front["area"] * shield_thickness,
    shield_inner["area"]*shield_inner["layers"]*shield_thickness,
    battery_mass,
    hydrazine_mass,
    metis_mass,
    cdm_mass
]


view_factors = [
    [0, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0, 0, 0, 0, 0.7*0.5],
    [0, 0, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0, 0, 0, 0, 0],
    [0, 0, 0, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0, 0, 0.1*0.25, 0, 0.01*0.25],
    [0, 0, 0, 0, 0.167*0.5, 0.167*0.5, 0, 0, 0.1*0.25, 0, 0.01*0.25],
    [0, 0, 0, 0, 0, 0.167*0.5, 0, 0, 0.1*0.25, 0, 0.01*0.25],
    [0, 0, 0, 0, 0, 0, 0, 0, 0.1*0.25, 0, 0.01*0.25],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.35]
]

view_factors_deployed = [
    [0, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0, 0, 0, 0.7*0.5],
    [0, 0, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0, 0, 0, 0],
    [0, 0, 0, 0.167*0.5, 0.167*0.5, 0.167*0.5, 0.0424413*0.25, 0.00089*0.25, 0.1*0.25, 0.01*0.25],
    [0, 0, 0, 0, 0.167*0.5, 0.167*0.5, 0.0424413*0.25, 0.00089*0.25, 0.1*0.25, 0.01*0.25],
    [0, 0, 0, 0, 0, 0.167*0.5, 0.0424413*0.25, 0.00089*0.25, 0.1*0.25, 0.01*0.25],
    [0, 0, 0, 0, 0, 0, 0.0424413*0.25, 0.00089*0.25, 0.1*0.25, 0.01*0.25],
    [0, 0, 0, 0, 0, 0, 0, 0.002, 0, 0.0001],
    [0, 0, 0, 0, 0, 0, 0, 0, 0.00178714, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

view_factors_deployed = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.35],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5],
    [0, 0, 0, 0, 0, 0, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5],
    [0, 0, 0, 0, 0, 0, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5],
    [0, 0, 0, 0, 0, 0, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5],
    [0, 0, 0, 0, 0, 0, 0, 0.002, 0.005, 0, 0, 0.0001],
    [0, 0, 0, 0, 0, 0, 0, 0, 0.0017814, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.35],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

if internal:
    view_factors_deployed = [
        [0, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0, 0, 0, 0, 0, 0.35, 0.0167*0.5, 0, 0.0167*0.5, 0.0167*0.5],
        [0, 0, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0, 0, 0, 0, 0, 0, 0, 0.0167*0.5, 0, 0],
        [0, 0, 0, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5, 0, 0.167*0.5, 0.0167*0.5, 0.0167*0.5],
        [0, 0, 0, 0, 0.0167*0.5, 0.0167*0.5, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5, 0.0167*0.5, 0.167*0.5, 0.0167*0.5, 0.0167*0.5],
        [0, 0, 0, 0, 0, 0.0167*0.5, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5, 0.0167*0.5, 0.167*0.5, 0.0167*0.5, 0.0167*0.5],
        [0, 0, 0, 0, 0, 0, 0.0424413*0.5, 0.00085*0.5, 0.01, 0.01, 0, 0.01*0.5, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5, 0.0167*0.5],
        [0, 0, 0, 0, 0, 0, 0, 0.002, 0.005, 0, 0, 0.0001, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0.0017814, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.35, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0002, 0.0002, 0.0002],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0002, 0.0002],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0002],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
        0.5*node_1["area"]*0.5*node_2["area"],
        0.5*node_1["area"]*0.5*node_3["area"],
        0.5*node_1["area"]*0.5*node_4["area"],
        0.5*node_1["area"]*0.5*node_5["area"],
        0.5*node_1["area"]*0.5*node_6["area"],
        0.5*node_1["area"]*0.125*sail["area"],
        0.5*node_1["area"]*0.125*booms["area"],
        0.5*node_1["area"]*solar_panels["area"],
        0.5*node_1["area"]*antenna["area"],
        0.5*node_1["area"]*0.5*shield_front["area"],
        0.5*node_1["area"]*0.5*shield_inner["area"],
    ],
    [
        0,
        0,
        0.5*node_2["area"]*0.5*node_3["area"],
        0.5*node_2["area"]*0.5*node_4["area"],
        0.5*node_2["area"]*0.5*node_5["area"],
        0.5*node_2["area"]*0.5*node_6["area"],
        0.5*node_2["area"]*0.125*sail["area"],
        0.5*node_2["area"]*0.125*booms["area"],
        0.5*node_2["area"]*solar_panels["area"],
        0.5*node_2["area"]*antenna["area"],
        0.5*node_2["area"]*0.5*shield_front["area"],
        0.5*node_2["area"]*0.5*shield_inner["area"]
    ],
    [
        0,
        0,
        0,
        0.5*node_3["area"]*0.5*node_4["area"],
        0.5*node_3["area"]*0.5*node_5["area"],
        0.5*node_3["area"]*0.5*node_6["area"],
        0.5*node_3["area"]*0.125*sail["area"],
        0.5*node_3["area"]*0.125*booms["area"],
        0.5*node_3["area"]*solar_panels["area"],
        0.5*node_3["area"]*antenna["area"],
        0.5*node_3["area"]*0.5*shield_front["area"],
        0.5*node_3["area"]*0.5*shield_inner["area"],
    ],
    [
        0,
        0,
        0,
        0,
        0.5*node_4["area"]*0.5*node_5["area"],
        0.5*node_4["area"]*0.5*node_6["area"],
        0.5*node_4["area"]*0.125*sail["area"],
        0.5*node_4["area"]*0.125*booms["area"],
        0.5*node_4["area"]*solar_panels["area"],
        0.5*node_4["area"]*antenna["area"],
        0.5*node_4["area"]*0.5*shield_front["area"],
        0.5*node_4["area"]*0.5*shield_inner["area"]
    ],
    [
        0,
        0,
        0,
        0,
        0,
        0.5*node_5["area"]*0.5*node_6["area"],
        0.5*node_5["area"]*0.125*sail["area"],
        0.5*node_5["area"]*0.125*booms["area"],
        0.5*node_5["area"]*solar_panels["area"],
        0.5*node_5["area"]*antenna["area"],
        0.5*node_5["area"]*0.5*shield_front["area"],
        0.5*node_5["area"]*0.5*shield_inner["area"]
    ],
    [0, 0, 0, 0, 0, 0, 0.5*node_6["area"]*0.125*sail["area"], 0.5*node_6["area"]*0.125*booms["area"], 0.5*node_6["area"]*solar_panels["area"], 0.5*node_6["area"]*antenna["area"], 0.5*node_6["area"]*0.5*shield_front["area"], 0.5*node_6["area"]*0.5*shield_inner["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0.5*sail["area"]*0.125*booms["area"], 0.5*sail["area"]*0.5*solar_panels["area"], 0.125*sail["area"]*antenna["area"], 0.125*sail["area"]*0.5*shield_front["area"], 0.125*sail["area"]*0.5*shield_inner["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0.125*booms["area"]*solar_panels["area"], 0.125*booms["area"]*antenna["area"], 0.125*booms["area"]*0.5*shield_front["area"],0.125*booms["area"]*0.5*shield_inner["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, solar_panels["area"]*antenna["area"], solar_panels["area"]*0.5*shield_front["area"], solar_panels["area"]*0.5*shield_inner["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, antenna["area"]*0.5*shield_front["area"], antenna["area"]*0.5*shield_front["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5*shield_front["area"]*0.5*shield_inner["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

if internal:
    radiative_area = [
    [
        0,
        0.5*node_1["area"]*0.5*node_2["area"],
        0.5*node_1["area"]*0.5*node_3["area"],
        0.5*node_1["area"]*0.5*node_4["area"],
        0.5*node_1["area"]*0.5*node_5["area"],
        0.5*node_1["area"]*0.5*node_6["area"],
        0.5*node_1["area"]*0.125*sail["area"],
        0.5*node_1["area"]*0.125*booms["area"],
        0.5*node_1["area"]*solar_panels["area"],
        0.5*node_1["area"]*antenna["area"],
        0.5*node_1["area"]*0.5*shield_front["area"],
        0.5*node_1["area"]*0.5*shield_inner["area"],
        0.5*node_1["area"]*0.167*batteries["area"],
        0.5*node_1["area"]*0.5*hydrazine["area"],
        0.5*node_1["area"]*0.083*metis["area"],
        0.5*node_1["area"]*0.083*cdm["area"]
    ],
    [
        0,
        0,
        0.5*node_2["area"]*0.5*node_3["area"],
        0.5*node_2["area"]*0.5*node_4["area"],
        0.5*node_2["area"]*0.5*node_5["area"],
        0.5*node_2["area"]*0.5*node_6["area"],
        0.5*node_2["area"]*0.125*sail["area"],
        0.5*node_2["area"]*0.125*booms["area"],
        0.5*node_2["area"]*solar_panels["area"],
        0.5*node_2["area"]*antenna["area"],
        0.5*node_2["area"]*0.5*shield_front["area"],
        0.5*node_2["area"]*0.5*shield_inner["area"],
        0.5*node_2["area"]*0.167*batteries["area"],
        0.5*node_2["area"]*0.5*hydrazine["area"],
        0.5*node_2["area"]*0.083*metis["area"],
        0.5*node_2["area"]*0.083*cdm["area"]
    ],
    [
        0,
        0,
        0,
        0.5*node_3["area"]*0.5*node_4["area"],
        0.5*node_3["area"]*0.5*node_5["area"],
        0.5*node_3["area"]*0.5*node_6["area"],
        0.5*node_3["area"]*0.125*sail["area"],
        0.5*node_3["area"]*0.125*booms["area"],
        0.5*node_3["area"]*solar_panels["area"],
        0.5*node_3["area"]*antenna["area"],
        0.5*node_3["area"]*0.5*shield_front["area"],
        0.5*node_3["area"]*0.5*shield_inner["area"],
        0.5*node_3["area"]*0.167*batteries["area"],
        0.5*node_3["area"]*0.5*hydrazine["area"],
        0.5*node_3["area"]*0.167*metis["area"],
        0.5*node_3["area"]*0.167*cdm["area"]
    ],
    [
        0,
        0,
        0,
        0,
        0.5*node_4["area"]*0.5*node_5["area"],
        0.5*node_4["area"]*0.5*node_6["area"],
        0.5*node_4["area"]*0.125*sail["area"],
        0.5*node_4["area"]*0.125*booms["area"],
        0.5*node_4["area"]*solar_panels["area"],
        0.5*node_4["area"]*antenna["area"],
        0.5*node_4["area"]*0.5*shield_front["area"],
        0.5*node_4["area"]*0.5*shield_inner["area"],
        0.5*node_4["area"]*0.167*batteries["area"],
        0.5*node_4["area"]*0.5*hydrazine["area"],
        0.5*node_4["area"]*0.167*metis["area"],
        0.5*node_4["area"]*0.167*cdm["area"]
    ],
    [
        0,
        0,
        0,
        0,
        0,
        0.5*node_5["area"]*0.5*node_6["area"],
        0.5*node_5["area"]*0.125*sail["area"],
        0.5*node_5["area"]*0.125*booms["area"],
        0.5*node_5["area"]*solar_panels["area"],
        0.5*node_5["area"]*antenna["area"],
        0.5*node_5["area"]*0.5*shield_front["area"],
        0.5*node_5["area"]*0.5*shield_inner["area"],
        0.5*node_5["area"]*0.167*batteries["area"],
        0.5*node_5["area"]*0.5*hydrazine["area"],
        0.5*node_5["area"]*0.167*metis["area"],
        0.5*node_5["area"]*0.167*cdm["area"]
    ],
    [0, 0, 0, 0, 0, 0, 0.5*node_6["area"]*0.125*sail["area"], 0.5*node_6["area"]*0.125*booms["area"], 0.5*node_6["area"]*solar_panels["area"], 0.5*node_6["area"]*antenna["area"], 0.5*node_6["area"]*0.5*shield_front["area"], 0.5*node_6["area"]*0.5*shield_inner["area"], 0.5*node_6["area"]*0.167*batteries["area"], 0.5*node_6["area"]*0.5*hydrazine["area"], 0.5*node_6["area"]*0.167*metis["area"], 0.5*node_6["area"]*0.167*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0.5*sail["area"]*0.125*booms["area"], 0.5*sail["area"]*0.5*solar_panels["area"], 0.125*sail["area"]*antenna["area"], 0.125*sail["area"]*0.5*shield_front["area"], 0.125*sail["area"]*0.5*shield_inner["area"], sail["area"]*0.167*batteries["area"], sail["area"]*0.5*hydrazine["area"], sail["area"]*0.167*metis["area"], sail["area"]*0.167*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0.125*booms["area"]*solar_panels["area"], 0.125*booms["area"]*antenna["area"], 0.125*booms["area"]*0.5*shield_front["area"],0.125*booms["area"]*0.5*shield_inner["area"], booms["area"]*0.167*batteries["area"], booms["area"]*0.5*hydrazine["area"], booms["area"]*0.167*metis["area"], booms["area"]*0.167*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, solar_panels["area"]*antenna["area"], solar_panels["area"]*0.5*shield_front["area"], solar_panels["area"]*0.5*shield_inner["area"], solar_panels["area"]*batteries["area"], solar_panels["area"]*hydrazine["area"], solar_panels["area"]*metis["area"], solar_panels["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, antenna["area"]*0.5*shield_front["area"], antenna["area"]*0.5*shield_front["area"], antenna["area"]*batteries["area"], antenna["area"]*hydrazine["area"], antenna["area"]*metis["area"], antenna["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5*shield_front["area"]*0.5*shield_inner["area"], shield_front["area"]*batteries["area"], shield_front["area"]*hydrazine["area"], shield_front["area"]*metis["area"], shield_front["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, shield_inner["area"]*batteries["area"], shield_inner["area"]*hydrazine["area"], shield_inner["area"]*metis["area"], shield_inner["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, batteries["area"]*hydrazine["area"], batteries["area"]*metis["area"], batteries["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, hydrazine["area"]*metis["area"], hydrazine["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, metis["area"]*cdm["area"]],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
                           materials.bus_material(node_1["internal"])["absorptivity"]*emissivities[0]*emissivities[6],
                           materials.bus_material(node_1["internal"])["absorptivity"]*emissivities[0]*emissivities[7],
                           emissivities[0]*emissivities[8],
                           emissivities[0]*emissivities[9],
                           emissivities[0]*emissivities[10],
                           emissivities[0]*emissivities[11],
                           emissivities[0]*emissivities[12],
                           emissivities[0]*emissivities[13],
                           emissivities[0]*emissivities[14],
                           emissivities[0]*emissivities[15]], 
                           [0, 
                            0,
                           emissivities[1]*emissivities[2],
                           emissivities[1]*emissivities[3],
                           emissivities[1]*emissivities[4],
                           emissivities[1]*emissivities[5],
                           materials.bus_material(node_2["internal"])["absorptivity"]*emissivities[1]*emissivities[6],
                           materials.bus_material(node_2["internal"])["absorptivity"]*emissivities[1]*emissivities[7],
                           emissivities[1]*emissivities[8],
                           emissivities[1]*emissivities[9],
                           emissivities[1]*emissivities[10],
                           emissivities[1]*emissivities[11],
                           emissivities[1]*emissivities[12],
                           emissivities[1]*emissivities[13],
                           emissivities[1]*emissivities[14],
                           emissivities[1]*emissivities[15]],
                           [0, 
                            0,
                            0,
                           emissivities[2]*emissivities[3],
                           emissivities[2]*emissivities[4],
                           emissivities[2]*emissivities[5],
                           materials.bus_material(node_3["internal"])["absorptivity"]*emissivities[2]*emissivities[6],
                           materials.bus_material(node_3["internal"])["absorptivity"]*emissivities[2]*emissivities[7],
                           emissivities[2]*emissivities[8],
                           emissivities[2]*emissivities[9],
                           emissivities[2]*emissivities[10],
                           emissivities[2]*emissivities[11],
                           emissivities[2]*emissivities[12],
                           emissivities[2]*emissivities[13],
                           emissivities[2]*emissivities[14],
                           emissivities[2]*emissivities[15]],
                           [0, 
                            0,
                            0,
                            0,
                           emissivities[3]*emissivities[4],
                           emissivities[3]*emissivities[5],
                           materials.bus_material(node_4["internal"])["absorptivity"]*emissivities[2]*emissivities[6],
                           materials.bus_material(node_4["internal"])["absorptivity"]*emissivities[2]*emissivities[7],
                           emissivities[3]*emissivities[8],
                           emissivities[3]*emissivities[9],
                           emissivities[3]*emissivities[10],
                           emissivities[3]*emissivities[11],
                           emissivities[3]*emissivities[12],
                           emissivities[3]*emissivities[13],
                           emissivities[3]*emissivities[14],
                           emissivities[3]*emissivities[15]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                           emissivities[4]*emissivities[5],
                           materials.bus_material(node_5["internal"])["absorptivity"]*emissivities[2]*emissivities[6],
                           materials.bus_material(node_5["internal"])["absorptivity"]*emissivities[2]*emissivities[7],
                           emissivities[4]*emissivities[8],
                           emissivities[4]*emissivities[9],
                           emissivities[4]*emissivities[10],
                           emissivities[4]*emissivities[11],
                           emissivities[4]*emissivities[12],
                           emissivities[4]*emissivities[13],
                           emissivities[4]*emissivities[14],
                           emissivities[4]*emissivities[15]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                            0,
                           materials.bus_material(node_6["internal"])["absorptivity"]*emissivities[2]*emissivities[6],
                           materials.bus_material(node_6["internal"])["absorptivity"]*emissivities[2]*emissivities[7],
                           emissivities[5]*emissivities[8],
                           emissivities[5]*emissivities[9],
                           emissivities[5]*emissivities[10],
                           emissivities[5]*emissivities[11],
                           emissivities[5]*emissivities[12],
                           emissivities[5]*emissivities[13],
                           emissivities[5]*emissivities[14],
                           emissivities[5]*emissivities[15]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            emissivities[6]*emissivities[7],
                            emissivities[6]*emissivities[8],
                           emissivities[6]*emissivities[9],
                           emissivities[6]*emissivities[10],
                           emissivities[6]*emissivities[11],
                           emissivities[6]*emissivities[12],
                           emissivities[6]*emissivities[13],
                           emissivities[6]*emissivities[14],
                           emissivities[6]*emissivities[15]],
                           [0, 
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            emissivities[7]*emissivities[8],
                            emissivities[7]*emissivities[9],
                            emissivities[7]*emissivities[10],
                            emissivities[7]*emissivities[11],
                           emissivities[7]*emissivities[12],
                           emissivities[7]*emissivities[13],
                           emissivities[7]*emissivities[14],
                           emissivities[7]*emissivities[15]],
                            [0,
                             0, 
                             0, 
                             0, 
                             0, 
                             0, 
                             0, 
                             0,
                             0, 
                             emissivities[8]*emissivities[9],
                             emissivities[8]*emissivities[10],
                             emissivities[8]*emissivities[11],
                           emissivities[8]*emissivities[12],
                           emissivities[8]*emissivities[13],
                           emissivities[8]*emissivities[14],
                           emissivities[8]*emissivities[15]],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, emissivities[9]*emissivities[10], emissivities[9]*emissivities[11],
                           emissivities[9]*emissivities[12],
                           emissivities[9]*emissivities[13],
                           emissivities[9]*emissivities[14],
                           emissivities[9]*emissivities[15]],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(shield_inner["layers"] + 1) *emissivities[10]*emissivities[11],
                           emissivities[10]*emissivities[12],
                           emissivities[10]*emissivities[13],
                           emissivities[10]*emissivities[14],
                           emissivities[10]*emissivities[15]],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           emissivities[11]*emissivities[12],
                           emissivities[11]*emissivities[13],
                           emissivities[11]*emissivities[14],
                           emissivities[11]*emissivities[15]],
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           emissivities[12]*emissivities[13],
                           emissivities[12]*emissivities[14],
                           emissivities[12]*emissivities[15]],
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           emissivities[13]*emissivities[14],
                           emissivities[13]*emissivities[15]],
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           emissivities[14]*emissivities[15]],
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                           ]

radiative_relationship_deployed = SB * np.multiply(
    np.multiply(view_factors_deployed, radiative_area), effective_emissivities
)
# radiative_relationship = SB * np.multiply(
#     np.multiply(view_factors, radiative_area), effective_emissivities
# )
conductive_coeff = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [130, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [130, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [130, 130, 130, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [130, 130, 130, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 130, 130, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 130, 130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [130, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

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


conductive_area = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [shell_thickness*spacecraft_width, shell_thickness*spacecraft_width, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [shell_thickness*spacecraft_width, shell_thickness*spacecraft_width,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0, 0, 0, 0, 0, 0, 0],
    [shell_thickness*spacecraft_width,
     shell_thickness*spacecraft_width,
     shell_thickness*spacecraft_length,
     shell_thickness*spacecraft_length, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [shell_thickness*spacecraft_width,
     shell_thickness*spacecraft_width,
     shell_thickness*spacecraft_length,
     shell_thickness*spacecraft_length, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4*np.pi*0.005**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 4*np.pi*0.005**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 4*np.pi*0.005**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 4*np.pi*0.005**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 4*np.pi*0.005**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

conductive_length = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [
        1 / (2 * shell_thickness),
        1 / (2 * shell_thickness),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0
    ],
    [1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 1 / (2 * shell_thickness), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

conductive_relationship = np.multiply(
    np.multiply(conductive_coeff, conductive_area), conductive_length
)
conductive_relationship[conductive_relationship == np.nan] = 0

capacities_matrix = np.diag(np.multiply(masses, capacities))

# node_relationship = conductive_relationship + capacities_matrix + radiative_relationship

node_relationship_deployed = conductive_relationship + capacities_matrix + radiative_relationship_deployed
