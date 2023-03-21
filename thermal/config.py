import numpy as np
import pandas as pd

'''
config.py

This file contains the inputs for the thermal model. The inputs will be automatically logged.

TODO:
    - Input correct values for configuration matrix.
'''

### General Parameters
spacecraft_radius = 1.5
spacecraft_length = 5
boom_length = 70
boom_diameter = 0.029
sail_area = 10000
panel_area = 18

### Structure Nodes
node_1 = {'name': 'Spacecraft +X', 'area': np.pi*spacecraft_radius*spacecraft_length, 'sun_vf': 1,
          'space_vf': 0, 'external': 'multi-layer insulation', 'temp_range': [-100, 100]}
node_2 = {'name': 'Spacecraft -X', 'area': np.pi*spacecraft_radius*spacecraft_length, 'sun_vf': 0,
          'space_vf': 1, 'external': 'az93 white paint', 'temp_range': [-100, 100]}
node_3 = {'name': 'Spacecraft +Z', 'area': np.pi*spacecraft_radius**2, 'sun_vf': 0.5,
          'space_vf': 0.5, 'external': 'multi-layer insulation', 'temp_range': [-100, 100]}
node_4 = {'name': 'Spacecraft -Z', 'area': np.pi*spacecraft_radius**2, 'sun_vf': 0.5,
          'space_vf': 0.5, 'external': 'multi-layer insulation', 'temp_range': [-100, 100]}
nodes = [node_1, node_2, node_3, node_4]

### Sail Nodes
node_front = {'name': 'Sail Front', 'area': sail_area, 'sun_vf': 1, 'space_vf': 0, 'external': 'aluminum', 'temp_range': [-100, 100]}
node_back = {'name': 'Sail Back', 'area': sail_area, 'sun_vf': 0, 'space_vf': 1, 'external': 'chromium', 'temp_range': [-100, 100]}
sail_nodes = [node_front, node_back]
### External Nodes
booms = {'name': 'Booms', 'area': 4*np.pi*(boom_diameter / 2)*boom_length, 'sun_vf': 0.5, 'space_vf': 0.5,
         'external': 'carbon fibre', 'temp_range': [-100, 100]}
solar_panels = {'name': 'Solar Panels', 'area': panel_area, 'sun_vf': 1, 'space_vf': 0,
                'external': 'osr', 'temp_range': [-100, 100]}
shield = {'name': 'Heat Shield', 'area': node_1['area'], 'sun_vf': 1, 'space_vf': 0,
          'external': 'ceramic cloth', 'temp_range': [-100, 100], 'layers': 7}

### Node Relationships Matrix.
'''
    The setup for this matrix is as follows.
    i/j     1       2       3       4   
    1     c1*m1    F_12    F_13    F_14
    2      G_21    c2*m2   F_23    F_24
    3      G_31    G_32    c3*m3   F_34
    4      G_41    G_42    G_43   c4*m4  

    Here, c_i*m_i are the thermal capacity of each node (the specific heat capacity c times the mass m.)
    Then, G_ij are the conductive coefficients, equal to k_ij * A_ij / L_ij
    Then, R_ij are the radiative coefficients, equal to e_ij * A_ij * F_ij.
'''
node_relationship = [[0, 0, 0, 0, 0.1, 0.1, 1, 1, 1],
                    [0, 0, 0, 0, 0.1, 0.1, 0, 1, 1],
                    [0, 0, 0, 0.25, 0.25, 0.25, 0.25, 1, 1],
                    [0, 0, 0, 0.25, 0.25, 0.25, 0.25, 1, 1],
                    [0, 0.25, 0.25, 0, 0.4, 0.7, 0, 1, 1],
                    [0, 0.25, 0.25, 0.4, 0, 0.4, 0.1, 1, 1],
                    [0.1, 0.25, 0.25, 0.7, 0.4, 0, 0.05, 1, 1],
                    [0, 0.25, 0.25, 0, 0.1, 0.05, 0, 1, 1],
                    [0, 0.25, 0.25, 0, 0.1, 0.05, 0, 1, 1]]
