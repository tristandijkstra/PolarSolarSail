import numpy as np
import pandas as pd

'''
config.py

This file contains the inputs for the thermal model. The inputs will be automatically logged.

TODO:
    - View factor to sun?
    - View factor to space?
    - Conductive contact with other nodes?
'''

### General Parameters
spacecraft_radius = 1.5
spacecraft_length = 5
boom_length = 70
boom_diameter = 0.029
sail_area = 10000
panel_area = 18

### Structure Nodes
node_1 = {'name': 'Spacecraft +X', 'area': np.pi*spacecraft_radius*spacecraft_length, 
          'external': 'multi-layer insulation'}
node_2 = {'name': 'Spacecraft -X', 'area': np.pi*spacecraft_radius*spacecraft_length, 
          'external': 'az93 white paint'}
node_3 = {'name': 'Spacecraft +Z', 'area': np.pi*spacecraft_radius**2, 
          'external': 'multi-layer insulation'}
node_4 = {'name': 'Spacecraft -Z', 'area': np.pi*spacecraft_radius**2, 
          'external': 'multi-layer insulation'}
nodes = [node_1, node_2, node_3, node_4]

### Sail Nodes
node_front = {'name': 'Sail Front', 'area': sail_area, 'external': 'aluminum'}
node_back = {'name': 'Sail Back', 'area': sail_area, 'external': 'chromium'}
sail_nodes = [node_front, node_back]
### External Nodes
booms = {'name': 'Booms', 'area': 4*np.pi*(boom_diameter / 2)*boom_length, 
         'external': 'carbon fibre'}
solar_panels = {'name': 'Solar Panels', 'area': panel_area,
                'external': 'osr'}
shield = {'name': 'Heat Shield', 'area': node_1['area'], 
          'external': 'ceramic cloth', 'layers': 7}

### View Factors
view_factors = [[0, 0, 0, 0, 0, 0.1, 0.1, 1],
                [0, 0, 0, 0, 0, 0.1, 0.1, 0],
                [0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25],
                [0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25],
                [0, 0, 0.25, 0.25, 0, 0.4, 0.7, 0],
                [0.1, 0, 0.25, 0.25, 0.4, 0, 0.4, 0.1],
                [0.1, 0.1, 0.25, 0.25, 0.7, 0.4, 0, 0.05],
                [1, 0, 0.25, 0.25, 0, 0.1, 0.05, 0]]