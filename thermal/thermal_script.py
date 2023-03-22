'''
thermal_script.py

This code runs the single-node thermal analysis at different time steps throughout the orbit, independently
from the orbital simulation.

-TODO:

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

AU = 149.6e9
yearInSeconds = 365 * 24 * 3600
solar_luminosity = 3.83e26

################################################################################################################
#############################################  EDITABLE PARAMETERS #############################################
################################################################################################################
dt = 1000                                   # How much to timeskip the orbital sim (to save computational time)
heaters_power = 0                           # How much heat power is coming from heaters [W]
coolers_power = 0                           # How much heat is being rejected by coolers [W]
electrical_heat = 800                       # How much electrical waste heat is being generated [W]
data_file = 'data/mass500_area10000.dat'    # Main data file
data_dep_file = 'data/mass500_area10000_dep.dat' # Dependent data file
log_file = 'data/thermal.log'               # Logging file path

################################################################################################################
#########################################  END OF EDITABLE PARAMETERS ##########################################
################################################################################################################

logging.basicConfig(filename=log_file, encoding='utf-8', level=logging.INFO, format='%(asctime)s | %(message)s')

# By default, the solar sail, spacecraft structure, booms, heat shield, and solar panels are written as nodes.
# The main spacecraft structure can then be split into additional nodes (up to 4).
node_structure = len(config.nodes)          # The number of nodes you would like placed on the spacecraft structure

shield_layers = config.shield['layers']     # How many layers of shielding would you like to analyze?
logging.info("**********************************")
logging.info("==================================")
logging.info("New Configuration")
logging.info("==================================")
logging.info("**********************************")
# Splitting the spacecraft into nodes
bus_materials = []
str_properties = []
sail_materials = []
sail_properties = []
total_properties = []
logging.info("Configuration Order: [emissivity, reflectivity, absorptivity, density, area]")
for i, node in enumerate(config.nodes):
    bus_materials.append(materials.bus_material(node['external']))
    str_properties.append([bus_materials[i]['emissivity'], bus_materials[i]['reflectivity'], 
                          bus_materials[i]['absorptivity'], bus_materials[i]['density'], 
                          node['area'], node['sun_vf'], node['space_vf']])
    logging.info("==================================")
    logging.info("Spacecraft Structure Configuration")
    logging.info("==================================")
    logging.info(str_properties[i])
total_properties.append(str_properties)
for i, node in enumerate(config.sail_nodes):
    sail_materials.append(materials.sail_material(node['external']))
    sail_properties.append([sail_materials[i]['emissivity'], sail_materials[i]['reflectivity'], 
                            sail_materials[i]['absorptivity'], sail_materials[i]['density'], 
                            node['area'], node['sun_vf'], node['space_vf']])
    logging.info("==================================")
    logging.info("Solar Sail Configuration")
    logging.info("==================================")
    logging.info(sail_properties[i])
total_properties.append(sail_properties)
panel = materials.panel_material(config.solar_panels['external'])
panel_properties = [panel['emissivity'], panel['reflectivity'], panel['absorptivity'], 
                    panel['density'], config.solar_panels['area'], config.solar_panels['sun_vf'],
                    config.solar_panels['space_vf']]
logging.info("==================================")
logging.info("Panel Configuration")
logging.info("==================================")
logging.info(panel_properties)
total_properties.append(panel_properties)
boom = materials.boom_material(config.booms['external'])
boom_properties = [boom['emissivity'], boom['reflectivity'], boom['absorptivity'], boom['density'], 
                   config.booms['area'], config.booms['sun_vf'], config.booms['space_vf']]
logging.info("==================================")
logging.info("Boom Configuration")
logging.info("==================================")
logging.info(boom_properties)
total_properties.append(boom_properties)
shield = materials.shield_material(config.shield['external'])
shield_properties = [shield['emissivity'], shield['reflectivity'], shield['absorptivity'],
                     shield['density'], config.shield['area'], config.shield['sun_vf'], 
                     config.shield['space_vf']]
logging.info("==================================")
logging.info("Shield Configuration")
logging.info("==================================")
logging.info(shield_properties)
total_properties.append(shield_properties)

relationships = np.asarray(config.node_relationship)
logging.info("=============================================")
logging.info("Relationship Matrix")
logging.info("=============================================")
logging.info(relationships)

logging.info("**********************************")


total_nodes = len(relationships)

# if shield == None:
#     heat_shield_A = 0
#     shield_layers = 0

# if shield_layers == 0:
#     heat_shield_A = 0

data = pd.read_csv(
    data_file, delimiter="	", names=["time", "x", "y", "z", "vx", "vy", "vz"], header=None,
    skiprows=lambda i: i % dt,
).assign(altitude=lambda x: np.sqrt(x.x**2 + x.y**2 + x.z**2) / AU)
depVars = [
    "time",
    "ThrustX",
    "ThrustY",
    "ThrustZ",
    "ThrustMagnitude",
    "a",
    "e",
    "i",
    "omega",
    "RAAN",
    "theta",
    "cone",
    "clock",
]
data2 = pd.read_csv(
        data_dep_file,
        delimiter="	",
        names=depVars,
        header=None,
        skiprows=lambda i: i % dt,
    )

time = (data.time - data.time.iloc[0])/yearInSeconds
time_interp = time
sun_dist = data.altitude
angle = data2.cone

node_temperatures = np.zeros((len(time_interp), total_nodes))

for idx, time_step in enumerate(tqdm(time_interp)):
    thermal_case = [sun_dist[idx], angle[idx]]
    spacecraft = []
    spacecraft_bus = []
    solar_sail = []  
    for node, node_properties in enumerate(str_properties):
        spacecraft_bus.append(nd.Node("Spacecraft Bus", node_properties, thermal_case))
    spacecraft.append(spacecraft_bus)
    for node, node_properties in enumerate(sail_properties):
        solar_sail.append(nd.Node("Solar Sail", node_properties, thermal_case))
    spacecraft.append(solar_sail)
    spacecraft = [item for sublist in spacecraft for item in sublist]
    solar_panel = nd.Node("Solar Panel", panel_properties, thermal_case)
    spacecraft.append(solar_panel)
    boom = nd.Node("Boom", boom_properties, thermal_case)
    spacecraft.append(boom)
    if shield_layers > 0:
        heat_shield = nd.Node("Heat Shield", shield_properties, thermal_case, shield_layers)
    else:
        heat_shield = None
    spacecraft.append(heat_shield)
    for node in spacecraft:
        if node != None:
            flux_in = node.solar_heat_in()
    node_temp_step = nd.temperatures(spacecraft, relationships, dt)
    node_temperatures[idx, :] = node_temp_step

# Plotting Results
fig, ax = plt.subplots(3, 2)
ax[0][0].set_title(f"Spacecraft Bus Temperatures")
ax[0][0].plot(time_interp, node_temperatures[:, 0], label='Spacecraft +X')
ax[0][0].plot(time_interp, node_temperatures[:, 1], label='Spacecraft -X')
ax[0][0].plot(time_interp, node_temperatures[:, 2], label='Spacecraft +Z')
ax[0][0].plot(time_interp, node_temperatures[:, 3], label='Spacecraft -Z')
ax[0][0].plot(time_interp, config.node_1['temp_range'][0]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed', label='Temperature Limit')
ax[0][0].plot(time_interp, config.node_1['temp_range'][1]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed')
ax[0][0].set_ylabel("Temperature [K]")
ax[0][0].set_xlabel("Time [years]")
ax[0][0].legend()
ax[0][1].set_title(f"Solar Sail Temperatures")
ax[0][1].plot(time_interp, node_temperatures[:, 4], label='Sail Front')
ax[0][1].plot(time_interp, node_temperatures[:, 5], label='Sail Back')
ax[0][1].plot(time_interp, config.node_front['temp_range'][0]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed', label='Temperature Limit')
ax[0][1].plot(time_interp, config.node_back['temp_range'][1]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed')
ax[0][1].set_ylabel("Temperature [K]")
ax[0][1].set_xlabel("Time [years]")
ax[0][1].legend()
ax[1][0].set_title(f"Boom Temperatures")
ax[1][0].plot(time_interp, node_temperatures[:, 6], label='Booms')
ax[1][0].plot(time_interp, config.booms['temp_range'][0]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed', label='Temperature Limit')
ax[1][0].plot(time_interp, config.booms['temp_range'][1]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed')
ax[1][0].set_ylabel("Temperature [K]")
ax[1][0].set_xlabel("Time [years]")
ax[1][0].legend()
ax[1][1].set_title(f"Solar Panel Temperatures")
ax[1][1].plot(time_interp, node_temperatures[:, 7], label='Solar Panels')
ax[1][1].plot(time_interp, config.solar_panels['temp_range'][0]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed', label='Temperature Limit')
ax[1][1].plot(time_interp, config.solar_panels['temp_range'][1]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed')
ax[1][1].set_ylabel("Temperature [K]")
ax[1][1].set_xlabel("Time [years]")
ax[1][1].legend()
ax[2][1].set_title(f"Heat Shield Temperatures")
ax[2][1].plot(time_interp, node_temperatures[:, 8], label='Heat Shield')
ax[2][1].plot(time_interp, config.shield['temp_range'][0]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed', label='Temperature Limit')
ax[2][1].plot(time_interp, config.shield['temp_range'][1]*np.ones(np.shape(time_interp)), 'r', linestyle='dashed')
ax[2][1].set_ylabel("Temperature [K]")
ax[2][1].set_xlabel("Time [years]")
ax[2][1].legend()
fig.set_tight_layout(True)

plt.savefig(f"data/thermal_{data_file[5:-4]}.png")

plt.show()

# Mass Estimation
orbiter = Properties()
passive_mass_estimate = orbiter.passive_thermal_mass(shield_layers)
cost_estimate = orbiter.cost()

# Log Information
logging.info("===========================================================")
logging.info("======================= SUMMARY ===========================")
logging.info("===========================================================")
log_input = f"Shielding Layers: {shield_layers}"
logging.info(log_input)
log_mass = f"Total Thermal Passive Mass: {str(round(passive_mass_estimate, 2))} kg"
log_cost = f"Total Materials Cost: € {str(round(cost_estimate, 2))} M"
logging.info(log_mass)
logging.info(log_cost)

print("===========================================================")
print("======================= SUMMARY ===========================")
print("===========================================================")
print(log_mass)
print(log_cost)