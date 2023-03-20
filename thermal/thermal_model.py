'''
thermal_model.py

This code runs the single-node thermal analysis at different time steps throughout the orbit.

-TODO:
    - Figure out how to do multiple nodes for the spacecraft bus (in a modular fashion).

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
max_temp = 343.0                            # maximum allowable temperature [K]
min_temp = 243.0                            # minimum allowable temperature [K]
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
logging.info("Configuration Order: [emissivity, reflectivity, absorptivity, density, area]")
for i, node in enumerate(config.nodes):
    bus_materials.append(materials.bus_material(node['external']))
    str_properties.append([bus_materials[i]['emissivity'], bus_materials[i]['reflectivity'], 
                          bus_materials[i]['absorptivity'], bus_materials[i]['density'], 
                          node['area']])
    logging.info("==================================")
    logging.info("Spacecraft Structure Configuration")
    logging.info("==================================")
    logging.info(str_properties[i])
for i, node in enumerate(config.sail_nodes):
    sail_materials.append(materials.sail_material(node['external']))
    sail_properties.append([sail_materials[i]['emissivity'], sail_materials[i]['reflectivity'], 
                            sail_materials[i]['absorptivity'], sail_materials[i]['density'], 
                            node['area']])
    logging.info("==================================")
    logging.info("Solar Sail Configuration")
    logging.info("==================================")
    logging.info(sail_properties[i])

shield = materials.shield_material(config.shield['external'])
shield_properties = [shield['emissivity'], shield['reflectivity'], shield['absorptivity'],
                     shield['density'], config.shield['area']]
logging.info("==================================")
logging.info("Shield Configuration")
logging.info("==================================")
logging.info(shield_properties)
panel = materials.panel_material(config.solar_panels['external'])
panel_properties = [panel['emissivity'], panel['reflectivity'], panel['absorptivity'], 
                    panel['density'], config.solar_panels['area']]
logging.info("==================================")
logging.info("Panel Configuration")
logging.info("==================================")
logging.info(panel_properties)
boom = materials.boom_material(config.booms['external'])
boom_properties = [boom['emissivity'], boom['reflectivity'], boom['absorptivity'], boom['density'], 
                   config.booms['area']]
logging.info("==================================")
logging.info("Boom Configuration")
logging.info("==================================")
logging.info(boom_properties)

vfs = config.view_factors
logging.info("=============================================")
logging.info("View Factor Matrix (in order of node listing)")
logging.info("=============================================")
logging.info(vfs)

logging.info("**********************************")

if shield == None:
    heat_shield_A = 0
    shield_layers = 0

if shield_layers == 0:
    heat_shield_A = 0

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

temps = []
spacecraft_bus = []

for idx, time_step in enumerate(tqdm(time_interp)):
    thermal_case = [sun_dist[idx], angle[idx]]
    for node, node_properties in enumerate(str_properties):
        spacecraft_bus.append(nd.Node("Spacecraft Bus", node_properties, thermal_case))
    solar_sail = nd.Node("Solar Sail", sail_properties, thermal_case)
    solar_panel = nd.Node("Solar Panel", panel_properties, thermal_case)
    boom = nd.Node("Boom", boom_properties, thermal_case)
    if shield_layers > 0:
        heat_shield = nd.Node("Heat Shield", shield_properties, thermal_case, shield_layers)
    else:
        heat_shield = None
    spacecraft = [spacecraft_bus, solar_sail, solar_panel, boom, heat_shield]
    for node in spacecraft:
        if node != None:
            flux_in = node.solar_heat_in()
            q_internal = node.internal_heat(heaters_power, electrical_heat)
    temp_bus, temp_shield, temp_sail = nd.temperatures(spacecraft_bus, solar_sail, heat_shield)
    temp_spacecraft = [temp_bus, temp_sail, temp_shield]
    temps.append(temp_bus)
    for idx, node in enumerate(spacecraft):
        if node != None:
            heat_added = node.heat_absorbed(temp_spacecraft[idx])
            if node.name == "Spacecraft Bus":
                heat_lost = node.heat_radiated(temp_spacecraft[idx])
                thermal_balance = node.thermal_balance()

lower_temp_limit = min_temp*np.ones(np.shape(time_interp))
upper_temp_limit = max_temp*np.ones(np.shape(time_interp))
# Plotting Results
fig, ax = plt.subplots(2, 1)
ax[0].set_title(f"Spacecraft Surface Temperature - Shielding Layers: {shield_layers}")
ax[0].plot(time_interp, temps, label='Spacecraft Surface Temperature')
limit = ax[0].plot(time_interp, lower_temp_limit, 'r', linestyle='dashed', label='Temperature Limit')
ax[0].plot(time_interp, upper_temp_limit, 'r', linestyle='dashed')
ax[0].set_ylabel("Temperature [K]")
ax[0].set_xlabel("Time [years]")
ax[0].legend()
ax[1].plot(time_interp, sun_dist)
ax[1].set_title("Altitude [AU]")
ax[1].set_ylabel("Altitude [AU]")
ax[1].set_xlabel("Time [years]")
ax[1].plot()
fig.set_tight_layout(True)

plt.savefig(f"data/thermal_{data_file[5:-4]}.png")

plt.show()

# Mass Estimation
orbiter = Properties()
passive_mass_estimate = orbiter.passive_thermal_mass(shield_layers)
cost_estimate = orbiter.cost()

# Distance at maximum temperature
distance_max = sun_dist[temps.index(max(temps))*dt]

# Percent Time Estimation
percent_outside_limits = ((np.asarray(temps) > max_temp).sum() + (np.asarray(temps) < min_temp).sum()) / len(temps) * 100

# Log Information
logging.info("===========================================================")
logging.info("======================= SUMMARY ===========================")
logging.info("===========================================================")
log_input = f"Shielding Layers: {shield_layers}"
logging.info(log_input)
log_temp = f"Maximum Temperature Achieved: {str(round(np.max(temps), 2))} K"
log_dist = f"Distance of Maximum Temperature Achieved: {str(round(distance_max, 2))} AU"
log_percent = f"Percent of Transit Outside Survivable Limits: {str(round(percent_outside_limits, 2))} %"
log_mass = f"Total Thermal Passive Mass: {str(round(passive_mass_estimate, 2))} kg"
log_cost = f"Total Materials Cost: € {str(round(cost_estimate, 2))} M"
logging.info(log_temp)
logging.info(log_dist)
logging.info(log_percent)
logging.info(log_mass)
logging.info(log_cost)

print("===========================================================")
print("======================= SUMMARY ===========================")
print("===========================================================")
print(log_temp)
print(log_dist)
print(log_percent)
print(log_mass)
print(log_cost)