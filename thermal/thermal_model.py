'''
thermal_model.py

This code runs the single-node thermal analysis at different time steps throughout the orbit.

-TODO:
    - Find a way to update the angle to the Sun.

'''
import node as nd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm
from budget_properties import Properties

################################################################################################################
#############################################  EDITABLE PARAMETERS #############################################
################################################################################################################
AU = 149.6e9
dt = 1000                                   # How much to timeskip the orbital sim (to save computational time)
yearInSeconds = 365 * 24 * 3600
solar_luminosity = 3.83e26
max_temp = 343.0                            # maximum allowable temperature [K]
min_temp = 243.0                            # minimum allowable temperature [K]
heaters_power = 0                           # How much heat power is coming from heaters in watts
coolers_power = 0                           # How much heat is being rejected by coolers in watts
electrical_heat = 800                       # How much electrical waste heat is being generated (taken from power budget)
data_file = 'data/mass680_area10000.dat'    # Main data file
data_dep_file = 'data/mass680_area10000_dep.dat' # Dependent data file
data = pd.read_csv(
    data_file, delimiter="	", names=["time", "x", "y", "z", "vx", "vy", "vz"], header=None
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
]
data2 = pd.read_csv(data_dep_file, delimiter="	", names=depVars, 
                    header=None)

logging.basicConfig(filename=f'data/thermal.log', encoding='utf-8', level=logging.INFO, format='%(asctime)s | %(message)s')


# pointing_angle = pd.DataFrame(index=range(len(data.axes[0])),columns=range(len(data.axes[1])))
# for index, row in data.iterrows():
#     for index2, row2 in data2.iterrows():
#         if index < 657:
#             index = index*1000
#             index2 = index*1000
#             pointing_angle[index] = np.arccos(np.dot([row2.ThrustX, row2.ThrustY, row2.ThrustZ], [-row.x, -row.y, -row.z]))
#         else:
#             break

time = (data.time - data.time.iloc[0])/yearInSeconds
time_interp = time[::dt]
sun_dist = data.altitude.iloc[::dt]



no_coatings = [[0.4, 0.30, 0.04, 0.95], [7, 7, 0]]
coatings = [[0.90, 0.30, 0.04, 0.95], [7, 0.5, 6.5]] # Sample inputs

shield_layers = 7                      # How many layers of shielding would you like to analyze?
temps = []

for idx, time_step in enumerate(tqdm(time_interp)):
    # Simplistic format for now
    if sun_dist[idx*dt] > 0.75:
        angle = np.pi / 2
    else:
        angle = 0 
    thermal_case = [sun_dist[idx*dt], angle]
    # Please read Node() class __init__ for details on how you could change these numbers.
    spacecraft_bus = nd.Node("Spacecraft Bus", no_coatings[0], no_coatings[1], thermal_case)
    solar_sail = nd.Node("Solar Sail", [0.63, 0.12, 0.88, 0.00], [5000, 5000, 0], thermal_case)
    if shield_layers > 0:
        heat_shield = nd.Node("Heat Shield", [0.44, 0.96, 0.16, 0.00], [6, 6, 0], thermal_case, shield_layers)
    else:
        heat_shield = None
    ################################################################################################################
    #########################################  END OF EDITABLE PARAMETERS ##########################################
    ################################################################################################################
    spacecraft = [spacecraft_bus, solar_sail, heat_shield]
    for node in spacecraft:
        if node != None:
            flux_in = node.solar_heat_in()
            q_internal = node.internal_heat(heaters_power, electrical_heat)
    temp_bus, temp_shield, temp_sail = nd.temperatures(spacecraft_bus, solar_sail, heat_shield)
    temp_spacecraft = [temp_bus, temp_sail, temp_shield]
    temps.append(temp_bus)
    for idx, node in enumerate(spacecraft):
        if node != None:
            # print("Solar Flux: " + str(flux_in) + " W")
            # print("Internal Heat: " + str(q_internal) + " W")
            # print("Surface Temperature: " + str(temp_spacecraft[idx]) + " K")
            heat_added = node.heat_absorbed(temp_spacecraft[idx])
            if node.name == "Spacecraft Bus":
                # print("Heat Added - " + node.name + " " + str(heat_added) + " W")
                heat_lost = node.heat_radiated(temp_spacecraft[idx])
                # print("Heat Lost - " + node.name + " " + str(heat_lost) + " W")
                thermal_balance = node.thermal_balance()
                # print("Thermal Balance - " + node.name + " " + str(thermal_balance) + " W")

# Plotting Results
fig, ax = plt.subplots(2, 1)
ax[0].set_title(f"Spacecraft Surface Temperature - Shielding Layers: {shield_layers}")
ax[0].plot(time_interp, temps)
ax[0].set_ylabel("Temperature [K]")
ax[0].set_xlabel("Time [years]")
ax[1].plot(time_interp, sun_dist)
ax[1].set_title("Altitude [AU]")
ax[1].set_ylabel("Altitude [AU]")
ax[1].set_xlabel("Time [years]")
ax[1].plot()

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