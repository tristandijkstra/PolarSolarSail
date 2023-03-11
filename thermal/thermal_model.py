import node as nd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

################################################################################################################
#############################################  EDITABLE PARAMETERS #############################################
################################################################################################################
AU = 149.6e9
yearInSeconds = 365 * 24 * 3600
solar_luminosity = 3.83e26
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
time_interp = time[::1000]
sun_dist = data.altitude.iloc[::1000]



no_coatings = [[0.4, 0.30, 0.04, 0.95], [7, 7, 0]]
coatings = [[0.90, 0.30, 0.04, 0.95], [7, 0.5, 6.5]] # Sample inputs

shield_layers = 0                      # How many layers of shielding would you like to analyze?
temps = []

for idx, time_step in enumerate(tqdm(time_interp)):
    thermal_case = [sun_dist[idx*1000], 0]
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


fig, ax = plt.subplots(2, 1)
ax[0].set_title("Spacecraft Surface Temperature - Shielding Layers: 0")
ax[0].plot(time_interp, temps)
ax[0].set_ylabel("Temperature [K]")
ax[0].set_xlabel("Time [years]")
ax[1].plot(time_interp, sun_dist)
ax[1].set_title("Altitude [AU]")
ax[1].set_ylabel("Altitude [AU]")
ax[1].set_xlabel("Time [years]")
ax[1].plot()

plt.show()

print("===========================================================")
print("======================= SUMMARY ===========================")
print("===========================================================")
print("Maximum Temperature Achieved: " + str(np.max(temps)) + " K")