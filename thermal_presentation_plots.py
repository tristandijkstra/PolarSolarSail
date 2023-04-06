import matplotlib.pyplot as plt
import numpy as np

from thermal.thermal_model import Thermal

step_size = 3600

thermal = Thermal(step_size)

interest_nodes = [thermal.metis, thermal.cdm, thermal.batteries, thermal.solar_sail, thermal.solar_panel, thermal.heat_shield]
interest_wch = [-14.41, 7.7, 10.23, 152.71, 149.22, 616.9]
for node in interest_nodes:
    node.temp_range.insert(0, node.temp_range[0] - (abs(node.temp_range[0])*0.5))
    node.temp_range.insert(2, (node.temp_range[1] + node.temp_range[2])/2)
    node.temp_range.append(node.temp_range[-1] + (abs(node.temp_range[0])*0.5))

labels = ["METIS", "CDM", "Batteries", "Solar Sails", "Solar Panels", "Heat Shield"]

c_heatshield = "#ebe2d3"
c_solarpanel = "#d8e6ec"
c_metis = "#e6e6fa"
c_cdm = "#d2f8d2"
c_battery = "#fbbf77"
c_solarsail = "#d3d3d3"

colors = [c_metis, c_cdm, c_battery, c_solarsail, c_solarpanel, c_heatshield]

plt.rcParams.update({'font.size': 19})
fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
bp = ax.boxplot([node.temp_range for node in interest_nodes], vert=0, patch_artist=True, 
                usermedians=interest_wch, showfliers=False, whis=0, showmeans=0, showcaps=0)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
ax.set_yticklabels(labels)
ax.set_xlabel("Temperature (°C)")
ax.xaxis.grid(True)
plt.show()