from simulation import simulationV1 as sim
from solarsail.sailPhysicalV2 import SolarSailGuidance
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import logging
import os

from itertools import product
from tqdm import tqdm


saveDirectory = "data"
if not os.path.exists(saveDirectory):
    os.mkdir(saveDirectory)

logging.basicConfig(filename=f'{saveDirectory}/simV0.log', encoding='utf-8', level=logging.INFO, format='%(asctime)s | %(message)s')

masses = [100, 200, 300, 400, 500, 600, 700, 800]
sailAreas = [5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500]

# masses = [100, 200, 300]
# sailAreas = [17500, 20000, 22500]

stepSize = 3600
stepSize = 36000 * 3

paramNames = ["mass", "area"]

yearsToRun = 40
yearInSeconds = 365 * 24 * 3600
targetAltitude = 0.5

combinations = list(product(masses, sailAreas))

saveFiles = []
resultslst = []

logging.info(f"=== Starting run with {len(combinations)} combinations ===")
for combination in tqdm(combinations):
    spacecraftName = f"{paramNames[0]}{combination[0]}_{paramNames[1]}{combination[1]}"

    print(f"Running Combination: mass = {combination[0]} | area = {combination[1]}")
    guidanceObject = SolarSailGuidance(
        None, sailName=spacecraftName, mass=combination[0], sailArea=combination[1], targetAltitude=targetAltitude
    )

    finalGuidanceObj, save, saveDep = sim.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=spacecraftName,
        yearsToRun=yearsToRun,
        simStepSize=stepSize,
    )

    (
        characteristicacceleration,
        spiralDuration,
        inclinationChangeDuration,
        finalInclination,
    ) = finalGuidanceObj.getInclinationChangeDuration()
    dur = round(inclinationChangeDuration/yearInSeconds, 3)
    print(f"Final inclination = {round(finalInclination, 3)} deg")
    print(f"Inclination change duration = {dur} years")

    extraTxt = f"\nFinal inclin. = {round(finalInclination, 3)} | Inclin. change duration = {dur} years"
    saveFiles.append([spacecraftName, save, saveDep, extraTxt])

    logStr = f"run: {spacecraftName} | mass = {combination[0]} | area = {combination[1]} | Final inclin. = {round(finalInclination, 3)} | Inclin. change duration = {dur} years"
    logging.info(logStr)

    resultslst.append([combination[0], combination[1], spiralDuration/yearInSeconds, inclinationChangeDuration/yearInSeconds])

# for u in saveFiles:
#     sim.plotSimulation(*u, quiverEvery=200)


# plt.show()
plt.style.use('Solarize_Light2')
things = np.array(resultslst)
levelss = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 14, 15, 16, 18, 20, 21, 22.5, 25, 30, 35]
fig, ax = plt.subplots(1,1, figsize=(18, 12))
ax.set_title("Spiral Transfer Time (years)")
ax.set_ylabel("Sail Area [m^2]")
ax.set_xlabel("Spacecraft Mass [kg]")
norm1=colors.LogNorm(vmin=things[:, 2].min(), vmax=things[:, 2].max())
im = ax.tricontour(things[:, 0], things[:, 1], things[:, 2], levels=levelss, norm=norm1, cmap="winter")
ax.clabel(im, inline=True, fontsize=10)
ax.grid(b=True, which='major', color='gray', linestyle='-')
# cbar = plt.colorbar(im, )
# cbar.set_label("Spiral Transfer Time [years]")
fig.set_tight_layout(True)
fig.savefig("data/_spiraltime.png")


fig2, ax2 = plt.subplots(1,1, figsize=(18, 12))
norm2=colors.LogNorm(vmin=things[:, 3].min(), vmax=things[:, 3].max())
ax2.set_title("Inclination Change Time (years)")
ax2.set_ylabel("Sail Area [m^2]")
ax2.set_xlabel("Spacecraft Mass [kg]")
im2 = ax2.tricontour(things[:, 0], things[:, 1], things[:, 3], levels=levelss, norm=norm2, cmap="winter")
ax2.clabel(im2, inline=True, fontsize=10)
ax2.grid(b=True, which='major', color='gray', linestyle='-')
# cbar2 = plt.colorbar(im2, )
# cbar2.set_label("Inclination Change Time [years]")
fig2.set_tight_layout(True)
fig2.savefig("data/_inclinationtime.png")


fig3, ax3 = plt.subplots(1,1, figsize=(18, 12))
ax3.set_title("Total Transfer Time (years)")
ax3.set_ylabel("Sail Area [m^2]")
ax3.set_xlabel("Spacecraft Mass [kg]")
norm3=colors.LogNorm(vmin=(things[:, 2] + things[:, 3]).min(), vmax=(things[:, 2] + things[:, 3]).max())
im3 = ax3.tricontour(things[:, 0], things[:, 1], things[:, 2] + things[:, 3], levels=levelss, norm=norm2, cmap="winter")
ax3.clabel(im3, inline=True, fontsize=10)
ax3.grid(b=True, which='major', color='gray', linestyle='-')
# cbar3 = plt.colorbar(im3, )
# cbar3.set_label("Total Transfer Time [years]")
fig3.set_tight_layout(True)
fig3.savefig("data/_totaltransfertime.png")

plt.show()