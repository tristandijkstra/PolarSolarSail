from simulation import simulationV1 as sim
from solarsail.sailPhysicalV3 import SolarSailGuidance
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os

from itertools import product
from tqdm import tqdm


saveDirectory = "data"
if not os.path.exists(saveDirectory):
    os.mkdir(saveDirectory)

logging.basicConfig(
    filename=f"{saveDirectory}/simV0.log",
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
)


masses = [500]
sailAreas = [10000, 12500, 15000, 17500, 20000, 22500]
sailAreas = [10000]
# sailAreas = [7000,10000]
sailAreas = sailAreas[::-1]

stepSize = 3600
stepSize = 36000
stepSize = 18000

paramNames = ["mass", "area"]

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600
targetAltitude = 0.48

combinations = list(product(masses, sailAreas))

saveFiles = []
resultslst = []

logging.info(f"=== Starting run with {len(combinations)} combinations ===")
for combination in tqdm(combinations):
    spacecraftName = f"{paramNames[0]}{combination[0]}_{paramNames[1]}{combination[1]}"

    print(f"Running Combination: mass = {combination[0]} | area = {combination[1]}")
    guidanceObject = SolarSailGuidance(
        None,
        sailName=spacecraftName,
        mass=combination[0],
        sailArea=combination[1],
        targetAltitude=targetAltitude,
        targetInclination=90,
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
    dur = round(inclinationChangeDuration / yearInSeconds, 3)
    totdur = round((spiralDuration + inclinationChangeDuration) / yearInSeconds, 3)
    print(f"Final inclination = {round(finalInclination, 3)} deg")
    print(f"Inclination change duration = {dur} years")

    extraTxt = f"\nFinal inclin. = {round(finalInclination, 3)} | Characteristic Acceleration = {round(characteristicacceleration*1000, 4)} mm/s^2\
                 \nInclin. change duration = {dur} years | Total duration = {totdur} years"
    saveFiles.append([spacecraftName, save, saveDep, extraTxt])

    logStr = f"run: {spacecraftName} | mass = {combination[0]} | area = {combination[1]} | Final inclin. = {round(finalInclination, 3)} | Inclin. change duration = {dur} years"
    logging.info(logStr)

    resultslst.append(
        [
            combination[0],
            combination[1],
            spiralDuration / yearInSeconds,
            inclinationChangeDuration / yearInSeconds,
        ]
    )

for u in saveFiles:
    sim.plotSimulation(*u, quiverEvery=1000)


plt.show()

# things = np.array(resultslst)
# levelss = 100

# fig, ax = plt.subplots(1,1)
# ax.set_title("Spiral Transfer Time")
# ax.set_ylabel("Sail Area [m^2]")
# ax.set_xlabel("Spacecraft Mass [kg]")
# im = ax.tricontour(things[:, 0], things[:, 1], things[:, 2], levels=levelss)
# cbar = plt.colorbar(im, )
# cbar.set_label("Spiral Transfer Time [years]")
# fig.savefig("doc/spiraltime2.png")


# fig2, ax2 = plt.subplots(1,1)
# ax2.set_title("Inclination Change Time")
# ax2.set_ylabel("Sail Area [m^2]")
# ax2.set_xlabel("Spacecraft Mass [kg]")
# im2 = ax2.tricontour(things[:, 0], things[:, 1], things[:, 3], levels=levelss)
# cbar2 = plt.colorbar(im2, )
# cbar2.set_label("Inclination Change Time [years]")
# fig2.savefig("doc/inclinationtime2.png")


# fig3, ax3 = plt.subplots(1,1)
# ax3.set_title("Total Transfer Time")
# ax3.set_ylabel("Sail Area [m^2]")
# ax3.set_xlabel("Spacecraft Mass [kg]")
# im3 = ax3.tricontour(things[:, 0], things[:, 1], things[:, 2] + things[:, 3], levels=levelss)
# cbar3 = plt.colorbar(im3, )
# cbar3.set_label("Total Transfer Time [years]")
# fig3.savefig("doc/totaltransfertime2.png")

# plt.show()
