from simulation import simulationV1 as sim
from solarsail.sailPhysicalV3 import SolarSailGuidance
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os

from itertools import product
from tqdm import tqdm
from thermal.thermal_model import Thermal

saveDirectory = "data"
if not os.path.exists(saveDirectory):
    os.mkdir(saveDirectory)

logging.basicConfig(
    filename=f"{saveDirectory}/simV0.log",
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
)


masses = [590]
sailAreas = [10000]

stepSize = 3600
stepSize = 36000
# stepSize = 100000

paramNames = ["mass", "area"]

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600
targetAltitude = 0.40
deepestAltitude = 0.37

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
        targetInclination=68,
        deepestAltitude=deepestAltitude,
        fastTransferOptimiseParameter=0.04,
        thermalModel=Thermal(10*3600)
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
    saveFiles.append([spacecraftName, save, saveDep, extraTxt, finalGuidanceObj.thermalModel.node_keys])

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
    sim.plotSimulation(*u[0:4], quiverEvery=1000, thermalOn=True)
    sim.plotThermal(u[0], u[2], u[4])


plt.show()
