from simulation import simulationV0, simulationV1
from solarsail.sailPhysicalV1 import SolarSailGuidance
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


masses = [590]
sailAreas = [10000]
stepSize = 36000

paramNames = ["mass", "area"]

yearsToRun = 20
yearInSeconds = 365 * 24 * 3600

combinations = list(product(masses, sailAreas))

saveFiles = []

logging.info(f"=== Starting run with {len(combinations)} combinations ===")
for combination in tqdm(combinations):
    spacecraftName = f"{paramNames[0]}{combination[0]}_{paramNames[1]}{combination[1]}"

    print(f"Running Combination: mass = {combination[0]} | area = {combination[1]}")
    guidanceObject = SolarSailGuidance(
        None,
        sailName=spacecraftName,
        mass=combination[0],
        sailArea=combination[1],
        targetAltitude=0.48,
        targetInclination=70
        # characteristicAcceleration=0.3027,
    )

    finalGuidanceObj, save, saveDep = simulationV0.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=spacecraftName,
        yearsToRun=yearsToRun,
        simStepSize=stepSize,
    )

    (
        inclinationChangeDuration,
        finalInclination,
    ) = finalGuidanceObj.getInclinationChangeDuration()
    dur = round(inclinationChangeDuration / yearInSeconds, 3)
    print(f"Final inclination = {round(finalInclination, 3)} deg")
    print(f"Inclination change duration = {dur} years")

    extraTxt = f"\nFinal inclin. = {round(finalInclination, 3)} | Inclin. change duration = {dur} years"
    saveFiles.append([spacecraftName, save, saveDep, extraTxt])

    logStr = f"run: {spacecraftName} | mass = {combination[0]} | area = {combination[1]} | Final inclin. = {round(finalInclination, 3)} | Inclin. change duration = {dur} years"
    logging.info(logStr)

for u in saveFiles:
    simulationV1.plotSimulation(*u, quiverEvery=0)

# plt.axis(False)
plt.show()
