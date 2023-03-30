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


masses = [760]
# sailAreas = [6000, 6500, 7000, 7500, 8000, 9000, 9500]
sailAreas = [10000]
# sailAreas = [120]
# sailAreas = [7000,10000]
sailAreas = sailAreas[::-1]

stepSize = 3600
stepSize = 14400
# stepSize = 18000

paramNames = ["mass", "area"]

yearsToRun = 50
yearInSeconds = 365 * 24 * 3600
targetAltitude = 0.48
deepestAltitude = 0.4
FTOP = 0.05727209780778128
initialEpoch = 1117800000
escapeBurnDV = np.array([0.084108927717811,0.240321913434112,4.946831358585164]) * 1000
targetInclination = 52.75

combinations = list(product(masses, sailAreas))

saveFiles = []
resultslst = []
finalEpoch = 0

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
        targetInclination=targetInclination,
        deepestAltitude=deepestAltitude,
        fastTransferOptimiseParameter=FTOP
    )

    finalGuidanceObj, save, saveDep = sim.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=spacecraftName,
        yearsToRun=yearsToRun,
        simStepSize=stepSize,
        initialEpoch=initialEpoch,
        C3BurnVector=escapeBurnDV
    )

    (
        characteristicacceleration,
        spiralDuration,
        inclinationChangeDuration,
        finalInclination,
    ) = finalGuidanceObj.getInclinationChangeDuration()

    finalEpoch = 1117886400 + spiralDuration + inclinationChangeDuration

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
    sim.plotSimulation(*u, quiverEvery=0)

planets = ["Earth", "Venus", "Mercury"]

# for planet in planets:
#     print(f"Running {planet}")
#     sim.simulatePlanet(planet, 1117886400, finalEpoch, stepSize)

plt.show()
