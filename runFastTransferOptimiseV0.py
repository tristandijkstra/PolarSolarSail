from simulation import simulationV1 as sim
from solarsail.sailPhysicalV3 import SolarSailGuidance
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os

from itertools import product
from tqdm import tqdm

import scipy
import time


saveDirectory = "data"
if not os.path.exists(saveDirectory):
    os.mkdir(saveDirectory)

logging.basicConfig(
    filename=f"{saveDirectory}/_optimiseV1.log",
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
)

targetAltitude = 0.48
deepestAltitude = 0.2
sailArea = 10000
mass = 500
timesOutwardMax = 1

stepSize = 72000

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600

targetInclination = 90

saveFiles = []
# resultslst = []
resultslst = {}


def runSim(
    w,
    timesOutwardMax=timesOutwardMax,
    mass=mass,
    sailArea=sailArea,
    targetInclination=targetInclination,
    yearsToRun=25,
    stepSize=stepSize,
):
    # print(w)
    w = w[0]
    spacecraftName = "sc_w=" + str(w)


    guidanceObject = SolarSailGuidance(
        None,
        sailName=spacecraftName,
        mass=mass,
        sailArea=sailArea,
        targetAltitude=targetAltitude,
        targetInclination=targetInclination,
        deepestAltitude=deepestAltitude,
        fastTransferOptimiseParameter=w,
        verbose=False
    )

    start = time.perf_counter()
    finalGuidanceObj, save, saveDep = sim.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=None,
        yearsToRun=yearsToRun,
        simStepSize=stepSize,
        verbose=False,
    )

    (
        inclinationChangeDuration,
        lastInclination,
        timesOutward,
    ) = finalGuidanceObj.getOptimiseOutput()

    incldur = inclinationChangeDuration / yearInSeconds
    # totdur = round((spiralDuration + inclinationChangeDuration) / yearInSeconds, 3)
    end = time.perf_counter()
    runtime = round(end - start, 2)
    logStr = f"duration = {runtime} s | run: {spacecraftName} | Final inclin. = {round(lastInclination, 3)} | Inclin. change duration = {incldur} years"
    logging.info(logStr)
    print(logStr)

    # extraTxt = f"\nFinal inclin. = {round(finalInclination, 3)} | Characteristic Acceleration = {round(characteristicacceleration*1000, 4)} mm/s^2\
    #             \nInclin. change duration = {dur} years | Total duration = {totdur} years"

    resultslst[w] = [spacecraftName, inclinationChangeDuration / yearInSeconds, save, saveDep]
    # resultslst.append({str(w): [spacecraftName, inclinationChangeDuration / yearInSeconds, save, saveDep]})

    if lastInclination < 90:
        return incldur + 100
    elif timesOutward != timesOutwardMax:
        return incldur + 100
    else:
        return incldur


res = scipy.optimize.minimize(
    runSim, np.array([0.001]), method="Powell", bounds=[(0.001, 0.2)]
)
# res = scipy.optimize.minimize(
#     runSim, np.array([0.1]), method="Nelder-Mead", bounds=[(0.001, 0.2)]
# )

print(res)

# bestValue = str(res.x)
bestValue = res.x[0]

print(resultslst[bestValue])
saveFiel = "w=" + str(bestValue)
guidanceObject = SolarSailGuidance(
    None,
    sailName="best",
    mass=mass,
    sailArea=sailArea,
    targetAltitude=targetAltitude,
    targetInclination=targetInclination,
    deepestAltitude=deepestAltitude,
    fastTransferOptimiseParameter=bestValue,
    verbose=True
)

finalGuidanceObj, save, saveDep = sim.simulate(
    spacecraftName="best",
    sailGuidanceObject=guidanceObject,
    saveFile="w=" + str(bestValue),
    yearsToRun=yearsToRun,
    simStepSize=36000,
    verbose=True,
)
sim.plotSimulation(satName=bestValue, dataFile=save, dataDepFile=saveDep, quiverEvery=1000)


plt.show()
