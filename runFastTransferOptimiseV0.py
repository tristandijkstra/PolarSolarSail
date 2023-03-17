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
sailArea = 10000
mass = 500
timesOutwardMax = 1

stepSize = 18000

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600

targetInclination = 90

saveFiles = []
resultslst = []


def runSim(
    w,
    timesOutwardMax=timesOutwardMax,
    mass=mass,
    sailArea=sailArea,
    targetInclination=targetInclination,
    yearsToRun=25,
    stepSize=stepSize
):
    
    # print(w)
    spacecraftName = "sc_w=" + str(w)

    guidanceObject = SolarSailGuidance(
        None,
        sailName=spacecraftName,
        mass=mass,
        sailArea=sailArea,
        targetAltitude=targetAltitude,
        targetInclination=targetInclination,
        deepestAltitude=0.2,
        fastTransferOptimiseParameter=w[0]
    )

    finalGuidanceObj, save, saveDep = sim.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=spacecraftName,
        yearsToRun=yearsToRun,
        simStepSize=stepSize,
    )

    (
        inclinationChangeDuration,
        lastInclination,
        timesOutward,
    ) = finalGuidanceObj.getOptimiseOutput()

    incldur = round(inclinationChangeDuration / yearInSeconds, 3)
    # totdur = round((spiralDuration + inclinationChangeDuration) / yearInSeconds, 3)

    logStr = f"run: {spacecraftName} | Final inclin. = {round(lastInclination, 3)} | Inclin. change duration = {incldur} years"
    logging.info(logStr)
    print(logStr)

    resultslst.append(
        [
            w,
            inclinationChangeDuration / yearInSeconds,
            save,
            saveDep
        ]
    )

    if lastInclination < 90:
        return incldur + 10
    elif timesOutward != timesOutwardMax:
        return incldur + 10
    else:
        return incldur


res = scipy.optimize.minimize(runSim, np.array([0.04]))#, bounds=[(0.001, 0.2)])

print(res)

for u in resultslst:
    sim.plotSimulation(satName=u[0], dataFile=u[2], dataDepFile=u[3], quiverEvery=1000)


plt.show()
