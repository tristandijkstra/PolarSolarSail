# Load standard modules
import math
import pygmo
import matplotlib
import numpy as np
import pandas as pd
from numpy import random
import time

# Load pygmo library
import pygmo as pg

# from ..solarsail.base import SolarSailGuidanceBase

from typing import Callable, Tuple, Union

# pd.DataFrame(a ba )
# np.array.shape


class SailOptimise:
    yearInSeconds = 31536000

    def __init__(
        self,
        FTOP_min: float,
        FTOP_max: float,
        deepestAltitude_min: float,
        deepestAltitude_max: float,
        solarSailGuidanceObject,
        simuFunction: Callable,
        timesOutwardMax,
        thermalModelObject=None,
        mass=500,
        sailArea=10000,
        targetInclination=90.0,
        yearsToRun=25,
        stepSize=72000,
        targetAltitude=0.48,
        initialEpoch=1117886400,
        C3BurnVector=np.array([0, 0, 0]),
        verbose=False,
        endPrecision=0.001,
    ):
        # Set input arguments as attributes, representaing the problem bounds for both design variables
        self.FTOP_min = FTOP_min
        self.FTOP_max = FTOP_max
        self.deepestAltitude_min = deepestAltitude_min
        self.deepestAltitude_max = deepestAltitude_max
        self.timesOutwardMax = timesOutwardMax
        self.mass = mass
        self.sailArea = sailArea
        self.targetInclination = targetInclination
        self.yearsToRun = yearsToRun
        self.stepSize = stepSize
        self.targetAltitude = targetAltitude
        self.solarSailGuidanceObject = solarSailGuidanceObject
        self.simuFunction = simuFunction
        self.thermalModelObject = thermalModelObject

        self.initialEpoch = initialEpoch
        self.C3BurnVector = C3BurnVector

        self.resultsDict = {}

        self.step = 0
        self.verbose = verbose

        if thermalModelObject is not None:
            self.thermalModel = thermalModelObject(10*3600)
        else:
            self.thermalModel = None

        self.endPrecision = endPrecision

        self.minSofar = 1000

    def __repr__(self) -> str:
        if self.thermalModel is not None:
            return f"SailOptimise V1 | {self.solarSailGuidanceObject} | {self.thermalModel}"
        else:
            return f"SailOptimise V1 | {self.solarSailGuidanceObject}"

    def __str__(self) -> str:
        if self.thermalModel is not None:
            return f"SailOptimise V1 | {self.solarSailGuidanceObject} | {self.thermalModel}"
        else:
            return f"SailOptimise V1 | {self.solarSailGuidanceObject}"

    def get_bounds(self):
        return (
            [self.FTOP_min, self.deepestAltitude_min],
            [self.FTOP_max, self.deepestAltitude_max],
        )

    def fitness(self, x):
        FTOP = x[0]
        deepestAltitude = x[1]
        spacecraftName = "sc_w=" + str(FTOP)

        guidanceObject = self.solarSailGuidanceObject(
            None,
            sailName=spacecraftName,
            mass=self.mass,
            sailArea=self.sailArea,
            targetAltitude=self.targetAltitude,
            targetInclination=self.targetInclination,
            deepestAltitude=deepestAltitude,
            fastTransferOptimiseParameter=FTOP,
            thermalModel=self.thermalModel,
            verbose=False,
            endPrecision=self.endPrecision,
        )

        start = time.perf_counter()
        finalGuidanceObj, save, saveDep = self.simuFunction(
            spacecraftName=spacecraftName,
            sailGuidanceObject=guidanceObject,
            saveFile=None,
            yearsToRun=self.yearsToRun,
            simStepSize=self.stepSize,
            verbose=False,
            initialEpoch=self.initialEpoch,
            C3BurnVector=self.C3BurnVector,
        )

        (
            inclinationChangeDuration,
            lastInclination,
            spiralInclPrecision,
        ) = finalGuidanceObj.getOptimiseOutput()

        incldur = inclinationChangeDuration / SailOptimise.yearInSeconds
        # totdur = round((spiralDuration + inclinationChangeDuration) / yearInSeconds, 3)
        end = time.perf_counter()
        runtime = round(end - start, 2)

        # extraTxt = f"\nFinal inclin. = {round(finalInclination, 3)} | Characteristic Acceleration = {round(characteristicacceleration*1000, 4)} mm/s^2\
        #             \nInclin. change duration = {dur} years | Total duration = {totdur} years"

        self.resultsDict[FTOP] = [
            spacecraftName,
            inclinationChangeDuration / SailOptimise.yearInSeconds,
            save,
            saveDep,
        ]

        expp = 1
            
        precPen = (min(self.endPrecision / spiralInclPrecision, 1))**expp
        precc = (0.5/precPen) - 0.5
        fun = (incldur + max(0, abs(self.targetInclination - lastInclination))**4) + precc
        if (self.thermalModel is not None) and self.thermalModel.thermalFailure:
            fun += 10000
        # if lastInclination < self.targetInclination:
        #     fun = incldur + 1e9
        # elif spiralInclPrecision != 1:
        #     fun = incldur / (min(self.endPrecision / spiralInclPrecision, 1))**expp
        # else:
        #     fun = incldur

        # logStr = f"duration = {runtime} s | run: {spacecraftName} | Final inclin. = {round(lastInclination, 3)} | Inclin. change duration = {incldur} years"
        self.step += 1
        self.minSofar = min(fun, self.minSofar)
        if self.verbose:
            logStr = f"Ev {self.step} | runtime = {runtime} s | FTOP = {round(FTOP,3)} | DA = {round(deepestAltitude, 2)} | prec = {precc} | lastIncl = {round(lastInclination, 1)} | fun = {round(fun, 5)} | min = {round(self.minSofar, 5)}"
            print(logStr)
            if self.thermalModel is not None:
                print()

        return [fun]
