from simulation import simulationV0
from solarsail.sailIntermediate import SolarSailGuidance
import numpy as np
import pandas as pd

from itertools import product
from tqdm import tqdm

masses = [500, 700]
sailAreas = [15000, 22500]
masses = [500]
sailAreas = [22500/4]

paramNames = ["mass", "area"]

yearsToRun = 7
yearInSeconds = 365 * 24 * 3600

combinations = list(product(masses, sailAreas))

saveFiles = []

for combination in tqdm(combinations):
    spacecraftName = f"{paramNames[0]}{combination[0]}_{paramNames[1]}{combination[1]}"

    print(f"Running Combination: mass = {combination[0]} | area = {combination[1]}")
    guidanceObject = SolarSailGuidance(
        None, sailName=spacecraftName, mass=combination[0], sailArea=combination[1]
    )

    finalGuidanceObj, save, saveDep = simulationV0.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=spacecraftName,
        yearsToRun=yearsToRun,
    )

    inclinationChangeDuration, finalInclination = finalGuidanceObj.getInclinationChangeDuration()
    dur = round(inclinationChangeDuration/yearInSeconds, 3)
    print(f"Final inclination = {round(finalInclination, 3)} deg")
    print(f"Inclination change duration = {dur} years")
    saveFiles.append([save, saveDep])

