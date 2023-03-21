from simulation import simulationV1 as sim
from solarsail.sailPhysicalV2_2 import SolarSailGuidance
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



stepSize = 3600
stepSize = 36000
stepSizes = [360, 720, 3600, 7200, 36000, 72000,]
stepSizes = [2250, 4500, 9000, 18000, 36000, 72000, 144000, 288000]
# skipSizes = [int(x/min(stepSizes)) for x in stepSizes]
skipSizes = [int(x/min(stepSizes)) for x in stepSizes][::-1]
# stepSize = 18000

paramNames = ["mass", "area"]

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600
targetAltitude = 0.48
deepestAltitude = targetAltitude


saveFiles = []
resultslst = []

for step in tqdm(stepSizes):
    spacecraftName = f"rk4_{step}"

    print(f"Running step: {step}")
    guidanceObject = SolarSailGuidance(
        None,
        sailName=spacecraftName,
        mass=500,
        sailArea=10000,
        targetAltitude=targetAltitude,
        targetInclination=68,
        deepestAltitude=deepestAltitude,
        # fastTransferOptimiseParameter=0.04
    )

    finalGuidanceObj, save, saveDep = sim.simulate(
        spacecraftName=spacecraftName,
        sailGuidanceObject=guidanceObject,
        saveFile=spacecraftName,
        yearsToRun=yearsToRun,
        simStepSize=step,
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

    resultslst.append(
        [
            step,
            spiralDuration / yearInSeconds,
            inclinationChangeDuration / yearInSeconds,
        ]
    )

AU = 149.6e9  # m
yearInSeconds = 365 * 24 * 3600

fig, ax = plt.subplots(3, 1)
fig3, ax3 = plt.subplots(1, 1, subplot_kw=dict(projection="3d"))
cols = ["time", "x", "y", "z", "vx", "vy", "vz"]
depVars = [
    "time",
    "ThrustX",
    "ThrustY",
    "ThrustZ",
    "ThrustMagnitude",
    "a",
    "e",
    "i",
    "omega",
    "RAAN",
    "theta",
    "cone",
    "clock",
]

spacecraftName, save, saveDep, extraTxt = saveFiles[0]
datasmallest = pd.read_csv(
    save,
    delimiter="	",
    names=cols,
    header=None,
    skiprows=lambda i: i % skipSizes[0],
    ).iloc[1:520]
for idx, u in enumerate(saveFiles):

    spacecraftName, save, saveDep, extraTxt = u
    data = pd.read_csv(
        save,
        delimiter="	",
        names=cols,
        header=None,
        skiprows=lambda i: i % skipSizes[idx],
    )
    ax[0].plot(data.time.iloc[1:520], (data.x - datasmallest.x).iloc[1:520].abs(), label=spacecraftName)
    ax[1].plot(data.time.iloc[1:520], (data.y - datasmallest.y).iloc[1:520].abs(), label=spacecraftName)
    ax[2].plot(data.time.iloc[1:520], (data.z - datasmallest.z).iloc[1:520].abs(), label=spacecraftName)

    ax3.plot(data.x, data.y, data.z)
    print(len(data))

print(skipSizes)

ax[0].grid()
ax[1].grid()
ax[2].grid()
ax[0].legend()
ax[1].legend()
ax[2].legend()
ax[0].set_yscale("log")
ax[1].set_yscale("log")
ax[2].set_yscale("log")
fig.set_tight_layout(True)

fig2, ax2 = plt.subplots(1,1)
print(np.array(resultslst)[:,0])
print(np.array(resultslst)[:,2])
ax2.bar(x=["dt="+str(x) for x in stepSizes], height=(np.array(resultslst)[:,2]-np.array(resultslst)[0,2])*365)
ax2.set_ylabel("Transfer time error (days)")
ax2.set_xlabel("Timestep (seconds)")
ax2.grid()
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
