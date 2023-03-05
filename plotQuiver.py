import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

AU = 149.6e9
yearInSeconds = 365 * 24 * 3600

fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection="3d"))

file = "SOLARSAILPropagationHistory_Q1.dat"
file2 = "SOLARSAILPropagationHistory_DependentVariables_Q1.dat"
file = "SOLARSAIL_Basic.dat"
file2 = "SOLARSAIL_Dep_Basic.dat"
file = "data/mass500_area45000.dat"
file2 = "data/mass500_area45000_dep.dat"
file = "data/mass500_area5625.0.dat"
file2 = "data/mass500_area5625.0_dep.dat"

data = pd.read_csv(
    file, delimiter="	", names=["time", "x", "y", "z", "vx", "vy", "vz"], header=None
).assign(altitude=lambda x: np.sqrt(x.x**2 + x.y**2 + x.z**2) / AU)

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
]
data2 = pd.read_csv(file2, delimiter="	", names=depVars, header=None)

time = (data2.time - data.time.iloc[0])/yearInSeconds
ax.plot(data.iloc[:, 1], data.iloc[:, 2], data.iloc[:, 3])

quiverEvery = 6000
ax.quiver(
    data.x[::quiverEvery],
    data.y[::quiverEvery],
    data.z[::quiverEvery],
    data2.ThrustX[::quiverEvery],
    data2.ThrustY[::quiverEvery],
    data2.ThrustZ[::quiverEvery],
    length=0.1 * AU,
    normalize=True,
    color="tab:orange",
)
ax.set_ylim(-0.5*AU, 0.5*AU)
ax.set_xlim(-0.5*AU, 0.5*AU)
# ax.set_zlim(-AU, AU)
ax.set_aspect("equal")
# fig.set_tight_layout()

fig2, ax2 = plt.subplots(3, 1, sharex=True)
ax2[0].plot(time, data2.iloc[:, 1:4], label=["x", "y", "z"])
ax2[0].plot(time, data2.iloc[:, 4], label="norm", linestyle="--")
ax2[0].legend()
ax2[0].set_title("Sail Acceleration")

ax2[1].set_title("Altitude (AU)")
ax2[1].plot(time, data.altitude)
ax2[2].set_title("Inclination (deg)")
ax2[2].plot(time, np.degrees(data2.i))
# fig2.set_tight_layout(True)
ax2[2].set_xlabel("Time [years]")
fig.set_tight_layout(True)
fig2.set_tight_layout(True)

plt.show()
