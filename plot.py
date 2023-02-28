import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

AU = 149.6e9

fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection="3d"))

file = "SOLARSAILPropagationHistory_Q1.dat"
file2 = "SOLARSAILPropagationHistory_DependentVariables_Q1.dat"

data = pd.read_csv(file, delimiter="	", names=["t", "x", "y", "z", "vx", "vy", "vz"], header=None).assign(
    altitude=lambda x: np.sqrt(x.x**2 + x.y**2 + x.z**2) / AU
)
data2 = pd.read_csv(file2, delimiter="	")


ax.plot(data.iloc[:, 1], data.iloc[:, 2], data.iloc[:, 3])
ax.set_ylim(-AU, AU)
ax.set_xlim(-AU, AU)
# ax.set_zlim(-AU, AU)
ax.set_aspect("equal")
# fig.set_tight_layout()

fig2, ax2 = plt.subplots(3,1)
# ax2[0].plot(data2.iloc[:, 1:5], label=["x", "y", "z", "norm"])
# ax2[0].legend()
# ax2[0].set_title("Sail Acceleration")
ax2[1].set_title("Altitude (AU)")
ax2[1].plot(data.altitude)
# ax2[2].set_title("Heading (deg)")
# ax2[2].plot(np.degrees(data2.iloc[:, 5]))
# fig2.set_tight_layout(True)

fig.set_tight_layout(True)

plt.show()
