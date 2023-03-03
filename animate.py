import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import os.path

AU = 149.6e9


def animate(speed, dataFile, bodyColor="cyan", plotLimits=AU, savename = "sailAnimation"):
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    # if Body is not None:
    #     assert (type(Body) == msp.Planet), "Incorrect Body type, should be main.Planet"
    #     # Sphere:
    #     u = np.linspace(0, 2 * np.pi, 100)
    #     v = np.linspace(0, np.pi, 100)
    #     x = Body.r * np.outer(np.cos(u), np.sin(v))
    #     y = Body.r * np.outer(np.sin(u), np.sin(v))
    #     z = Body.r * np.outer(np.ones(np.size(u)), np.cos(v))
    #     # Plot the surface
    #     ax.plot_surface(x, y, z, color=str(bodyColor))

    (line,) = ax.plot([], [], lw=1)
    (scatt,) = ax.plot([], [], linestyle="", marker="o", color="white")

    ax.set_ylim(-plotLimits, plotLimits)
    ax.set_xlim(-plotLimits, plotLimits)
    ax.set_zlim(-plotLimits, plotLimits)
    plt.axis("off")

    # SOLARSAILPropagationHistoryAnimate
    data = (
        pd.read_csv(
            dataFile,
            delimiter="	",
            names=["t", "x", "y", "z", "vx", "vy", "vz"],
            header=None,
        )
        .assign(altitude=lambda x: np.sqrt(x.x**2 + x.y**2 + x.z**2) / AU)
        .loc[:, ["x", "y", "z"]]
    )

    print(f"Original amount of data points: {len(data)}", end="")

    droppedPoints = []
    for u in range(len(data)):
        if u % speed != 0:
            droppedPoints.append(u)

    data = data.drop(droppedPoints, axis=0).to_numpy()
    print(f", reduced to: {len(data)}")

    # savename = "sailAnimation"
    headerText = savename + " " + str(speed*100) + "x speed"
    title = ax.set_title(headerText)
    steps = len(data)

    def animate(i):
        # if manoeuvreFileAvailable:
        #     # MANOEUVRE
        #     manoeuverPosList = []
        #     manoeuverTimeList = []
        #     for d, manoeuver in manoeuvreData.iterrows():
        #         if int((manoeuver["clock"] - clock.to_numpy()[0][0])/speed) < i:
        #             manoeuverPos = manoeuver["r"][1:-1].split(" ")
        #             manoeuverPos = [float(x) for x in manoeuverPos if x]
        #             if manoeuver["clock"] not in manoeuverTimeList:
        #                 manoeuverPosList.append(manoeuverPos)
        #                 manoeuverTimeList.append(manoeuver["clock"])

        #     if len(manoeuverPosList) > 0:
        #         mandisplaydata = np.array(manoeuverPosList).T
        #         scatt.set_data(mandisplaydata[0], mandisplaydata[1])
        #         scatt.set_3d_properties(mandisplaydata[2])

        #     if i == 0:
        #         scatt.set_data(100000, 100000)
        #         scatt.set_3d_properties(100000)

        # TRAJECTORY
        displaydata = np.array(data[0 : i + 1]).T
        line.set_data(displaydata[0], displaydata[1])
        line.set_3d_properties(displaydata[2])

        # ax.view_init(elev=20., azim=60)
        ax.view_init(elev=20., azim=i/4)
        # title.set_text(headerText + ' radius = {}'.format(round(np.sqrt(displaydata[0].dot(displaydata[0])),0)))
        return (
            title,
            line,
            scatt,
        )

    anim = FuncAnimation(fig, animate, frames=steps, interval=1, blit=True)
    if not os.path.exists("animations"):
        os.makedirs("animations")
    fullSaveName = "animations/" + savename + ".mp4"
    print("Saving to " + fullSaveName)
    anim.save(fullSaveName, fps=120, extra_args=["-vcodec", "libx264"], dpi=100)

    # plt.show()


if __name__ == "__main__":
    SPEED = 1000
    DATAFILE = r"SOLARSAILPropagationHistoryAnimate.dat"
    savename = "sailAnimationRot"

    animate(SPEED, DATAFILE, savename=savename)
