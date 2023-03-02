import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm


AU = 149.6e9
muSun = 1.327124400e20


def Vc(radiusAU: float, mu: float = muSun) -> float:
    return np.sqrt(mu / (radiusAU * AU))


def dV(radiusAU: float, inclChangeDeg: float, mu: float = muSun):
    v1 = Vc(radiusAU, mu=mu)

    dv = 2 * v1 * np.sin(np.radians(inclChangeDeg) / 2)

    return dv


if __name__ == "__main__":
    defaultVal = (dV(0.48, 90))

    radii = np.linspace(0.25, 0.75, 100)
    inclinations = np.linspace(60, 90, 30)
    rs, incls = np.meshgrid(radii, inclinations)

    dvs = dV(rs, incls)

    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})

    surf = ax.plot_surface(
        rs, incls, dvs, cmap=cm.coolwarm, linewidth=0, antialiased=False
    )
    ax.set_ylabel(r"Inclination change [$\deg$]")
    ax.set_xlabel(r"Altitude $AU$")
    ax.set_zlabel(r"delta V $ms^{-1}$")
    fig.colorbar(surf, shrink=0.5, aspect=5)

    fig.suptitle(f"For 0.48 AU, 90 deg => dV = {round(defaultVal, 0)} m/s")
    plt.show()
