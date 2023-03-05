import numpy as np

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion

from typing import Callable, Tuple


class SolarSailGuidanceBase:
    AU = constants.ASTRONOMICAL_UNIT  # m
    I_AU = 1366  # W/m2 irradiance at 1 AU
    c = constants.SPEED_OF_LIGHT

    MUsun = 1.32712440042e20  #

    gAU = MUsun / AU**2

    def __init__(
        self,
        bodies: environment.SystemOfBodies,
        sailName: str = "Sail",
        mass: float = 100,
        sailArea: float = 22500,
        targetAltitude: float = 0.48,
        deepestAltitude: float = 0.48,
        targetInclination: float = 90,
        # Deprecated: here for backwards compatibility
        maximum_thrust: float = 0,
    ):
        self.bodies = bodies
        self.sailName = sailName
        self.targetAltitude = targetAltitude
        self.deepestAltitude = deepestAltitude

        self.mass = mass
        self.sigma = mass / sailArea
        self.sigmaC = 2 * (
            SolarSailGuidanceBase.I_AU
            / (SolarSailGuidanceBase.c * SolarSailGuidanceBase.gAU)
        )

        print(f"sigma = {self.sigma} | sigmaC = {self.sigmaC}")

        self.targetInclination = np.radians(targetInclination)

        # 0 = Transfer, 1 = Pause, 2 = Inclination change
        self.currentPhase = 0

        self.lastTimeMeasured = 0
        self.lastAccelVector: np.ndarray

        self.inclinationChangeStart = 0
        self.inclinationChangeEnd = 0
        self.lastInclination = 0

    def norm(self, vec: np.ndarray):
        return np.sqrt(np.square(vec).sum())

    def computeSail(self, current_time) -> np.ndarray:
        return np.zeros([3, 1])

    def compute_thrust_direction(self, current_time: float) -> np.ndarray:
        # Check if computation is to be done
        if current_time == current_time:
            if self.lastTimeMeasured == current_time:
                acc = self.lastAccelVector
                return acc / np.sqrt(np.square(acc).sum())
            else:
                acc = self.computeSail(current_time)

                self.lastAccelVector = acc
                self.lastTimeMeasured = current_time

                accNorm = self.norm(acc)
                if accNorm == 0:
                    return np.zeros([3, 1])
                else:
                    direction = acc / accNorm

                return direction

        # If no computation is to be done, return zeros
        else:
            return np.zeros([3, 1])

    def compute_thrust_magnitude(self, current_time: float) -> float:
        # Check if computation is to be done
        if current_time == current_time:
            if self.lastTimeMeasured == current_time:
                acc = self.lastAccelVector
                return np.sqrt(np.square(acc).sum()) * self.mass
            else:
                acc = self.computeSail(current_time)

                self.lastAccelVector = acc
                self.lastTimeMeasured = current_time

                magnitude = np.sqrt(np.square(acc).sum()) * self.mass
                return magnitude

        # If no computation is to be done, return zeros
        else:
            return 0.0

    def getInclinationChangeDuration(self):
        inclinationChangeDuration = (
            self.inclinationChangeEnd - self.inclinationChangeStart
        )
        return inclinationChangeDuration, self.lastInclination
