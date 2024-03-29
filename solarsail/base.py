import numpy as np

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion

from typing import Callable, Tuple, Union


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
        characteristicAcceleration = None,
        verbose=True,
    ):
        """Base sail model object

        Args:
            bodies (environment.SystemOfBodies): tudat bodies object
            sailName (str, optional): name of the sailcraft for tudat. Defaults to "Sail".
            mass (float, optional): mass in kg. Defaults to 100.
            sailArea (float, optional): sail area in m^2. Defaults to 22500.
            targetAltitude (float, optional): target altitude above the sun in AU. Defaults to 0.48.
            deepestAltitude (float, optional): deepest altitude to spiral to. Defaults to 0.48.
            targetInclination (float, optional): target final inclination. Defaults to 90.
            characteristicAcceleration (_type_, optional): optional set characteristic acceleration isntead of area. Defaults to None.
            verbose (bool, optional): print outputs, disable to improve performance. Defaults to True.
        """
        self.bodies = bodies
        self.sailName = sailName
        self.targetAltitude = targetAltitude
        self.deepestAltitude = deepestAltitude
        self.verbose = verbose

        self.mass = mass
        # self.sigmaC = 2 * (
        #     SolarSailGuidanceBase.I_AU
        #     / (SolarSailGuidanceBase.c * SolarSailGuidanceBase.gAU)
        # )

        reflectivity = 0.9 # TODO set a realistic value for this
        if characteristicAcceleration is None:
            self.sigma = mass / sailArea
            self.charAccel = ((9.08 * reflectivity) / (self.sigma * 1000)) / 1000
        else:
            if verbose:
                print(f"Characteristic acceleration overide: {characteristicAcceleration}")
            self.charAccel = characteristicAcceleration
            self.sigma = ((9.08 * reflectivity) / characteristicAcceleration) / 1e6
       
        if verbose:
            print(f"Characteristic acceleration = {round(self.charAccel*1000, 4)}")
            print(f"Mass = {self.mass} | sigma = {round(self.sigma,4)}")
            #    | sigmaC = {round(self.sigmaC, 4)}")

        self.targetInclination = np.radians(targetInclination)

        # 0 = Transfer, 1 = Pause, 2 = Inclination change
        self.currentPhase:int = 0
        self.alpha = 0
        self.delta = 0
        self.extraDependentVariables = 2
        
        self.lastTimeMeasured = 0
        self.lastAccelVector: np.ndarray

        self.inclinationChangeStart = 0
        self.inclinationChangeEnd = 0
        self.lastInclination = 0

        self.thermalAvailable = False

    def dependantVariables(self) -> np.ndarray:
        return np.array([self.alpha, self.delta])

    def norm(self, vec: np.ndarray) -> float:
        return np.sqrt(np.square(vec).sum())

    def computeSail(self, current_time) -> np.ndarray:
        return np.zeros([3, 1])

    def stopPropagation(self, time) -> bool:
        return False

    def compute_thrust_direction(self, current_time: float) -> np.ndarray:
        """Uses the computesail function to calculate the thrust direction

        Args:
            current_time (float): current time from env

        Returns:
            np.ndarray: normalised thrust direction
        """
        # Check if computation is to be done
        if current_time == current_time:
            # This if/else is to avoid duplicate compute sail calls
            if self.lastTimeMeasured == current_time:
                acc = self.lastAccelVector
                return acc / np.sqrt(np.square(acc).sum())
            else:
                acc = self.computeSail(current_time)

                self.lastAccelVector = acc
                self.lastTimeMeasured = current_time

                accNorm = self.norm(acc)
                if accNorm == 0:
                    direction = np.zeros([3, 1])
                else:
                    direction = acc / accNorm

                return direction

        # If no computation is to be done, return zeros
        else:
            return np.zeros([3, 1])

    def compute_thrust_magnitude(self, current_time: float) -> float:
        """Uses the computesail function to calculate the thrust direction

        Args:
            current_time (float): current time from env

        Returns:
            float: thrust magnitude in newton
        """
        # Check if computation is to be done
        if current_time == current_time:
            # This if/else is to avoid duplicate compute sail calls
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
