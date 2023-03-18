from .base import SolarSailGuidanceBase


import numpy as np

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion
from scipy.optimize import fsolve

from typing import Union


def lambdFunc(alpha, lamd=0):
    return (2 - 4 * np.tan(alpha) ** 2) / (
        np.cos(alpha) * (2 - np.tan(alpha) ** 2)
    ) - lamd


class SolarSailGuidance(SolarSailGuidanceBase):
    def __init__(
        self,
        bodies: environment.SystemOfBodies,
        sailName: str = "Sail",
        mass: float = 100,
        sailArea: float = 22500,
        targetAltitude: float = 0.48,
        # TODO: add fast transfer stuff
        deepestAltitude: float = 0.2,
        targetInclination: float = 90,
        fastTransferOptimiseParameter: float = 0.0,
        characteristicAcceleration: Union[None, float] = None,
        verbose=True,
    ):
        super().__init__(
            bodies,
            sailName,
            mass,
            sailArea,
            targetAltitude,
            deepestAltitude,
            targetInclination,
            characteristicAcceleration,
            verbose
        )

        lambd = 168.6284 * self.charAccel
        self.spiralAlpha = fsolve(lambdFunc, 30 * np.pi / 180, args=(lambd))[0]
        self.inclinationChangeAlpha = np.arctan(1 / np.sqrt(2))
        self.fastTransferOptimiseParameter = fastTransferOptimiseParameter
        
        if self.verbose:
            print("=== Cone angles ===")
            print(f"Spiral Alpha = {round(self.spiralAlpha, 4)} rad = {round(np.degrees(self.spiralAlpha), 4)} deg")
            print(f"Inclination Change Alpha = {round(self.inclinationChangeAlpha, 4)} rad = {round(np.degrees(self.inclinationChangeAlpha), 4)} deg")
            print(f"Fast transfer optimisation parameter: {self.fastTransferOptimiseParameter}")
        
        self.currentPhase:int = -1
        self.startTime = 0
        self.timesOutwards = 0

        
    def stopPropagation(self, time):
        if self.currentPhase == 10:
            if self.verbose:
                print(f"Stopping Propagation.")
                print(f"Times outwards: {self.timesOutwards}")
            return True
        else:
            return False
        
    def getInclinationChangeDuration(self):
        inclinationChangeDuration = (
            self.inclinationChangeEnd - self.inclinationChangeStart
        )
        spiralDuration = (
            self.inclinationChangeStart - self.startTime
        )
        return self.charAccel, spiralDuration, inclinationChangeDuration, self.lastInclination
    

    def getOptimiseOutput(self):
        inclinationChangeDuration = (
            self.inclinationChangeEnd - self.inclinationChangeStart
        )
        return inclinationChangeDuration, self.lastInclination, self.timesOutwards


    def computeSail(self, current_time) -> np.ndarray:
        current_cartesian_state = (
            self.bodies.get(self.sailName).state - self.bodies.get("Sun").state
        )
        current_cartesian_position = (
            self.bodies.get(self.sailName).position - self.bodies.get("Sun").position
        )

        current_alt = np.sqrt(np.square(current_cartesian_position).sum())

        mu = self.bodies.get("Sun").gravitational_parameter

        current_keplerian_state = element_conversion.cartesian_to_keplerian(
            current_cartesian_state, mu
        )

        inclination = current_keplerian_state[2]
        semimajoraxis = current_keplerian_state[0]

        trueAnomaly = current_keplerian_state[5]
        argPeriapsis = current_keplerian_state[3]
        RAAN = current_keplerian_state[4]

        current_altAU = current_alt / SolarSailGuidance.AU

        ###########################################
        ############## Flight Phases ##############
        ###########################################
        if self.currentPhase == -1:
            self.startTime = current_time
            self.currentPhase = 0
            if self.verbose:
                print("Spiraling")

        # Transfer
        elif self.currentPhase == 0:
            if current_altAU < self.deepestAltitude:
                if self.verbose:
                    print("Altitude reached -> Inclination Change")
                self.inclinationChangeStart = current_time
                self.currentPhase = 2
                self.timesOutwards += 1

            self.delta = 1.5 * np.pi
            self.alpha = self.spiralAlpha

        # # wait for correct orbit
        # elif self.currentPhase == 1:
        #     # TODO
        #     self.alpha = 0
        #     self.delta = 0

        # inclination change + spiral back out
        elif self.currentPhase == 2:
            self.inclinationChangeEnd = current_time
            self.lastInclination = np.degrees(inclination)

            if inclination > self.targetInclination:
                if self.verbose:
                    print("Inclination Change complete -> Science")
                self.currentPhase = 9

            if current_altAU > self.targetAltitude:
                if self.verbose:
                    print("Spiraling In")
                self.currentPhase = 3
                self.timesOutwards += 1

            self.alpha = self.inclinationChangeAlpha

            if np.cos(argPeriapsis + trueAnomaly) >= 0:
                self.delta = np.pi * (self.fastTransferOptimiseParameter)
            else:
                self.delta = np.pi * (1 - self.fastTransferOptimiseParameter)

        # spiral back in
        elif self.currentPhase == 3:
            self.inclinationChangeEnd = current_time
            self.lastInclination = np.degrees(inclination)

            if inclination > self.targetInclination:
                if self.verbose:
                    print("Inclination Change complete -> Science")
                self.currentPhase = 9

            if current_altAU < self.deepestAltitude:
                if self.verbose:
                    print("Spiraling out")
                self.currentPhase = 2

            self.alpha = self.inclinationChangeAlpha

            if np.cos(argPeriapsis + trueAnomaly) >= 0:
                self.delta = np.pi * (2 - self.fastTransferOptimiseParameter)
            else:
                self.delta = np.pi * (1 + self.fastTransferOptimiseParameter)

        # Transfer
        elif self.currentPhase == 9:
            if current_altAU > self.targetAltitude:
                # print("Inclination reached -> spiraling back out")
                # self.inclinationChangeStart = current_time
                if self.verbose:
                    print("finished spiraling")
                self.currentPhase = 10

            self.delta = 0.5 * np.pi
            self.alpha = self.spiralAlpha


        # post-inclination change
        elif self.currentPhase == 10:
            # TODO
            self.alpha = 0
            self.delta = 0
        else:
            print(f"Error: current phase unbounded = {self.currentPhase}")
            self.alpha = 0
            self.delta = 0


        nVecAlt = np.array(
            [
                np.cos(self.alpha),
                np.sin(self.alpha) * np.sin(self.delta),
                np.sin(self.alpha) * np.cos(self.delta),
            ]
        )

        F0 = (SolarSailGuidanceBase.AU / current_alt) ** 2 * self.charAccel

        ForceRTN = F0 * (np.cos(self.alpha) ** 2) * nVecAlt

        A1 = np.array(
            [
                [
                    np.cos(argPeriapsis + trueAnomaly),
                    np.sin(argPeriapsis + trueAnomaly),
                    0,
                ],
                [
                    -np.sin(argPeriapsis + trueAnomaly),
                    np.cos(argPeriapsis + trueAnomaly),
                    0,
                ],
                [
                    0,
                    0,
                    1
                ],
            ]
        )

        A2 = np.array(
            [
                [
                    1,
                    0,
                    0,
                ],
                [
                    0,
                    np.cos(inclination),
                    np.sin(inclination),
                ],
                [
                    0,
                    -np.sin(inclination),
                    np.cos(inclination),
                ],
            ]
        )
        A3 = np.array(
            [
                [
                    np.cos(RAAN),
                    np.sin(RAAN),
                    0,
                ],
                [
                    -np.sin(RAAN),
                    np.cos(RAAN),
                    0,
                ],
                [
                    0,
                    0,
                    1,
                ],
            ]
        )

        AtotalInv = np.linalg.inv(A1 @ A2 @ A3)

        b = AtotalInv @ ForceRTN

        return b
