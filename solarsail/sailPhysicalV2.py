from .base import SolarSailGuidanceBase


import numpy as np

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion
from scipy.optimize import fsolve


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
        deepestAltitude: float = 0.48,
        targetInclination: float = 90,
        # Deprecated: here for backwards compatibility
        maximum_thrust: float = 0,
        # ac=0.0005,
        reflectivity: float = 1,
    ):
        super().__init__(
            bodies,
            sailName,
            mass,
            sailArea,
            targetAltitude,
            deepestAltitude,
            targetInclination,
            maximum_thrust,
        )

        lambd = 168.6284 * self.charAccel
        self.spiralAlpha = fsolve(lambdFunc, 30 * np.pi / 180, args=(lambd))[0]
        
        print(f"Spiral Alpha = {self.spiralAlpha}")

        
    def stopPropagation(self):
        if self.currentPhase == 3:
            return True
        else:
            return False

    def computeSail(self, current_time) -> np.ndarray:
        current_cartesian_state = (
            self.bodies.get(self.sailName).state - self.bodies.get("Sun").state
        )
        # current_cartesian_velocity = (
        #     self.bodies.get(self.sailName).velocity
        #     # - self.bodies.get("Sun").velocity
        # )

        current_cartesian_position = (
            self.bodies.get(self.sailName).position - self.bodies.get("Sun").position
        )

        # progradeDirection = current_cartesian_velocity / np.sqrt(
        #     np.square(current_cartesian_velocity).sum()
        # )

        current_alt = np.sqrt(np.square(current_cartesian_position).sum())

        # radialDirection = current_cartesian_position / np.sqrt(
        #     np.square(current_cartesian_position).sum()
        # )

        mu = self.bodies.get("Sun").gravitational_parameter

        current_keplerian_state = element_conversion.cartesian_to_keplerian(
            current_cartesian_state, mu
        )
        # H = np.cross(current_cartesian_velocity, current_cartesian_position)
        # Hdirection = H / self.norm(H)

        inclination = current_keplerian_state[2]

        trueAnomaly = current_keplerian_state[5]
        argPeriapsis = current_keplerian_state[3]
        RAAN = current_keplerian_state[4]

        ###########################################
        ############## Flight Phases ##############
        ###########################################
        # Transfer
        if self.currentPhase == 0:
            if current_alt / SolarSailGuidance.AU < self.targetAltitude:
                print("Altitude reached -> Inclination Change")
                self.inclinationChangeStart = current_time
                self.currentPhase = 2

            delta = 1.5 * np.pi
            alpha = self.spiralAlpha

        # wait for correct orbit
        elif self.currentPhase == 1:
            # TODO
            alpha = 0
            delta = 0

        # inclination change
        elif self.currentPhase == 2:
            self.inclinationChangeEnd = current_time
            self.lastInclination = np.degrees(inclination)

            if inclination > self.targetInclination:
                print("Inclination Change complete -> Science")
                self.currentPhase = 3

            alpha = np.arctan(1 / np.sqrt(2))

            if np.cos(argPeriapsis + trueAnomaly) >= 0:
                delta = 0
            else:
                delta = np.pi

        # post-inclination change
        elif self.currentPhase == 3:
            # TODO
            alpha = 0
            delta = 0
        else:
            alpha = 0
            delta = 0

        # B = (mu / (current_alt * current_alt)) * (0.5 * (self.sigmaC / self.sigma))

        nVecAlt = np.array(
            [
                np.cos(alpha),
                np.sin(alpha) * np.sin(delta),
                np.sin(alpha) * np.cos(delta),
            ]
        )

        F0 = (SolarSailGuidanceBase.AU / current_alt) ** 2 * self.charAccel

        ForceRTN = F0 * (np.cos(alpha) ** 2) * nVecAlt

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
