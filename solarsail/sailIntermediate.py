from .base import SolarSailGuidanceBase


import numpy as np

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion


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

    def computeSail(self, current_time) -> np.ndarray:
        T = 100 + 273

        current_cartesian_state = (
            self.bodies.get(self.sailName).state - self.bodies.get("Sun").state
        )
        current_cartesian_velocity = (
            self.bodies.get(self.sailName).velocity
            # - self.bodies.get("Sun").velocity
        )

        current_cartesian_position = (
            self.bodies.get(self.sailName).position - self.bodies.get("Sun").position
        )

        progradeDirection = current_cartesian_velocity / np.sqrt(
            np.square(current_cartesian_velocity).sum()
        )

        current_alt = np.sqrt(np.square(current_cartesian_position).sum())

        radialDirection = current_cartesian_position / np.sqrt(
            np.square(current_cartesian_position).sum()
        )

        mu = self.bodies.get("Sun").gravitational_parameter

        current_keplerian_state = element_conversion.cartesian_to_keplerian(
            current_cartesian_state, mu
        )
        H = np.cross(current_cartesian_velocity, current_cartesian_position)
        # Hnorm = np.cross(radialDirection, progradeDirection)
        Hdirection = H / self.norm(H)

        inclination = current_keplerian_state[2]

        ###########################################
        ############## Flight Phases ##############
        ###########################################
        # Transfer
        if self.currentPhase == 0:
            if current_alt / SolarSailGuidance.AU < self.targetAltitude:
                print("Altitude reached -> Inclination Change")
                self.inclinationChangeStart = current_time
                self.currentPhase = 2
            b = -progradeDirection

        # wait for correct orbit
        elif self.currentPhase == 1:
            b = np.zeros([3, 1])

        # inclination change
        elif self.currentPhase == 2:
            self.inclinationChangeEnd = current_time
            self.lastInclination = np.degrees(inclination)
            if inclination > self.targetInclination:
                print("Inclination Change complete -> Science")
                self.currentPhase = 3
            
            thetaNode = current_keplerian_state[3]
            # thetaNode = (2 * np.pi) - current_keplerian_state[3]
            trueAnomaly = current_keplerian_state[5]

            nearNode1 = abs(trueAnomaly - thetaNode) < 0.5 * np.pi
            nearNode2 = abs(trueAnomaly - (thetaNode + np.pi)) < 0.5 * np.pi
            
            if nearNode1:
                b = Hdirection
            elif nearNode2:
                b = -Hdirection
            else:
                print((trueAnomaly - thetaNode), trueAnomaly)
                b = Hdirection

        # post-inclination change
        elif self.currentPhase == 3:
            b = np.zeros([3, 1])

        else:
            b = np.zeros([3, 1])

        B = (mu / (current_alt * current_alt)) * (0.5 * (self.sigmaC / self.sigma))

        return B * b
