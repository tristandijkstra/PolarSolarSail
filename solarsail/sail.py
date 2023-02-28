import os

import numpy as np
from matplotlib import pyplot as plt

# from tudatpy.io import save2txt
# from tudatpy.kernel import constants
# from tudatpy.kernel.interface import spice
# from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment

# from tudatpy.kernel.numerical_simulation import environment_setup
# from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.astro import element_conversion


class SolarSailGuidance:
    AU = 149.6e9

    def __init__(
        self,
        maximum_thrust: float,
        bodies: environment.SystemOfBodies,
        sailName: str = "Sail",
        targetAltitude=0.48,
    ):
        self.maximum_thrust = maximum_thrust
        self.bodies = bodies
        self.sailName = sailName
        self.targetAltitude = targetAltitude

        # 0 = Transfer, 1 = Pause, 2 = Inclination change
        self.currentPhase = 0

    def compute_thrust_direction(self, current_time: float) -> np.ndarray:
        # Check if computation is to be done
        if current_time == current_time:
            # Transfer
            if self.currentPhase == 0:
                # Retrieve current JUICE Cartesian state w.r.t. Ganymede from environment
                current_cartesian_state = (
                    self.bodies.get(self.sailName).state - self.bodies.get("Sun").state
                )
                current_cartesian_velocity = (
                    self.bodies.get(self.sailName).velocity
                    # - self.bodies.get("Sun").velocity
                )

                progradeDirection = current_cartesian_velocity / np.sqrt(
                    np.square(current_cartesian_velocity).sum()
                )
                current_alt = self.bodies.get(self.sailName).flight_conditions.altitude

                if current_alt / SolarSailGuidance.AU < self.targetAltitude:
                    print("Altitude reached -> Inclination Change")
                    self.currentPhase = 2

                # Compute and return current thrust direction (3x1 vector)
                # NOTE THIS CAN BE BETTER:
                return -progradeDirection
                # return np.array([0, 0, 0])

            elif self.currentPhase == 1:
                return np.zeros([3, 1])

            # inclination change
            elif self.currentPhase == 2:
                current_cartesian_state = (
                    self.bodies.get(self.sailName).state
                    - self.bodies.get("Sun").state
                )
                current_cartesian_position = (
                    self.bodies.get(self.sailName).position
                    - self.bodies.get("Sun").position
                )
                mu = self.bodies.get("Sun").gravitational_parameter

                current_keplerian_state = element_conversion.cartesian_to_keplerian(
                    current_cartesian_state, mu
                )
                current_cartesian_velocity = (
                    self.bodies.get(self.sailName).velocity
                    # - self.bodies.get("Sun").velocity
                )

                progradeDirection = current_cartesian_velocity / np.sqrt(
                    np.square(current_cartesian_velocity).sum()
                )
                radialDirection = current_cartesian_position / np.sqrt(
                    np.square(current_cartesian_position).sum()
                )

                Hdirection = np.cross(radialDirection, progradeDirection)

                inclination = current_keplerian_state[2]

                if inclination > np.pi / 2:
                    print("Inclination Change complete -> Science")
                    self.currentPhase = 3
                return Hdirection

            # post-inclination change
            elif self.currentPhase == 3:
                return np.zeros([3, 1])

        # If no computation is to be done, return zeros
        else:
            return np.zeros([3, 1])

    def compute_thrust_magnitude(self, current_time: float):
        # Check if computation is to be done
        if current_time == current_time:
            # # Retrieve current JUICE Cartesian  and Keplerian state w.r.t. Ganymede from environment
            # current_cartesian_state = (
            #     self.bodies.get("JUICE").state
            #     - self.bodies.get("Ganymede").state
            # )
            # gravitational_parameter = self.bodies.get(
            #     "Ganymede"
            # ).gravitational_parameter
            # current_keplerian_state = element_conversion.cartesian_to_keplerian(
            #     current_cartesian_state, gravitational_parameter
            # )

            # Compute and return current thrust magnitude (scalar)

            if (self.currentPhase == 0):
                thrust_magnitude = self.maximum_thrust/2
            elif (self.currentPhase == 2):
                thrust_magnitude = self.maximum_thrust*4
            else:
                thrust_magnitude = 0

            return thrust_magnitude
        # If no computation is to be done, return zeros
        else:
            return 0.0
