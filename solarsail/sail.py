import os

import numpy as np
from matplotlib import pyplot as plt

# from tudatpy.io import save2txt
from tudatpy.kernel import constants

# from tudatpy.kernel.interface import spice
# from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment

# from tudatpy.kernel.numerical_simulation import environment_setup
# from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.astro import element_conversion
from typing import Callable


def eps_f_func(T: float) -> float:
    return 0.73


def eps_b_func(T: float) -> float:
    return 0.3


def kappaf(
    T: float,
    Xf: float = 1,
    Xb: float = 1,
    epsF: Callable[[float], float] = eps_f_func,
    epsB: Callable[[float], float] = eps_b_func,
):
    kappa = (Xf * epsF(T) - Xb * epsB(T)) / (epsF(T) + epsB(T))
    return kappa, Xf, Xb


class SolarSailGuidance:
    AU = constants.ASTRONOMICAL_UNIT  # m
    I_AU = 1366  # W/m2 irradiance at 1 AU
    c = constants.SPEED_OF_LIGHT

    MUsun = 1.32712440042e20  #

    gAU = MUsun / AU**2

    def __init__(
        self,
        bodies: environment.SystemOfBodies,
        sailName: str = "Sail",
        mass=500,
        sailArea=22500,
        a: float = 1,
        rspec: float = 1,
        rdiff: float = 1,
        targetAltitude=0.48,
        targetInclination=90,
    ):
        self.bodies = bodies
        self.sailName = sailName
        self.targetAltitude = targetAltitude

        self.mass = mass
        self.sigma = mass / sailArea

        self.a = a
        self.rspec = rspec
        self.rdiff = rdiff

        self.targetInclination = np.radians(targetInclination)

        self.sigmaC = 2 * (
            SolarSailGuidance.I_AU / (SolarSailGuidance.c * SolarSailGuidance.gAU)
        )

        # 0 = Transfer, 1 = Pause, 2 = Inclination change
        self.currentPhase = 0
        self.lastTimeMeasured = 0
        self.lastAccelVector: np.ndarray

    def computeSail(self):
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

        Hdirection = np.cross(radialDirection, progradeDirection)

        inclination = current_keplerian_state[2]

        ###########################################
        ############## Flight Phases ##############
        ###########################################
        # Transfer
        if self.currentPhase == 0:
            if current_alt / SolarSailGuidance.AU < self.targetAltitude:
                print("Altitude reached -> Inclination Change")
                self.currentPhase = 2
            # return -progradeDirection
            # alpha = np.pi/2
            delta = 0
            alpha = np.arctan(1 / np.sqrt(2))
            # delta = np.pi/2
        # wait for correct orbit
        elif self.currentPhase == 1:
            # return np.zeros([3, 1])
            alpha = 0.1
            delta = 0.1
        # inclination change
        elif self.currentPhase == 2:
            if inclination > self.targetInclination:
                print("Inclination Change complete -> Science")
                self.currentPhase = 3
            # return Hdirection
            alpha = 0.1
            delta = 0.1
        # post-inclination change
        elif self.currentPhase == 3:
            # return np.zeros([3, 1])
            alpha = 0.1
            delta = 0.1
        else:
            alpha = 0.1
            delta = 0.1

        ##########################################
        ############### SAIL STUFF ###############
        ##########################################

        # TODO make this nice
        ##########################################

        nvec = np.array(
            [
                (np.cos(alpha) * np.cos(delta)),
                (np.sin(alpha) * np.cos(delta)),
                np.sin(delta),
            ]
        )
        ##########################################

        nx = nvec[0]

        B = (
            (self.bodies.get("Sun").gravitational_parameter / current_alt**2)
            * (0.5 * (self.sigmaC / self.sigma))
            * nx
        )

        kappa, Xf, Xb = kappaf(T)

        # TODO nvec uvec
        bVec = (((2 * self.rspec * nx) + (Xf * self.rdiff) + kappa * self.a) * nvec) + (
            (self.a + self.rdiff) * radialDirection
        )

        # Transform to inertial
        # TODO

        accelerationSC = B * bVec

        # self.lastAccelVector = accelerationSC

        return accelerationSC

    def compute_thrust_direction(self, current_time: float) -> np.ndarray:
        # Check if computation is to be done
        if self.lastTimeMeasured == current_time:
            acc = self.lastAccelVector
            return acc / np.sqrt(np.square(acc).sum())
        elif current_time == current_time:
            acc = self.computeSail()

            self.lastAccelVector = acc
            self.lastTimeMeasured = current_time

            direction = acc / np.sqrt(np.square(acc).sum())

            return direction
        # If no computation is to be done, return zeros
        else:
            return np.zeros([3, 1])

    def compute_thrust_magnitude(self, current_time: float):
        # Check if computation is to be done
        if self.lastTimeMeasured == current_time:
            acc = self.lastAccelVector
            return np.sqrt(np.square(acc).sum())
        elif current_time == current_time:
            acc = self.computeSail()

            self.lastAccelVector = acc
            self.lastTimeMeasured = current_time

            magnitude = np.sqrt(np.square(acc).sum())
            return magnitude

        # If no computation is to be done, return zeros
        else:
            return 0.0
