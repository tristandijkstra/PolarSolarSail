import numpy as np

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion

from typing import Callable, Tuple


def norm(vec:np.ndarray):
    return np.sqrt(np.square(vec).sum())

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
) -> Tuple[float, float, float]:
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
        mass: float = 100,
        sailArea: float = 22500,
        targetAltitude: float = 0.48,
        # TODO: add fast transfer stuff
        deepestAltitude: float = 0.48,
        targetInclination: float = 90,
        # TODO: get realistic values for these:
        a: float = 1,
        rspec: float = 1,
        rdiff: float = 1,
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
            SolarSailGuidance.I_AU / (SolarSailGuidance.c * SolarSailGuidance.gAU)
        )

        print(f"sigma = {self.sigma} | sigmaC = {self.sigmaC}")

        self.a = a
        self.rspec = rspec
        self.rdiff = rdiff

        self.targetInclination = np.radians(targetInclination)

        # 0 = Transfer, 1 = Pause, 2 = Inclination change
        self.currentPhase = 0

        self.lastTimeMeasured = 0
        self.lastAccelVector: np.ndarray

    def computeSail(self) -> np.ndarray:
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
        Hnorm = np.cross(radialDirection, progradeDirection)

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
            # alpha = 0 * np.pi / 4
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

        B = (mu / (current_alt * current_alt)) * (0.5 * (self.sigmaC / self.sigma)) * nx

        kappa, Xf, Xb = kappaf(T)

        # uvec = np.array([1, 0, 0])

        # TODO nvec uvec
        bVec = (
            ((2 * self.rspec * nx) + (Xf * self.rdiff) + (kappa * self.a)) * nvec
        ) + ((self.a + self.rdiff) * radialDirection)

        # TODO Transform to inertial
        # rh = radialDirection * Hdirection
        n_temp = H / norm(H)

        ncrossr = np.cross(n_temp, current_cartesian_position)
        t_temp = ncrossr / norm(ncrossr)
        r_temp = current_cartesian_position / norm(current_cartesian_position)

        transformMatrix = np.array([n_temp, t_temp, r_temp])
        # transform = np.cross(rh, rh)

        accelerationSC = transformMatrix @ (B * bVec)
        # accelerationSC = (B * bVec)
        # print(accelerationSC)

        return accelerationSC

    def compute_thrust_direction(self, current_time: float) -> np.ndarray:
        if self.lastTimeMeasured == current_time:
            acc = self.lastAccelVector
            return acc / np.sqrt(np.square(acc).sum())
        elif current_time == current_time:
            acc = self.computeSail()

            self.lastAccelVector = acc
            self.lastTimeMeasured = current_time

            # print(f"thrust direction {acc}")
            direction = acc / np.sqrt(np.square(acc).sum())

            return direction
        # If no computation is to be done, return zeros
        else:
            return np.zeros([3, 1])

    def compute_thrust_magnitude(self, current_time: float) -> float:
        # Check if computation is to be done
        if self.lastTimeMeasured == current_time:
            acc = self.lastAccelVector
            return np.sqrt(np.square(acc).sum()) * self.mass
        elif current_time == current_time:
            # print("thrust magnitude")
            acc = self.computeSail()

            self.lastAccelVector = acc
            self.lastTimeMeasured = current_time

            magnitude = np.sqrt(np.square(acc).sum()) * self.mass
            return magnitude

        # If no computation is to be done, return zeros
        else:
            return 0.0
