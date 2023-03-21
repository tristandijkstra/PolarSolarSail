from .base import SolarSailGuidanceBase
from .transform import simplifiedSailToInertial


import numpy as np
import math

from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.astro import element_conversion
from scipy.optimize import fsolve

from typing import Union, Tuple


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
        deepestAltitude: float = 0.2,
        targetInclination: float = 90,
        fastTransferOptimiseParameter: float = 0.0,
        characteristicAcceleration: Union[None, float] = None,
        verbose=True,
    ):
        """V3 of the sail guidance model:
            - Semi-Physical sail model
            - Fast Transfers
            - Compatible with optimisers
            - Outputs alpha and delta
            - Faster computation time

        Args:
            bodies (environment.SystemOfBodies): tudat bodies object
            sailName (str, optional): name of the sailcraft for tudat. Defaults to "Sail".
            mass (float, optional): mass in kg. Defaults to 100.
            sailArea (float, optional): sail area in m^2. Defaults to 22500.
            targetAltitude (float, optional): target final altitude above the sun in AU. Defaults to 0.48.
            deepestAltitude (float, optional): deepest altitude to spiral to. Defaults to 0.48.
            targetInclination (float, optional): target final inclination. Defaults to 90.
            fastTransferOptimiseParameter (float, optional): mixing value. Defaults to 0.0.
            characteristicAcceleration (_type_, optional): optional set characteristic acceleration isntead of area. Defaults to None.
            verbose (bool, optional): print outputs, disable to improve performance. Defaults to True.
        """

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
        self.timesOutwardsCompleted = 0
        self.goingOutwards = False

        
    def stopPropagation(self, time:float) -> bool:
        """Stop Tudat propagation

        Args:
            time (float): environment time

        Returns:
            bool: True if stopping, false otherwise
        """
        if self.currentPhase == 10:
            if self.verbose:
                print(f"Stopping Propagation.")
                print(f"Times outwards: {self.timesOutwardsCompleted}")
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
    

    def getOptimiseOutput(self) -> Tuple[float, float, float]:
        """Returns a set of parameters to the the optimiser

        Returns:
            Tuple[float, float, float]: inclination change duration (yrs), final inclination (deg), times going outwards
        """
        inclinationChangeDuration = (
            self.inclinationChangeEnd - self.inclinationChangeStart
        )
        return inclinationChangeDuration, self.lastInclination, self.timesOutwardsCompleted


    def computeSail(self, current_time:float) -> np.ndarray:
        """The core of the sail object, gives guidance based on environment and calculates thrust vector.
            The flight is split into multiple phases:
                1. TODO
            The result is served to the thrust direction and thrust magnitude function.

        Args:
            current_time (float): current time from integrator

        Returns:
            np.ndarray: Thrust vector
        """
        current_cartesian_state = (
            self.bodies.get(self.sailName).state
        )
        current_cartesian_position = (
            self.bodies.get(self.sailName).position
        )

        current_alt = np.sqrt(np.square(current_cartesian_position).sum())

        mu = self.bodies.get("Sun").gravitational_parameter

        current_keplerian_state = element_conversion.cartesian_to_keplerian(
            current_cartesian_state, mu
        )

        inclination = current_keplerian_state[2]
        # semimajoraxis = current_keplerian_state[0]

        trueAnomaly = current_keplerian_state[5]
        argPeriapsis = current_keplerian_state[3]
        RAAN = current_keplerian_state[4]

        current_altAU = current_alt / SolarSailGuidance.AU

        ###########################################
        ############## Flight Phases ##############
        ###########################################

        # This is just to record the start time
        if self.currentPhase == -1:
            self.startTime = current_time
            self.currentPhase = 0

            if self.verbose:
                print("Spiraling to deepest altitude")


        # Phase 0: Spiraling inwards
        elif self.currentPhase == 0:
            # if we're just past our target altitude start spiraling
            if current_altAU < self.deepestAltitude:
                self.inclinationChangeStart = current_time
                self.currentPhase = 2
                self.goingOutwards = True
            
                if self.verbose:
                    print("Deepest altitude reached -> Inclination Change")

            self.delta = 1.5 * np.pi
            self.alpha = self.spiralAlpha

        # # wait for correct orbit
        # elif self.currentPhase == 1:
        #     # TODO
        #     self.alpha = 0
        #     self.delta = 0


        # Phase 2: Inclination change + spiral back out
        elif self.currentPhase == 2:
            self.inclinationChangeEnd = current_time
            self.lastInclination = np.degrees(inclination)

            # if we reach target inclination go to end
            if inclination > self.targetInclination:
                if current_altAU > self.targetAltitude:
                    self.currentPhase = 10
                else:
                    self.currentPhase = 9

                if self.verbose:
                    print("Inclination Change complete -> Spiraling to target")

            # If we are past outmost altitude, spiral back in
            if current_altAU > self.targetAltitude:
                self.currentPhase = 3
                if self.goingOutwards:
                    self.timesOutwardsCompleted += 1
                self.goingOutwards = False

                if self.verbose:
                    print("Inclination change + Spiraling In")

            self.alpha = self.inclinationChangeAlpha

            # Switching cone angle based depending on ascending/descending
            if np.cos(argPeriapsis + trueAnomaly) >= 0:
                self.delta = np.pi * (self.fastTransferOptimiseParameter)
            else:
                self.delta = np.pi * (1 - self.fastTransferOptimiseParameter)


        # Phase 3: Inclination change + spiral back in
        elif self.currentPhase == 3:
            self.inclinationChangeEnd = current_time
            self.lastInclination = np.degrees(inclination)

            # if we reach target inclination go to end
            if inclination > self.targetInclination:
                self.currentPhase = 9

                if self.verbose:
                    print("Inclination Change complete -> Spiraling to target")

            # If we are past innermost altitude, spiral back out
            if current_altAU < self.deepestAltitude:
                self.currentPhase = 2
                self.goingOutwards = True

                if self.verbose:
                    print("Inclination change + Spiraling out")

            self.alpha = self.inclinationChangeAlpha

            # Switching cone angle based depending on ascending/descending
            if np.cos(argPeriapsis + trueAnomaly) >= 0:
                self.delta = np.pi * (2 - self.fastTransferOptimiseParameter)
            else:
                self.delta = np.pi * (1 + self.fastTransferOptimiseParameter)


        # Phase 9: Spiral back out to target altitude.
        elif self.currentPhase == 9:
            if current_altAU > self.targetAltitude:
                self.currentPhase = 10
                # self.timesOutwardsCompleted = 100

                if self.verbose:
                    print("Starting Science Phase")

            self.delta = 0.5 * np.pi
            self.alpha = self.spiralAlpha


        # Phase 10: Ending, set everything to 0 thrust, 
        # this also calls the termination function
        elif self.currentPhase == 10:
            # TODO
            self.alpha = 0.5 * np.pi
            self.delta = 0


        # If the phase is unbounded, print a bunch of errors
        else:
            print(f"Error: current phase unbounded = {self.currentPhase}")
            self.alpha = 0.5 * np.pi
            self.delta = 0


        # Normal vector.
        # Simplified physics -> normal vec = thrust vec
        nVecAlt = np.array(
            [
                math.cos(self.alpha),
                math.sin(self.alpha) * math.sin(self.delta),
                math.sin(self.alpha) * math.cos(self.delta),
            ]
        )

        # Force Magnitude scaled on altitude
        # F0 = ((SolarSailGuidanceBase.AU / current_alt) ** 2) * self.charAccel
        F0 = (1/(current_altAU * current_altAU)) * self.charAccel

        # Force Vector in radial tangential normal reference frame
        forceRTN = F0 * (math.cos(self.alpha) * math.cos(self.alpha)) * nVecAlt

        # Transform to inertial reference frame
        forceINERTIAl = simplifiedSailToInertial(inclination, argPeriapsis, trueAnomaly, RAAN) @ forceRTN

        return forceINERTIAl
