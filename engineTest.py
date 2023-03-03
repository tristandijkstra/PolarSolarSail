###########################################################################
#
# # Numerical Astrodynamics 2022/2023
#
# # Assignment 1, Question 5 - Propagation Settings
#
###########################################################################


''' 
Copyright (c) 2010-2020, Delft University of Technology
All rights reserved
This file is part of Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

import os

import numpy as np
from matplotlib import pyplot as plt

from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.astro import element_conversion

from solarsail.sail import SolarSailGuidance



# Retrieve current directory
current_directory = os.getcwd()


###########################################################################
# INPUTS ##################################################################
###########################################################################

daysToRun = 365 * 1
fixed_step_size = 100.0
spacecraftMass = 400 # kg
# solar sail length
solarSailDim = 150 #m
# fixed_step_size = 10.0

###########################################################################

# Retrieve current directory
current_directory = os.getcwd()
print(current_directory)

simulation_start_epoch = (
    30 * constants.JULIAN_YEAR
    + 0 * 7.0 * constants.JULIAN_DAY
    - 0.5 * constants.JULIAN_DAY
)
simulation_end_epoch = simulation_start_epoch + (daysToRun * constants.JULIAN_DAY)

###########################################################################
# CREATE ENVIRONMENT ######################################################
###########################################################################


# Load spice kernels.
spice.load_standard_kernels()

# Create settings for celestial bodies
bodies_to_create = ["Sun"]
global_frame_origin = "Sun"
global_frame_orientation = "ECLIPJ2000"

body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation
)

# Create environment
bodies = environment_setup.create_system_of_bodies(body_settings)

###########################################################################
# CREATE VEHICLE ##########################################################
###########################################################################

# Create vehicle object
spacecraftName = "SOLARSAIL"
bodies.create_empty_body(spacecraftName)
bodies.get(spacecraftName).mass = spacecraftMass

###########################################################################
# CREATE THRUST MODEL #####################################################
###########################################################################

# Create thrust guidance object (e.g. object that calculates direction/magnitude of thrust)
thrust_magnitude = 0.2
solar_sail_object = SolarSailGuidance(bodies, sailName=spacecraftName)

# Create engine model (default JUICE-fixed pointing direction) with custom thrust magnitude calculation
constant_specific_impulse = 0
thrust_magnitude_settings = (
    propagation_setup.thrust.custom_thrust_magnitude_fixed_isp(
        solar_sail_object.compute_thrust_magnitude,
        constant_specific_impulse ) )
environment_setup.add_engine_model(
    'SOLARSAIL', 'SAILENGINE', thrust_magnitude_settings, bodies )

# Create vehicle rotation model such that thrust points in required direction in inertial frame
thrust_direction_function = solar_sail_object.compute_thrust_direction
rotation_model_settings = environment_setup.rotation_model.custom_inertial_direction_based(
    thrust_direction_function,
    "SOLARSAIL-fixed",
    "ECLIPJ2000" )

environment_setup.add_rotation_model( bodies, spacecraftName, rotation_model_settings)


###########################################################################
# CREATE ACCELERATIONS ####################################################
###########################################################################

# Define bodies that are propagated, and their central bodies of propagation.
bodies_to_propagate = [spacecraftName]
central_bodies = ["Sun"]

# Define accelerations acting on vehicle.
acceleration_settings_on_vehicle = dict(
        SOLARSAIL=[
        # Define the thrust acceleration from its direction and magnitude
        propagation_setup.acceleration.thrust_from_engine('SAILENGINE')
    ],
    Sun=[propagation_setup.acceleration.point_mass_gravity()]
)

# Create global accelerations dictionary.
acceleration_settings = {spacecraftName: acceleration_settings_on_vehicle}

# Create acceleration models.
acceleration_models = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings, bodies_to_propagate, central_bodies
)

###########################################################################
# CREATE PROPAGATION SETTINGS #############################################
###########################################################################

# Define initial state.
system_initial_state = spice.get_body_cartesian_state_at_epoch(
    target_body_name="Venus",  # NOTE start at earth's position (CHANGE LATER)
    observer_body_name="Sun",
    reference_frame_name="ECLIPJ2000",
    aberration_corrections="NONE",
    ephemeris_time=simulation_start_epoch,
)

# Define required outputs  panelled_radiation_pressure_acceleration_type
# acctype = propagation_setup.acceleration.panelled_radiation_pressure_acceleration_type
acctype = propagation_setup.acceleration.thrust_acceleration_type 
dependent_variables_to_save = [
    propagation_setup.dependent_variable.single_acceleration(acctype, spacecraftName, spacecraftName),
    propagation_setup.dependent_variable.single_acceleration_norm(acctype, spacecraftName, spacecraftName),
    # propagation_setup.dependent_variable.single_acceleration_norm(acctype, "SOLARSAIL", "Sun"),
    # propagation_setup.dependent_variable.heading_angle("SOLARSAIL", "Sun")
    # propagation_setup.dependent_variable.keplerian_state("SOLARSAIL", "Sun")
]

# Create numerical integrator settings.
integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)

# Create propagation settings.
termination_settings = propagation_setup.propagator.time_termination(
    simulation_end_epoch
)
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    system_initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_settings,
    # propagation_setup.propagator.cowell,
    output_variables=dependent_variables_to_save,
)

propagator_settings.print_settings.print_initial_and_final_conditions = True


###########################################################################
# PROPAGATE ORBIT #########################################################
###########################################################################

# Create simulation object and propagate dynamics.
dynamics_simulator = numerical_simulation.create_dynamics_simulator(
    bodies, propagator_settings
)

# Retrieve all data produced by simulation
propagation_results = dynamics_simulator.propagation_results

# Extract numerical solution for states and dependent variables
state_history = propagation_results.state_history
dependent_variables = propagation_results.dependent_variable_history

###########################################################################
# SAVE RESULTS ############################################################
###########################################################################

save2txt(
    solution=state_history,
    filename="SOLARSAILPropagationHistory_Q1.dat",
    directory="./",
)

save2txt(
    solution=dependent_variables,
    filename="SOLARSAILPropagationHistory_DependentVariables_Q1.dat",
    directory="./",
)

###########################################################################
# PLOT RESULTS ############################################################
###########################################################################

# Extract time and Kepler elements from dependent variables
kepler_elements = np.vstack(list(dependent_variables.values()))
time = dependent_variables.keys()
time_days = [
    t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
    for t in time
]
