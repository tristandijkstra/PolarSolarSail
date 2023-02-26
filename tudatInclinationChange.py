###########################################################################
#
# # Numerical Astrodynamics 2022/2023
#
# # Assignment 1 - Propagation Settings
#
###########################################################################


""" 
Copyright (c) 2010-2020, Delft University of Technology
All rights reserved

This file is part of Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup


###########################################################################
# INPUTS ##################################################################
###########################################################################

daysToRun = 365 * 3
fixed_step_size = 100.0
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
# spice.load_kernel(current_directory + "/solarSail_mat_crema_5_1_150lb_v01.bsp")

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
bodies.create_empty_body("SOLARSAIL")
bodies.get("SOLARSAIL").mass = 400.0

# solar sail len
solarSailDim = 300
# define parameters of the panelled model
emissivities = [0.73]  # emissivity of each panel
areas = [solarSailDim * solarSailDim]  # area of each panel
diffusion_coefficients = [0.4]  # diffusion coefficient of each panel
panel_surface_normals = (
    [  # normals of each panel surfaces in body-fixed reference frame
        # [0.2, 0.8, 0.0],
        # [1.0, 0.0, 0.0],
        # [0.0, 1.0, 0.0],
        # [0.5, 0.5, 0.0],
        [0.0, 0.6, -0.4],
    ]
)
# define parameter for occulting body
# create radiation pressure interface settings
radiation_pressure_settings = environment_setup.radiation_pressure.panelled(
    "Sun", emissivities, areas, diffusion_coefficients, panel_surface_normals
)

rotation_model_settings = (
    environment_setup.rotation_model.orbital_state_direction_based(
        "Sun", True, False, "J2000", "VehicleFixed"
    )
)
# add radiation pressure interface to "Spacecraft" body
environment_setup.add_rotation_model(bodies, "SOLARSAIL", rotation_model_settings)
environment_setup.add_radiation_pressure_interface(
    bodies, "SOLARSAIL", radiation_pressure_settings
)

###########################################################################
# CREATE ACCELERATIONS ####################################################
###########################################################################

# Define bodies that are propagated, and their central bodies of propagation.
bodies_to_propagate = ["SOLARSAIL"]
central_bodies = ["Sun"]

# Define accelerations acting on vehicle.
acceleration_settings_on_vehicle = dict(
    Sun=[propagation_setup.acceleration.point_mass_gravity(),propagation_setup.acceleration.panelled_radiation_pressure(),]
)

# Create global accelerations dictionary.
acceleration_settings = {"SOLARSAIL": acceleration_settings_on_vehicle}

# Create acceleration models.
acceleration_models = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings, bodies_to_propagate, central_bodies
)

###########################################################################
# CREATE PROPAGATION SETTINGS #############################################
###########################################################################

# Define initial state.
system_initial_state = spice.get_body_cartesian_state_at_epoch(
    target_body_name="Mercury",  # NOTE start at earth's position (CHANGE LATER)
    observer_body_name="Sun",
    reference_frame_name="ECLIPJ2000",
    aberration_corrections="NONE",
    ephemeris_time=simulation_start_epoch,
)

# Define required outputs  panelled_radiation_pressure_acceleration_type
acctype = propagation_setup.acceleration.panelled_radiation_pressure_acceleration_type
dependent_variables_to_save = [
    propagation_setup.dependent_variable.single_acceleration(acctype, "SOLARSAIL", "Sun"),
    propagation_setup.dependent_variable.keplerian_state("SOLARSAIL", "Sun")
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
