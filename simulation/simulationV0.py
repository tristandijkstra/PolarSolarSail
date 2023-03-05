import os

import numpy as np
from matplotlib import pyplot as plt

from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup

# from solarsail.sailBasic import SolarSailGuidance


def simulate(
    spacecraftName,
    sailGuidanceObject,
    saveFile: str,
    yearsToRun: float = 7,
    simStepSize: float = 100.0,
):
    ###########################################################################
    # INPUTS ##################################################################
    ###########################################################################

    daysToRun = 365 * yearsToRun
    spacecraftMass = sailGuidanceObject.mass  # kg

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
    bodies.create_empty_body(spacecraftName)
    bodies.get(spacecraftName).mass = spacecraftMass

    ###########################################################################
    # CREATE THRUST MODEL #####################################################
    ###########################################################################
    sailGuidanceObject.bodies = bodies

    # Create engine model (default JUICE-fixed pointing direction) with custom thrust magnitude calculation
    constant_specific_impulse = 0
    thrust_magnitude_settings = (
        propagation_setup.thrust.custom_thrust_magnitude_fixed_isp(
            sailGuidanceObject.compute_thrust_magnitude, constant_specific_impulse
        )
    )
    environment_setup.add_engine_model(
        spacecraftName, "SAILENGINE", thrust_magnitude_settings, bodies
    )

    # Create vehicle rotation model such that thrust points in required direction in inertial frame
    thrust_direction_function = sailGuidanceObject.compute_thrust_direction
    rotation_model_settings = (
        environment_setup.rotation_model.custom_inertial_direction_based(
            thrust_direction_function, "SOLARSAIL-fixed", "ECLIPJ2000"
        )
    )

    environment_setup.add_rotation_model(
        bodies, spacecraftName, rotation_model_settings
    )

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated, and their central bodies of propagation.
    bodies_to_propagate = [spacecraftName]
    central_bodies = ["Sun"]

    acceleration_settings_on_vehicle = {
        spacecraftName: [propagation_setup.acceleration.thrust_from_engine("SAILENGINE")],
        "Sun": [propagation_setup.acceleration.point_mass_gravity()],
    }

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
        propagation_setup.dependent_variable.single_acceleration(
            acctype, spacecraftName, spacecraftName
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            acctype, spacecraftName, spacecraftName
        ),
        propagation_setup.dependent_variable.keplerian_state(spacecraftName, "Sun"),
    ]

    # Create numerical integrator settings.
    integrator_settings = propagation_setup.integrator.runge_kutta_4(simStepSize)

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

    # fileName = current_directory + "data/" + saveFile + ".dat"
    # fileNameDep = current_directory + "data/" + saveFile + "_dep" + ".dat"
    fileName = "data/" + saveFile + ".dat"
    fileNameDep = "data/" + saveFile + "_dep" + ".dat"
    save2txt(
        solution=state_history,
        filename=fileName,
        directory="./",
    )

    save2txt(
        solution=dependent_variables,
        filename=fileNameDep,
        directory="./",
    )

    return sailGuidanceObject, fileName, fileNameDep
