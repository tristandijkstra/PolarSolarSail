from simulation import simulationV1 as sim
from solarsail.sailPhysicalV3 import SolarSailGuidance
from simulation.pygmoOptV1 import SailOptimise
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os
import pygmo
from tqdm import tqdm

import time

from thermal.thermal_model import Thermal
from multiprocessing import Process, freeze_support

saveDirectory = "data"
if not os.path.exists(saveDirectory):
    os.mkdir(saveDirectory)

logging.basicConfig(
    filename=f"{saveDirectory}/_optimiseV1.log",
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
)


def runMP(
    udp, evolutions=5, islands=2, generations=10, popsize=10, seed=42, algo=pygmo.gwo
):
    prob = pygmo.problem(udp)
    # print(prob)
    # Create Differential Evolution object by passing the number of generations as input
    de_algo = algo(gen=generations, seed=seed)
    algo = pygmo.algorithm(de_algo)

    # Create population
    # pop = pygmo.population(prob, size=pop_size, seed=current_seed)
    archi = pygmo.archipelago(
        n=islands, algo=algo, prob=prob, pop_size=popsize, seed=seed
    )
    print(archi)

    start = time.perf_counter()
    fits_log, vectors_log = [], []
    for i in tqdm(range(evolutions)):
        archi.evolve()
        archi.wait()
        # individuals_list.append(pop.get_x()[pop.best_idx()])
        # fitness_list.append(pop.get_f()[pop.best_idx()])

        vectors = [isl.get_population().get_x() for isl in archi]
        vectors_log.append(vectors)
        fits = [isl.get_population().get_f() for isl in archi]
        fits_log.append(fits)

    
    # totalEvals = 0
    # for isl in archi:
    #     isl.
    print(f"Evals= {prob.get_fevals()}")

    champsf = [x[0] for x in archi.get_champions_f()]
    champ_i = np.argmin(champsf)

    bestVector = archi.get_champions_x()[champ_i]
    end = time.perf_counter()
    duration = round((end - start) / 60, 2)

    return bestVector, duration, fits_log, vectors_log


if __name__ == "__main__":
    initialEpoch = 1117886400 - (86400 * 4)
    C3BurnVec = np.array([0, 0, 4000])

    targetAltitude = 0.42
    deepestAltitude = 0.2
    sailArea = 10000
    mass = 700
    timesOutwardMax = 1

    stepSize = 144000

    yearsToRun = 25
    yearInSeconds = 365 * 24 * 3600

    targetInclination = 65

    FTOPmin = 0.01
    FTOPmax = 0.08
    deepestAltitude_min = 0.2
    deepestAltitude_max = 0.35

    # Instantiation of the UDP problem
    udp = SailOptimise(
        FTOP_min=FTOPmin,
        FTOP_max=FTOPmax,
        deepestAltitude_min=deepestAltitude_min,
        deepestAltitude_max=deepestAltitude_max,
        solarSailGuidanceObject=SolarSailGuidance,
        mass=mass,
        sailArea=sailArea,
        thermalModelObject=Thermal,
        simuFunction=sim.simulate,
        timesOutwardMax=timesOutwardMax,
        stepSize=stepSize,
        initialEpoch=initialEpoch,
        C3BurnVector=C3BurnVec,
        verbose=False,
    )

    # Run optimisation:
    optInput = dict(
        udp=udp,
        evolutions=10,
        islands=16,
        generations=10,
        popsize=2,
        seed=42,
        algo=pygmo.gwo,
    )
    freeze_support()
    # Process(target=runMP, args=optInput).start()
    bestVector, duration, fits_log, vectors_log = runMP(**optInput)

    print(f"Duration = {duration} mins")


    best_FTOP, best_DEEPESTALT = bestVector

    # Plot
    saveFiel = "w=" + str(best_FTOP)
    guidanceObject = SolarSailGuidance(
        None,
        sailName="best",
        mass=mass,
        sailArea=sailArea,
        targetAltitude=targetAltitude,
        targetInclination=targetInclination,
        deepestAltitude=best_DEEPESTALT,
        fastTransferOptimiseParameter=best_FTOP,
        verbose=True,
    )
    namee = "FTOP=" + str(best_FTOP) + "deepestAlt=" + str(best_DEEPESTALT)
    _, save, saveDep = sim.simulate(
        spacecraftName="best",
        sailGuidanceObject=guidanceObject,
        saveFile=namee,
        yearsToRun=yearsToRun,
        simStepSize=3600,
        verbose=True,
        initialEpoch=initialEpoch,
        C3BurnVector=C3BurnVec,
    )

    (
        characteristicacceleration,
        spiralDuration,
        inclinationChangeDuration,
        finalInclination,
    ) = guidanceObject.getInclinationChangeDuration()

    dur = round(inclinationChangeDuration / yearInSeconds, 3)
    totdur = round((spiralDuration + inclinationChangeDuration) / yearInSeconds, 3)

    extraTxt = f"\nMass = {int(guidanceObject.mass)} kg | Area = {int(guidanceObject.mass/guidanceObject.sigma)} m^2"
    extraTxt2 = (
        f"\nFinal inclin. = {round(finalInclination, 3)} |"
        + f"Characteristic Acceleration = {round(characteristicacceleration*1000, 4)} mm/s^2"
        + f"\nInclin. change duration = {dur} years | Total duration = {totdur} years"
    )
    sim.plotSimulation(
        satName=namee,
        extraText=extraTxt + extraTxt2,
        dataFile=save,
        dataDepFile=saveDep,
        quiverEvery=1000,
    )

    plt.show()


# # Set number of evolutions
# number_of_evolutions = 2

# # Initialize empty containers
# individuals_list = []
# fitness_list = []

# freeze_support()
# # archi.get_champions_f()

# # Extract the best individual
# print("\n########### PRINTING CHAMPION INDIVIDUALS ###########\n")
# print("Fitness (= function) value: ", pop.champion_f)
# print("Decision variable vector: ", pop.champion_x)
# print("Number of function evaluations: ", pop.problem.get_fevals())
# print("Difference wrt the minimum: ", pop.champion_x - np.array([3, 2]))
# print("\n########### RUN TIME ###########\n")
# print(f"Run time = {round((end-start), 2)} seconds = {duration} minutes")
# print()

# print(fitness_list)
# # Extract best individuals for each generation
# best_x = [ind[0] for ind in individuals_list]
# best_y = [ind[1] for ind in individuals_list]

# # Extract problem bounds
# (x_min, y_min), (x_max, y_max) = udp.get_bounds()

# # Plot fitness over generations
# fig, ax2 = plt.subplots(figsize=(9, 5))
# ax2.plot(np.arange(0, number_of_evolutions), fitness_list, label="Function value")
# champion_n = np.argmin(np.array(fitness_list))


# print(champion_n)
# best_FTOP = [ind[0] for ind in individuals_list][0]
# best_DEEPESTALT = [ind[1] for ind in individuals_list][0]

# print(best_FTOP, best_DEEPESTALT)
# saveFiel = "w=" + str(best_FTOP)
# guidanceObject = SolarSailGuidance(
#     None,
#     sailName="best",
#     mass=mass,
#     sailArea=sailArea,
#     targetAltitude=targetAltitude,
#     targetInclination=targetInclination,
#     deepestAltitude=best_DEEPESTALT,
#     fastTransferOptimiseParameter=best_FTOP,
#     verbose=True,
# )
# namee = "FTOP=" + str(best_FTOP) + "deepestAlt=" + str(best_DEEPESTALT)
# _, save, saveDep = sim.simulate(
#     spacecraftName="best",
#     sailGuidanceObject=guidanceObject,
#     saveFile=namee,
#     yearsToRun=yearsToRun,
#     simStepSize=3600,
#     verbose=True,
#     initialEpoch=initialEpoch,
#     C3BurnVector=C3BurnVec,
# )

# (
#     characteristicacceleration,
#     spiralDuration,
#     inclinationChangeDuration,
#     finalInclination,
# ) = guidanceObject.getInclinationChangeDuration()

# dur = round(inclinationChangeDuration / yearInSeconds, 3)
# totdur = round((spiralDuration + inclinationChangeDuration) / yearInSeconds, 3)

# extraTxt = f"\nMass = {int(guidanceObject.mass)} kg | Area = {int(guidanceObject.mass/guidanceObject.sigma)} m^2"
# extraTxt2 = (
#     f"\nFinal inclin. = {round(finalInclination, 3)} |"
#     + f"Characteristic Acceleration = {round(characteristicacceleration*1000, 4)} mm/s^2"
#     + f"\nInclin. change duration = {dur} years | Total duration = {totdur} years"
# )
# sim.plotSimulation(
#     satName=namee,
#     extraText=extraTxt + extraTxt2,
#     dataFile=save,
#     dataDepFile=saveDep,
#     quiverEvery=1000,
# )


# plt.show()
