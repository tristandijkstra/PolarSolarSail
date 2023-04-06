from simulation import simulationV1 as sim
from solarsail.sailPhysicalV3_opt import SolarSailGuidance
from simulation.pygmoOptV1 import SailOptimise
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import os
import pygmo

import time

from thermal.thermal_model import Thermal

saveDirectory = "data"
if not os.path.exists(saveDirectory):
    os.mkdir(saveDirectory)

logging.basicConfig(
    filename=f"{saveDirectory}/_optimiseV1.log",
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
)

initialEpoch = 1117800000
C3BurnVec = np.array([0.084108927717811,0.240321913434112,4.946831358585164]) * 1000

targetAltitude = 0.48
# deepestAltitude = 0.4
sailArea = 10000
mass = 590
timesOutwardMax = 1

stepSize = 144000

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600

# targetInclination = 75-7.25
targetInclination = 60-7.25


FTOPmin = 0.016
FTOPmax = 0.03
deepestAltitude_min = 0.372
deepestAltitude_max = 0.372
endPrecision = 0.01



# Instantiation of the UDP problem
udp = SailOptimise(
    FTOP_min=FTOPmin,
    FTOP_max=FTOPmax,
    deepestAltitude_min=deepestAltitude_min,
    deepestAltitude_max=deepestAltitude_max,
    solarSailGuidanceObject=SolarSailGuidance,
    mass=mass,
    sailArea=sailArea,
    targetInclination=targetInclination,
    # thermalModelObject=Thermal,
    simuFunction=sim.simulate,
    timesOutwardMax=timesOutwardMax,
    stepSize=stepSize,
    initialEpoch=initialEpoch,
    C3BurnVector=C3BurnVec,
    verbose=True,
    endPrecision=endPrecision
)

# Creation of the pygmo problem object
prob = pygmo.problem(udp)

# Print the problem's information
print(prob)


# Define number of generations
number_of_generations = 7

# Fix seed
current_seed = 420

# Create Differential Evolution object by passing the number of generations as input
# de_algo = pygmo.gwo(gen=number_of_generations, seed=current_seed)
de_algo = pygmo.de(gen=number_of_generations, seed=current_seed,CR=0.8, xtol=1e-4)
# de_algo = pygmo.sga(gen=number_of_generations, seed=current_seed)
# de_algo = pygmo.bee_colony(gen=number_of_generations, seed=current_seed)
# de_algo = pygmo.gaco(gen=number_of_generations, seed=current_seed)

# Create pygmo algorithm object
algo = pygmo.algorithm(de_algo)

# Print the algorithm's information
print(algo)


# Set population size
pop_size = 5

# Create population
pop = pygmo.population(prob, size=pop_size, seed=current_seed)

# Inspect population (this is going to be long, uncomment if desired)
inspect_pop = False
if inspect_pop:
    print(pop)

# Set number of evolutions
number_of_evolutions = 5

# Initialize empty containers
individuals_list = []
fitness_list = []

start = time.perf_counter()
# Evolve population multiple times
for i in range(number_of_evolutions):
    pop = algo.evolve(pop)
    individuals_list.append(pop.get_x()[pop.best_idx()])
    fitness_list.append(pop.get_f()[pop.best_idx()])

end = time.perf_counter()
# Extract the best individual
print("\n########### PRINTING CHAMPION INDIVIDUALS ###########\n")
print("Fitness (= function) value: ", pop.champion_f)
print("Decision variable vector: ", pop.champion_x)
print("Number of function evaluations: ", pop.problem.get_fevals())
print("Difference wrt the minimum: ", pop.champion_x - np.array([3, 2]))
print("\n########### RUN TIME ###########\n")
print(
    f"Run time = {round((end-start), 2)} seconds = {round((end-start)/60, 2)} minutes"
)
print()

print(fitness_list)
# Extract best individuals for each generation
best_x = [ind[0] for ind in individuals_list]
best_y = [ind[1] for ind in individuals_list]

# Extract problem bounds
(x_min, y_min), (x_max, y_max) = udp.get_bounds()

# Plot fitness over generations
fig, ax2 = plt.subplots(figsize=(9, 5))
ax2.plot(np.arange(0, number_of_evolutions), np.minimum(np.array(fitness_list),25), label="Function value")
# plt.show()

# Plot champion
champion_n = np.argmin(np.array(fitness_list))


# print(champion_n)
best_FTOP = pop.champion_x[0]
best_DEEPESTALT = pop.champion_x[1]

print("FTOP:", best_FTOP, "deepestAlt:",  best_DEEPESTALT)
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
    thermalModel=Thermal(24*3600, verbose=True),
    verbose=True,
)
namee = "FTOP=" + str(best_FTOP) + " deepestAlt=" + str(best_DEEPESTALT)
_, save, saveDep = sim.simulate(
    spacecraftName="best",
    sailGuidanceObject=guidanceObject,
    saveFile=namee,
    yearsToRun=yearsToRun,
    simStepSize=7200,
    verbose=True,
    initialEpoch=initialEpoch,
    C3BurnVector=C3BurnVec
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
extraTxt2 = f"\nFinal inclin. = {round(finalInclination, 3)} |" + \
            f"Characteristic Acceleration = {round(characteristicacceleration*1000, 4)} mm/s^2" + \
            f"\nInclin. change duration = {dur} years | Total duration = {totdur} years"
sim.plotSimulation(satName=namee, extraText=extraTxt+extraTxt2, dataFile=save, dataDepFile=saveDep, quiverEvery=1000, thermalOn=True)
sim.plotThermal(satName=namee, dataDepFile=saveDep, thermalNodeNames=guidanceObject.thermalModel.node_keys)


planets = ["Earth", "Venus", "Mercury"]
finalEpoch = initialEpoch + spiralDuration + inclinationChangeDuration



for planet in planets:
    print(f"Running {planet}")
    sim.simulatePlanet(planet, initialEpoch, finalEpoch, 7200)


plt.show()
