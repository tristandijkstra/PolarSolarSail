from simulation import simulationV1 as sim
from solarsail.sailPhysicalV3 import SolarSailGuidance
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

initialEpoch = 1117886400 - (86400 * 4)
C3BurnVec = np.array([0,0,4000])

targetAltitude = 0.42
deepestAltitude = 0.2
sailArea = 12000
mass = 400
timesOutwardMax = 1

stepSize = 144000

yearsToRun = 25
yearInSeconds = 365 * 24 * 3600

targetInclination = 65


FTOPmin = 0.01
FTOPmax = 0.2
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
    # thermalModelObject=Thermal,
    simuFunction=sim.simulate,
    timesOutwardMax=timesOutwardMax,
    stepSize=stepSize,
    initialEpoch=initialEpoch,
    C3BurnVector=C3BurnVec
)

# Creation of the pygmo problem object
prob = pygmo.problem(udp)

# Print the problem's information
print(prob)


# Define number of generations
number_of_generations = 10

# Fix seed
current_seed = 171015

# Create Differential Evolution object by passing the number of generations as input
# de_algo = pygmo.gwo(gen=number_of_generations, seed=current_seed)
de_algo = pygmo.gaco(gen=number_of_generations, seed=current_seed)

# Create pygmo algorithm object
algo = pygmo.algorithm(de_algo)

# Print the algorithm's information
print(algo)


# Set population size
pop_size = 100

# Create population
pop = pygmo.population(prob, size=pop_size, seed=current_seed)

# Inspect population (this is going to be long, uncomment if desired)
inspect_pop = False
if inspect_pop:
    print(pop)

# Set number of evolutions
number_of_evolutions = 7

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
ax2.plot(np.arange(0, number_of_evolutions), fitness_list, label="Function value")
# plt.show()

# Plot champion
champion_n = np.argmin(np.array(fitness_list))


print(champion_n)
best_FTOP = [ind[0] for ind in individuals_list][0]
best_DEEPESTALT = [ind[1] for ind in individuals_list][0]

print(best_FTOP, best_DEEPESTALT)
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
sim.plotSimulation(satName=namee, extraText=extraTxt+extraTxt2, dataFile=save, dataDepFile=saveDep, quiverEvery=1000)


plt.show()
