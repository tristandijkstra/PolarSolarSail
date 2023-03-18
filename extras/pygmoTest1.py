# Load standard modules
import math
import pygmo
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
import numpy as np
from numpy import random

# Load pygmo library
import pygmo as pg

class HimmelblauOptimization:

    def __init__(self,
                 x_min: float,
                 x_max: float,
                 y_min: float,
                 y_max: float):

        # Set input arguments as attributes, representaing the problem bounds for both design variables
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

    def get_bounds(self):
        return ([self.x_min, self.y_min], [self.x_max, self.y_max])

    def fitness(self, x):
        # Compute Himmelblau function value
        function_value = math.pow(x[0] * x[0] + x[1] - 11.0, 2.0) + math.pow(x[0] + x[1] * x[1] - 7.0, 2.0)

        # Return list
        return [function_value]
    

# Instantiation of the UDP problem
udp = HimmelblauOptimization(-5.0, 5.0, -5.0, 5.0)

# Creation of the pygmo problem object
prob = pygmo.problem(udp)

# Print the problem's information
print(prob)



# Define number of generations
number_of_generations = 1

# Fix seed
current_seed = 171015

# Create Differential Evolution object by passing the number of generations as input
de_algo = pygmo.de(gen=number_of_generations, seed=current_seed)

# Create pygmo algorithm object
algo = pygmo.algorithm(de_algo)

# Print the algorithm's information
print(algo)



# Set population size
pop_size = 1000

# Create population
pop = pygmo.population(prob, size=pop_size, seed=current_seed)

# Inspect population (this is going to be long, uncomment if desired)
inspect_pop = False
if inspect_pop:
    print(pop)

# Set number of evolutions
number_of_evolutions = 100

# Initialize empty containers
individuals_list = []
fitness_list = []

# Evolve population multiple times
for i in range(number_of_evolutions):
    pop = algo.evolve(pop)
    individuals_list.append(pop.get_x()[pop.best_idx()])
    fitness_list.append(pop.get_f()[pop.best_idx()])

# Extract the best individual
print('\n########### PRINTING CHAMPION INDIVIDUALS ###########\n')
print('Fitness (= function) value: ', pop.champion_f)
print('Decision variable vector: ', pop.champion_x)
print('Number of function evaluations: ', pop.problem.get_fevals())
print('Difference wrt the minimum: ', pop.champion_x - np.array([3,2]))




# Extract best individuals for each generation
best_x = [ind[0] for ind in individuals_list]
best_y = [ind[1] for ind in individuals_list]

# Extract problem bounds
(x_min, y_min), (x_max, y_max) = udp.get_bounds()

# Plot fitness over generations
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(np.arange(0, number_of_evolutions), fitness_list, label='Function value')

# Plot champion
champion_n = np.argmin(np.array(fitness_list))
ax.scatter(champion_n, np.min(fitness_list), marker='x', color='r', label='All-time champion')

# Prettify
ax.set_xlim((0, number_of_evolutions))
ax.grid('major')
ax.set_title('Best individual of each generation', fontweight='bold')
ax.set_xlabel('Number of generation')
ax.set_ylabel(r'Himmelblau function value $f(x,y)$')
ax.legend(loc='upper right')
ax.set_yscale('log')
plt.tight_layout()

# Show the figure
plt.show()



# Plot Himmelblau function
grid_points = 100
x_vector = np.linspace(x_min, x_max, grid_points)
y_vector = np.linspace(y_min, y_max, grid_points)
x_grid, y_grid = np.meshgrid(x_vector, y_vector)
z_grid = np.zeros((grid_points, grid_points))
for i in range(x_grid.shape[1]):
    for j in range(x_grid.shape[0]):
        z_grid[i, j] = udp.fitness([x_grid[i, j], y_grid[i, j]])[0]

# Create figure
fig, ax = plt.subplots(figsize=(9,5))
cs = ax.contour(x_grid, y_grid, z_grid, 50)

# Plot best individuals of each generation
ax.scatter(best_x, best_y, marker='x', color='r')

# Prettify
ax.set_xlim((x_min, x_max))
ax.set_ylim((y_min, y_max))
ax.set_title('Himmelblau function', fontweight='bold')
ax.set_xlabel('X-coordinate')
ax.set_ylabel('Y-coordinate')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel(r'Himmelblau function value $f(x,y)$')
plt.tight_layout()

# Show the plot
plt.show()