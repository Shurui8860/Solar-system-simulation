# experiments
# author: Shurui Zhang
# UUN: s2193407

import copy
import numpy as np
import matplotlib.pyplot as plt
from body import Body
import dataProvider
from simulation import Simulation


# This method generates animation of the solar system,
# and plot the movement track of each celestial body.
def solar_system_animation():
    simulation1 = dataProvider.read_init_v('demo.txt')
    simulation1.run("animation of solar system")
    simulation1.init_simulation()
    simulation1.plot_orbital_graph()


# Two graphs will be plotted:
# the first one is to observe how the total energy fluctuates in a small range.
# the other graph is plotted to observe the minor fluctuation.
def energy_conservation_experiment():
    simulation2 = dataProvider.read_init_v('demo.txt')
    simulation2.energy_conservation_test()


# The average orbital period of each body will be written in the file average_period.txt
def period_experiment():
    simulation3 = dataProvider.read_init_v('demo.txt')
    simulation3.period_test()


# an ideal x_velocity will be returned.
def satellite_launching_searching():
    simulation4 = dataProvider.read_init_v('satellite_simulation.txt')
    ideal_x_velocity, min_distance, time_taken, returned\
        = simulation4.search_suitable_velocity(4e6, 1e4, 1.5e4)
    print("\nthe suitable initial velocity of the satellite is " + str(ideal_x_velocity) + " m/s."
          + "\nThe minimal distance to Mars is " + str(min_distance) + "m."
          + "\nThe total time taken to Mars is " + str(time_taken) + " days."
          + "\nEver returned to earth: " + str(returned) + "\n")
    return ideal_x_velocity


# simulation of satellite launching with the ideal x_velocity,
# and the movement track graph will be plotted.
def satellite_launching_experiment():
    ideal_velocity = satellite_launching_searching()
    simulation5 = dataProvider.read_init_v('demo.txt')
    satellite = Body("Satellite", 3500, np.array([1.496664e11, 0]),
                     np.array([ideal_velocity, 29800]), 1, 'k')
    simulation5.add_body(satellite)
    simulation6 = copy.deepcopy(simulation5)
    simulation5.run("animation of satellite launching")
    simulation6.plot_orbital_graph()

# plot the variance graph of the total energy.
def energy_conservation_variance():
    timesteps = np.linspace(1, 10000, 100)
    variances = []
    for timestep in timesteps:
        simulation = dataProvider.read_init_v('energy_conservation_variance.txt')
        simulation.set_timestep(timestep)
        variances.append(simulation.get_energy_data_variance())
    plt.plot(timesteps, variances)
    plt.title('Energy variance of the System')
    plt.xlabel("length of each timestep")
    plt.ylabel("total energy variance")
    plt.show()

