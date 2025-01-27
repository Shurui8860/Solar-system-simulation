# dataProvider
# author: Shurui Zhang
# UUN: s2193407

# This python file is used to read the input data from the file,
# and initialise the velocity of every body in this system if necessary.

from body import Body
from simulation import Simulation
import numpy as np


# read_input_data is responsible for read the simulation information from a txt file,
# and convert the information to an instance of class Simulation
def read_input_data(filename):
    # several lists are created to remember the information of these bodies
    # including their names, masses, initial positions, initial velocities, radii, and colours.
    name_list = []
    mass_list = []
    position_list = []
    velocity_list = []
    colour_list = []
    radius_list = []
    body_list = []

    with open(filename) as file:
        for line in file:
            # blank spaces are removed from each line
            newline = line.replace(" ", "")
            # determining what kind of information this line stores by the starting characters of each line
            # and then remembering the information in the corresponding list.
            if newline.startswith('timestep:'):
                timestep = float(newline[9:].strip("\n"))

            elif newline.startswith('numberofiterations:'):
                num_of_iteration = int(newline[19:].strip("\n"))

            elif newline.startswith('name:'):
                name_list.append(newline[5:].strip("\n"))

            elif newline.startswith('mass:'):
                mass_list.append(float(newline[5:].strip("\n")))

            elif newline.startswith('initialposition:'):
                # from the file formate: initial position is stored as x,y
                # splitting x,y to (x,y), converting both parts from string to float,
                # and converting the data type to numpy array.
                position = list(map(float, newline[16:].strip("\n").split(',')))
                position_list.append(np.array(position))

            elif newline.startswith('initialvelocity:'):
                # from the file formate: initial velocity is stored as x,y
                # splitting x,y to (x,y), converting both parts to float,
                # and converting the data type to numpy array.
                velocity = list(map(float, newline[16:].strip("\n").split(',')))
                velocity_list.append(np.array(velocity))

            elif newline.startswith('radius:'):
                radius_list.append(float(newline[7:].strip("\n")))

            elif newline.startswith('colour:'):
                colour_list.append(newline[7:].strip("\n"))

    # using the information copied from the file to create several body instances.
    for i in range(len(name_list)):
        body_i = Body(name_list[i], mass_list[i], position_list[i],
                      velocity_list[i], radius_list[i], colour_list[i])
        body_list.append(body_i)
    # then create and return the instance of class simulation
    return timestep, num_of_iteration, body_list


# In this simulation, data for the initial velocity is not given,
# initialise the velocities of every celestial body in this system
# by assuming uniform circular motion around a central body.
def read_init_v(filename):
    timestep, num_of_iteration, body_list = read_input_data(filename)
    # the 0th element in the body_list is assumed to be the central body in this simulation.
    # attention should be paid when writing the txt file.
    # central body should always be added to the body_list first
    central_body = body_list[0]
    for i in range(1, len(body_list)):
        (body_list[i]).init_velocity(central_body)
    simulation = Simulation(timestep, num_of_iteration, body_list)
    return simulation
