# simulation
# author: Shurui Zhang
# UUN: s2193407
import math
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from body import Body
import copy


class Simulation(object):

    # --------------------constructor and init methods--------------------

    def __init__(self, timestep, num_of_iterations, body_list):
        self.timestep = timestep
        self.num_of_iterations = num_of_iterations
        self.body_list = copy.deepcopy(body_list)
        self.num_of_bodies = len(self.body_list)
        self.original_body_list = body_list

        # the acceleration of each celestial body in this system is initialised.
        self.init_accs()

        # period_count is used to count how much time is passed before a celestial body finishes a period.
        # each iteration, all the count will be added by 1, and when a body finishes one period,
        # the corresponding count will be set to 0. (refer to method period_compute)
        self.period_count = np.zeros(self.num_of_bodies)

        # patches will be used in the animation to add a visual representation of each body.
        self.patches = []

    def init_accs(self):
        initial_accs = [self.calc_acc(body) for body in self.body_list]
        for (body, new_acc) in zip(self.body_list, initial_accs):
            body.init_acc(new_acc)

    def init_simulation(self):
        self.num_of_bodies = len(self.original_body_list)
        self.body_list = [copy.deepcopy(body) for body in self.original_body_list]
        # the acceleration of each celestial body in this system is initialised.
        self.init_accs()
        self.period_count = np.zeros(self.num_of_bodies)
        self.patches = []

    # --------------------getter and setter methods--------------------
    # The following are some get methods to get access to the accelerations of each body,
    # or get the position, velocity, or colour of a given body.
    # These methods are used to make the code readable and easy to follow.

    # get the acceleration of a given body.
    def calc_acc(self, body):
        acc = [body.calc_acc(other_body)
               for other_body in self.body_list if body != other_body]
        return sum(acc)

    # get the acceleration of each body in this system.
    def get_accs(self):
        return [self.calc_acc(body) for body in self.body_list]

    # get the names of all bodies in this system.
    def get_names(self):
        return [body.name for body in self.body_list]

    # get the position of the body with index i in body_list
    def get_position(self, i):
        return (self.body_list[i]).position

    # get the velocity of the body with index i in body_list
    def get_velocity(self, i):
        return (self.body_list[i]).velocity

    # get the colour of the body with index i in body_list
    def get_colour(self, i):
        return (self.body_list[i]).colour

    # get the radius of the body with index i in body_list
    def get_radius(self, i):
        return (self.body_list[i]).radius

    # get the index of the body with the given name.
    def get_index(self, name):
        for i in range(0, self.num_of_bodies):
            if self.body_list[i].name == name:
                return i

    # reset the number of iterations.
    def set_num_of_iterations(self, num_of_iterations):
        self.num_of_iterations = num_of_iterations

    # --------------------energy methods--------------------
    # The following are some methods to calculate the total kinetic energy,
    # total potential energy, and the total energy, and write the total energy in a txt file.

    # This method is used to calculate the total kinetic energy in this system.
    def cal_tot_KE(self):
        # every body's kinetic energy in this system is computed, and then its sum is returned.
        kinetic_energy = [body.calc_KE() for body in self.body_list]
        return sum(kinetic_energy)

    # This method is to calculate a given body's potential energy in this system.
    def calc_PE(self, body):
        potential_energy = [body.calc_PE(other_body)
                            for other_body in self.body_list if other_body != body]
        return sum(potential_energy)

    # This method is to calculate the total potential energy in this system.
    def calc_tot_PE(self):
        tot_PE = [self.calc_PE(body) for body in self.body_list]
        return sum(tot_PE) / 2

    # This method is used to calculate the total energy in this system.
    def calc_tot_energy(self):
        return self.cal_tot_KE() + self.calc_tot_PE()

    # This method is used to write the total energy in a txt file called "total_energy.txt".
    def write_energy(self):
        with open("total_energy.txt", mode='a') as f:
            f.write("The total energy in this solar system is " + str(self.calc_tot_energy()) + " j.\n")
            f.close()

    # --------------------animation methods-------------------

    # This method is to update the state of each body for one timestep.
    def step_forward(self):
        # update the period_count.
        self.period_count = self.period_count + 1

        # Update positions: loop over all bodies.
        for i in range(self.num_of_bodies):
            (self.body_list[i]).update_position(self.timestep)

        # Calculate new accelerations: loop over all bodies.
        new_accs = self.get_accs()

        # Update velocities: loop over all bodies.
        for (body, new_acc) in zip(self.body_list, new_accs):
            body.update_velocity(new_acc, self.timestep)

        # Update accelerations (for next timestep): loop over all bodies.
        for (body, new_acc) in zip(self.body_list, new_accs):
            body.update_acc(new_acc)

    # print the period in earth year
    # if some celestial body (except the central body) has finished one period.
    def period_compute(self):
        for i in range(1, self.num_of_bodies):
            body = self.body_list[i]
            #  check whether some celestial body has finished one period
            if body.check_orbital_period():
                # change the period unit from second to earth year.
                period = self.period_count[i] * self.timestep / (3600 * 24 * 365)
                print("The orbital period of " + body.name + " is " + str(period) + " earth years.")
                # set the corresponding count to its initial value 0.
                self.period_count[i] = 0

    def crash_test(self):
        for i in range(self.num_of_bodies - 1):
            for j in range(i, self.num_of_bodies):
                distance = norm(self.get_position(i) - self.get_position(j))
                sum_of_radii = self.get_radius(i) + self.get_radius(j)
                if distance < sum_of_radii:
                    return True
        return False

    # update the visual representation of the animation.
    def update_display(self, frame):
        # update each celestial body's state.
        self.step_forward()

        # print the period if necessary.
        self.period_compute()

        # calculate the kinetic energy at suitable intervals.
        # here the interval is set to be 100 iterations.
        if frame % 100 == 0:
            self.write_energy()
        # update the position of these circles.
        for i in range(self.num_of_bodies):
            self.patches[i].center = self.get_position(i)
        return self.patches

    # This method is used to generate the animation.
    def run(self, title):
        # create plot elements
        fig = plt.figure()
        ax = plt.axes()
        # loop over all bodies to initialise the list patches.
        for i in range(self.num_of_bodies):
            self.patches.append(plt.Circle(self.get_position(i), 8e9,
                                           color=self.get_colour(i), animated=True))
        # add circles to axes
        for i in range(self.num_of_bodies):
            ax.add_patch(self.patches[i])

        # set up the axes
        ax.set_title(title)
        ax.axis('scaled')
        ax.set_xlim(- 2.5e11, 2.5e11)
        ax.set_ylim(- 2.5e11, 2.5e11)
        ax.set_xlabel('x axis')
        ax.set_ylabel('y axis')

        # create the animator
        self.anim = FuncAnimation(fig, self.update_display, repeat=False,
                             frames=self.num_of_iterations, interval=0.01, blit=True)
        plt.show()

    # This method is to plot the movement track of each body in this solar system.
    def plot_orbital_graph(self):
        # body_movement_track remembers each body's position in each iteration.
        body_movement_track = np.zeros(shape=(self.num_of_bodies, self.num_of_iterations, 2))
        # go through the simulation, and update the body_movement_track.
        for i in range(self.num_of_iterations):
            self.step_forward()
            for j in range(self.num_of_bodies):
                # update the position of jth body in the ith iteration.
                body_movement_track[j][i] = self.get_position(j)
        # each celestial body's movement track is plotted.
        for i in range(self.num_of_bodies):
            x = [x for x, y in body_movement_track[i]]
            y = [y for x, y in body_movement_track[i]]
            plt.plot(x, y, self.get_colour(i), label=self.body_list[i].name)

        leg = plt.legend(loc='best', ncol=2, mode="expand", shadow=True, fancybox=True)
        leg.get_frame().set_alpha(0.3)
        plt.title('Movement Track of Celestial Bodies in Solar System')
        plt.show()

    # --------------------experiment methods-------------------

    # --------------------experiment 1-------------------
    # --------------------periods test-------------------
    def period_test(self):
        # element of period_list has the format [total time taken, number of periods finished]
        period_list = np.zeros(shape=(self.num_of_bodies - 1, 2))

        for i in range(self.num_of_iterations):
            self.step_forward()  # each time, the state of each body, and the period_list is updated.
            # loop over all bodies except the central body.
            for j in range(1, self.num_of_bodies):
                body = self.body_list[j]
                # check whether the jth body has just finished one period.
                if body.check_orbital_period():
                    # the number of periods finished is updated.
                    period_list[j - 1][1] = period_list[j - 1][1] + 1
                    # change the unit of the total time taken from second to day.
                    period_list[j - 1][0] = i * self.timestep / (3600 * 24)
        # This if statement is to avoid division by zero,
        # since some celestial body may ont have finished one period yet.
        if all(pair[1] != 0 for pair in period_list):
            # pair = [total_time_taken, number of periods finishes]
            periods = [pair[0] / pair[1] for pair in period_list]
        else:
            print("Please enter a larger number of iterations.")
            return

        names = self.get_names()

        # the names of non-central bodies and their periods are zipped.
        data = zip(names[1:self.num_of_bodies], periods)

        # write the average orbital period in this simulation in a file.
        # and compares them with the actual ones.
        with open("average_periods.txt", mode='a') as file:
            for (name, period) in data:
                file.write("The average orbital period of " + name +
                           " in this simulation is " + str(period) + " days.\n")
        file.close()

    # --------------------    experiment  two     -------------------
    # --------------------energy conservation test-------------------

    # provide the total energy in each iteration.
    def get_energy_data(self):
        # the energy_list remembers the total energy in this system for each iteration.
        list_of_tot_energy = []
        for i in range(self.num_of_iterations):
            self.step_forward()
            # for each iteration, add the total energy to the energy list
            list_of_tot_energy.append(self.calc_tot_energy())
        return list_of_tot_energy

    def get_energy_data_variance(self):
        list_of_tot_energy = self.get_energy_data()
        num_of_data = len(list_of_tot_energy)
        average = sum(list_of_tot_energy)/len(list_of_tot_energy)
        variance = sum([(energy-average)**2/num_of_data for energy in list_of_tot_energy])
        return variance

    def energy_conservation_test(self):
        list_of_tot_energy = self.get_energy_data()
        # Observe how the total energy fluctuates in a small range.
        # In this simulation, the range is set to be (1+-1e-8)*(the initial total energy).
        range_min = (1 - 1e-8) * list_of_tot_energy[1]
        range_max = (1 + 1e-8) * list_of_tot_energy[1]
        plt.ylim([range_min, range_max])
        plt.plot(list_of_tot_energy)
        plt.title('Total Energy of the System')
        plt.xlabel("Number of iterations")
        plt.ylabel("Total energy")
        plt.show()

        # another graph is plotted to observe the minor fluctuation.
        # zoom in to observe the fluctuation of the total energy graph.
        plt.plot(list_of_tot_energy)
        plt.title('Total Energy of the System')
        plt.xlabel("Number of iterations")
        plt.ylabel("Total energy")
        plt.show()

    # -------------------- experiment 3 -------------------
    # --------------------satellite test-------------------

    # add a body to the simulation.
    # In this simulation, the method is mainly used to launch a satellite,
    # or generate Asteroids.
    def add_body(self, satellite):
        # the body is added to both original list and current body_list
        self.original_body_list.append(copy.deepcopy(satellite))
        self.body_list.append(copy.deepcopy(satellite))
        # num_of_bodies and period_count is updated to avoid index_out_of_range error
        self.num_of_bodies = self.num_of_bodies + 1
        self.period_count = np.append(self.period_count, 0)

    # This method is to launch a satellite
    # return the total time taken for the satellite to get near Mars,
    # and the smallest distance to Mars.
    def start_mission(self, satellite):
        simulation = copy.deepcopy(self)
        simulation.add_body(satellite)  # add a satellite to the system
        returned = False
        min_to_mars = math.inf  # smallest distance to Mars
        time_taken = 0  # total time taken to get to the nearest place
        # get the index of Mars, Earth, and Satellite
        index_of_mars = simulation.get_index("Mars")
        index_of_earth = simulation.get_index("Earth")
        index_of_satellite = simulation.get_index("Satellite")

        for i in range(simulation.num_of_iterations):
            simulation.step_forward()
            distance_to_mars = \
                norm(simulation.get_position(index_of_satellite) - simulation.get_position(index_of_mars))
            distance_to_earth = \
                norm(simulation.get_position(index_of_satellite) - simulation.get_position(index_of_earth))
            # if the distance is shorter, update the min_to_mars, and time_taken
            if min_to_mars > distance_to_mars:
                min_to_mars = distance_to_mars
                time_taken = i
            # if the satellite has already landed on Mars, stop the for loop.
            if distance_to_mars <= simulation.get_radius(index_of_mars):
                break
            # if the satellite ever come back to Earth,
            # stop the loop, and change the value of returned
            if distance_to_earth <= simulation.get_radius(index_of_earth):
                returned = True
                break
        time_taken = time_taken * simulation.timestep / (3600 * 24)  # change the unit of time to day
        return min_to_mars, time_taken, returned

    # search in a range of initial x_velocity, [range_min, range_max]
    # and return the ideal x_velocity, which will take the satellite close to Mars
    def satellite_to_mars(self, range_min, range_max, search_density):
        velocity_range = np.linspace(range_min, range_max, search_density)  # the range of velocities to check
        ideal_x_velocity = 0
        ideal_launching_time_taken = 0
        returned = False
        # the initial position is the [(radius of earth + orbital radius of earth), 0]
        initial_position = np.array([1.496664e11, 0])
        min_distance_to_mars = math.inf  # the current minimal distance to mars

        for v in velocity_range:
            # 2.98e4 is the orbital speed od earth
            # assume the fuel is all used to get the x_velocity
            initial_velocity = np.array([v, 2.98e4])

            # create the satellite
            satellite = Body("Satellite", 3500, initial_position, initial_velocity, 1, 'y')

            # get the data of this launching, including minimal distance of this launching, and time taken
            (min_distance_this_time, time_taken, returned) = self.start_mission(satellite)
            print("\ninitial x-velocity: " + str(v) + " m/s." +
                  "\nsmallest distance to Mars: " + str(min_distance_this_time) + " m." +
                  "\ntime taken: " + str(time_taken) + " days." +
                  "\never returned to earth: " + str(returned))

            # update the minimal distance to Mars and the ideal x_velocity and time taken.
            if min_distance_this_time < min_distance_to_mars:
                min_distance_to_mars = min_distance_this_time
                ideal_x_velocity = v
                ideal_launching_time_taken = time_taken
        return ideal_x_velocity, min_distance_to_mars, ideal_launching_time_taken, returned

    # This method is used to search for satellite launching until an ideal one has been found.
    # i.e. until the satellite can get really close to Mars (the distance is less than the max_distance)
    def search_suitable_velocity(self, max_distance, range_min, range_max):
        # get the ideal data from this searching
        # (current minimal distance to Mars, the ideal x_velocity so far, the time taken of the ideal launching)
        ideal_velocity, current_min_distance, time_taken, returned \
            = self.satellite_to_mars(range_min, range_max, 20)
        range_length = range_max - range_min
        range_min = ideal_velocity - range_length / 10
        range_max = ideal_velocity + range_length / 10

        while current_min_distance > max_distance:
            # get the ideal data from this searching
            ideal_velocity, current_min_distance, time_taken, returned \
                = self.satellite_to_mars(range_min, range_max, 5)

            # calculate the length of the range, and update the range.
            range_length = range_max - range_min
            range_min = ideal_velocity - range_length / 10
            range_max = ideal_velocity + range_length / 10

        return ideal_velocity, current_min_distance, time_taken, returned
