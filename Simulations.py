"""Simulation is a parent class of differnt simulation models. Actuall simulations are run through other files
which import this file to get global functions and classes"""

import hoomd
import hoomd.md
import math
import matplotlib.pyplot as plt
import numpy
import time
import matplotlib.animation as animation
from Equations import *

class Simulation:

    def __init__(self, dt=.0001):

        self.k1 = 0
        self.k2 = 0
        self.k3 = 0
        self.static_forces = []
        self.dynamic_forces = []
        self.force_eq = []
        self.system = None
        self.disabled = False
        self.first_run = True
        self.Eq1 = None
        self.Eq2 = None
        self.Eq3 = None
        self.harmonic = None
        self.integrate_group = None
        self.xml = None
        self.p_v_recorded = False
        self.p_p_recorded = False
        self.p_k_recorded = False
        self.dt = dt
        self.m = 1
        hoomd.context.initialize("")


    def create_lattice(self):
        raise NotImplementedError("Need a way to create the lattice")

    def graph_p_pos(self):
        """Create a graph of the system

        > simulation.graph_p_pos()
        """

        positions_x = []
        positions_y = []

        #record the positions of the systems particles=
        for p in self.system.particles:
            positions_x.append(p.position[0])
            positions_y.append(p.position[1])

        # find largest particle position value
        largest_value = 0
        for x in positions_x:
            if abs(x) > largest_value:
                largest_value = abs(x)
        for y in positions_y:
            if abs(y) > largest_value:
                largest_value = abs(y)

        # graph particle positions
        plt.figure(figsize=(4, 4), dpi=140)
        plt.plot(positions_x, positions_y, 'ro')
        plt.axis([-int(largest_value) - 1, int(largest_value) + 1, -int(largest_value) - 1, int(largest_value) + 1])
        for i in range(len(positions_x)):
            x = positions_x[i]
            y = positions_y[i]
            plt.text(x + 0.1, y + 0.1, i, fontsize=9)

    def dimension_constrain(self, constraint_vec: list):
        """Constrian the system to move in a single dimension
        constraint_vec is a list of three numbers [x,y,z]. 0 means particles are pinned from moving in that dimension
        1 means particles are free to move.

        #1D in x direction
        > dimension_constrian([1,0,0])
        """

        all = hoomd.group.all()
        hoomd.md.constrain.oneD(group=all, constraint_vector=constraint_vec)

    def particle_dimension_constrian(self, constraint_vec: list, p_list: list):
        """Constrian a set of particles to move in a single dimension
        constraint_vec is a list of three numbers [x,y,z]. 0 means particles are pinned from moving in that dimension
        1 means particles are free to move.

        #1D in x direction
        > particle_dimension_constrian([1,0,0], [0])
        """

        group = hoomd.group.tag_list(name="group", tags=p_list)
        hoomd.md.constrain.oneD(group=group, constraint_vector=constraint_vec)

    def pin_particles(self, tags: list):
        """Pins particles in place by updating the group of particles that is integrated to one the does not include
        the points with tags equal to the tags given in the tag list

        > simulation.pin_particles([0,10,15])
        """

        non_excluded = []
        for p in self.system.particles:
            if p.tag not in tags:
                non_excluded.append(p.tag)
        self.integrate_group = hoomd.group.tag_list(name="group", tags=non_excluded)

    def unpin_particles(self):
        """Resets the integrated group to all particles in the system

        > simulation.unpin_particles()
        """

        self.integrate_group = hoomd.group.all()

    def apply_static_force(self, tag_list: list, fx: int, fy: int, fz: int):
        """Sets a force to act on the lattice when system is moved forward in time.
        tag_list is a list of the particle tags of the particles on which the force should be applied.
        fx, fy, and fz are negative magnitudes of the applied force

        > simulation.apply_force([0,10,15], 10, 10, 0)
        """

        # group together the points to which the force will be applied
        force_pt = hoomd.group.tag_list(name='group', tags=tag_list)

        # set the force
        self.static_forces.append(hoomd.md.force.constant(fvec=(fx,fy,fz), fx=fx, fy=fy, fz=fz, group=force_pt))

    def disable_static_force(self, force_i: int):
        """Command to tell the simulation to ignore an applied static force when calculating future particle
        positions. Force will remain disabled until enabled"""
        force = self.static_forces[i]
        hoomd.md.force.constant.disable(force)

    def enable_static_force(self, force_i: int):
        """Command to re-enable a static force that has been disabled"""
        force = self.static_forces[i]
        hoomd.md.force.constant.disable(force)

    def log_PE(self, period=500):
        """Record the potential energy and temperature of the system in a file, base period is 500 time steps
        """

        hoomd.analyze.log(filename="PE_and_Temp.log",
                          quantities=['potential_energy', 'temperature'],
                          period=period,
                          overwrite=True)

    def log_velocity(self, timestep):
        """Record all particle velocities to a file so that they can be latter be analysed"""
        if self.p_v_recorded:
            with open("velocities.txt", "a") as f:
                line = "\n" + str(timestep) + "   "
                for p in self.system.particles:
                    velocity = p.velocity
                    line += str(velocity)
                    line += "   "
                f.write(line)
        else:
            with open("velocities.txt", "w") as f:
                N = len(self.system.particles)
                dt = self.dt
                m = self.m
                line = f"Simulation particle velocites, run dt= {dt} mass= {m} N= {N} \n"
                f.write(line)
                line = "timestep   "
                for p in self.system.particles:
                    line += f"p{p.tag}   "
                f.write(line)
                self.p_v_recorded = True

    def log_positions(self, timestep):
        """Record all particle positions to a file so that they can be latter be analysed"""
        if self.p_p_recorded:
            with open("positions.txt", "a") as f:
                line = "\n" + str(timestep) + "   "
                for p in self.system.particles:
                    velocity = p.position
                    line += str(velocity)
                    line += "   "
                f.write(line)
        else:
            with open("positions.txt", "w") as f:
                N = len(self.system.particles)
                dt = self.dt
                m = self.m
                a = self.a
                line = f"Simulation particle positions, run dt= {dt} mass= {m} a= {a} N= {N} \n"
                f.write(line)
                line = "timestep   "
                for p in self.system.particles:
                    line += f"p{p.tag}   "
                f.write(line)
                self.p_p_recorded = True

    def log_p_info(self, timestep):
        """Callback function for log both particle velocities and positions when a simulation is run"""
        self.log_positions(timestep)
        self.log_velocity(timestep)

    def log_system_kinetic_energy(self, fname: str, period=100):
        """Record the kinetic energy of the system into a file that can be latter analysed"""
        hoomd.analyze.log(filename=fname,
                          quantities=["kinetic_energy"],
                          period=period,
                          overwrite=True)

    def log_stiffnesses(self, timestep):
        """Record all spring stiffnesses to a file so that they can be latter be analysed"""
        if self.p_k_recorded:
            with open("stiffnesses.txt", "a") as f:
                line = "\n" + str(timestep) + "   "
                #for p in self.system.particles:
                    #velocity = p.position
                line += str(self.k1)
                line += "   "
                line += str(self.k2)
                line += "   "
                line += str(self.k3)
                line += "   "
                f.write(line)
        else:
            with open("stiffnesses.txt", "w") as f:
                N = len(self.system.particles)
                dt = self.dt
                m = self.m
                a = self.a
                line = f"Simulation spring stiffnesses, run dt= {dt} mass= {m} a= {a} N= {N} \n"
                f.write(line)
                line = "timestep   "
                #for p in self.system.particles:
                line += "k1   k2   k3   "
                f.write(line)
                self.p_k_recorded = True

    def change_spring_eq(self, num: int, equation: "Spring_Eq"):
        """change the base spring equations into ones that will cause the spring constant to vary in time"""
        if num == 0:
            self.Eq1 = equation
        if num == 1:
            self.Eq2 = equation
        if num == 2:
            self.Eq3 = equation

    def add_dynamic_force(self, equation: "Force_Eq", tag_list: list, Negative=False):
        """Create a dynamic force and apply it to the system. The dynamic is really a static force with an
        equation relating to it so the magnitude can be updated every timestep"""

        if Negative:
            f_mag = -equation.calc_mag(0)
        else:
            f_mag = equation.calc_mag(0)

        # group together the points to which the force will be applied
        force_pt = hoomd.group.tag_list(name='group', tags=tag_list)

        # store the dynamic force equation and the group of particles it is to be applied to
        self.force_eq.append((equation, force_pt))

        # set the force
        self.dynamic_forces.append(hoomd.md.force.constant(fx=f_mag, fy=0, fz=0, group=force_pt))

        with open("force.txt", "w") as f:
            line = f"Simulation dynamic force \n"
            f.write(line)
            line = "timestep   force_mag \n"
            f.write(line)
            f.write("0    " + str(f_mag))


    def update_dynamic_force(self, time: float, index: int, record: bool):
        """Updates the dynamic force by recalculating the force magnitude and creating a new force
        with the new magnitude to replace the old one

        index = the index of the dynamic force that you want to update the magnitude of

        """

        #pull out the stored force equation
        equation = self.force_eq[index][0]

        #calculate a new force magnitude for the moment in time
        f_mag = equation.calc_mag(time)

        #pull out the stored group that the force is to be applied to
        force_pt = self.force_eq[index][1]

        #disable the previous force
        force = self.dynamic_forces[index]
        hoomd.md.force.constant.disable(force) #disable the previous force to insure it no longer effects the system

        #set a new force with the new magnitude in place of the previous one
        self.dynamic_forces[index] = hoomd.md.force.constant(fx=f_mag, fy=0, fz=0, group=force_pt)

        with open("force.txt", "a") as f:
            line = "\n" + str(time) + "    " + str(f_mag)
            f.write(line)


    def force_frequency(self, timestep):
        """Callback function for turning the static forces on a system on and off as a simulation is run"""

        if self.disabled:
            for force in self.static_forces:
                hoomd.md.force.constant.enable(force)
            self.disabled = False
        else:
            for force in self.static_forces:
                hoomd.md.force.constant.disable(force)
            self.disabled = True
        return 1

    def create_animation(self, fname: str, period=500):
        """fname is what you want the name of the animation to be. GSD file will be created in the project folder"""
        all = hoomd.group.all()
        hoomd.dump.gsd(f"{fname}.gsd", period=period, group=all, overwrite=True)

    def run(self, run_time: int, kt=.0001, dt=.0001, callback=None, callback_period=0, quiet=True, dynamic_f_stop=0, record_force=False):
        """ run_time = the number of time steps the simulation should move forward
            kt = temperature in the system (The higher the temature the more energy the system has)
            dt = change of time per time step
            callback = a function to be called periodicly during the run
            callback_period = the period at chich the callback function is triggered

        > sim.run(50000,callback=sim.force_frequency,callback_period=25000)
        """

        # set up the integrator for the run
        if self.first_run:
            nl = hoomd.md.nlist.cell()
            dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=kt, seed=1)
            dpd.pair_coeff.set('A', 'A', A=25.0, gamma=1.0)
            nl.reset_exclusions(exclusions=[])
            hoomd.md.integrate.mode_standard(dt=dt)
            if self.integrate_group is None:
                self.integrate_group = hoomd.group.all()
            integrator = hoomd.md.integrate.nve(group=self.integrate_group)
            self.first_run = False

        start = time.time()
        for i in range(run_time):
        #for i in range(1000):

            if callback is None:
                hoomd.run((run_time // 1000), quiet=quiet)
            else:
                hoomd.run((1), quiet=quiet, callback=callback, callback_period=callback_period)
                #hoomd.run((run_time // 1000), quiet=quiet, callback=callback, callback_period=callback_period)
            time_step = hoomd.get_step()

            #update the stiffnesses of the springs if they are dynamic
            if self.Eq1 is not None:
                self.k1 = self.Eq1.calc(time_step)
            if self.Eq2 is not None:
                self.k2 = self.Eq2.calc(time_step)
            if self.Eq3 is not None:
                self.k3 = self.Eq3.calc(time_step)
            self.harmonic.bond_coeff.set('polymer1', k=self.k1)
            self.harmonic.bond_coeff.set('polymer2', k=self.k2)
            self.harmonic.bond_coeff.set('polymer3', k=self.k3)

            #update the dynamic force magnitudes untill stop point
            if i <= dynamic_f_stop:
                for force_i in range (len(self.dynamic_forces)):
                    Time = time_step * dt
                    self.update_dynamic_force(Time, force_i, record_force)
            if i == dynamic_f_stop+1:
                for force in self.dynamic_forces:
                    hoomd.md.force.constant.disable(force)
                #print(time_step)

            end = time.time()
        print(f"\nsimulation took {round((end - start), 5)} secs to run")

    def run_langevin(self, run_time: int, kt=.0001, dt=.0001, callback=None, callback_period=0, quiet=True, dynamic_f_stop=0, record_force=False):
        """ run_time = the number of time steps the simulation should move forward
            kt = temperature in the system (The higher the temature the more energy the system has)
            dt = change of time per time step
            callback = a function to be called periodicly during the run
            callback_period = the period at chich the callback function is triggered

        > sim.run(50000,callback=sim.force_frequency,callback_period=25000)
        """

        # set up the integrator for the run
        if self.first_run:
            nl = hoomd.md.nlist.cell()
            dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=kt, seed=1)
            dpd.pair_coeff.set('A', 'A', A=25.0, gamma=1.0)
            nl.reset_exclusions(exclusions=[])
            hoomd.md.integrate.mode_standard(dt=dt)
            if self.integrate_group is None:
                self.integrate_group = hoomd.group.all()
            integrator = hoomd.md.integrate.langevin(group=self.integrate_group,kT=kt, seed=5)
            integrator.set_gamma('A', gamma=1)
            self.first_run = False

        start = time.time()
        for i in range(run_time):
        #for i in range(1000):

            if callback is None:
                hoomd.run((run_time // 1000), quiet=quiet)
            else:
                hoomd.run((1), quiet=quiet, callback=callback, callback_period=callback_period)
                #hoomd.run((run_time // 1000), quiet=quiet, callback=callback, callback_period=callback_period)
            time_step = hoomd.get_step()

            #update the stiffnesses of the springs if they are dynamic
            if self.Eq1 is not None:
                self.k1 = self.Eq1.calc(time_step)
            if self.Eq2 is not None:
                self.k2 = self.Eq2.calc(time_step)
            if self.Eq3 is not None:
                self.k3 = self.Eq3.calc(time_step)
            self.harmonic.bond_coeff.set('polymer1', k=self.k1)
            self.harmonic.bond_coeff.set('polymer2', k=self.k2)
            self.harmonic.bond_coeff.set('polymer3', k=self.k3)

            #update the dynamic force magnitudes untill stop point
            if i <= dynamic_f_stop:
                for force_i in range (len(self.dynamic_forces)):
                    Time = time_step * dt
                    self.update_dynamic_force(Time, force_i, record_force)
            if i == dynamic_f_stop+1:
                for force in self.dynamic_forces:
                    hoomd.md.force.constant.disable(force)
                #print(time_step)

            end = time.time()
        print(f"\nsimulation took {round((end - start), 5)} secs to run")
