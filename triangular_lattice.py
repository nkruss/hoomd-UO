"""Molecular Dynamic simulation of a triangular lattice

how to get to simulations in terminal get to directory where triangle_lattice.py is stored
    from triangle_lattice import *
example:
        simulation = Tri_Lattice()
        simulation.create_lattice(5, 5, 10, 10, 10)
        simulation.apply_force([0],10,0,0)
        simulation.print_bonds()
        simulation.graph_p_pos(graph_bonds=True)
        simulation.create_animation("dump.gsd")
        simulation.run(50000,0)
        simulation.graph_p_pos()

        #changing spring constants
        eq1 = Spring_Eq(10, 5, math.pi)
        simulation.change_spring_eq(1, eq1)
        eq2 = eq1.change_theta((math.pi / 2))
        simulation.change_spring_eq(2, eq2)
        eq3 = eq1.change_theta(0)
        simulation.change_spring_eq(3, eq3)
        simulation.run(50000,0)

        #option two for changing spring constants
        equations = Spring_Eq.get_equations(10, 5, [math.pi, (math.pi / 2), 0])
        for equation_i in range(len(equations)):
            simulation.change_spring_eq(equation_i, equations[equation_i])
        simulation.run(50000, 0)

The simulation runs with an NVE integrator
"""

from Simulations import *
import matplotlib.animation as animation
import random
import numpy as np

class Tri_Lattice(Simulation):

    def __init__(self):
        Simulation.__init__(self)
        self.thirty_degree_bonds = []
        self.horizontal_bonds = []
        self.one_fifty_degree_bonds = []

    def create_lattice(self, x: int, y: int, k1:int, k2:int, k3:int, add_periodic_bonds=False, a=1.0, randomized_pos=False, two_D=True):
        """Create a snapshot of a triangular lattice from a hexagonal unitcell, with three different types of bonds
        connecting the particles of the system

        x = number of times unitcell is replicated in x direction
        y = number of times unitcell is replicated in y direction
        k1 = initial spring constant value of the 30 degree bonds
        k2 = initial spring constant value of the 150 degree bonds
        k3 = initial spring constant value of the horizontal bonds
        set 'add_periodic_bonds' to True to bond edge particles for periodic boundary conditions
        """

        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.a = a
        self.x = x
        self.y = y

        uc = hoomd.lattice.hex(a=a)

        snap = uc.get_snapshot()

        snap.replicate(x, y, 1)
        snap.bonds.types = ["polymer1", "polymer2", "polymer3"]

        snap.box= hoomd.data.boxdim(Lx=25, Ly=25, Lz=3, dimensions=2)

        if randomized_pos:
            print(snap.box)
            pos_list = []
            x_rand = 1000000000
            y_rand = 1000000000
            for p_i in range(snap.particles.N):
                pos_x = snap.particles.position[p_i][0]
                pos_y = snap.particles.position[p_i][1]
                x_new = x_rand + pos_x
                y_new = y_rand + pos_x

                while x_new < -(x / 2) or x_new > (x / 2):
                    x_rand = random.uniform(-(self.a/4),(self.a/4))
                    x_new = x_rand + pos_x

                while y_new < -(y * math.sqrt(3) / 2) or y_new > (y * math.sqrt(3) / 2):
                    y_rand = random.uniform(-(self.a/4),(self.a/4))
                    y_new = y_rand + pos_y


                if x_new < -(x / 2) or x_new > (x / 2):
                    print("error")
                if y_new < -(y * math.sqrt(3) / 2) or y_new > (y * math.sqrt(3) / 2):
                    print("error")

                x_rand = 1000000000
                y_rand = 1000000000

                pos_list.append([x_new, y_new, 0])

            snap.particles.position[:] = pos_list

        '''pos_list = []
        for p_i in range(snap.particles.N):
            pos_x = snap.particles.position[p_i][0]
            pos_y = snap.particles.position[p_i][1]
            if p_i == 88:
                pos_y += -.3
            pos_list.append([pos_x,pos_y,0])
        snap.particles.position[:] = pos_list'''

        self.system = hoomd.init.read_snapshot(snap)

        if two_D:
            hoomd.md.update.enforce2d()

        self.harmonic = hoomd.md.bond.harmonic()
        self.harmonic.bond_coeff.set("polymer1", k=self.k1, r0=1)
        self.harmonic.bond_coeff.set("polymer2", k=self.k2, r0=1)
        self.harmonic.bond_coeff.set("polymer3", k=self.k3, r0=1)

        # add in the 30 degree bonds
        for particle_i in range(0, len(self.system.particles), 2):
            self.system.bonds.add("polymer1", self.system.particles[particle_i].tag,
                                  self.system.particles[(particle_i + 1)].tag)

            self.thirty_degree_bonds.append((self.system.particles[particle_i].tag,
                                             self.system.particles[(particle_i + 1)].tag))

        # part 2 of adding in the 30 degree bonds, counter needed to prevent possible errors from going from
        # the top of the lattice to the next column
        counter = 1
        for particle_i in range(((2 * y) + 3), len(self.system.particles), 2):
            if (counter % y) != 0:  # counter could b an issue
                self.system.bonds.add("polymer1", self.system.particles[particle_i].tag,
                                 self.system.particles[(particle_i - (2 * y) - 3)].tag)

                self.thirty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                 self.system.particles[(particle_i - (2 * y) - 3)].tag))
            counter += 1

        # part 3 of adding in the 30 degree bonds, adding in the bonds for periodic boundary conditions
        if add_periodic_bonds:
            for particle_i in range((2 * y) - 2, len(self.system.particles) - 2, (2 * y)):
                self.system.bonds.add("polymer1", self.system.particles[particle_i].tag,
                                      self.system.particles[(particle_i + 3)].tag)

                self.thirty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                 self.system.particles[(particle_i + 3)].tag))

        # add in the horizontal bonds between particles
        for particle_i in range((2 * y), len(self.system.particles)):
            self.system.bonds.add("polymer3", self.system.particles[particle_i].tag,
                                  self.system.particles[(particle_i - (2 * y))].tag)

            self.horizontal_bonds.append((self.system.particles[particle_i].tag,
                                          self.system.particles[(particle_i - (2 * y))].tag))

        # add in periodic boundary horizontal bonds
        if add_periodic_bonds:
            for particle_i in range(len(self.system.particles)):
                if particle_i < (2 * y):
                    self.system.bonds.add("polymer3", self.system.particles[particle_i].tag,
                                     self.system.particles[particle_i + (2 * x * y) - (2 * y)].tag)

                    self.horizontal_bonds.append((self.system.particles[particle_i].tag,
                                                  self.system.particles[particle_i + (2 * x * y) - (2 * y)].tag))

        # add in the 150 degree bonds
        counter2 = 1
        for particle_i in range(0, len(self.system.particles)):
            if ((counter2 % y) != 0) or ((particle_i % 2) != 0):
                # check for particles at the left edge of the lattice
                if particle_i < (2 * y):
                    if (particle_i % 2) == 0:  # check if particle tag is even
                        self.system.bonds.add("polymer2", self.system.particles[particle_i].tag,
                                         self.system.particles[(particle_i + 3)].tag)

                        self.one_fifty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                            self.system.particles[(particle_i + 3)].tag))
                else:
                    if (particle_i % 2) == 0:  # check if particle tag is even
                        self.system.bonds.add("polymer2", self.system.particles[particle_i].tag,
                                         self.system.particles[(particle_i + 3)].tag)

                        self.one_fifty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                            self.system.particles[(particle_i + 3)].tag))
                    else:
                        self.system.bonds.add("polymer2", self.system.particles[particle_i].tag,
                                         self.system.particles[(particle_i - (2 * y) - 1)].tag)

                        self.one_fifty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                            self.system.particles[(particle_i - (2 * y) - 1)].tag))
            if (particle_i % 2) != 0:
                counter2 += 1

        # part 2 for 150 degree, adding in bonds for periodic boundary condition bond
        if add_periodic_bonds:
            for particle_i in range((2 * y) - 2, len(self.system.particles), (2 * y)):
                self.system.bonds.add("polymer2", self.system.particles[particle_i].tag,
                                 self.system.particles[particle_i - ((2 * y) - 3)].tag)

                self.one_fifty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                    self.system.particles[particle_i - ((2 * y) - 3)].tag))

            for particle_i in range(1, (2 * y), (2)):
                self.system.bonds.add("polymer2", self.system.particles[particle_i].tag,
                                 self.system.particles[(particle_i + ((2 * x * y) - (2 * y) - 1))].tag)

                self.one_fifty_degree_bonds.append((self.system.particles[particle_i].tag,
                                                    self.system.particles[(particle_i + ((2 * x * y) - (2 * y) - 1))].tag))

    def print_bonds(self):
        """Print out lists of the different bonds applied though-out the lattice. Bonds are represented by a tuple of
        the two bonded particle tags
        """

        print(f"\n30 degree bonds: {self.thirty_degree_bonds} \n")
        print(f"150 degree bonds: {self.one_fifty_degree_bonds} \n")
        print("horizontal bonds:", self.horizontal_bonds, "\n")

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
                x = self.x
                y = self.y
                line = f"Simulation particle positions, run dt= {dt} mass= {m} a= {a} x= {x} y= {y} N= {N} \n"
                f.write(line)
                line = "timestep   "
                for p in self.system.particles:
                    line += f"p{p.tag}   "
                f.write(line)
                self.p_p_recorded = True

    def graph_p_pos(self, graph_bonds=False):
        """Create a graph of the system, set graph_bonds to True to show the system bonds on the graph
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
        if graph_bonds:
            for i in range(len(self.thirty_degree_bonds)):
                plt.plot([positions_x[self.thirty_degree_bonds[i][0]], positions_x[self.thirty_degree_bonds[i][1]]],
                        [positions_y[self.thirty_degree_bonds[i][0]], positions_y[self.thirty_degree_bonds[i][1]]], 'b-')
            for i in range(len(self.one_fifty_degree_bonds)):
                plt.plot([positions_x[self.one_fifty_degree_bonds[i][0]], positions_x[self.one_fifty_degree_bonds[i][1]]],
                        [positions_y[self.one_fifty_degree_bonds[i][0]], positions_y[self.one_fifty_degree_bonds[i][1]]], 'g-')
            for i in range(len(self.horizontal_bonds)):
                plt.plot([positions_x[self.horizontal_bonds[i][0]], positions_x[self.horizontal_bonds[i][1]]],
                        [positions_y[self.horizontal_bonds[i][0]], positions_y[self.horizontal_bonds[i][1]]], 'y-')
        plt.show()
