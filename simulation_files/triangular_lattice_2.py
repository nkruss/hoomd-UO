#Tri_Lattice different bond patterning

from Simulations import *
import matplotlib.animation as animation
import random
import numpy as np

class Tri_Lattice_diff(Simulation):

    def __init__(self):
        Simulation.__init__(self)
        self.thirty_degree_bonds = []
        self.horizontal_bonds = []
        self.one_fifty_degree_bonds = []

    def create_lattice(self, x: int, y: int, k1:int, k2:int, k3:int, add_periodic_bonds=False, a=1.0, randomized_pos=False, two_D=True, no_boundery=False):
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

        if no_boundery:
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

        counter = 0
        A = 1
        B = 3
        C = 2
        for p_i in range(0, len(self.system.particles), 2):

            bondA = f"polymer{A}"
            bondB = f"polymer{B}"
            bondC = f"polymer{C}"

            try:
                self.system.bonds.add(bondA, self.system.particles[p_i].tag, self.system.particles[p_i + 1].tag)
            except:
                pass

            try:
                if counter == (y - 1):
                    pass
                else:
                    self.system.bonds.add(bondA, self.system.particles[p_i].tag, self.system.particles[p_i + 3].tag)
            except:
                pass

            try:
                self.system.bonds.add(bondA, self.system.particles[p_i].tag, self.system.particles[p_i + (2*y)].tag)
            except:
                pass

            try:
                self.system.bonds.add(bondB, self.system.particles[p_i].tag, self.system.particles[p_i - (2*y)].tag)
            except:
                pass

            try:
                self.system.bonds.add(bondB, self.system.particles[p_i].tag, self.system.particles[p_i + (2*y) + 1].tag)
            except:
                pass

            try:
                if counter == (y - 1):
                    pass
                else:
                    self.system.bonds.add(bondB, self.system.particles[p_i].tag, self.system.particles[p_i + (2*y) + 3].tag)
            except:
                pass

            try:
                self.system.bonds.add(bondC, self.system.particles[p_i+1].tag, self.system.particles[p_i + (2*y) + 1].tag)
            except:
                pass

            counter += 1

            if counter == y:
                counter = 0
                A += 1
                B += 1
                C += 1

                if A == 4:
                    A = 1
                if B == 4:
                    B = 1
                if C == 4:
                    C = 1

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
