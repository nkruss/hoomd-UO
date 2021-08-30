"""
Author - Noah Kruss

File containing python class for initilizing a Molecular Dynamic simulation
of a line of n particles bounded by (num_bonds) number of different springs

Example:
        simulation = Line()
        simulation.create_lattice(10, 10, 10, p_mass=.001)
        simulation.apply_force([0],-10,0,0)
        simulation.log
        simulation.run(50000,0)
        simulation.graph_p_pos()
"""

from Simulations import *

class Line(Simulation):

    def create_lattice(self,
                        k: int,
                        num_bonds = 3,
                        add_periodic_bonds = False,
                        p_mass = 1,
                        N = 3,
                        a = 1,
                        pos_offsets = [],
                        vel_init = []):

        """
        Create a snapshot of N particles arranged in a line with bonds
        connecting each particle to the ones next to it in the line
        """

        self.m = p_mass
        self.a = a
        self.force_counter = True

        self.k_list = []
        bond_types = []
        for i in range(num_bonds):
            self.k_list.append(k)
            bond_types.append(f"polymer_{i}")

        #create a hoomd snapshot framework for the initial condition
        snapshot = hoomd.data.make_snapshot(N=N,
                                            box=hoomd.data.boxdim(Lx=(N), Ly=3, Lz=3),
                                            particle_types=['A'],
                                            bond_types=bond_types)

        #set particle mass
        mass_list = [ ]
        for i in range(N):
            mass_list.append(p_mass)
        snapshot.particles.mass[:] = mass_list

        #set particle positions
        pos = -(float(N) / 2)
        pos_list = []
        for p_i in range(N):
            pos_list.append([pos, 0, 0])
            pos += 1
        for offset_i in range(len(pos_offsets)):
            offset = pos_offsets[offset_i]
            pos_list[offset_i][0] += offset
        snapshot.particles.position[:] = pos_list

        #set initial velocities
        snapshot.particles.velocity[:] = vel_init

        #initialize created snapshot as the initial condition
        self.system = hoomd.init.read_snapshot(snapshot)

        # create lattice bonds
        for p_i in range(N - 1):
            bond_type = (p_i) % num_bonds
            self.system.bonds.add(f"polymer_{bond_type}", p_i, p_i + 1)

        #add in periodic bond
        if add_periodic_bonds:
            bond_type = (N - 1) % num_bonds
            self.system.bonds.add(f"polymer_{bond_type}", N - 1, 0)


        self.harmonic = hoomd.md.bond.harmonic()
        for i in range(num_bonds):
            self.harmonic.bond_coeff.set(f"polymer_{i}", k=k, r0=a)
