"""Molecular Dynamic simulation of a line of three particles bounded by different springs

how to get to simulations in terminal get to directory where triangle_lattice.py is stored
    from spring_line import *

example:
        simulation = Line()
        simulation.create_lattice(10, 10, 10, p_mass=.001)
        simulation.apply_force([0],-10,0,0)
        simulation.log
        simulation.run(50000,0)
        simulation.graph_p_pos()
"""

from Simulations import *

class Line(Simulation):

    def create_lattice(self, k: int, num_bonds=3, add_periodic_bonds=False, p_mass=1, N=3, a=1):
        """Create a snapshot of N particles arranged in a line with bonds connecting each particle to the
        ones next to it in the line """

        self.k1 = k
        self.k2 = k
        self.k3 = k
        self.m = p_mass
        self.a = a
        self.force_counter = True

        snapshot = hoomd.data.make_snapshot(N=N,
                                            box=hoomd.data.boxdim(Lx=(N), Ly=3, Lz=3),
                                            particle_types=['A'],
                                            bond_types=['polymer1', 'polymer2', 'polymer3'])

        #set particle mass
        mass_list = [ ]
        for i in range(N):
            mass_list.append(p_mass)
        snapshot.particles.mass[:] = mass_list

        #set particle positions
        pos = -(float(N) / 2)
        pos_list = []
        for p_i in range(N):
            pos_list.append([pos,0,0])
            pos += 1
        snapshot.particles.position[:] = pos_list

        self.system = hoomd.init.read_snapshot(snapshot)

        # create lattice bonds
        counter = 1
        for p_i in range(N - 1):
            if num_bonds == 3:
                if counter == 1:
                    self.system.bonds.add('polymer1', p_i, (p_i + 1))
                elif counter == 2:
                    self.system.bonds.add('polymer2', p_i, (p_i + 1))
                elif counter == 3:
                    self.system.bonds.add('polymer3', p_i, (p_i + 1))
                    counter = 0
            elif num_bonds == 2:
                if counter == 1:
                    self.system.bonds.add('polymer1', p_i, (p_i + 1))
                elif counter == 2:
                    self.system.bonds.add('polymer2', p_i, (p_i + 1))
                    counter = 0
            elif num_bonds == 1:
                self.system.bonds.add('polymer1', p_i, (p_i + 1))
            counter += 1

        #add in periodic bond
        if add_periodic_bonds:
            if num_bonds == 3:
                if counter == 1:
                    self.system.bonds.add('polymer1', 0, (N-1))
                elif counter == 2:
                    self.system.bonds.add('polymer2', 0, (N-1))
                elif counter == 3:
                    self.system.bonds.add('polymer3', 0, (N-1))
                    counter = 0
            elif num_bonds == 2:
                if counter == 1:
                    self.system.bonds.add('polymer1', 0, (N-1))
                elif counter == 2:
                    self.system.bonds.add('polymer2', 0, (N-1))
                    counter = 0
            elif num_bonds == 1:
                self.system.bonds.add('polymer1', 0, (N-1))

        self.harmonic = hoomd.md.bond.harmonic()
        self.harmonic.bond_coeff.set('polymer1', k=k, r0=a)
        self.harmonic.bond_coeff.set('polymer2', k=k, r0=a)
        self.harmonic.bond_coeff.set('polymer3', k=k, r0=a)

    def force_call(self, timestep):
        """callback function for recording particle velocities and fluctuating the force on the system
        (specific to simulation of calculating SED with a single force acting on the system """
        if self.force_counter:
            self.disable_force(0)
            self.enable_force(1)
            self.force_counter = False
        else:
            self.disable_force(1)
            self.enable_force(0)
            self.force_counter = True
        self.log_velocities(timestep)
