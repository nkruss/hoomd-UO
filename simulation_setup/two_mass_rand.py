from Simulations import *
import random
import numpy as np

class Multi_Mass_Rand_Line(Simulation):

    def create_lattice(self, k1: int, k2: int, p_mass1: int, p_mass2: int, N=78, a=1):

        self.k1 = k1
        self.k2 = k2
        #self.k3 = k
        self.m = [p_mass1,p_mass2]
        self.a = a

        snapshot = hoomd.data.make_snapshot(N=N,
                                            box=hoomd.data.boxdim(Lx=N, Ly=3, Lz=3),
                                            particle_types=['A'],
                                            bond_types=['polymer1', 'polymer2', 'polymer3'])
        #set particle mass
        mass_list = [ ]
        for i in range(N):
            if i % 2 == 0:
              mass_list.append(p_mass1)
            else:
              mass_list.append(p_mass2)
        snapshot.particles.mass[:] = mass_list

        #set randamized particle positions
        pos = -(float(N) / 2)
        pos_list = []
        for p_i in range(N):
            if p_i == 0 or p_i == (N-1):
                pos_list.append([pos,0,0])
            else:
                    #rand = random.uniform(-(a/2),(a/2))
                    rand = 0.5*np.sin(1*p_i*a) # use this to probe system response to sinusoidal initial displacements

                    pos_list.append([rand+pos,0,0])
            pos += 1
        snapshot.particles.position[:] = pos_list

        self.system = hoomd.init.read_snapshot(snapshot)

        # create lattice bonds
        counter = 1
        for p_i in range(N - 1):
            if counter == 1:
                self.system.bonds.add('polymer1', p_i, (p_i + 1))
            elif counter == 2:
                self.system.bonds.add('polymer2', p_i, (p_i + 1))
            #elif counter == 3:
            #    self.system.bonds.add('polymer3', p_i, (p_i + 1))
                counter = 0
            counter += 1
        self.system.bonds.add('polymer2', 0, (N - 1))


        self.harmonic = hoomd.md.bond.harmonic()
        self.harmonic.bond_coeff.set('polymer1', k=self.k1, r0=a)
        self.harmonic.bond_coeff.set('polymer2', k=self.k2, r0=a)
        self.harmonic.bond_coeff.set('polymer3', k=self.k1, r0=a)
