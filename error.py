import hoomd
import hoomd.md

from triangular_lattice import *
from spring_line import *
from random_line import *
from two_mass_rand import *
from analysis import *
from Equations import *
import math

def x_y_displacement():
    """why is particle 0 moving in the y dimension when constrianed to move only in x as done in line 76"""

    '''#create a callback for loging particle position
    def log_positions(timestep):
        """Record all particle positions to a file so that they can be latter be analysed"""
        if timestep > 0:
            with open("positions.txt", "a") as f:
                line = "\n" + str(timestep) + "   "
                for p in system.particles:
                    velocity = p.position
                    line += str(velocity)
                    line += "   "
                f.write(line)
        else:
            with open("positions.txt", "w") as f:
                N = len(system.particles)
                line = f"Simulation particle positions, N= {N} \n"
                f.write(line)
                line = "timestep   "
                for p in system.particles:
                    line += f"p{p.tag}   "
                f.write(line)
                p_p_recorded = True

    #initialize the system
    hoomd.context.initialize("")

    snapshot = hoomd.data.make_snapshot(N=3,
                                        box=hoomd.data.boxdim(Lx=(3), Ly=3, Lz=3),
                                        particle_types=['A'],
                                        bond_types=['polymer1', 'polymer2', 'polymer3'])
    #set particle positions
    pos = -(float(3) / 2)
    pos_list = []
    for p_i in range(3):
        pos_list.append([pos,0,0])
        pos += 1.0
    snapshot.particles.position[:] = pos_list

    system = hoomd.init.read_snapshot(snapshot)

    # create lattice bonds
    counter = 1
    for p_i in range(3 - 1):
        if counter == 1:
            system.bonds.add('polymer1', p_i, (p_i + 1))
        elif counter == 2:
            system.bonds.add('polymer2', p_i, (p_i + 1))
        elif counter == 3:
            system.bonds.add('polymer3', p_i, (p_i + 1))
            counter = 0
        counter += 1
    system.bonds.add('polymer3', 0, (3 - 1))

    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set('polymer1', k=500, r0=1)
    harmonic.bond_coeff.set('polymer2', k=500, r0=1)
    harmonic.bond_coeff.set('polymer3', k=500, r0=1)

    #constrian lattice to only move in x and y dimensions
    all = hoomd.group.all()
    hoomd.md.constrain.oneD(group=all, constraint_vector=[1,1,0])

    #create an animation file
    hoomd.dump.gsd(f"test.gsd", period=100, group=all, overwrite=True)

    # Apply a force to the system
    force_pt = hoomd.group.tag_list(name='group', tags=[0])
    hoomd.md.force.constant(fvec=(-2000,0,0), fx=-2000, fy=0, fz=0, group=force_pt)

    #constrain the leftmost particle to only move in the x dimension
    group = hoomd.group.tag_list(name="group", tags=[0])
    hoomd.md.constrain.oneD(group=group, constraint_vector=[1,0,0])

    #run the simulation
    nl = hoomd.md.nlist.cell()
    dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=.0001, seed=1)
    dpd.pair_coeff.set('A', 'A', A=25.0, gamma=1.0)
    nl.reset_exclusions(exclusions=[])
    hoomd.md.integrate.mode_standard(dt=.0001)
    integrator = hoomd.md.integrate.nve(group=all)
    hoomd.run(100000, quiet=True, callback=log_positions, callback_period=1)'''

    data = Analysis()
    sim = Line()

    sim.create_lattice(500, 500, 500)
    sim.dimension_constrain([1,1,0])

    sim.create_animation(f"error_test", period=100)
    sim.log_system_kinetic_energy("kinetic.txt")

    sim.run(1000,callback=sim.log_p_info,callback_period=1, quiet=True, dynamic_f_stop=1000)
    sim.apply_static_force([0], -5000, 0, 0)
    #sim.apply_static_force(tag_list_right, 500, 0, 0)
    sim.particle_dimension_constrian([1,0,0], [0])

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run(100000,callback=sim.log_p_info,callback_period=10, quiet=True, dynamic_f_stop=100000)

    data.read_pos("positions.txt")
    data.graph_mult_pos_3D([0,1,2])


def force_error():
    """why is particle 0 moving in the y dimension when constrianed to move only in x as done in line __ , also why
    are the other particles moving?"""

    #create a callback for loging particle position
    def log_positions(timestep):
        '''Record all particle positions to a file so that they can be latter be analysed'''
        if timestep > 0:
            with open("positions.txt", "a") as f:
                line = "\n" + str(timestep) + "   "
                for p in system.particles:
                    velocity = p.position
                    line += str(velocity)
                    line += "   "
                f.write(line)
        else:
            with open("positions.txt", "w") as f:
                N = len(system.particles)
                line = f"Simulation particle positions, N= {N} \n"
                f.write(line)
                line = "timestep   "
                for p in system.particles:
                    line += f"p{p.tag}   "
                f.write(line)
                p_p_recorded = True

    data = Analysis()

    #initialize the system
    hoomd.context.initialize("")


    snapshot = hoomd.data.make_snapshot(N=1,
                                        box=hoomd.data.boxdim(Lx=(3), Ly=3, Lz=3, dimensions=2),
                                        particle_types=['A'],
                                        bond_types=['polymer1', 'polymer2', 'polymer3'])
    #set particle positions
    pos = -(float(1) / 2)
    pos_list = []
    for p_i in range(1):
        pos_list.append([pos,0,0])
        pos += 1.0
    snapshot.particles.position[:] = pos_list


    system = hoomd.init.read_snapshot(snapshot)
    '''
    # create lattice bonds
    counter = 1
    for p_i in range(3 - 1):
        if counter == 1:
            system.bonds.add('polymer1', p_i, (p_i + 1))
        elif counter == 2:
            system.bonds.add('polymer2', p_i, (p_i + 1))
        elif counter == 3:
            system.bonds.add('polymer3', p_i, (p_i + 1))
            counter = 0
        counter += 1
    system.bonds.add('polymer3', 0, (3 - 1))'''


    #harmonic = hoomd.md.bond.harmonic()
    #harmonic.bond_coeff.set('polymer1', k=500, r0=1)
    #harmonic.bond_coeff.set('polymer2', k=500, r0=1)
    #harmonic.bond_coeff.set('polymer3', k=500, r0=1)

    #constrian lattice to only move in x and y dimensions
    all = hoomd.group.all()
    #hoomd.md.constrain.oneD(group=all, constraint_vector=[1,0,0])
    #hoomd.md.constrain.oneD(group=all, constraint_vector=[0,1,0])

    hoomd.md.update.enforce2d()

    #create an animation file
    hoomd.dump.gsd(f"test.gsd", period=100, group=all, overwrite=True)

    # Apply a force to the system
    force_pt = hoomd.group.tag_list(name='group', tags=[0])
    hoomd.md.force.constant(fvec=(100,0,0), fx=100, fy=0, fz=0, group=force_pt)

    #constrain the leftmost particle to only move in the x dimension
    group = hoomd.group.tag_list(name="group", tags=[0])
    #hoomd.md.constrain.oneD(group=group, constraint_vector=[1,0,0])

    #run the simulation
    nl = hoomd.md.nlist.cell()
    dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=.0001, seed=1)
    dpd.pair_coeff.set('A', 'A', A=25.0, gamma=1.0)
    nl.reset_exclusions(exclusions=[])
    hoomd.md.integrate.mode_standard(dt=.0001)
    integrator = hoomd.md.integrate.nve(group=all)
    hoomd.run(100000, quiet=True, callback=log_positions, callback_period=1)

    data.read_pos("positions.txt")
    data.graph_mult_pos_3D([0])

    """
    data = Analysis()
    sim = Line()

    sim.create_lattice(500, 500, 500, N=1)
    sim.dimension_constrain([1,1,0])

    sim.create_animation(f"error_test", period=100)
    sim.log_system_kinetic_energy("kinetic.txt")

    sim.run(1000,callback=sim.log_p_info,callback_period=1, quiet=True, dynamic_f_stop=1000)
    sim.apply_static_force([0], -5000, 0, 0)
    #sim.apply_static_force(tag_list_right, 500, 0, 0)
    sim.particle_dimension_constrian([1,0,0], [0])

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run(100000,callback=sim.log_p_info,callback_period=10, quiet=True, dynamic_f_stop=100000)

    data.read_pos("positions.txt")
    data.graph_mult_pos_3D([0,1,2])"""
