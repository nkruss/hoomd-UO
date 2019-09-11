from triangular_lattice import *
from spring_line import *
from random_line import *
from line_force import *
from two_mass_rand import *
from three_mass_rand import *
from analysis import *
from Equations import *

def tri_lattice_force_in_z():
    sim = Tri_Lattice()
    sim.create_lattice(5,5,10,10,10, add_periodic_bonds=True)
    sim.apply_static_force([10, 12, 14], 0, 0, 20)
    sim.create_animation("dump1")
    sim.run(50000,callback=sim.force_frequency,callback_period=10000)

def spring_line_force_in_x():
    sim = Line()
    sim.create_lattice(100,100,100)
    sim.apply_static_force([2], 100, 0, 0)
    #sim.pin_particles([0])
    #sim.create_animation("spring100e")
    #sim.log_velocity()
    equations = get_spring_equations(100, 5, [math.pi, (math.pi / 2), 0])
    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run(500000,callback=sim.log_p_info,callback_period=1000)

    data = Analysis()
    data.read_velocity("velocities.txt")
    data.graph_p_v(1)


def random_line(n: int, Nk=15, Nw=100, spring_const=4000):
    SED = np.zeros((Nk,Nw,n))
    cb = 100
    for i in range(n):
        print('\n','\n')
        print(i+1, ' out of ' , n, '\n')
        sim = Rand_Line()
        sim.create_lattice(spring_const)
        sim.create_animation("randline",period=cb)
        #equations = Spring_Eq.get_equations(100, 5, [math.pi, (math.pi / 2), 0])
        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run(10000, callback=sim.log_p_info, callback_period=cb)
        data = Analysis()
        data.read_velocity("velocities.txt")
        SED[:,:,i] = data.SED(Nk=Nk,Nw=Nw,cb=cb, spring_const=spring_const)
    print('The SED array has shape ' , SED.shape)
    data.Plot_SED1(SED, spring_const=spring_const)

def line_force(n: int, Nk=15, Nw=100, spring_const=4000):
    SED = np.zeros((Nk,Nw,n))
    cb = 100
    frequency = 100
    for i in range(n):
        print('\n','\n')
        print(i+1, ' out of ' , n, '\n')
        sim = Line_Force()
        sim.create_lattice(spring_const)
        sim.create_animation("randline",period=cb)
        f_equation = Force_Eq(0, 1000, frequency, 0)
        tag_list = []
        for particle_i in range(len(sim.system.particles)):
            if particle_i is 39:
                tag_list.append(sim.system.particles[particle_i].tag)

        sim.add_dynamic_force(f_equation, tag_list)

        sim.run(50000, callback=sim.log_p_info, callback_period=cb, dynamic_f_stop=5000)
        data = Analysis()
        data.read_velocity("velocities.txt")
        SED[:,:,i] = data.SED(Nk=Nk,Nw=Nw,cb=cb, spring_const=spring_const)
    print('The SED array has shape ' , SED.shape)
    data.Plot_SED1(SED, spring_const=spring_const)

def random_line_2mass(n: int, Nk=15, Nw=100, mass1=1, mass2=2, k1=4000, k2=4000):
    SED = np.zeros((Nk,Nw,n))
    cb = 100
    for i in range(n):
        print('\n','\n')
        print(i+1, ' out of ' , n, '\n')
        sim = Multi_Mass_Rand_Line()
        sim.create_lattice(k1, k2, p_mass1=mass1, p_mass2=mass2)
        sim.create_animation("randline")
        #equations = Spring_Eq.get_equations(100, 5, [math.pi, (math.pi / 2), 0])
        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run(10000, callback=sim.log_p_info, callback_period=cb)
        data = Analysis()
        data.read_velocity_mult("velocities.txt")
        SED[:,:,i] = data.SED_2mass(Nk=Nk,Nw=Nw,cb=cb,p_mass1=mass1,p_mass2=mass2,k1=k1,k2=k2)
    print('The SED array has shape ' , SED.shape)
    data.Plot_SED2(SED,p_mass1=mass1,p_mass2=mass2,k1=k1,k2=k2)

def random_line_3mass(n: int, Nk=15, Nw=100, mass1=1, mass2=2, mass3=3, k1=4000, k2=4000, k3=4000):
    SED = np.zeros((Nk,Nw,n))
    cb = 100
    for i in range(n):
        print('\n','\n')
        print(i+1, ' out of ' , n, '\n')
        sim = Three_Mass_Rand_Line()
        sim.create_lattice(k1, k2, k3, p_mass1=mass1, p_mass2=mass2, p_mass3=mass3)
        sim.create_animation("randline")
        #equations = Spring_Eq.get_equations(100, 5, [math.pi, (math.pi / 2), 0])
        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run(10000, callback=sim.log_p_info, callback_period=cb)
        data = Analysis()
        data.read_velocity_mult("velocities.txt")
        SED[:,:,i] = data.SED_3mass(Nk=Nk,Nw=Nw,cb=cb,p_mass1=mass1,p_mass2=mass2,p_mass3=mass3,k1=k1,k2=k2,k3=k3)
    print('The SED array has shape ' , SED.shape)
    data.Plot_SED3(SED,p_mass1=mass1,p_mass2=mass2,p_mass3=mass3,k1=k1, k2=k2,k3=k3)

def tri_lattice_changing_force():
    equations = get_spring_equations(1, .5, [0, (2 * math.pi / 3), (4 * math.pi / 3)])
    w = equations[0].w
    frequencies = [(w/10), (w/2), w, (2*w), (4*w)]

    for frequency in frequencies:
        sim = Tri_Lattice()
        sim.create_lattice(100,100,10,10,10, add_periodic_bonds=True)
        f_equation = Force_Eq(1, frequency, 0)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(len(sim.system.particles)):
            if particle_i < (2 * 100):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 100 * 100) - (2 * 100)].tag)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.create_animation(f"tri_force_freq {w}")
        sim.run(50000,callback=sim.log_p_info,callback_period=1000, forces=True)

def test():

    equations = get_spring_equations(1, .5, [0, (2 * math.pi / 3), (4 * math.pi / 3)])
    w = equations[0].w
    frequencies = [w, (4*w)]

    for frequency in frequencies:
        print(f"runing sim with frequency {frequency}")
        sim = Tri_Lattice()
        sim.create_lattice(100,100,10,10,10, add_periodic_bonds=True)
        f_equation = Force_Eq(1, frequency, 0)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(len(sim.system.particles)):
            if particle_i < (2 * 100):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 100 * 100) - (2 * 100)].tag)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.create_animation(f"tri_force_freq {w}")
        sim.run(50000,callback=sim.log_p_info,callback_period=1000, quiet=True, forces=True)

def test_lan():

    equations = get_spring_equations(1, .5, [0, (2 * math.pi / 3), (4 * math.pi / 3)])
    w = equations[0].w
    frequencies = w

    sim = Tri_Lattice()
    sim.create_lattice(20,20,10,10,10, add_periodic_bonds=True)
    f_equation = Force_Eq(1, w, 0)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(0,len(sim.system.particles),2):
        if particle_i < (2 * 20):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 20 * 20) - (2 * 20)].tag)

    sim.add_dynamic_force(f_equation, tag_list_left)
    sim.add_dynamic_force(f_equation, tag_list_right, Negative=True)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.create_animation(f"tri_force_freq test", period=100)
    sim.run_langevin(10000, quiet=True, forces=True)

def test_dynamic():

    sim = Line()
    sim.create_lattice(20,20,20,N=3)
    sim.create_animation("dynamictest", period=100)
    sim.dimension_constrain([1,0,0])
    f_equation = Force_Eq(0,1000,(2 * math.pi * .01),0)
    sim.add_dynamic_force(f_equation, [0])
    #sim.apply_static_force([0], 20, 0, 0)
    #equations = get_spring_equations(100, 5, [math.pi, (math.pi / 2), 0])
    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run(100000, callback=sim.log_p_info, callback_period=100, dynamic_f_stop=100000)
    data = Analysis()
    data.read_velocity("velocities.txt")
    data.read_pos("positions.txt")
    data.read_force("force.txt")
    data.graph_pos_v_force(0)

def test_damp():
    """leave test unchange proof that the damping of the langevin run works """
    sim = Rand_Line()
    sim.create_lattice(20,N=30)
    sim.create_animation("randlinetest", period=100)
    sim.dimension_constrain([1,0,0])
    #equations = get_spring_equations(100, 5, [math.pi, (math.pi / 2), 0])
    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(100000, callback=sim.log_p_info, callback_period=100)
    data = Analysis()
    data.read_velocity("velocities.txt")
    data.graph_p_v(1)


def dynamic_3dof(n: int, Nk=15, Nw=100, mass1=1, mass2=1, mass3=1, k1=4000, k2=4000, k3=4000):
    sim = Three_Mass_Rand_Line()
    SED = np.zeros((Nk,Nw,n))
    cb = 100
    k = k1
    dk = 0.5 * k
    thetas = [np.pi, np.pi/2, 0]
    equations = sim.get_equations(k,dk,thetas)
    k_init = [equations[i].calc(0) for i in range (3)]
    for i in range(n):
        print('\n','\n')
        print(i+1, ' out of ' , n, '\n')
        #sim = Three_Mass_Rand_Line()
        sim.create_lattice(k_init[0], k_init[1], k_init[2], p_mass1=mass1, p_mass2=mass2, p_mass3=mass3)
        sim.create_animation("randline",period=cb)
        for j in range (3):
            sim.change_spring_eq(j+1,equations[j])
        sim.run(50000, callback=sim.log_p_k_info, callback_period=cb)
        data = Analysis()
        data.read_velocity_mult("velocities.txt")
        #SED[:,:,i] = data.SED_3mass(Nk=Nk,Nw=Nw,cb=cb,p_mass1=mass1,p_mass2=mass2,p_mass3=mass3,k1=k1,k2=k2,k3=k3)
    print('The SED array has shape ' , SED.shape)
    #data.Plot_SED3(SED,p_mass1=mass1,p_mass2=mass2,p_mass3=mass3,k1=k1, k2=k2,k3=k3)
    data.read_stiffnesses('stiffnesses.txt')
    data.graph_stiffness(0)

def test_dynamic_spring(dynamic=False):
    """A test to confirm that making the springs dynamic actually has an effect on the system. Test contains a system
    of three particles 1 distance unit apart conected by springs with an equilibriam length of .5 of a distance unit.

    When test is run without dynamic springs the system is stationary with the force from the springs cancelling
    out with one another.

    When test is run with dynamic springs the particles in the system move back and forth untill the damping of the
    integrator takes effect.
    """

    data = Analysis()
    sim = Line()
    sim.create_lattice(1,1,1, N=3, a=.5)
    sim.dimension_constrain([1,0,0])

    equations = get_spring_equations(10, 5, 3.3, [0, (2 * math.pi / 3.0), (-4 * math.pi / 3.0)])

    if dynamic:
        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])

    sim.create_animation(f"dynamic_spring_test")
    sim.log_system_kinetic_energy(f"dynamic_spring_test.txt")
    sim.run_langevin(100000,callback=sim.log_p_info,callback_period=1000, quiet=True)
    data.read_kinetic("dynamic_spring_test.txt")
    data.graph_k_energy(0)

def kinetic_test(compare=False):
    dynamic_data = Analysis()
    equations = get_spring_equations(5000, 2500, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    dynamic_data.spring_frequency = w
    #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
    frequencies = [w/2,  w, (w*2)]
    for frequency in range(20,150,7):

        dynamic_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 500, frequency, 0)
        neg_f_equation = Force_Eq(0, 500, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        #sim.create_animation(f"kinetic {frequency}")
        sim.log_system_kinetic_energy("kinetic.txt", period=100)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)
        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        dynamic_data.read_kinetic("kinetic.txt")

    for frequency in frequencies:

        dynamic_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 500, frequency, 0)
        neg_f_equation = Force_Eq(0, 500, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        #sim.create_animation(f"kinetic {frequency}")
        sim.log_system_kinetic_energy("kinetic.txt", period=100)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)
        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        dynamic_data.read_kinetic("kinetic.txt")

    #data.graph_avg_k_energy()

    if compare:
        static_data = Analysis()
        equations = get_spring_equations(5000, 2500, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
        w = equations[0].w
        static_data.spring_frequency = w
        #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
        frequencies = [w/4, w/2,  w, (w*2), w*4]
        for frequency in range(20,150,7):

            static_data.dynamic_frequencies.append(frequency)

            sim = Tri_Lattice()
            sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)
            sim.dimension_constrain([1,1,0])

            f_equation = Force_Eq(0, 500, frequency, 0)
            neg_f_equation = Force_Eq(0, 500, frequency, math.pi)

            tag_list_left = []
            tag_list_right = []
            for particle_i in range(1,len(sim.system.particles),):
                if particle_i < (2 * 10):
                    tag_list_left.append(sim.system.particles[particle_i].tag)
                    tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

            #sim.create_animation(f"kinetic {frequency}")
            sim.log_system_kinetic_energy("kinetic.txt", period=100)

            sim.add_dynamic_force(f_equation, tag_list_left)
            sim.add_dynamic_force(neg_f_equation, tag_list_right)
            #for equation_i in range(len(equations)):
                #sim.change_spring_eq(equation_i, equations[equation_i])
            sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
            static_data.read_kinetic("kinetic.txt")

        for frequency in frequencies:

            static_data.dynamic_frequencies.append(frequency)

            sim = Tri_Lattice()
            sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

            f_equation = Force_Eq(0, 500, frequency, 0)
            neg_f_equation = Force_Eq(0, 500, frequency, math.pi)

            tag_list_left = []
            tag_list_right = []
            for particle_i in range(1,len(sim.system.particles),):
                if particle_i < (2 * 10):
                    tag_list_left.append(sim.system.particles[particle_i].tag)
                    tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

            #sim.create_animation(f"kinetic {frequency}")
            sim.log_system_kinetic_energy("kinetic.txt", period=100)

            sim.add_dynamic_force(f_equation, tag_list_left)
            sim.add_dynamic_force(neg_f_equation, tag_list_right)
            #for equation_i in range(len(equations)):
                #sim.change_spring_eq(equation_i, equations[equation_i])
            sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
            static_data.read_kinetic("kinetic.txt")

        kinetic_compare(static_data, dynamic_data)

        #data.graph_avg_k_energy()

def test_sim():
    '''data = Analysis()
    equations = get_spring_equations(5000, 2500, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w
    for i in range(2):
        data.dynamic_frequencies.append(w)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 500, w, 0)
        neg_f_equation = Force_Eq(0, 500, w, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"kinetic {i}")
        sim.log_system_kinetic_energy(f"kinetic{i}.txt")

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)
        if i == 0:
            for equation_i in range(len(equations)):
                sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(50000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=50000)
        data.read_kinetic("kinetic.txt")
    data.graph_avg_k_energy()'''

    data = Analysis()
    sim = Line()
    sim.create_lattice(5000,5000,5000, N=3)
    sim.dimension_constrain([1,0,0])

    equations = get_spring_equations(5000, 5000, [0, (2 * math.pi / 3.0), (-4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w
    data.dynamic_frequencies.append(w)

    f_equation = Force_Eq(0, 500, w, 0)
    neg_f_equation = Force_Eq(0, 500, w, math.pi)

    sim.create_animation(f"kinetic line")
    sim.log_system_kinetic_energy(f"kinetic_line.txt")

    sim.add_dynamic_force(f_equation, [0])
    sim.add_dynamic_force(neg_f_equation, [2])

    sim.run_langevin(100000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=75000)
    #data.read_kinetic("kinetic_line.txt")
    #data.graph_k_energy()

    sim.run_langevin(75000,callback=sim.log_p_info,callback_period=1000, quiet=True)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])

    sim.run_langevin(100000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=75000)
    data.read_kinetic("kinetic_line.txt")
    data.graph_k_energy(0)

def test_SED():
    SED = np.zeros((15,100,1))
    data = Analysis()
    equations = get_spring_equations(100, 50, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w

    data.dynamic_frequencies.append(w)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,100,100,100, randomized_pos=True)

    f_equation = Force_Eq(0, 50, w, 0)
    neg_f_equation = Force_Eq(0, 50, w, math.pi)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"randamized_tri_3_spring", period=100)
    sim.log_system_kinetic_energy("random_tri_3_spring.txt")

    #sim.add_dynamic_force(f_equation, tag_list_left)
    #sim.add_dynamic_force(neg_f_equation, tag_list_right)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run(100000,callback=sim.log_p_info,callback_period=100, quiet=True, dynamic_f_stop=100000)
    data.read_kinetic("random_tri_3_spring.txt")
    #data.graph_k_energy(0)
    data.read_velocity("velocities.txt")
    data.read_pos("positions.txt")

    center = 88
    p_list = []
    p_list = [center]
    p_list.append(center + 1)
    p_list.append(center + 3)
    p_list.append(center - (2 * data.y))
    p_list.append(center + (2 * data.y))
    p_list.append(center + (2 * data.y) + 3)
    p_list.append(center + (2 * data.y) + 1)

    data.graph_mult_pos_3D(p_list)
    #SED[:,:,0] = data.SED_3D(k1=10,k2=10,k3=10)
    #data.Plot_SED1(SED, spring_const=10)

def avg_kinetic():
    data = Analysis()
    equations = get_spring_equations(5000, 2500, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w

    data.dynamic_frequencies.append(w)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

    f_equation = Force_Eq(0, 500, w, 0)
    neg_f_equation = Force_Eq(0, 500, w, math.pi)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    #sim.create_animation(f"kinetic {frequency}")
    sim.log_system_kinetic_energy("kinetic.txt")

    sim.add_dynamic_force(f_equation, tag_list_left)
    sim.add_dynamic_force(neg_f_equation, tag_list_right)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
    data.read_kinetic("kinetic.txt")
    data.graph_k_energy(0)


def multi_pos_2D_test(center: int):
    data = Analysis()

    equations = get_spring_equations(5000, 500, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w

    data.dynamic_frequencies.append(w)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,5000,5000,5000)

    f_equation = Force_Eq(0, 1000, 100, 0)
    neg_f_equation = Force_Eq(0, 1000, 100, math.pi)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles)):
        if particle_i < (2 * 10) and particle_i > (10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"muli_D")
    sim.log_system_kinetic_energy("kinetic.txt")

    sim.add_dynamic_force(f_equation, tag_list_left)
    sim.add_dynamic_force(neg_f_equation, tag_list_right)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(100000,callback=sim.log_p_info,callback_period=100, quiet=True, dynamic_f_stop=100000)

    data.read_pos("positions.txt")
    p_list = []

    p_list.append(center + 3)
    p_list.append(center - (2 * data.y))
    p_list.append(center + 1)
    p_list.append(center)
    p_list.append(center + (2 * data.y) + 1)
    p_list.append(center + (2 * data.y))
    p_list.append(center + (2 * data.y) + 3)

    #data.graph_mult_pos_3D(p_list)


def compression_test():
    center = 90
    data = Analysis()

    equations = get_spring_equations(500, 500, math.sqrt(500), [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w

    data.dynamic_frequencies.append(w)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,500,500,500, a=1)

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(0,len(sim.system.particles)):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"compression_test")
    sim.log_system_kinetic_energy("kinetic.txt")

    #tag_list_left = [7,9,11,13]
    #tag_list_right = [186, 188, 190, 192]

    #sim.run(1000,callback=sim.log_p_info,callback_period=10, quiet=True, dynamic_f_stop=1000)
    sim.apply_static_force(tag_list_left, -100, 0, 0)
    sim.apply_static_force(tag_list_right, 100, 0, 0)

    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=10, quiet=True, dynamic_f_stop=150000)

    data.read_pos("positions.txt")
    p_list = [center]
    p_list.append(center + 1)
    p_list.append(center + 3)
    p_list.append(center - (2 * data.y))
    p_list.append(center + (2 * data.y))
    p_list.append(center + (2 * data.y) + 3)
    p_list.append(center + (2 * data.y) + 1)
    par_list = [7,9,11,13,187,189,191,193]


    #data.read_stiffnesses("stiffnesses.txt")
    #data.graph_stiffness(2)

    data.contraction_expansion_tri()
    #data.read_kinetic("kinetic.txt")
    #data.graph_k_energy(0)

    #data.graph_mult_pos_3D(par_list)

def static_dynamic_kinetic_compair():
    dynamic_data = Analysis()
    static_data = Analysis()
    equations = get_spring_equations(5000, 2500, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    dynamic_data.spring_frequency = w
    #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
    frequencies = [w/4, w/2,  w, (w*2), w*4]

    for w in range(20,150,7):

        dynamic_data.dynamic_frequencies.append(w)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 750, w, 0)
        neg_f_equation = Force_Eq(0, 750, w, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        #sim.create_animation(f"kinetic {frequency}")
        sim.log_system_kinetic_energy("kinetic.txt")

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        dynamic_data.read_kinetic("kinetic.txt")

    for w in frequencies:
        dynamic_data.dynamic_frequencies.append(w)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 750, w, 0)
        neg_f_equation = Force_Eq(0, 750, w, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        #sim.create_animation(f"kinetic {frequency}")
        sim.log_system_kinetic_energy("kinetic.txt")

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        dynamic_data.read_kinetic("kinetic.txt")

    for w in range(20,150,7):

        static_data.dynamic_frequencies.append(w)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 750, w, 0)
        neg_f_equation = Force_Eq(0, 750, w, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        #sim.create_animation(f"kinetic {frequency}")
        sim.log_system_kinetic_energy("kinetic.txt")

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        static_data.read_kinetic("kinetic.txt")

    print("runing static")

    for w in frequencies:

        static_data.dynamic_frequencies.append(w)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=True)

        f_equation = Force_Eq(0, 750, w, 0)
        neg_f_equation = Force_Eq(0, 750, w, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        #sim.create_animation(f"kinetic {frequency}")
        sim.log_system_kinetic_energy("kinetic.txt")

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        static_data.read_kinetic("kinetic.txt")

    kinetic_compare(static_data, dynamic_data)
