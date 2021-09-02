#---------------IMPORT STATEMENTS------------------------------
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import shutil
from analysis import *

sys.path.append('./../simulation_files')
from Equations import *
from triangular_lattice import *
from spring_line import *
from random_line import *
from line_force import *
from two_mass_rand import *
from three_mass_rand import *
from triangular_lattice_2 import *
sys.path.append('./../sandbox')

#---------------Function to Run Simulation---------------------
def tri_lattice_force_in_y():
    """
    Function to run a simulation of a triangular lattice of 3 springs with
    opposing forces in the center of the lattice pushing out in the y dimension
    """

    dim_lattice = 10
    sim = Tri_Lattice()
    sim.create_lattice(dim_lattice,dim_lattice,1,1,1, add_periodic_bonds=False)

    #pinch particles in the middle of the lattice
    f_equation = Force_Eq(0, 1, .75, 0)
    neg_f_equation = Force_Eq(0, 1, .75, math.pi)
    # sim.add_dynamic_force(f_equation, [86], dimension="y")
    # sim.add_dynamic_force(neg_f_equation, [88], dimension="y")

    #sim.apply_static_force([90], 0, 10, 0)
    sim.create_animation("dump1_dynamic")
    sim.log_system_kinetic_energy("kinetic.txt")

    # #code for identifing the indexes of the edge particles of the lattice
    # edge_particle_list = []
    # for i in range(1, (2 * dim_lattice), 2):
    #     edge_particle_list.append(i)
    # for i in range((2 * dim_lattice * dim_lattice) - 2, (2 * dim_lattice * dim_lattice) - (2 * dim_lattice), -2):
    #     edge_particle_list.append(i)
    #
    # for particle_i in range((2 * dim_lattice) - 2, (2 * dim_lattice * dim_lattice), (2 * dim_lattice)):
    #     edge_particle_list.append(particle_i)
    # for particle_i in range(0, (2 * dim_lattice * dim_lattice) - 2, (2 * dim_lattice)):
    #     edge_particle_list.append(particle_i)
    #
    # print(edge_particle_list)

    equations = get_spring_equations(1, .5, .75, [(2 * math.pi / 3.0), (4 * math.pi / 3.0), 0])
    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])

    #sim.run_langevin(1000000,callback=sim.log_p_info,gamma=1.0, callback_period=500, quiet=True)
    sim.run_langevin(1000000,callback=sim.log_p_info, gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=1000000)

    data = Analysis()
    data.read_pos("positions.txt")
    data.read_velocity("velocities.txt")
    data.read_kinetic("kinetic.txt")
    title = "test dynamic, runtime=1000000"
    data.heat_map(title)


def random_line(n: int, Nk=15, Nw=100, spring_const=4000):
    """
    Function for project with Nathan Villager. Function creates a simulation of
    a 1D line with the particles offset by their equilibrium positions by
    """
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

def only_data(force_mag: int):

    dynamic_position_file = "positions"
    dynamic_velocity_file = "velocities"

    static_position_file = "positions"
    static_velocity_file = "velocities"

    if force_mag == .1:
        dynamic_sim_prop = "_dynamic_point1.txt"
        dynamic_title = "Dynamic bond trial, runtime=1000000, force mag = .1"

        static_sim_prop = "_static_point1.txt"
        static_title = "Static bond trial, runtime=1000000, force mag = .1"
        compair_title = "Kinetic energy percent difference between static and dynamic systems \nforce mag = .1"

    elif force_mag == 1:
        dynamic_sim_prop = "_dynamic_1.txt"
        dynamic_title = "Dynamic bond trial, runtime=1000000, force mag = 1"

        static_sim_prop = "_static_1.txt"
        static_title = "Static bond trial, runtime=1000000, force mag = 1"

        compair_title = "Kinetic energy percent difference between static and dynamic systems \nforce mag = 1"

    elif force_mag == 0:
        dynamic_sim_prop = "_dynamic_0.txt"
        dynamic_title = "Dynamic bond trial, runtime=1000000, force mag = 0"

        static_sim_prop = "_static_0.txt"
        static_title = "Static bond trial, runtime=1000000, force mag = 0"
        compair_title = "Kinetic energy percent difference between static and dynamic systems \nforce mag = 0"

        # dynamic_sim_prop = "_rotated_dynamic_0.txt"
        # dynamic_title = "Dynamic bond trial, runtime=1000000, force mag = 0"
        #
        # static_sim_prop = "_static_0.txt"
        # static_title = "Static bond trial, runtime=1000000, force mag = 0"
        # compair_title = "Kinetic energy percent difference between static and rotated dynamic systems \nforce mag = 0"

    dynamic_position_file += dynamic_sim_prop
    dynamic_velocity_file += dynamic_sim_prop

    static_position_file += static_sim_prop
    static_velocity_file += static_sim_prop


    static_data = Analysis()
    static_data.read_pos(static_position_file)
    static_data.read_velocity(static_velocity_file)
    static_data.read_kinetic("kinetic.txt")
    static_energies = static_data.heat_map(static_title)

    dynamic_data = Analysis()
    dynamic_data.read_pos(dynamic_position_file)
    dynamic_data.read_velocity(dynamic_velocity_file)
    dynamic_data.read_kinetic("kinetic.txt")
    dynamic_energies = dynamic_data.heat_map(dynamic_title)

    #compairison of static and dynamic trial
    E_k_diffs = []
    for energy_i in range(len(static_energies[2])):
        #diff = abs(dynamic_energies[2][energy_i] - static_energies[2][energy_i])
        diff = (dynamic_energies[2][energy_i] - static_energies[2][energy_i]) / static_energies[2][energy_i] * 100
        E_k_diffs.append(diff)

    H, xe, ye = np.histogram2d(static_energies[0], static_energies[1], [75,75], weights=E_k_diffs)
    H = H.T
    ax1 = plt.gca()
    plot = ax1.imshow(H, interpolation='nearest', origin='low', extent=[xe[0]*2.5, xe[-1]*2.5, ye[0], ye[-1]], cmap='Blues')
    ax1.set_xticks([])
    ax1.set_yticks([])

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    ax1.figure.colorbar(plot, cax=cax)

    plt.suptitle(compair_title)
    plt.show()

def test_width(dynamic=False, lattice_width=10, gamma=1.0):

    equations = get_spring_equations(1, .5, 1, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    frequencies = [w]

    data = Analysis()

    for frequency in frequencies:
        print(f"runing sim with frequency {frequency}")

        sim = Tri_Lattice()
        sim.create_lattice(lattice_width, lattice_width,1,1,1, add_periodic_bonds_in_y=True)
        f_equation = Force_Eq(0, .1, frequency, 0)
        neg_f_equation = Force_Eq(0, .1, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(len(sim.system.particles)):
            if particle_i < (2 * lattice_width):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i +
                                      (2 * lattice_width * lattice_width) -
                                      (2 * lattice_width)].tag
                                      )

        # sim.add_dynamic_force(f_equation, tag_list_left)
        # sim.add_dynamic_force(f_equation, tag_list_right)

        sim.apply_static_force(tag_list_left, -.1, 0, 0)
        sim.apply_static_force(tag_list_right, .1, 0, 0)
        sim.log_system_kinetic_energy("kinetic.txt", period=500)


        if dynamic:
            for equation_i in range(len(equations)):
                sim.change_spring_eq(equation_i, equations[equation_i])
            sim.create_animation(f"width_strech_test_dynamic")
        else:
            sim.create_animation(f"width_strech_test_static")
        sim.create_animation(f"tri_force_freq {w}")
        sim.run_langevin(3000000,callback=sim.log_p_info,gamma=gamma, callback_period=500, quiet=True, dynamic_f_stop=0)

        data.read_pos("positions.txt")
        data.read_kinetic("kinetic.txt")
        #data.graph_k_energy(0)
        if dynamic:
            data.width_of_tri(0, "dynamic", gamma, w)
        else:
            data.width_of_tri(0, "static", gamma, w)

def test_ploting():
    width_analysis("static_dynamic_compair")

def static_dynamic_compair():
    #create folders to store graphs and gsd files
    path = os.getcwd() + "/animations1"
    bin_path = os.getcwd() + "/bin_graphs1"
    sim_data_path = os.getcwd() + "/sim_data"
    try:
        shutil.rmtree(path)
        shutil.rmtree(bin_path)
        shutil.rmtree(sim_data_path)
    except:
        pass
    os.mkdir(path, 0o755)
    os.mkdir(bin_path, 0o755)
    os.mkdir(sim_data_path, 0o755)

    #add a file with the simulations run conditions
    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 1 ± .5, dynamic force of 0 ± .1 \nsimulation runtime of 3000000, data_recorded every 500 timesteps, periodic bonds in the y direction \ngamma =.05")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", "bin_graphs1")
    shutil.copy("run_conditions.txt", "animations1")
    shutil.copy("run_conditions.txt", "sim_data")


    #get a list of different spring frequencies to test over
    # spring_frequencies = [0.0, 1]
    # force_frequencies = [0.0, 1]
    spring_frequencies = [1, 0.0]
    force_frequencies = [1, 0.0]
    w = .01
    while(w<=4):
        spring_frequencies.append(w)
        force_frequencies.append(w)
        w = w*4

    results = []
    kinetic_differences = []
    for spring_freq in spring_frequencies:

        equations = get_spring_equations(1, .5, spring_freq, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])

        for force_freq in force_frequencies:

            dynamic_data = Analysis()
            static_data = Analysis()
            print("\n-------------------------------------------------")
            print(f"running force_freq{force_freq}, spring_freq{spring_freq} \n")

            ##dynamic_data.spring_frequency = spring_freq

            #deal with static force
            if force_freq == 0.0:

            #dynamic_data.spring_frequency = spring_freq

                ##static_force dynamic bonds

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

                tag_list_left = []
                tag_list_right = []
                for particle_i in range(1,len(sim.system.particles),):
                    if particle_i < (2 * 10):
                        tag_list_left.append(sim.system.particles[particle_i].tag)
                        tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

                sim.create_animation(f"dynamic_bonds_static_force")

                sim.log_system_kinetic_energy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", period=500)

                sim.apply_static_force(tag_list_left, -.1, 0, 0)
                sim.apply_static_force(tag_list_right, .1, 0, 0)

                for equation_i in range(len(equations)):
                    sim.change_spring_eq(equation_i, equations[equation_i])
                #sim.run_langevin(300000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=300000)
                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                #dynamic_data.read_kinetic("kinetic.txt")

                # data.read_pos("positions.txt")
                # data.width_of_tri(force_freq, "dynamic", .05, spring_freq)
                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data")
                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data")


                shutil.copy("dynamic_bonds_static_force.gsd", "animations1")

                # ------------------------------------------

                ##static_force static bonds

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

                tag_list_left = []
                tag_list_right = []
                for particle_i in range(1,len(sim.system.particles),):
                    if particle_i < (2 * 10):
                        tag_list_left.append(sim.system.particles[particle_i].tag)
                        tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

                sim.create_animation(f"static_bonds_static_force")

                sim.log_system_kinetic_energy(f"kinetic_static_{spring_freq}_{force_freq}.txt", period=500)

                sim.apply_static_force(tag_list_left, -.1, 0, 0)
                sim.apply_static_force(tag_list_right, .1, 0, 0)

                # for equation_i in range(len(equations)):
                #     sim.change_spring_eq(equation_i, equations[equation_i])
                #sim.run_langevin(300000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=300000)
                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                #static_data.read_kinetic("kinetic.txt")

                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data")

                shutil.copy("static_bonds_static_force.gsd", "animations1")

            else:

                f_equation = Force_Eq(0, .1, force_freq, 0)
                neg_f_equation = Force_Eq(0, .1, force_freq, math.pi)

                #dynamic test

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

                tag_list_left = []
                tag_list_right = []
                for particle_i in range(1,len(sim.system.particles),):
                    if particle_i < (2 * 10):
                        tag_list_left.append(sim.system.particles[particle_i].tag)
                        tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

                sim.create_animation(f"dynamic_bonds_force_at_frequency{force_freq}")

                sim.log_system_kinetic_energy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", period=500)

                sim.add_dynamic_force(f_equation, tag_list_left)
                sim.add_dynamic_force(neg_f_equation, tag_list_right)

                for equation_i in range(len(equations)):
                    sim.change_spring_eq(equation_i, equations[equation_i])
                #sim.run_langevin(300000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=300000)
                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)
                #dynamic_data.read_kinetic("kinetic.txt")
                # dynamic_data.read_velocity("velocities.txt")
                # dynamic_data.my_kinetic()
                # dynamic_data.kinetic_bin()

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data")
                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data")

                shutil.copy(f"dynamic_bonds_force_at_frequency{force_freq}.gsd", "animations1")

                #------------------------------------------------

                #static test

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

                tag_list_left = []
                tag_list_right = []
                for particle_i in range(1,len(sim.system.particles),):
                    if particle_i < (2 * 10):
                        tag_list_left.append(sim.system.particles[particle_i].tag)
                        tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

                sim.create_animation(f"static_bonds_force_at_frequency{force_freq}")

                sim.log_system_kinetic_energy(f"kinetic_static_{spring_freq}_{force_freq}.txt", period=500)

                sim.add_dynamic_force(f_equation, tag_list_left)
                sim.add_dynamic_force(neg_f_equation, tag_list_right)

                # for equation_i in range(len(equations)):
                #     sim.change_spring_eq(equation_i, equations[equation_i])
                #sim.run_langevin(300000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=300000)
                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)
                # static_data.read_kinetic("kinetic.txt")
                # static_data.read_velocity("velocities.txt")
                # static_data.my_kinetic()

                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data")

                shutil.copy(f"static_bonds_force_at_frequency{force_freq}.gsd", "animations1")

            #----------
            #create a tuple of the result infomation
            # static_avg = static_data.average_kinetic[-1]
            # print(f"static_avg ={static_avg}")
            # dynamic_avg = dynamic_data.average_kinetic[-1]
            # print(f"dynamic_avg = {dynamic_avg}")
            # ratio = dynamic_avg / static_avg
            #
            # results.append((spring_freq, force_freq, ratio))
            # kinetic_differences.append(ratio)

    print("saving results")
    np.savetxt("static_dynamic_compair.npy",results,header="tuples of data points in the format (spring_freq, force_freq, ratio)")
    ## need to graph results. Treat each result tuple as a three dimensional point
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    #
    # ax.scatter(spring_frequencies, force_frequencies, kinetic_differences, c='r', marker='o')
    #
    # ax.set_xlabel('Spring frequency')
    # ax.set_ylabel('Force frequency')
    # ax.set_zlabel('Kinetic energy ratio')
    #
    # plt.savefig(f"static_dynamic_compair.pdf")
    # plt.clf()

def compairison_anaylisis():
    data_files_locations = os.getcwd() + "/sim_data"

    spring_frequencies = []
    force_frequencies = []
    dynamic_files = []
    static_files = []
    results = []

    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)
        if filename.startswith("kinetic_dynamic"):
             # print(os.path.join(directory, filename))
             fname_parts = filename.split("_")
             spring_freq = float(fname_parts[2])
             force_freq = fname_parts[3].split('t')
             force_freq = float(force_freq[0][:-1])

             if spring_freq not in spring_frequencies:
                 spring_frequencies.append(spring_freq)
             if force_freq not in force_frequencies:
                 force_frequencies.append(force_freq)

             dynamic_files.append(filename)

             continue
        elif filename.startswith("kinetic_static"):
            static_files.append(filename)
        else:
            continue

    #no freqs less then .5
    # for file in sorted(os.listdir(data_files_locations)):
    #     filename = os.fsdecode(file)
    #     #print(filename)
    #     if filename.startswith("kinetic_dynamic"):
    #          # print(os.path.join(directory, filename))
    #          fname_parts = filename.split("_")
    #          spring_freq = float(fname_parts[2])
    #          force_freq = fname_parts[3].split('t')
    #          force_freq = float(force_freq[0][:-1])
    #          if force_freq < .5 or spring_freq < .5:
    #              pass
    #          else:
    #              if spring_freq not in spring_frequencies:
    #                  spring_frequencies.append(spring_freq)
    #              if force_freq not in force_frequencies:
    #                  force_frequencies.append(force_freq)
    #
    #              dynamic_files.append(filename)
    #
    #          continue
    #     elif filename.startswith("kinetic_static"):
    #         #static_files.append(filename)
    #         # print(os.path.join(directory, filename))
    #         fname_parts = filename.split("_")
    #         spring_freq = float(fname_parts[2])
    #         force_freq = fname_parts[3].split('t')
    #         force_freq = float(force_freq[0][:-1])
    #         if force_freq < .5 or spring_freq < .5:
    #             pass
    #         else:
    #             static_files.append(filename)
    #     else:
    #         continue

    print(spring_frequencies)
    print(force_frequencies)
    #print(dynamic_files)
    #print(static_files)

    for i in range(len(dynamic_files)):
        print()
        static_data = Analysis()
        dynamic_data = Analysis()
        print(static_files[i])
        print(dynamic_files[i])

        fname_parts = static_files[i].split("_")
        spring_freq = float(fname_parts[2])
        force_freq = fname_parts[3].split('t')
        force_freq = float(force_freq[0][:-1])

        #read corresponding kinetic energy files
        static_data.read_kinetic(data_files_locations + '/' + static_files[i])
        dynamic_data.read_kinetic(data_files_locations + '/' + dynamic_files[i])

        static_avg = static_data.average_kinetic[-1]
        print(f"static_avg ={static_avg}")
        dynamic_avg = dynamic_data.average_kinetic[-1]
        print(f"dynamic_avg = {dynamic_avg}")
        #(dyn-statoc) / static

        diff = abs(dynamic_avg - static_avg)
        ratio = abs((dynamic_avg - static_avg)) / static_avg

        #ratio
        results.append((spring_freq, force_freq, ratio))

        #absolute difference
        # results.append((spring_freq, force_freq, diff))

    print("\n saving results")
    np.savetxt("static_dynamic_compair.npy",results,header="tuples of data points in the format (spring_freq, force_freq, ratio)")


    static_dynamic_analysis()

def test_lan():

    equations = get_spring_equations(1, .5, [0, (2 * math.pi / 3), (4 * math.pi / 3)])
    w = equations[0].w
    frequencies = w

    sim = Tri_Lattice()
    sim.create_lattice(10,10,1,1,1, add_periodic_bonds=True)
    f_equation = Force_Eq(.1, w, 0)

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

def test_dynamic_spring(dynamic=True):
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
    equations = get_spring_equations(5000, 2500, math.sqrt(5000), [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = equations[0].w
    data.spring_frequency = w

    data.dynamic_frequencies.append(w)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,5000,5000,5000)

    f_equation = Force_Eq(0, 500, w, 0)
    neg_f_equation = Force_Eq(0, 500, w, math.pi)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"kinetic")
    sim.log_system_kinetic_energy("kinetic.txt")

    sim.add_dynamic_force(f_equation, tag_list_left)
    sim.add_dynamic_force(neg_f_equation, tag_list_right)

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])

    #sim.cite()
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
    data.read_kinetic("kinetic.txt")
    data.graph_k_energy(0)

    data2 = Analysis()
    """leave test unchange proof that the damping of the langevin run works """
    sim = Rand_Line()
    sim.create_lattice(20,N=30)
    sim.create_animation("randlinetest", period=100)
    sim.dimension_constrain([1,0,0])
    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.log_system_kinetic_energy("kinetic.txt")
    sim.run_langevin(100000, callback=sim.log_p_info, callback_period=100)
    data2.read_kinetic("kinetic.txt")
    data2.graph_k_energy(0)

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

    sim = Tri_Lattice_diff()
    sim.create_lattice(10,10,500,500,500, a=1, no_boundery=True)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(0,len(sim.system.particles)):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"compression_test_new")
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

def static_force_test():
    sim = Line()
    sim.create_lattice(1,1,1,N=3)
    sim.create_animation("static_force_test", period=100)
    sim.dimension_constrain([1,0,0])
    sim.apply_static_force([0], -5, 0, 0)
    sim.apply_static_force([2], 5, 0, 0)

    equations = get_spring_equations(1,.5, math.sqrt(1), [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    # for equation_i in range(len(equations)):
    #     sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000, callback=sim.log_p_info, callback_period=100, dynamic_f_stop=100000)
    # data = Analysis()
    # data.read_velocity("velocities.txt")
    # data.read_pos("positions.txt")

def static_dynamic_kinetic_compair(w_new: float):
    #create folders to store graphs and gsd files
    path = os.getcwd() + "/animations1"
    bin_path = os.getcwd() + "/bin_graphs1"
    try:
        shutil.rmtree(path)
        shutil.rmtree(bin_path)
    except:
        pass
    os.mkdir(path, 0o755)
    os.mkdir(bin_path, 0o755)

    #add a file with the simulations run conditions
    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 5000 ± 2500, dynamic force of 0 ± 250, simulation runtime of 150000, spring ocilation frquency of {w_new}")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", "bin_graphs1")
    shutil.copy("run_conditions.txt", "animations1")

    #run simulation with static springs
    dynamic_data = Analysis()
    static_data = Analysis()

    equations = get_spring_equations(1, .5, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(1)
    dynamic_data.spring_frequency = w

    ##static_force dynamic bonds
    dynamic_data.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"dynamic_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic.txt", period=500)

    sim.apply_static_force(tag_list_left, -.1, 0, 0)
    sim.apply_static_force(tag_list_right, .1, 0, 0)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(200000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    dynamic_data.read_kinetic("kinetic.txt")
    dynamic_data.read_pos("positions.txt")
    #dynamic_data.width_of_tri(0, "dynamic")

    ##static_force dynamic bonds
    static_data.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"static_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic.txt", period=500)

    sim.apply_static_force(tag_list_left, -.1, 0, 0)
    sim.apply_static_force(tag_list_right, .1, 0, 0)

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    static_data.read_kinetic("kinetic.txt")
    static_data.read_pos("positions.txt")
    static_data.width_of_tri(0, "static")

    shutil.copy("dynamic_bonds_force_at_0.gsd", "animations1")
    shutil.copy("static_bonds_force_at_0.gsd", "animations1")


    #---------------

    #repeat simulation with dynamic force at different frequencies

    frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200]
    for frequency in frequencies:
        print(frequency)

        dynamic_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, .1, frequency, 0)
        neg_f_equation = Force_Eq(0, .1, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        dynamic_data.read_kinetic("kinetic.txt")
        dynamic_data.read_velocity("velocities.txt")
        dynamic_data.read_pos("positions.txt")
        dynamic_data.my_kinetic()
        dynamic_data.kinetic_bin()
        dynamic_data.width_of_tri(frequency, "dynamic")

        shutil.copy(f"dynamic_bonds_force_at_{frequency}.gsd", "animations1")
        #dynamic_data.graph_k_energy()

        #static
        static_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,1,1,1, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, .1, frequency, 0)
        neg_f_equation = Force_Eq(0, .1, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"static_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        static_data.read_kinetic("kinetic.txt")
        static_data.read_velocity("velocities.txt")
        static_data.read_pos("positions.txt")
        static_data.my_kinetic()
        static_data.kinetic_bin()
        static_data.width_of_tri(frequency, "static")

        shutil.copy(f"static_bonds_force_at_{frequency}.gsd", "animations1")

        kinetic_bin_compare(static_data, dynamic_data, frequency, save_location="bin_graphs1")

    kinetic_compare(static_data, dynamic_data)
    print(w_new)

def static_dynamic_kinetic_compair_2(w_new: float):
    path = os.getcwd() + "/animations2"
    bin_path = os.getcwd() + "/bin_graphs2"
    try:
        shutil.rmtree(path)
        shutil.rmtree(bin_path)
    except:
        pass
    os.mkdir(path, 0o755)
    os.mkdir(bin_path, 0o755)

    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 5000 ± 2500, dynamic force of 0 ± 250, simulation runtime of 150000, spring ocilation frquency of {w_new}")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", "bin_graphs2")
    shutil.copy("run_conditions.txt", "animations2")


    dynamic_data2 = Analysis()
    static_data2 = Analysis()

    equations = get_spring_equations(5000, 2500, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(5000)
    dynamic_data2.spring_frequency = w
    frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200]
    #frequencies = [w/4, w/2,  w, (w*2), w*4]
    #frequencies = [110]

    #static_force dynamic bonds
    dynamic_data2.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.velocity_fname ="velocities2.txt"
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"dynamic_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic2.txt", period=500)

    sim.apply_static_force(tag_list_left, 750, 0, 0)
    sim.apply_static_force(tag_list_right, -750, 0, 0)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    dynamic_data2.read_kinetic("kinetic2.txt")

    #static_force dynamic bonds
    static_data2.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"static_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic2.txt", period=500)

    sim.apply_static_force(tag_list_left, 750, 0, 0)
    sim.apply_static_force(tag_list_right, -750, 0, 0)

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    static_data2.read_kinetic("kinetic2.txt")

    shutil.copy("dynamic_bonds_force_at_0.gsd", "animations2")
    shutil.copy("static_bonds_force_at_0.gsd", "animations2")

    #---------------

    for frequency in frequencies:
        print(frequency)

        dynamic_data2.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.velocity_fname ="velocities2.txt"
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 750, frequency, 0)
        neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic2.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        dynamic_data2.read_kinetic("kinetic2.txt")
        dynamic_data2.read_velocity("velocities2.txt")
        dynamic_data2.my_kinetic()
        dynamic_data2.kinetic_bin()

        shutil.copy(f"dynamic_bonds_force_at_{frequency}.gsd", "animations2")
        #dynamic_data.graph_k_energy()

        #static
        static_data2.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.velocity_fname ="velocities2.txt"
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 750, frequency, 0)
        neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"static_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic2.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        static_data2.read_kinetic("kinetic2.txt")
        static_data2.read_velocity("velocities2.txt")
        static_data2.my_kinetic()
        static_data2.kinetic_bin()

        shutil.copy(f"static_bonds_force_at_{frequency}.gsd", "animations2")

        kinetic_bin_compare(static_data2, dynamic_data2, frequency, save_location="bin_graphs2")

    kinetic_compare(static_data2, dynamic_data2)
    print(w_new)

def static_dynamic_kinetic_compair_3(w_new: float):
    path = os.getcwd() + "/animations3"
    bin_path = os.getcwd() + "/bin_graphs3"
    try:
        shutil.rmtree(path)
        shutil.rmtree(bin_path)
    except:
        pass
    os.mkdir(path, 0o755)
    os.mkdir(bin_path, 0o755)

    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 5000 ± 2500, dynamic force of 0 ± 250, simulation runtime of 150000, spring ocilation frquency of {w_new}")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", "bin_graphs3")
    shutil.copy("run_conditions.txt", "animations3")

    dynamic_data3 = Analysis()
    static_data3 = Analysis()

    equations = get_spring_equations(5000, 2500, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(5000)
    dynamic_data3.spring_frequency = w
    frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200]
    #frequencies = [w/4, w/2,  w, (w*2), w*4]
    #frequencies = [110]

    #static_force dynamic bonds
    dynamic_data3.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"dynamic_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic3.txt", period=500)

    sim.apply_static_force(tag_list_left, 250, 0, 0)
    sim.apply_static_force(tag_list_right, -250, 0, 0)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    dynamic_data3.read_kinetic("kinetic3.txt")

    #static_force dynamic bonds
    static_data3.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"static_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic3.txt", period=500)

    sim.apply_static_force(tag_list_left, 250, 0, 0)
    sim.apply_static_force(tag_list_right, -250, 0, 0)

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    static_data3.read_kinetic("kinetic3.txt")

    shutil.copy("dynamic_bonds_force_at_0.gsd", "animations3")
    shutil.copy("static_bonds_force_at_0.gsd", "animations3")

    #---------------

    for frequency in frequencies:
        print(frequency)

        dynamic_data3.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.velocity_fname ="velocities3.txt"
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 250, frequency, 0)
        neg_f_equation = Force_Eq(0, 250, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic3.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        dynamic_data3.read_kinetic("kinetic3.txt")
        dynamic_data3.read_velocity("velocities3.txt")
        dynamic_data3.my_kinetic()
        dynamic_data3.kinetic_bin()

        shutil.copy(f"dynamic_bonds_force_at_{frequency}.gsd", "animations3")
        #dynamic_data.graph_k_energy()

        #static
        static_data3.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.velocity_fname ="velocities3.txt"
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 250, frequency, 0)
        neg_f_equation = Force_Eq(0, 250, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"static_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic3.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        static_data3.read_kinetic("kinetic3.txt")
        static_data3.read_velocity("velocities3.txt")
        static_data3.my_kinetic()
        static_data3.kinetic_bin()

        shutil.copy(f"static_bonds_force_at_{frequency}.gsd", "animations3")

        kinetic_bin_compare(static_data3, dynamic_data3, frequency, save_location="bin_graphs3")

    kinetic_compare(static_data3, dynamic_data3)
    print(w_new)

def static_dynamic_kinetic_compair_4(w_new: float):
    path = os.getcwd() + "/animations4"
    bin_path = os.getcwd() + "/bin_graphs4"
    try:
        shutil.rmtree(path)
        shutil.rmtree(bin_path)
    except:
        pass
    os.mkdir(path, 0o755)
    os.mkdir(bin_path, 0o755)

    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 5000 ± 2500, dynamic force of 0 ± 250, simulation runtime of 150000, spring ocilation frquency of {w_new}")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", "bin_graphs4")
    shutil.copy("run_conditions.txt", "animations4")

    dynamic_data4 = Analysis()
    static_data4 = Analysis()

    equations = get_spring_equations(5000, 2500, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(5000)
    dynamic_data4.spring_frequency = w
    frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200]
    #frequencies = [w/4, w/2,  w, (w*2), w*4]
    #frequencies = [110]

    #static_force dynamic bonds
    dynamic_data4.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.velocity_fname ="velocities4.txt"
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"dynamic_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic4.txt", period=500)

    sim.apply_static_force(tag_list_left, 250, 0, 0)
    sim.apply_static_force(tag_list_right, -250, 0, 0)

    for equation_i in range(len(equations)):
        sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    dynamic_data4.read_kinetic("kinetic4.txt")

    #static_force dynamic bonds
    static_data4.dynamic_frequencies.append(0)

    sim = Tri_Lattice()
    sim.velocity_fname ="velocities4.txt"
    sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

    tag_list_left = []
    tag_list_right = []
    for particle_i in range(1,len(sim.system.particles),):
        if particle_i < (2 * 10):
            tag_list_left.append(sim.system.particles[particle_i].tag)
            tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

    sim.create_animation(f"static_bonds_force_at_0")

    sim.log_system_kinetic_energy("kinetic4.txt", period=500)

    sim.apply_static_force(tag_list_left, 250, 0, 0)
    sim.apply_static_force(tag_list_right, -250, 0, 0)

    #for equation_i in range(len(equations)):
        #sim.change_spring_eq(equation_i, equations[equation_i])
    sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
    static_data4.read_kinetic("kinetic4.txt")

    shutil.copy("dynamic_bonds_force_at_0.gsd", "animations4")
    shutil.copy("static_bonds_force_at_0.gsd", "animations4")

    #---------------

    for frequency in frequencies:
        print(frequency)

        dynamic_data4.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.velocity_fname ="velocities4.txt"
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)
        f_equation = Force_Eq(0, 250, frequency, 0)
        neg_f_equation = Force_Eq(0, 250, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic4.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        dynamic_data4.read_kinetic("kinetic4.txt")
        dynamic_data4.read_velocity("velocities4.txt")
        dynamic_data4.my_kinetic()
        dynamic_data4.kinetic_bin()

        shutil.copy(f"dynamic_bonds_force_at_{frequency}.gsd", "animations4")
        #dynamic_data.graph_k_energy()

        #static
        static_data4.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.velocity_fname ="velocities4.txt"
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 250, frequency, 0)
        neg_f_equation = Force_Eq(0, 250, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"static_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic4.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        static_data4.read_kinetic("kinetic4.txt")
        static_data4.read_velocity("velocities4.txt")
        static_data4.my_kinetic()
        static_data4.kinetic_bin()

        shutil.copy(f"static_bonds_force_at_{frequency}.gsd", "animations4")

        kinetic_bin_compare(static_data4, dynamic_data4, frequency, save_location="bin_graphs4")

    kinetic_compare(static_data4, dynamic_data4)
    print(w_new)

def test_bins():
    data = Analysis()

    w = math.sqrt(5000)
    equations = get_spring_equations(5000, 2500, w, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(5000)
    data.spring_frequency = w
    #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
    #frequencies = [w/4, w/2,  w, (w*2), w*4]
    frequencies = [100]

    for frequency in frequencies:
        #print(frequency)

        data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=False)

        f_equation = Force_Eq(0, 750, frequency, 0)
        neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"static_100")

        #sim.log_system_kinetic_energy("kinetic.txt")

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
        data.read_kinetic("kinetic.txt")

    data.read_velocity("velocities.txt")
    data.kinetic_bin()

def test_bin_compare():

    dynamic_data = Analysis()
    static_data = Analysis()

    w = math.sqrt(5000)
    equations = get_spring_equations(5000, 2500, w, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    dynamic_data.spring_frequency = w
    #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
    #frequencies = [w/4, w/2,  w, (w*2), w*4]
    frequencies = [100]

    for frequency in frequencies:
        #print(frequency)

        dynamic_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 750, frequency, 0)
        neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic.txt", period = 500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        #dynamic_data.read_kinetic("kinetic.txt")

        dynamic_data.read_velocity("velocities.txt")
        dynamic_data.my_kinetic()
        dynamic_data.kinetic_bin()

    print("runing static")

    for frequency in frequencies:

        static_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 750, frequency, 0)
        neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 10):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

        sim.create_animation(f"static_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        #static_data.read_kinetic("kinetic.txt")

        static_data.read_velocity("velocities.txt")
        static_data.my_kinetic()
        static_data.kinetic_bin()

    kinetic_bin_compare(static_data, dynamic_data)

def static_dynamic_kinetic_compair_large(w_new: float):
    dynamic_data = Analysis()
    static_data = Analysis()

    equations = get_spring_equations(5000, 2500, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(5000)
    dynamic_data.spring_frequency = w
    frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200]
    #frequencies = [w/4, w/2,  w, (w*2), w*4]
    #frequencies = [w*2]

    for frequency in frequencies:
        #print(frequency)

        dynamic_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(30,30,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 250, frequency, 0)
        neg_f_equation = Force_Eq(0, 250, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 30):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 30 * 30) - (2 * 30)].tag)

        sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic.txt", period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        dynamic_data.read_kinetic("kinetic.txt")

    print("runing static")

    for frequency in frequencies:

        static_data.dynamic_frequencies.append(frequency)

        sim = Tri_Lattice()
        sim.create_lattice(30,30,5000,5000,5000, add_periodic_bonds_in_y=True)

        f_equation = Force_Eq(0, 250, frequency, 0)
        neg_f_equation = Force_Eq(0, 250, frequency, math.pi)

        tag_list_left = []
        tag_list_right = []
        for particle_i in range(1,len(sim.system.particles),):
            if particle_i < (2 * 30):
                tag_list_left.append(sim.system.particles[particle_i].tag)
                tag_list_right.append(sim.system.particles[particle_i + (2 * 30 * 30) - (2 * 30)].tag)

        sim.create_animation(f"static_bonds_force_at_{frequency}")

        sim.log_system_kinetic_energy("kinetic.txt",period=500)

        sim.add_dynamic_force(f_equation, tag_list_left)
        sim.add_dynamic_force(neg_f_equation, tag_list_right)

        #for equation_i in range(len(equations)):
            #sim.change_spring_eq(equation_i, equations[equation_i])
        sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
        static_data.read_kinetic("kinetic.txt")

    kinetic_compare(static_data, dynamic_data)
    print(w_new)

def static_dynamic_kinetic_compair_multi_size(w_new: float):
    dynamic_data = Analysis()
    static_data = Analysis()

    equations = get_spring_equations(5000, 2500, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    w = math.sqrt(5000)
    dynamic_data.spring_frequency = w
    #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
    frequencies = [w, w*2]
    sizes = [10,20,30]
    #frequencies = [w*2]

    for size in sizes:

        for frequency in frequencies:
            #print(frequency)

            dynamic_data.dynamic_frequencies.append(frequency)

            sim = Tri_Lattice()
            sim.create_lattice(30,30,5000,5000,5000, add_periodic_bonds=False)

            f_equation = Force_Eq(0, 750, frequency, 0)
            neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

            tag_list_left = []
            tag_list_right = []
            for particle_i in range(1,len(sim.system.particles),):
                if particle_i < (2 * size):
                    tag_list_left.append(sim.system.particles[particle_i].tag)
                    tag_list_right.append(sim.system.particles[particle_i + (2 * size * size) - (2 * size)].tag)

            sim.create_animation(f"#dynamic_bonds_force_at_{frequency}_size_{size}")

            sim.log_system_kinetic_energy("kinetic.txt")

            sim.add_dynamic_force(f_equation, tag_list_left)
            sim.add_dynamic_force(neg_f_equation, tag_list_right)

            for equation_i in range(len(equations)):
                sim.change_spring_eq(equation_i, equations[equation_i])
            sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
            dynamic_data.read_kinetic("kinetic.txt")

        print("runing static")

        for frequency in frequencies:

            static_data.dynamic_frequencies.append(frequency)

            sim = Tri_Lattice()
            sim.create_lattice(30,30,5000,5000,5000, add_periodic_bonds=False)

            f_equation = Force_Eq(0, 750, frequency, 0)
            neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

            tag_list_left = []
            tag_list_right = []
            for particle_i in range(1,len(sim.system.particles),):
                if particle_i < (2 * size):
                    tag_list_left.append(sim.system.particles[particle_i].tag)
                    tag_list_right.append(sim.system.particles[particle_i + (2 * size * size) - (2 * size)].tag)

            sim.create_animation(f"#static_bonds_force_at_{frequency}_size_{size}")

            sim.log_system_kinetic_energy("kinetic.txt")

            sim.add_dynamic_force(f_equation, tag_list_left)
            sim.add_dynamic_force(neg_f_equation, tag_list_right)

            #for equation_i in range(len(equations)):
                #sim.change_spring_eq(equation_i, equations[equation_i])
            sim.run_langevin(150000,callback=sim.log_p_info,callback_period=1000, quiet=True, dynamic_f_stop=150000)
            static_data.read_kinetic("kinetic.txt")

        kinetic_compare(static_data, dynamic_data)
        print(w_new, size)

def kinetic_energy_issue():
        dynamic_data = Analysis()
        static_data = Analysis()

        w_new = math.sqrt(5000)
        equations = get_spring_equations(5000, 2500, w_new, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
        w = math.sqrt(5000)
        dynamic_data.spring_frequency = w
        #frequencies = [10, (w/4), 20, 30, (w/2), 50, 60, w, 80, 90, 100, 110, 120, 130, (w*2), 170, 200, 230, 260, (w*4)]
        #frequencies = [w/4, w/2,  w, (w*2), w*4]
        frequencies = [w/2]

        for frequency in frequencies:
            #print(frequency)

            dynamic_data.dynamic_frequencies.append(frequency)

            sim = Tri_Lattice()
            sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=False)

            f_equation = Force_Eq(0, 750, frequency, 0)
            neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

            tag_list_left = []
            tag_list_right = []
            for particle_i in range(1,len(sim.system.particles),):
                if particle_i < (2 * 10):
                    tag_list_left.append(sim.system.particles[particle_i].tag)
                    tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

            sim.create_animation(f"dynamic_bonds_force_at_{frequency}")

            sim.log_system_kinetic_energy("kinetic.txt", period = 500)

            sim.add_dynamic_force(f_equation, tag_list_left)
            sim.add_dynamic_force(neg_f_equation, tag_list_right)

            for equation_i in range(len(equations)):
                sim.change_spring_eq(equation_i, equations[equation_i])
            sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
            dynamic_data.read_kinetic("kinetic.txt")

            dynamic_data.read_velocity("velocities.txt")
            dynamic_data.my_kinetic()
            dynamic_data.kinetic_bin()

        print("runing static")

        for frequency in frequencies:

            static_data.dynamic_frequencies.append(frequency)

            sim = Tri_Lattice()
            sim.create_lattice(10,10,5000,5000,5000, add_periodic_bonds=False)

            f_equation = Force_Eq(0, 750, frequency, 0)
            neg_f_equation = Force_Eq(0, 750, frequency, math.pi)

            tag_list_left = []
            tag_list_right = []
            for particle_i in range(1,len(sim.system.particles),):
                if particle_i < (2 * 10):
                    tag_list_left.append(sim.system.particles[particle_i].tag)
                    tag_list_right.append(sim.system.particles[particle_i + (2 * 10 * 10) - (2 * 10)].tag)

            sim.create_animation(f"static_bonds_force_at_{frequency}")

            sim.log_system_kinetic_energy("kinetic.txt", period = 500)

            sim.add_dynamic_force(f_equation, tag_list_left)
            sim.add_dynamic_force(neg_f_equation, tag_list_right)

            #for equation_i in range(len(equations)):
                #sim.change_spring_eq(equation_i, equations[equation_i])
            sim.run_langevin(150000,callback=sim.log_p_info,callback_period=500, quiet=True, dynamic_f_stop=150000)
            static_data.read_kinetic("kinetic.txt")

            static_data.read_velocity("velocities.txt")
            static_data.my_kinetic()
            static_data.kinetic_bin()

        kinetic_bin_compare(static_data, dynamic_data)
        kinetic_compare(static_data, dynamic_data)
        print(w_new)
