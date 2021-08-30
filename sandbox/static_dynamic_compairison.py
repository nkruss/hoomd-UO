from triangular_lattice import *
from spring_line import *
from random_line import *
from line_force import *
from two_mass_rand import *
from three_mass_rand import *
from analysis import *
from Equations import *
from triangular_lattice_2 import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import shutil

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
                                      (2 * lattice_width)].tag)

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
        if dynamic:
            data.width_of_tri(0, "dynamic", gamma, w)
        else:
            data.width_of_tri(0, "static", gamma, w)

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


            #deal with static force
            if force_freq == 0.0:

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)


                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data")

                shutil.copy(f"static_bonds_force_at_frequency{force_freq}.gsd", "animations1")

def static_dynamic_compair_2():
    #different lattice
    # #runs compairison on the new latice configuration
    # equations = get_spring_equations(1, .5, 1, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
    # w = equations[0].w
    # frequencies = [w]
    #
    # data = Analysis()
    #
    # for frequency in frequencies:
    #     print(f"runing sim with frequency {frequency}")
    #
    #     sim = Tri_Lattice_diff()
    #     sim.create_lattice(10, 10,1,1,1)
    #     f_equation = Force_Eq(0, .1, frequency, 0)
    #     neg_f_equation = Force_Eq(0, .1, frequency, math.pi)
    #
    #     tag_list_left = []
    #     tag_list_right = []
    #     for particle_i in range(len(sim.system.particles)):
    #         if particle_i < (2 * 10):
    #             tag_list_left.append(sim.system.particles[particle_i].tag)
    #             tag_list_right.append(sim.system.particles[particle_i +
    #                                   (2 * 10 * 10) -
    #                                   (2 * 10)].tag)
    #
    #     sim.apply_static_force(tag_list_left, -.1, 0, 0)
    #     sim.apply_static_force(tag_list_right, .1, 0, 0)
    #
    #     sim.create_animation(f"new_lattice")
    #     sim.run_langevin(30000,callback=sim.log_p_info,gamma=1, callback_period=500, quiet=True, dynamic_f_stop=0)

    #create folders to store graphs and gsd files
    path = os.getcwd() + "/animations2"
    sim_data_path = os.getcwd() + "/sim_data2"
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
    shutil.copy("run_conditions.txt", "animations2")
    shutil.copy("run_conditions.txt", "sim_data2")


    #get a list of different spring frequencies to test over
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


            #deal with static force
            if force_freq == 0.0:

                ##static_force dynamic bonds

                sim = Tri_Lattice_diff()
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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data2")
                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data2")


                shutil.copy("dynamic_bonds_static_force.gsd", "animations2")

                # ------------------------------------------

                ##static_force static bonds

                sim = Tri_Lattice_diff()
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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data2")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data2")

                shutil.copy("static_bonds_static_force.gsd", "animations2")

            else:

                f_equation = Force_Eq(0, .1, force_freq, 0)
                neg_f_equation = Force_Eq(0, .1, force_freq, math.pi)

                #dynamic test

                sim = Tri_Lattice_diff()
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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data2")
                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data2")

                shutil.copy(f"dynamic_bonds_force_at_frequency{force_freq}.gsd", "animations2")

                #------------------------------------------------

                #static test

                sim = Tri_Lattice_diff()
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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)


                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data2")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data2")

                shutil.copy(f"static_bonds_force_at_frequency{force_freq}.gsd", "animations2")


def static_dynamic_compair_around_1():
    #create folders to store graphs and gsd files
    path = os.getcwd() + "/animations_around1"
    bin_path = os.getcwd() + "/bin_graphs_around1"
    sim_data_path = os.getcwd() + "/sim_data_around1"
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
    shutil.copy("run__conditions.txt", "bin_graphs_around1")
    shutil.copy("run_conditions.txt", "animations_around1")
    shutil.copy("run_conditions.txt", "sim_data_around1")


    #get a list of different spring frequencies to test over
    spring_frequencies = []
    force_frequencies = []
    w = .75
    while(w<=1.25):
        spring_frequencies.append(w)
        force_frequencies.append(w)
        w = w+.10

    results = []
    kinetic_differences = []
    for spring_freq in spring_frequencies:

        equations = get_spring_equations(1, .5, spring_freq, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])

        for force_freq in force_frequencies:

            dynamic_data = Analysis()
            static_data = Analysis()
            print("\n-------------------------------------------------")
            print(f"running force_freq{force_freq}, spring_freq{spring_freq} \n")


            #deal with static force
            if force_freq == 0.0:

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_around1")
                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data_around1")


                shutil.copy("dynamic_bonds_static_force.gsd", "animations_around1")

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_around1")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data_around1")

                shutil.copy("static_bonds_static_force.gsd", "animations_around1")

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_around1")
                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data_around1")

                shutil.copy(f"dynamic_bonds_force_at_frequency{force_freq}.gsd", "animations_around1")

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

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)


                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_around1")
                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data_around1")

                shutil.copy(f"static_bonds_force_at_frequency{force_freq}.gsd", "animations_around1")

def static_dynamic_center_compair():
    #create folders to store graphs and gsd files
    path = os.getcwd() + "/animations_center"
    sim_data_path = os.getcwd() + "/sim_data_center"
    try:
        shutil.rmtree(path)
        shutil.rmtree(sim_data_path)
    except:
        pass
    os.mkdir(path, 0o755)
    os.mkdir(sim_data_path, 0o755)

    #add a file with the simulations run conditions
    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 1 ± .5, dynamic force of 0 ± .1 \nsimulation runtime of 3000000, data_recorded every 500 timesteps, periodic bonds in the y direction \ngamma =.05")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run_conditions.txt", "animations_center")
    shutil.copy("run_conditions.txt", "sim_data_center")


    #get a list of different spring frequencies to test over
    spring_frequencies = [1, 0.0]
    force_frequencies = [1, 0.0]
    w = .01
    # while(w<=4):
    #     spring_frequencies.append(w)
    #     force_frequencies.append(w)
    #     w = w*4

    results = []
    kinetic_differences = []
    for spring_freq in spring_frequencies:

        equations = get_spring_equations(1, .5, spring_freq, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])

        for force_freq in force_frequencies:

            dynamic_data = Analysis()
            static_data = Analysis()
            print("\n-------------------------------------------------")
            print(f"running force_freq{force_freq}, spring_freq{spring_freq} \n")


            #deal with static force
            if force_freq == 0.0:

                ##static_force dynamic bonds

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds=True)

                tag_list = [90]
                sim.apply_static_force(tag_list, 0, 1, 0)

                sim.create_animation(f"dynamic_bonds_static_force_center")

                sim.log_system_kinetic_energy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", period=500)

                for equation_i in range(len(equations)):
                    sim.change_spring_eq(equation_i, equations[equation_i])

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_center")

                velocity_file = f"velocity_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("velocities.txt", velocity_file)
                shutil.copy(velocity_file, "sim_data_center")

                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data_center")

                shutil.copy("dynamic_bonds_static_force.gsd", "animations_center")

                # ------------------------------------------

                ##static_force static bonds

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds=True)

                tag_list = [90]
                sim.apply_static_force(tag_list, 0, 1, 0)

                sim.create_animation(f"static_bonds_static_force_center")

                sim.log_system_kinetic_energy(f"kinetic_static_{spring_freq}_{force_freq}.txt", period=500)

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_center")

                velocity_file = f"velocity_static_{spring_freq}_{force_freq}.txt"
                os.rename("velocities.txt", velocity_file)
                shutil.copy(velocity_file, "sim_data_center")

                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data_center")

                shutil.copy("static_bonds_static_force.gsd", "animations_center")

            else:

                f_equation = Force_Eq(0, 1, force_freq, 0)
                neg_f_equation = Force_Eq(0, 1, force_freq, math.pi)

                #dynamic test

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds=True)

                tag_list = [90]
                sim.add_dynamic_force(f_equation, tag_list)

                sim.create_animation(f"dynamic_bonds_force_at_frequency{force_freq}_center")

                sim.log_system_kinetic_energy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", period=500)

                for equation_i in range(len(equations)):
                    sim.change_spring_eq(equation_i, equations[equation_i])

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_center")

                velocity_file = f"velocity_dynamic_{spring_freq}_{force_freq}.txt"
                os.rename("velocities.txt", velocity_file)
                shutil.copy(velocity_file, "sim_data_center")

                shutil.copy(f"kinetic_dynamic_{spring_freq}_{force_freq}.txt", "sim_data_center")

                shutil.copy(f"dynamic_bonds_force_at_frequency{force_freq}.gsd", "animations_center")

                #------------------------------------------------

                #static test

                sim = Tri_Lattice()
                sim.create_lattice(10,10,1,1,1, add_periodic_bonds=True)

                tag_list = [90]
                sim.add_dynamic_force(f_equation, tag_list)

                sim.create_animation(f"static_bonds_force_at_frequency{force_freq}_center")

                sim.log_system_kinetic_energy(f"kinetic_static_{spring_freq}_{force_freq}.txt", period=500)

                sim.run_langevin(3000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=3000000)

                position_file = f"positions_static_{spring_freq}_{force_freq}.txt"
                os.rename("positions.txt", position_file)
                shutil.copy(position_file, "sim_data_center")

                velocity_file = f"velocity_static_{spring_freq}_{force_freq}.txt"
                os.rename("velocities.txt", velocity_file)
                shutil.copy(velocity_file, "sim_data_center")

                shutil.copy(f"kinetic_static_{spring_freq}_{force_freq}.txt", "sim_data_center")

                shutil.copy(f"static_bonds_force_at_frequency{force_freq}.gsd", "animations_center")

def analysis(filename: str, graph_type: str, compair_type: str):
    """
    filename - name of folder where position files are stored
    graph type options - "skyscraper", "scatter"
    compair type options - "ratio", "log_ratio", "abs_difference"
    """
    #location of target files
    data_files_locations = os.getcwd() + f"/{filename}"

    spring_frequencies = []
    force_frequencies = []
    dynamic_files_kinetic = []
    static_files_kinetic = []
    dynamic_files_positions = []
    static_files_positions = []
    results = []

    #collect files to analyse
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)
        if filename.startswith("kinetic_dynamic"):
             fname_parts = filename.split("_")
             spring_freq = float(fname_parts[2])
             force_freq = fname_parts[3].split('t')
             force_freq = float(force_freq[0][:-1])

             if spring_freq not in spring_frequencies:
                 spring_frequencies.append(spring_freq)
             if force_freq not in force_frequencies:
                 force_frequencies.append(force_freq)
             dynamic_files_kinetic.append(filename)
             continue
        elif filename.startswith("kinetic_static"):
            static_files_kinetic.append(filename)
        else:
            continue

    #collect position files to analyse
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)
        if filename.startswith("positions_dynamic"):
             fname_parts = filename.split("_")
             spring_freq = float(fname_parts[2])
             force_freq = fname_parts[3].split('t')
             force_freq = float(force_freq[0][:-1])

             if spring_freq not in spring_frequencies:
                 spring_frequencies.append(spring_freq)
             if force_freq not in force_frequencies:
                 force_frequencies.append(force_freq)
             dynamic_files_positions.append(filename)
             continue
        elif filename.startswith("positions_static"):
            static_files_positions.append(filename)
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

    print(static_files_kinetic)
    print(dynamic_files_kinetic)

    #kinetic analysis
    for i in range(len(dynamic_files_kinetic)):
        print("in loop", i)
        print()
        static_data = Analysis()
        dynamic_data = Analysis()
        # print(static_files[i])
        # print(dynamic_files[i])

        fname_parts = static_files_kinetic[i].split("_")
        spring_freq = float(fname_parts[2])
        force_freq = fname_parts[3].split('t')
        force_freq = float(force_freq[0][:-1])

        #read corresponding kinetic energy files
        static_data.read_kinetic(data_files_locations + '/' + static_files_kinetic[i])
        dynamic_data.read_kinetic(data_files_locations + '/' + dynamic_files_kinetic[i])

        static_avg = static_data.average_kinetic[-1]
        print(f"static_avg ={static_avg}")
        dynamic_avg = dynamic_data.average_kinetic[-1]
        print(f"dynamic_avg = {dynamic_avg}")

        diff = abs(dynamic_avg - static_avg)
        ratio = abs((dynamic_avg - static_avg)) / static_avg
        log_ratio = math.log10(ratio)

        if compair_type == "ratio":
            results.append((spring_freq, force_freq, ratio))

        #log (base 10) ratio
        elif compair_type == "log_ratio":
            results.append((spring_freq, force_freq, log_ratio))

        elif comapir_type == "abs_difference":
            results.append((spring_freq, force_freq, diff))

    print("\n saving results")
    np.savetxt("static_dynamic_compair.npy",results,header="tuples of data points in the format (spring_freq, force_freq, ratio)")

    #static_dynamic_analysis(graph_type)

    ##########
    static_data = Analysis()
    fname_parts = static_files_positions[0].split("_")
    spring_freq = float(fname_parts[2])
    force_freq = fname_parts[3].split('t')
    force_freq = float(force_freq[0][:-1])

    title = f"static \nspring_freq={spring_freq}, force_freq={force_freq}"

    #read corresponding files
    static_data.read_pos(data_files_locations + '/' + static_files_positions[0])
    static_data.read_velocity(data_files_locations + '/' + static_files_positions[0])
    static_data.heat_map(title)

    ####
    dynamic_data = Analysis()
    fname_parts = dynamic_files_positions[0].split("_")
    spring_freq = float(fname_parts[2])
    force_freq = fname_parts[3].split('t')
    force_freq = float(force_freq[0][:-1])

    title = f"dynamic \nspring_freq={spring_freq}, force_freq={force_freq}"

    #read corresponding files
    dynamic_data.read_pos(data_files_locations + '/' + static_files_positions[0])
    dynamic_data.read_velocity(data_files_locations + '/' + static_files_positions[0])
    dynamic_data.heat_map(title)
