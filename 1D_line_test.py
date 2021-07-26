from triangular_lattice import *
from spring_line import *
from random_line import *
from standing_line import *
from line_force import *
from two_mass_rand import *
from three_mass_rand import *
from analysis import *
from Equations import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import shutil

def line_test(num_bonds: int):
    #known parameters that work are force (2.5), damp (.1), chain (210), time (1000000)
    #.1 works for damping

    # force_freq_list = []
    # freq = 0
    # while freq < 2.5:
    #     freq += .1
    #     force_freq_list.append(freq)

    force_freq_list = [.6]

    for freq in force_freq_list:

        #add a file with the simulations run conditions
        f = open("run_conditions.txt", "w+")
        f.write(f"springs constants of 1 ± .5, spring_freq = .75, dynamic force of 0 ± 2.5 \nsimulation runtime of 1000000, data_recorded every 500 timesteps \n")
        f.close()
        shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
        shutil.copy("run__conditions.txt", f"1Dline_{num_bonds}")

        #create and run simulation
        sim = Line()
        if num_bonds == 0:
            sim.create_lattice(1, num_bonds=1, add_periodic_bonds=False, N=210)
        else:
            sim.create_lattice(1, num_bonds=num_bonds, add_periodic_bonds=False, N=210)
        sim.dimension_constrain([1,0,0])
        sim.create_animation(f"line_test_{num_bonds}bonds_{freq}")
        sim.log_system_kinetic_energy(f"kinetic_{freq}.txt", period=500)

        #add force and pin one edge of the line
        f_equation = Force_Eq(0, 2.5, freq, 0)
        sim.add_dynamic_force(f_equation, [0], dimension="x")
        sim.pin_particles([209])

        #run with static bonds
        if num_bonds == 0:
            pass
        #change bonds to be dynamic
        else:
            if num_bonds == 3:
                equations = get_spring_equations(1, .5, 2, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
            elif num_bonds == 2:
                equations = get_spring_equations(1, .5 , 2, [0, math.pi])
            elif num_bonds == 1:
                equations = get_spring_equations(1, .5, 2, [0])

            for equation_i in range(len(equations)):
                sim.change_spring_eq(equation_i, equations[equation_i])

        sim.run_langevin(2500000, callback=sim.log_p_info, gamma=.005, callback_period=500, dynamic_f_stop=2500000)
        #sim.run(1000000, callback=sim.log_p_info, callback_period=500, dynamic_f_stop=1000000)

        #store data into folder
        position_file = f"positions_1Dline_{num_bonds}bonds_{freq}.txt"
        velocity_file = f"velocities_1Dline_{num_bonds}bonds_{freq}.txt"
        os.rename("positions.txt", position_file)
        os.rename("velocities.txt", velocity_file)
        shutil.copy(position_file, f"1Dline_{num_bonds}bonds")
        shutil.copy(velocity_file, f"1Dline_{num_bonds}bonds")
        shutil.copy(f"line_test_{num_bonds}bonds_{freq}.gsd", f"1Dline_{num_bonds}bonds")
        shutil.copy(f"kinetic_{freq}.txt", f"1Dline_{num_bonds}bonds")

def run_sims():
    #create folders to store data and gsd files
    main_path = os.getcwd() + f"/1Dline"
    try:
        shutil.rmtree(main_path)
    except:
        pass
    os.mkdir(main_path, 0o755)

    for i in range(1, 2):
        num_bonds = i

        #create folders to store data and gsd files
        path = os.getcwd() + f"/1Dline_{num_bonds}bonds"
        target = main_path + f"/1Dline_{num_bonds}bonds"
        try:
            shutil.rmtree(path)
        except:
            pass
        os.mkdir(path, 0o755)

        #run simulation for current number of bonds
        #line_test(num_bonds)

        test = [7]
        for p in test:
            standing_test(num_bonds, p)

        # for p in range(1,35):
        #     standing_test(num_bonds, p)

        shutil.move(path, target)

def standing_test(num_bonds: int, p: int):
    """
    Run function for standing wave trials. Creates a 1D_line with (num_bonds)
    number of bonds and a standing wave with p_value = p constriand in the x
    dimension. Then runs the simulation forward in time for 10000000 timesteps.
    """

    A_p = .2
    B_p = 0
    N = 36

    #add a file with the simulations run conditions
    f = open("run_conditions.txt", "w+")
    f.write(f"springs constants of 1 ± .5, spring_freq = 2, p = {p} \nsimulation runtime of 10000000, data_recorded every 500 timesteps \n")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", f"1Dline_{num_bonds}")

    #create and run simulation
    sim = Standing_Line()
    if num_bonds == 0:
        sim.create_lattice(1, p, num_bonds=1, add_periodic_bonds=False, N=N, A_p=A_p, B_p=B_p)
    else:
        sim.create_lattice(1, p, num_bonds=num_bonds, add_periodic_bonds=False, N=N, A_p=A_p, B_p=B_p)

    sim.dimension_constrain([1,0,0])
    sim.create_animation(f"line_test_{num_bonds}bonds_p_of_{p}")
    sim.log_system_kinetic_energy(f"kinetic_p_of_{p}.txt", period=500)

    #pin edges of the lattice
    sim.pin_particles([0, N-1])

    #run with static bonds
    if num_bonds == 0:
        pass
    #change bonds to be dynamic
    else:
        if num_bonds == 3:
            equations = get_spring_equations(1, .5, 2, [0, (2 * math.pi / 3.0), (4 * math.pi / 3.0)])
        elif num_bonds == 2:
            equations = get_spring_equations(1, .5 , 2, [0, math.pi])
        elif num_bonds == 1:
            equations = get_spring_equations(1, .5, .2995382286, [0])

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])

    sim.run(10000000, kt=0, dt=.0001, callback=sim.log_p_info, callback_period=500)
    #sim.run(50000000, kt=0, dt=.0001, callback=sim.log_p_info, callback_period=500)

    #store data into folder
    position_file = f"positions_1Dline_{num_bonds}bonds_p_of_{p}.txt"
    velocity_file = f"velocities_1Dline_{num_bonds}bonds_p_of_{p}.txt"
    os.rename("positions.txt", position_file)
    os.rename("velocities.txt", velocity_file)
    shutil.copy(position_file, f"1Dline_{num_bonds}bonds")
    shutil.copy(velocity_file, f"1Dline_{num_bonds}bonds")
    shutil.copy(f"line_test_{num_bonds}bonds_p_of_{p}.gsd", f"1Dline_{num_bonds}bonds")
    shutil.copy(f"kinetic_p_of_{p}.txt", f"1Dline_{num_bonds}bonds")

def mathieu():
    #Solved solution parameters

    # # a_2 solution
    # w = .1654983743
    # dk = .5
    # k_0 = 1
    # p = 2

    # # a_1 solution
    # w = .3044566161
    # dk = .5
    # k_0 = 1
    # p = 2

    # #b_2 solution
    # w = .1762997065
    # dk = .5
    # k_0 = 1
    # p = 2

    # #a_2 solution
    # w = 2.81339015
    # dk = .5
    # k_0 = 1
    # p = 34

    # #a_4 solution
    # w = 1.45124937
    # dk = .5
    # k_0 = 1
    # p = 34

    #p = 2 unstable
    w = .1711641876
    dk = .5
    k_0 = 1
    p = 2

    #p = 2 unstable
    # first w is precise but doesn't seem to fall into unstable
    w = .2949026294
    w = .2995382286
    dk = .5
    k_0 = 1
    p = 7


    #----------------------------------------

    ## code for plotting each p value point along the line a = q/gamma
    for p in range(36):
        b = math.pi * p / 36

        gamma = dk / (2 * k_0)
        a = 4 * (b ** 2) * k_0 / (w ** 2)
        q = gamma * a
        plt.scatter(q, a)
        plt.text(q+.03, a, p, fontsize=4)

    #----------------------------------------

    p = 2

    b = math.pi * p / 36

    gamma = dk / (2 * k_0)
    a = 4 * (b ** 2) * k_0 / (w ** 2)
    q = gamma * a
    print(q, a)
    plt.scatter(q, a)

    q_values = np.linspace(-10, 100, 1000)

    odd_dic = {}
    even_dic = {}
    m_values = []
    for i in range(20):
        m_values.append(i)
    legend_titles = []
    for m in m_values:
        odd_solutions = []
        even_solutions = []
        for q in q_values:
            odd = scipy.special.mathieu_b(m+1, q)
            even = scipy.special.mathieu_a(m, q)

            odd_solutions.append(odd)
            even_solutions.append(even)

        odd_dic[m] = odd_solutions
        even_dic[m] = even_solutions
        legend_titles.append(f"b_{m+1}")
        legend_titles.append(f"a_{m}")

        plt.plot(q_values, odd_solutions)
        plt.plot(q_values, even_solutions)

    plt.plot(q_values, q_values / abs(gamma))
    legend_titles.append(f"a = q / gamma")

    plt.xlabel("q")
    plt.ylabel("a")
    # plt.xlim([0,5])
    # plt.ylim([0,20])
    #plt.legend(legend_titles)
    plt.title(f"Mathiue stability plot: w = {w}")
    plt.show()

    # ----------------------------------
    # Individual mathieu theoretical plot
    # for p in range(36):
    #     t = np.linspace(0, 500, 500)
    #     z = w * t / 2 * (180 / math.pi)
    #     b = math.pi * p / 36
    #     gamma = -dk / (2 * k_0)
    #     a = 4 * (b ** 2) * k_0 / (w ** 2)
    #     q = gamma * a
    #
    #     plt.plot(z, scipy.special.mathieu_cem(m, q, z)[0] - (36 / 2) + 2)
    #     plt.show()

    t = np.linspace(0, 10000, 1000)
    z_t = w * t / 2
    b = math.pi * p / 36
    gamma = -dk / (2 * k_0)
    a = 4 * (b ** 2) * k_0 / (w ** 2)
    q = gamma * a

    plt.plot(z_t, scipy.special.mathieu_cem(1, q, z_t)[0] - (36 / 2) + 2)
    #plt.plot(t, scipy.special.mathieu_cem(1, -q, z_t)[0] - (36 / 2) + 2)
    #plt.plot(t, scipy.special.mathieu_sem(1, q, z_t)[0] - (36 / 2) + 2)
    plt.xlabel("z")
    plt.ylabel("position")
    plt.show()

    p_list = np.linspace(0, 36, 100)
    a_1 = []
    a_2 = []
    for p in p_list:
        a_1.append(((2 / w) ** 2) * 2 * (1 - math.cos(p * math.pi / 36)))
        a_2.append((4 * ((p * math.pi / 36) ** 2) / (w ** 2)))
    plt.plot(p_list, a_1)
    plt.plot(p_list, a_2)
    plt.legend(["a = (2 / w)^2 * 2 * (1 - cos(p * pi / 36))", "a = 4 * ((p * pi / 36)^2) / w^2)"])
    plt.xlabel("p")
    plt.ylabel("w")
    plt.show()


def data_analysis(data_folder, num_bonds):

    #create folders to store created graphs
    path = os.getcwd() + "/1Dline_analysis"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path, 0o755)

    #location of target files
    data_files_locations = os.getcwd() + f"/{data_folder}/1Dline_{num_bonds}bonds"

    force_frequencies = []
    files_kinetic = []
    files_positions = []
    files_velocities = []
    results = []

    #collect files to analyse
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)

        #collect kinetic energy files and record force frequencies
        if filename.startswith("kinetic"):
             fname_parts = filename.split("_")
             force_freq = fname_parts[1].split('t')
             force_freq = float(force_freq[0][:-1])
             if force_freq not in force_frequencies:
                 force_frequencies.append(force_freq)
             files_kinetic.append(filename)
             continue

        #collect position files
        elif filename.startswith("positions"):
            files_positions.append(filename)

        #collect velocity files
        elif filename.startswith("velocities"):
            files_velocities.append(filename)

        else:
            continue

    wavelengths = []
    for i in range(len(files_kinetic)):

        data = Analysis()

        #get the force freq for plot titles
        fname_parts = files_kinetic[i].split("_")
        force_freq = fname_parts[1].split('t')
        force_freq = force_freq[0][:-1]

        force_freq_str = force_freq.split('.')
        if len(force_freq_str) == 1:
            force_freq_str = force_freq_str[0]
        else:
            force_freq_str = force_freq_str[0] + "pt" + force_freq_str[1]

        force_freq = round(float(force_freq), 2)
        print(force_freq)

        if num_bonds == 0:
            kinetic_plot_title = f"1Dline_system_forcefreq_{force_freq_str}_static_bonds"
            title = f"1Dline test: num_bonds = {num_bonds}, force freq = {force_freq}"
        else:
            kinetic_plot_title = f"1Dline_system_forcefreq_{force_freq_str}_numbonds_{num_bonds}"
            title = f"1Dline test: num_bonds = {num_bonds}, force freq = {force_freq} \nspring freq = 2, spring mag = 1±.5"


        #read corresponding kinetic energy files
        data.read_kinetic(data_files_locations + '/' + files_kinetic[i])

        #read corresponding position files
        data.read_pos(data_files_locations + '/' + files_positions[i])

        #read corresponding velocity files
        data.read_velocity(data_files_locations + '/' + files_velocities[i])

        #create plots of data
        peak = data.graph_spring_line_particle_kinetic(title, path)
        data.graph_k_energy(0, kinetic_plot_title, path)

        freq = peak[0]
        wavelength = 1 / freq
        wavelengths.append(wavelength)

    print(wavelengths)
    plt.plot(force_frequencies, wavelengths)
    plt.xlabel('frequencies')
    plt.ylabel('wavelength')
    plt.show()
    plt.clf()

    for i in range(len(wavelengths)):
        wavelengths[i] = (2 * math.pi) / (wavelengths[i] * 2)
    plt.plot(wavelengths[3:], force_frequencies[3:])
    plt.ylabel('frequencies')
    plt.xlabel('q_a')
    plt.show()

def standing_analysis(data_folder, num_bonds):

    #location of target files
    data_files_locations = os.getcwd() + f"/{data_folder}/1Dline_{num_bonds}bonds"

    #create folders to store created graphs
    path = data_files_locations + "/Analysis_plots"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path, 0o755)
    save_loc = path

    p_values = []
    files_kinetic = []
    files_positions = []
    files_velocities = []
    results = []

    #collect files to analyse
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)

        #collect kinetic energy files and record force frequencies
        if filename.startswith("kinetic"):
             fname_parts = filename.split("_")
             p = fname_parts[3].split('t')
             p = float(p[0][:-1])
             if p not in p_values:
                 p_values.append(p)
             files_kinetic.append(filename)
             continue

        #collect position files
        elif filename.startswith("positions"):
            files_positions.append(filename)

        #collect velocity files
        elif filename.startswith("velocities"):
            files_velocities.append(filename)

        else:
            continue

    w_p_theoretical = []
    w_p_experimatal = []
    p_values = []

    #Analyse each run
    for i in range(len(files_kinetic)):

        data = Analysis()

        #get the force freq for plot titles
        fname_parts = files_kinetic[i].split("_")
        p = fname_parts[3].split('t')
        p = float(p[0][:-1])

        #Create apropriate save titles
        if num_bonds == 0:
            save_title = f"1Dline_static_bonds"
            title = f"1Dline test: Static bonds, spring mag = 1"
        else:
            save_title = f"1Dline_{num_bonds}_dynamic_bonds"
            title = f"1Dline test: num_bonds = {num_bonds}, spring mag = 1±.5"

        #read corresponding kinetic energy files
        data.read_kinetic(data_files_locations + '/' + files_kinetic[i])

        #read corresponding position files
        data.read_pos(data_files_locations + '/' + files_positions[i])

        #read corresponding velocity files
        data.read_velocity(data_files_locations + '/' + files_velocities[i])

        freq_info = data.standing_wave_1D(p, .2, 0, save_loc, num_bonds)

        w_p_theoretical.append(freq_info[0])
        w_p_experimatal.append(freq_info[1])
        p_values.append(p)

    #plot the data results
    plt.scatter(p_values, w_p_theoretical)
    plt.scatter(p_values, w_p_experimatal)
    plt.xlabel("p value")
    plt.ylabel("w_p")
    plt.title(title)
    plt.savefig(save_title)
    shutil.copy(save_title + ".png", os.getcwd() + f"/{data_folder}")

    plt.show()
