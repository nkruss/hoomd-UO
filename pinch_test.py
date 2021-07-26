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

def various_force_freq_test():
    #create folders to store data and gsd files
    path = os.getcwd() + "/pinch_animations"
    bin_path = os.getcwd() + "/pinch_bin_graphs"
    sim_data_path = os.getcwd() + "/pinch_sim_data"
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
    f.write(f"springs constants of 1 ± .5, spring_freq = .75, dynamic force of 0 ± 1, force_freq = various \nsimulation runtime of 1000000, data_recorded every 500 timesteps \ngamma =.05")
    f.close()
    shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run__conditions.txt")
    shutil.copy("run__conditions.txt", "pinch_animations")
    shutil.copy("run_conditions.txt", "pinch_bin_graphs")
    shutil.copy("run_conditions.txt", "pinch_sim_data")

    equations = get_spring_equations(1, .5, .75, [(2 * math.pi / 3.0), (4 * math.pi / 3.0), 0])

    force_frequencies = [1]
    w = .01
    while(w<=10):
        force_frequencies.append(w)
        w = w*4

    dim_lattice = 10

    for force_freq in force_frequencies:

        ##RUN DYNAMIC TRIAL
        sim = Tri_Lattice()
        sim.create_lattice(dim_lattice,dim_lattice,1,1,1, add_periodic_bonds=False)

        #pinch particles in the middle of the lattice
        f_equation = Force_Eq(0, 1, force_freq, 0)
        neg_f_equation = Force_Eq(0, 1, force_freq, math.pi)

        sim.add_dynamic_force(f_equation, [86], dimension="y")
        sim.add_dynamic_force(neg_f_equation, [88], dimension="y")

        sim.create_animation(f"dynamic_bonds_pinch_{force_freq}")
        sim.log_system_kinetic_energy(f"kinetic_dynamic_pinch_{force_freq}.txt")

        for equation_i in range(len(equations)):
            sim.change_spring_eq(equation_i, equations[equation_i])

        sim.run_langevin(1000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=1000000)

        #store data
        position_file = f"positions_dynamic_pinch_{force_freq}.txt"
        velocity_file = f"velocities_dynamic_pinch_{force_freq}.txt"
        os.rename("positions.txt", position_file)
        os.rename("velocities.txt", velocity_file)
        shutil.copy(position_file, "pinch_sim_data")
        shutil.copy(velocity_file, "pinch_sim_data")
        shutil.copy(f"kinetic_dynamic_pinch_{force_freq}.txt", "pinch_sim_data")
        shutil.copy(f"dynamic_bonds_pinch_{force_freq}.gsd", "pinch_animations")

        #-----------------------------------------------------------------------

        ##RUN STATIC TRIAL
        sim = Tri_Lattice()
        sim.create_lattice(dim_lattice,dim_lattice,1,1,1, add_periodic_bonds=False)

        #pinch particles in the middle of the lattice
        f_equation = Force_Eq(0, 1, force_freq, 0)
        neg_f_equation = Force_Eq(0, 1, force_freq, math.pi)

        sim.add_dynamic_force(f_equation, [86], dimension="y")
        sim.add_dynamic_force(neg_f_equation, [88], dimension="y")

        sim.create_animation(f"static_bonds_pinch_{force_freq}")
        sim.log_system_kinetic_energy(f"kinetic_static_pinch_{force_freq}.txt")

        # for equation_i in range(len(equations)):
        #     sim.change_spring_eq(equation_i, equations[equation_i])

        sim.run_langevin(1000000,callback=sim.log_p_info,gamma=.05, callback_period=500, quiet=True, dynamic_f_stop=1000000)

        #store data
        position_file = f"positions_static_pinch_{force_freq}.txt"
        velocity_file = f"velocities_static_pinch_{force_freq}.txt"
        os.rename("positions.txt", position_file)
        os.rename("velocities.txt", velocity_file)
        shutil.copy(position_file, "pinch_sim_data")
        shutil.copy(velocity_file, "pinch_sim_data")
        shutil.copy(f"kinetic_static_pinch_{force_freq}.txt", "pinch_sim_data")
        shutil.copy(f"static_bonds_pinch_{force_freq}.gsd", "pinch_animations")

def analysis_data(force_mag: int):

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

def pinch_data_analysis(data_folder_name: str):
    """
    filename - name of folder where position files are stored
    """
    #create folders to store created graphs
    path = os.getcwd() + "/pinch_analysis"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path, 0o755)

    #location of target files
    data_files_locations = os.getcwd() + f"/{data_folder_name}"

    force_frequencies = []
    dynamic_files_kinetic = []
    static_files_kinetic = []
    dynamic_files_positions = []
    static_files_positions = []
    dynamic_files_velocities = []
    static_files_velocities = []
    results = []

    #collect files to analyse
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)

        #collect kinetic energy files and record force frequencies
        if filename.startswith("kinetic_dynamic_pinch"):
             fname_parts = filename.split("_")
             force_freq = fname_parts[3].split('t')
             force_freq = float(force_freq[0][:-1])
             if force_freq not in force_frequencies:
                 force_frequencies.append(force_freq)
             dynamic_files_kinetic.append(filename)
             continue
        elif filename.startswith("kinetic_static_pinch"):
            static_files_kinetic.append(filename)

        #collect position files
        elif filename.startswith("positions_static_pinch"):
            static_files_positions.append(filename)
        elif filename.startswith("positions_dynamic_pinch"):
            dynamic_files_positions.append(filename)

        #collect velocity files
        elif filename.startswith("velocities_static_pinch"):
            static_files_velocities.append(filename)
        elif filename.startswith("velocities_dynamic_pinch"):
            dynamic_files_velocities.append(filename)

        else:
            continue

    print(force_frequencies)

    #kinetic analysis
    for i in range(len(dynamic_files_kinetic)):

        static_data = Analysis()
        dynamic_data = Analysis()

        fname_parts = static_files_kinetic[i].split("_")
        force_freq = fname_parts[3].split('t')
        force_freq = float(force_freq[0][:-1])
        print(force_freq)

        dynamic_title = f"Dynamic bond trial, runtime=1000000, force mag = {force_freq}"
        static_title = f"Static bond trial, runtime=1000000, force mag = {force_freq}"
        compair_title = "Kinetic energy percent difference between static and dynamic systems \nforce mag = {force_freq}"


        #read corresponding kinetic energy files
        static_data.read_kinetic(data_files_locations + '/' + static_files_kinetic[i])
        dynamic_data.read_kinetic(data_files_locations + '/' + dynamic_files_kinetic[i])
        # static_data.read_kinetic("kinetic.txt")
        # dynamic_data.read_kinetic("kinetic.txt")

        #read corresponding position files
        static_data.read_pos(data_files_locations + '/' + static_files_positions[i])
        dynamic_data.read_pos(data_files_locations + '/' + dynamic_files_positions[i])

        #read corresponding velocity files
        static_data.read_velocity(data_files_locations + '/' + static_files_velocities[i])
        dynamic_data.read_velocity(data_files_locations + '/' + dynamic_files_velocities[i])

        #create plots of data
        static_energies = static_data.heat_map(static_title, display=False)
        dynamic_energies = dynamic_data.heat_map(dynamic_title, display=False)

        compairison_heatmap(force_freq, static_energies, dynamic_energies)


def compairison_heatmap(force_freq, static_info, dynamic_info):

    compair_title = f"Kinetic energy percent difference between static and dynamic systems \nforce mag = {force_freq}"

    #compairison of static and dynamic trial
    E_k_diffs = []
    for energy_i in range(len(static_info[2])):
        #diff = abs(dynamic_energies[2][energy_i] - static_energies[2][energy_i])
        diff = (dynamic_info[2][energy_i] - static_info[2][energy_i]) / static_info[2][energy_i] * 100
        E_k_diffs.append(diff)

    plt.clf()
    H, xe, ye = np.histogram2d(static_info[0], static_info[1], [75,75], weights=E_k_diffs)
    H = H.T
    ax1 = plt.gca()
    plot = ax1.imshow(H, interpolation='nearest', origin='low', extent=[xe[0]*2.5, xe[-1]*2.5, ye[0], ye[-1]], cmap='Blues')
    ax1.set_xticks([])
    ax1.set_yticks([])

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    ax1.figure.colorbar(plot, cax=cax)

    plt.suptitle(compair_title)
    plt.savefig(f"Heatmap_compair_{force_freq}.png")
    shutil.copy(f"Heatmap_compair_{force_freq}.png", "pinch_analysis")
