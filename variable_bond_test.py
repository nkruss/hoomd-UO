from analysis import *
from Equations import *
from line_variable_bonds import *
from aural_analysis import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import shutil
import random

def run_sim(num_bonds: int, condition_name: str, alt_damp = False):
    directory = os.getcwd() + "/aural_metamaterial"
    main_path = directory + f"/{num_bonds}bonds"
    try:
        shutil.rmtree(main_path)
    except:
        pass
    os.mkdir(main_path, 0o755)

    #get initial position offsets and initial velocities
    pos_offsets = []
    vel_init = []

    pos_file = open(directory + "/initial_conditions" + f"/{condition_name}/InitPos.txt", "r")
    vel_file = open(directory + "/initial_conditions" + f"/{condition_name}/InitVel.txt", "r")

    for line in pos_file:
        offset = float(line.strip())
        pos_offsets.append(offset)
    for line in vel_file:
        vel = float(line.strip())
        vel_init.append([vel, 0, 0])

    pos_file.close()
    vel_file.close()

    #define run conditions variables
    N = len(vel_init)
    print(f"N = {N}");
    k_0 = 1
    m = 1
    a = 1
    # dk = .18
    dk = .09
    c_0 = math.sqrt(k_0 / m)
    w = 2 * math.pi * c_0 / num_bonds / a
    #damping = 0.01926
    dt = .0001
    damping = 0.01584
    runtime = int(10000 / .0001)
    #runtime = 100000000
    damp_cond = f"Damping - mag = {damping}, alternating damping = {alt_damp} \n"
    spring_cond = f"Springs - number of bonds = {num_bonds}, dk = {dk}, k_0 = {k_0}, w = {w}, a = {a} \n"
    time_cond = f"Simulation runtime of {runtime}, data_recorded every 500 timesteps \n"

    #add a file with the simulations run conditions
    f = open("run_conditions.txt", "w+")
    f.write(time_cond + spring_cond + damp_cond)
    f.close()
    #shutil.copyfile(os.getcwd() + "/run_conditions.txt", os.getcwd() + "/run_conditions.txt")
    shutil.copy("run_conditions.txt", main_path)

    #create and run simulation
    sim = Line()
    sim.create_lattice(1,
                        num_bonds = num_bonds,
                        add_periodic_bonds = True,
                        N = N,
                        pos_offsets = pos_offsets,
                        vel_init = vel_init)

    sim.dimension_constrain([1,0,0])
    sim.create_animation(f"line_test_{num_bonds}bonds")
    sim.log_system_kinetic_energy(f"kinetic_{num_bonds}bonds.txt", period=500)

    #change bonds to be dynamic
    if dk != 0:
        equations = get_spring_equations_2(k_0, dk, w, num_bonds)
        sim.Equation_list = equations
        for spring_i in range(len(sim.Equation_list)):
            k = sim.Equation_list[spring_i].calc(0)
            sim.k_list[spring_i] = k

    sim.run_variable_bonds(runtime,
                           kt=0,
                           dt=dt,
                           gamma=damping,
                           alt_damp = alt_damp)

    #store data into folder
    # position_file = f"positions_1Dline_{num_bonds}bonds.txt"
    # velocity_file = f"velocities_1Dline_{num_bonds}bonds.txt"
    # os.rename("positions.txt", position_file)
    # os.rename("velocities.txt", velocity_file)
    # shutil.copy(position_file, main_path)
    # shutil.copy(velocity_file, main_path)
    shutil.copy(f"line_test_{num_bonds}bonds.gsd", main_path)
    shutil.copy(f"kinetic_{num_bonds}bonds.txt", main_path)

def data_analysis(direct_name):

    #create folders to store created graphs
    data_files_locations = os.getcwd() + f"/aural_metamaterial/{direct_name}"
    path = data_files_locations + "/analysis"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path, 0o755)

    #get data files
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)

        #collect kinetic energy files and record force frequencies
        if filename.startswith("kinetic"):
             kinetic_file = filename

        #collect position files
        elif filename.startswith("positions"):
            pos_file = filename

        #collect velocity files
        elif filename.startswith("velocities"):
            vel_file = filename

        #collect run condition file
        elif filename.startswith("run"):
            cond_file = data_files_locations + '/' + filename

    with open(cond_file, "r") as f:
        time_cond = f.readline().strip()
        spring_cond = f.readline().strip()
        damp_cond = f.readline().strip()
        f.close()

    kinetic_plot_title = f"Kinetic Energy Plot\n{damp_cond}\n{spring_cond}\n"
    waterfall_plot_title = f"Gaussian Plot\n{damp_cond}\n{spring_cond}\n"

    data = Analysis()
    #read corresponding position files
    data.read_pos(data_files_locations + '/' + pos_file)
    print("Read pos file")

    #read corresponding velocity files
    data.read_velocity(data_files_locations + '/' + vel_file)
    print("Read vel file")

    #read corresponding kinetic energy files
    data.read_kinetic(data_files_locations + '/' + kinetic_file)
    print("Read kinetic file")

    data.graph_k_energy(store_loc = path, plot_title = kinetic_plot_title)
    data.wave_packet(store_loc = path, plot_title = waterfall_plot_title, num_samples = 10)


def gsd_multi_analysis():
    dir_list = ["8bonds_damp_pt01925", "8bonds_damp_pt01925_dt_pt00025", "8bonds_damp_pt01925_dt_pt0005",
                "8bonds_damp_pt019", "8bonds_damp_pt019_dt_pt00025", "8bonds_damp_pt019_dt_pt0005",
                "8bonds_damp_pt0195", "8bonds_damp_pt0195_dt_pt00025", "8bonds_damp_pt0195_dt_pt0005",
                ]

    count = 0
    for dir_name in dir_list:
        if count in [0, 3, 6]:
            gsd_analysis(dir_name, dt = 0.0001)
        elif count in [1, 4, 7]:
            gsd_analysis(dir_name, dt = 0.00025)
        else:
            gsd_analysis(dir_name, dt = 0.0005)

        count += 1


def gsd_analysis(dir_name, dt = 0.0001):

    #create folders to store created graphs
    data_files_locations = os.getcwd() + f"/aural_metamaterial/{dir_name}"
    path = data_files_locations + "/analysis"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path, 0o755)

    #get data files
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)

        #collect kinetic energy files and record force frequencies
        if filename.startswith("kinetic"):
             kinetic_file = filename

        #collect position files
        elif filename.startswith("positions"):
            pos_file = filename

        #collect velocity files
        elif filename.startswith("velocities"):
            vel_file = filename

        #collect run condition file
        elif filename.startswith("run"):
            cond_file = data_files_locations + '/' + filename

        #collect run condition file
        elif filename.startswith("line"):
            gsd_file = data_files_locations + '/' + filename

    with open(cond_file, "r") as f:
        time_cond = f.readline().strip()
        spring_cond = f.readline().strip()
        damp_cond = f.readline().strip()
        f.close()

    #define plot titles
    kinetic_plot_title = f"Kinetic Energy Plot\n{damp_cond}\n{spring_cond}\n"
    waterfall_plot_title = f"Gaussian Plot\n{damp_cond}\n{spring_cond}\n"
    mean_error_plot_title = f"RMS Error (Over Simulation)\n{damp_cond}\n{spring_cond}\n"
    normalized_error_plot_title = f"Normalised Fourier Amplitude Difference (Over Simulation)\n{damp_cond}\n{spring_cond}\n"

    data = Aural_Analysis()
    data.read_data(gsd_file)

    # #target times
    target_times = [0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 9999.95]
    # data.wave_packet(dt, store_loc = path, plot_title = waterfall_plot_title, num_samples = 10, target_times = target_times)
    # data.RMS_error(dt, store_loc = path, plot_title = mean_error_plot_title, num_samples = 10, target_times = target_times)
    # data.fractional_error(dt, store_loc = path, plot_title = fractional_error_plot_title, num_samples = 10, target_times = target_times)
    # data.peak_error(dt, store_loc = path, num_samples = 10, target_times = target_times)

    #target times
    #target_times = [0.0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400]
    #target_times = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35]
    data.wave_packet(dt, store_loc = path, plot_title = waterfall_plot_title, num_samples = 10, target_times = target_times)
    #data.gaussian_fitting(dt, store_loc = path, num_samples = 1000, target_times = target_times)
    #data.gaussian_fitting(dt, store_loc = path, num_samples = 1000)
    data.RMS_error(dt, store_loc = path, plot_title = mean_error_plot_title, num_samples = 1000)
    data.normalized_error(dt, store_loc = path, plot_title = normalized_error_plot_title, num_samples = 1000)
    print("---")
    #data.peak_error(dt, store_loc = path, num_samples = 1000)

def integrety_analysis(dir_name, dt = 0.0001):

    #create folders to store created graphs
    data_files_locations = os.getcwd() + f"/aural_metamaterial/{dir_name}"
    path = data_files_locations + "/analysis"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path, 0o755)

    #get data files
    for file in sorted(os.listdir(data_files_locations)):
        filename = os.fsdecode(file)

        #collect kinetic energy files and record force frequencies
        if filename.startswith("kinetic"):
             kinetic_file = filename

        #collect run condition file
        elif filename.startswith("run"):
            cond_file = data_files_locations + '/' + filename

        #collect run condition file
        elif filename.startswith("line"):
            gsd_file = data_files_locations + '/' + filename

    #get static gsd file
    static_files_locations = os.getcwd() + f"/aural_metamaterial/static_bonds"
    for file in sorted(os.listdir(static_files_locations)):
        filename = os.fsdecode(file)

        #collect run condition file
        if filename.startswith("line"):
            static_gsd_file = static_files_locations + '/' + filename

    with open(cond_file, "r") as f:
        time_cond = f.readline().strip()
        spring_cond = f.readline().strip()
        damp_cond = f.readline().strip()
        f.close()

    #target times
    target_times = [0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 9999.95]

    #set up analyis data objects
    data = Aural_Analysis()
    data.read_data(gsd_file)
    static_data = Aural_Analysis()
    static_data.read_data(static_gsd_file)

    #calculate integrety values
    # run_integrety_data = data.integrety_test(dt, num_samples = 10, target_times = target_times)
    # baseline_integrety_data = static_data.integrety_test(dt, num_samples = 10, target_times = target_times)

    run_integrety_data = data.normalized_error(dt, num_samples = 1000)
    baseline_integrety_data = static_data.normalized_error(dt, num_samples = 1000)
    integrety_value_list = []
    for i in range(len(run_integrety_data[0])):

        integrety_value = abs(run_integrety_data[0][i] - baseline_integrety_data[0][i])
        integrety_value_list.append(integrety_value)

    # print("integrety_value_list")
    # print(integrety_value_list)

    #plot integrety values
    plot_title = "Wave Packet Integrety - (over course of Simulation)"
    plt.plot(run_integrety_data[1], integrety_value_list)
    plt.title(plot_title, fontsize = 7)
    plt.savefig(plot_title)
    shutil.move(f"{plot_title}.png", path)
    plt.clf()

def Gaussian_center_fit_analysis(dt = 0.0001):
    try:
        os.getcwd() + "/aural_metamaterial/Gaussian_Fit_Center_All_Data.png"
        os.getcwd() + "/aural_metamaterial/Gaussian_Fit_Amp_All_Data.png"
    except:
        pass

    dir_list = ["static_bonds",
                "8bonds_damp_pt019",
                "8bonds_damp_pt019175",
                "8bonds_damp_pt0192",
                "8bonds_damp_pt01925",
                "8bonds_damp_pt02"]
    label_list = ["Static Bonds",
                "Damping 0.019000",
                "Damping 0.019175",
                "Damping 0.019200",
                "Damping 0.019250",
                "Damping 0.020000"]
    # fmt_list = ["k-",
    #             "m:",
    #             "g:",
    #             "b-",
    #             "g--",
    #             "m--"]

    gaussian_results = []

    for dir_name in dir_list:

        #create folders to store created graphs
        data_files_locations = os.getcwd() + f"/aural_metamaterial/{dir_name}"
        path = data_files_locations + "/analysis"
        try:
            shutil.rmtree(path)
        except:
            pass
        os.mkdir(path, 0o755)

        #get data files
        for file in sorted(os.listdir(data_files_locations)):
            filename = os.fsdecode(file)

            #collect kinetic energy files and record force frequencies
            if filename.startswith("kinetic"):
                 kinetic_file = filename

            #collect position files
            elif filename.startswith("positions"):
                pos_file = filename

            #collect velocity files
            elif filename.startswith("velocities"):
                vel_file = filename

            #collect run condition file
            elif filename.startswith("run"):
                cond_file = data_files_locations + '/' + filename

            #collect run condition file
            elif filename.startswith("line"):
                gsd_file = data_files_locations + '/' + filename

        with open(cond_file, "r") as f:
            time_cond = f.readline().strip()
            spring_cond = f.readline().strip()
            damp_cond = f.readline().strip()
            f.close()

        data = Aural_Analysis()
        data.read_data(gsd_file)

        target_times = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
        waterfall_plot_title = f"Gaussian Plot\n{damp_cond}\n{spring_cond}\n"
        data.wave_packet(dt, store_loc = path, plot_title = waterfall_plot_title, target_times = target_times)

        gaussian_data = data.gaussian_fitting(dt, store_loc = path, num_samples = 1000)
        gaussian_results.append(gaussian_data)
        #plt.plot(center_data[0], np.unwrap(center_data[1]))

    #plot all center data
    g = plt.figure(1)
    plt.clf()
    plt.xlabel('Dimensionless Time')
    plt.ylabel('Center Position')
    plt.title("Gaussian Fit Center Position")
    for i in range(len(gaussian_results)):
        data = gaussian_results[i]
        plt.plot(data[0], data[1], label=label_list[i])
        print(f"{dir_list[i]}, speed = {data[2]}")
    plt.legend()
    plt.savefig("Gaussian_Fit_Center_All_Data")
    shutil.move("Gaussian_Fit_Center_All_Data.png", os.getcwd() + f"/aural_metamaterial")

    #plot all amp data
    plt.clf()
    plt.xlabel('Dimensionless Time')
    plt.ylabel('Fit Amplitude Percent Error')
    plt.title("Gaussian Fitting Amplitude Error")
    #plot static trial
    data = gaussian_results[0]
    plt.plot(data[0], data[4], c = k, label=label_list[i])
    #plot rest of the data
    colors = pl.cm.jet(np.linspace(0,1,len(dir_list) - 1))
    for i in range(1, len(gaussian_results)):
        data = gaussian_results[i]
        # plt.plot(data[0], data[4], fmt_list[i], label=label_list[i])
        plt.plot(data[0], data[4], color=colors[i-1], label=label_list[i], cmap="coolwarm")
    plt.legend(loc = 'upper left', fontsize = 'small')
    plt.savefig("Gaussian_Fit_Amp_All_Data")
    shutil.move("Gaussian_Fit_Amp_All_Data.png", os.getcwd() + f"/aural_metamaterial")
