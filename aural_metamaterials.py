"""
Author - Noah Kruss

Functions for running the aural_metamaterial simulations
"""

#---------------IMPORT STATEMENTS------------------------------
import sys
import numpy as np
import os
import shutil

from Equations import *
from aural_analysis import *

sys.path.append('./simulation_setup')
from line_variable_bonds import *
sys.path.append('./..')


#---------------Function to Run Simulation---------------------
def run_sim(num_bonds: int, condition_name: str):
    """
    Function for running a simulation on a spring-mass chain with (num_bonds)
    number of dynamic bonds within the unit cell with a complex gaussian
    wavepacket initilized within the system.

    Inputs:
        num_bonds - number of bonds within the unit cell of the spring mass chain
        condition_name - name of directory within initial_condition folder in
                         aural_metamaterial directory with the inital position
                         and velocity values of the spring-mass chain

    """

    #initilize a directory for the simulation data
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
    dk = .36
    c_0 = math.sqrt(k_0 / m)
    w = 2 * math.pi * c_0 / num_bonds / a
    dt = .0001
    damping = 0.0225
    runtime = int(10000 / .0001)

    #add a file with the simulations run conditions
    f = open("run_conditions.txt", "w+")
    damp_cond = f"Damping - mag = {damping}\n"
    spring_cond = f"Springs - number of bonds = {num_bonds}, dk = {dk}, k_0 = {k_0}, w = {w}, a = {a} \n"
    time_cond = f"Simulation runtime of {runtime}, data_recorded every 500 timesteps \n"
    f.write(time_cond + spring_cond + damp_cond)
    f.close()
    shutil.copy("run_conditions.txt", main_path)

    #Initilize simulation system
    sim = Line()
    sim.create_lattice(1,
                        num_bonds = num_bonds,
                        add_periodic_bonds = True,
                        N = N,
                        pos_offsets = pos_offsets,
                        vel_init = vel_init)
    sim.dimension_constrain([1,0,0])

    #create logs for gsd file
    sim.create_animation(f"line_test_{num_bonds}bonds")
    sim.log_system_kinetic_energy(f"kinetic_{num_bonds}bonds.txt", period=500)

    #change bonds to be dynamic
    if dk != 0:
        equations = get_spring_equations_2(k_0, dk, w, num_bonds)
        sim.Equation_list = equations
        for spring_i in range(len(sim.Equation_list)):
            k = sim.Equation_list[spring_i].calc(0)
            sim.k_list[spring_i] = k

    #run simulation
    sim.run_variable_bonds(runtime,
                           kt=0,
                           dt=dt,
                           gamma=damping)

    #store data into simulation directory
    shutil.copy(f"line_test_{num_bonds}bonds.gsd", main_path)
    shutil.copy(f"kinetic_{num_bonds}bonds.txt", main_path)


def gsd_multi_analysis():
    #list of simulation folders to run analysis on
    dir_list = ["8bonds_damp_pt01925", "8bonds_damp_pt01925_dt_pt00025", "8bonds_damp_pt01925_dt_pt0005",
                "8bonds_damp_pt019", "8bonds_damp_pt019_dt_pt00025", "8bonds_damp_pt019_dt_pt0005",
                "8bonds_damp_pt0195", "8bonds_damp_pt0195_dt_pt00025", "8bonds_damp_pt0195_dt_pt0005",
                ]

    #
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
    data.gaussian_fitting(dt, store_loc = path, num_samples = 1000)
    #data.RMS_error(dt, store_loc = path, plot_title = mean_error_plot_title, num_samples = 1000)
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
        os.remove(os.getcwd() + "/aural_metamaterial/Gaussian_Fit_Center_All_Data.png")
        os.remove(os.getcwd() + "/aural_metamaterial/Gaussian_Fit_Amp_All_Data.png")
    except:
        pass

    dir_list = ["static_bonds",
                "8bonds_damp_pt019",
                "8bonds_damp_pt019175",
                "8bonds_damp_pt0192",
                "8bonds_damp_pt01925",
                "8bonds_damp_pt02"]
    label_list = ["Static",
                  "0.019000",
                  "0.019175",
                  "0.019200",
                  "0.019250",
                  "0.020000"]
    dash_list = [[1,0],
                [6,1],
                [3,1],
                [1,1],
                [3,4],
                [6,4]]

    gaussian_results = []
    xs = []
    ys = []
    yints = []

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

        xs.append(gaussian_data[0])
        ys.append(gaussian_data[4])
        yints.append(gaussian_data[4][500])

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
    plt.xlabel('Time')
    plt.ylabel('Fit Amplitude Factor')
    plt.title("Gaussian Fitting Amplitude Factor")

    #plot static trial
    # data = gaussian_results[0]
    # plt.plot(data[0], data[4], c = "k", label=label_list[i])

    cmap = plt.get_cmap('winter')
    yints = np.asarray(yints)
    for x, y, color, label, dash in zip(xs, ys, normalize(yints), label_list, dash_list):
        plt.plot(x, y, label=label, color=cmap(color), dashes=dash)

    plt.legend(loc = 'upper left', fontsize = 9)
    #plt.figure(figsize=fig_size)
    plt.save('Gaussian_Fit_Amp_All_Data.pdf')
    shutil.move("Gaussian_Fit_Amp_All_Data.pdf", os.getcwd() + f"/aural_metamaterial")
