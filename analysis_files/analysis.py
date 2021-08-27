import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import os
import shutil
import statistics as stats
import scipy.signal
import scipy.special

def avg(quantity: list, time: list, start=0, finish=-1):
    '''
    Calculates the mean value of the numbers in a list (quantity) over a some number of timesteps,
    by using the sum function to get a sum of the list and then
    deviding that by the number of timesteps

    '''
    if quantity == []: #checks for an empty list
        return 0.
    #print(start, finish)
    #print(quantity[start:finish])
    quantitySum = 0
    if finish == -1:
        finish = len(quantity)
        for energy in quantity[start:]:
            quantitySum += energy
    else:
        for energy in quantity[start:finish+1]:
            quantitySum += energy
    #print(time[finish], time[start])
    #timeSum = time[finish] - time[start]
    timeSum = finish - start
    #print(f"timeSum = {timeSum}, quantitySum={quantitySum}")
    mean = quantitySum/timeSum
    return mean

def kinetic_compare(static_data: "Analysis", dynamic_data: "Analysis"):
    """Comparison function for kinetic energy. Takes two different Analysis objects as parameters
    and then compairs the kinetic energy data stored in the objects and creates a plot of the ratio
    between the kinetic energies for each dynamic force frequency"""

    static_averages = static_data.average_kinetic
    dynamic_averages = dynamic_data.average_kinetic

    print(f"static averages = {static_averages}")
    print(f"dynamic averages = {dynamic_averages}")


    ratio_list = []
    for i in range(len(static_averages)):
        static_avg = static_averages[i]
        dynamic_avg = dynamic_averages[i]
        ratio = (dynamic_avg / static_avg)
        ratio_list.append(abs(ratio))

    colors = []

    g1 = plt.figure(1)
    plt.xlabel('angular frequency')
    plt.ylabel('average kinetic energy comparison')
    print(ratio_list)
    plt.scatter(static_data.dynamic_frequencies, ratio_list, c='b')
    plt.show()

def kinetic_bin_compare(data1: "Analysis", data2: "Analysis", force_frequency: int, save_location=" "):

    data1_bin_index = []
    data2_bin_index = []

    data1_bin_k = []
    data2_bin_k = []

    #data1
    for bin_num in range(len(data1.bins)):
        data1_bin_index.append(bin_num)
        data1_bin_k.append(data1.bins[bin_num][1])

    #data2
    for bin_num in range(len(data2.bins)):
        data2_bin_index.append(bin_num)
        data2_bin_k.append(data2.bins[bin_num][1])

    plot_title = f"Kinetic energy of triangular lattice at different bins, force frequency of {math.sqrt(5000)},\n spring frequency = {math.sqrt(5000)}, spring mag = 5000 ± 2500, dynamic force frequency = {force_frequency}"

    g1 = plt.figure(1)
    plt.xlabel('bin')
    plt.ylabel('kinetic energy')
    plt.scatter(data1_bin_index, data1_bin_k, c="b", label="static")
    plt.scatter(data2_bin_index, data2_bin_k, c="r", label="dynamic")
    plt.legend(loc='upper right')
    plt.title(plot_title, fontsize=10)
    #plt.show()
    plt.savefig(f"binplot_force_freq{force_frequency}.pdf")
    if save_location != " ":
        shutil.copy(f"binplot_force_freq{force_frequency}.pdf", save_location)
    plt.clf()

def static_dynamic_analysis(graph_type: str):
    data = np.loadtxt("static_dynamic_compair.npy")
    print(data)
    #print(data[:,2])
    # x = [0,0,1,1]
    # y = [0,0,1,1]
    # z = [1,1,1,1]

    if graph_type == "contour":
        spring_freq_list = []
        for freq in data[:,0]:
            if freq not in spring_freq_list:
                spring_freq_list.append(freq)

        force_freq_list = []
        for freq in data[:,1]:
            if freq not in force_freq_list:
                force_freq_list.append(freq)

        ratios = data[:,2]
        ratio_counter = 0
        array = np.zeros((len(spring_freq_list), len(force_freq_list)))

        for spring_freq_i in range(len(spring_freq_list)):
            for force_freq_i in range(len(force_freq_list)):
                array[spring_freq_i, force_freq_i] = ratios[ratio_counter]
                ratio_counter += 1

        X,Y = np.meshgrid(data[:,0],data[:,1])
        plt.contourf(spring_freq_list,force_freq_list, array)
        plt.colorbar()
        plt.xlabel('force_freq')
        plt.ylabel('spring_freq')
        plt.title('Kinetic energy difference')
        plt.show()

    if graph_type == "skyscraper":
        spring_freq_list = []
        for freq in data[:,0]:
            if freq not in spring_freq_list:
                spring_freq_list.append(freq)

        force_freq_list = []
        for freq in data[:,1]:
            if freq not in force_freq_list:
                force_freq_list.append(freq)

        ratios = data[:,2]
        ratio_counter = 0
        array = np.zeros((len(spring_freq_list), len(force_freq_list)))

        for spring_freq_i in range(len(spring_freq_list)):
            for force_freq_i in range(len(force_freq_list)):
                array[spring_freq_i, force_freq_i] = ratios[ratio_counter]
                ratio_counter += 1

        X,Y = np.meshgrid(force_freq_list,spring_freq_list)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,array)

        plt.xlabel('force_freq')
        plt.ylabel('spring_freq')
        plt.title('Kinetic energy difference')
        plt.show()

    elif graph_type == "scatter":

        spring_freq_list = []
        for freq in data[:,0]:
            spring_freq_list.append(freq)

        force_freq_list = []
        for freq in data[:,1]:
            force_freq_list.append(freq)

        ratio_list = []
        for ratio in data[:,2]:
            ratio_list.append(ratio)

        ratios = data[:,2]
        ratio_counter = 0
        array = np.zeros((len(spring_freq_list), len(force_freq_list)))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection= "3d")
        ax.scatter3D(spring_freq_list, force_freq_list, ratio_list, c='r', marker='o')

        ax.set_xlabel('Spring frequency')
        ax.set_ylabel('Force frequency')
        ax.set_zlabel('Kinetic energy ratio')

        # plt.savefig(f"static_dynamic_compair.pdf")
        # plt.clf()
        fig.show()


class Analysis():

    def __init__(self):
        #self.time is in timesteps
        self.time = []
        self.realTime = []
        self.particle_v = {}
        self.particle_p = {}
        self.stiffness = {}
        self.dt = None
        self.m = None
        self.N = None
        self.kspace_vec = None
        self.wspace_vec = None
        self.kinetic_energy = []
        self.dynamic_frequencies = []
        self.spring_frequencies = []
        self.average_kinetic = []
        self.spring_frequency = None
        self.stabalization_step = 0

# Read Functions
    def read_pos(self, fname: str):
        """Read through a simulation position file and record all of the particle
        positions into a dictionary of tuples and store the simulation paramaters as analysis object porperties"""

        #open the position file
        with open(fname, "r") as f:
            #recorded the simulation general infomation
            self.time = []
            self.realTime = []
            line = f.readline().strip().split()
            self.dt = float(line[5])
            self.m = float(line[7])
            self.N = int(line[-1])
            try:
                self.y = int(line[-3])
                self.x = int(line[-5])
            except:
                pass
            #create a key in a dictionary conected to a list for every particle in the simulation
            for i in range(self.N):
                self.particle_p[i] = []

            f.readline()

            #go through each line in the file
            for line in f:
                #split the line up into a list of the timestep and each particles position
                data = line.strip().split("(")
                timesteps = data[0].strip()
                #record the lines timestep
                self.time.append(int(timesteps))
                self.realTime.append(int(timesteps)*self.dt)

                #for each particle record it's position as a tuple to the list stored in the dictionary
                #of particles
                for item_i in range(1, len(data)):
                    p_data = data[item_i].strip()
                    p_data = p_data[:-1]
                    p_data = p_data.split(',')
                    self.particle_p[item_i - 1].append((float(p_data[0]),float(p_data[1]),float(p_data[2])))
        f.close()

    def read_velocity(self, fname: str):
        """Read through a simulation velocity file and record all of the particle
        multi-dimensional velociy components into a dictionary of tuples and
        store the simulation paramaters as analysis object porperties"""

        with open(fname, "r") as f:
            self.time = []
            self.realTime = []
            line = f.readline().strip().split()
            self.dt = float(line[5])
            self.m = float(line[7])
            self.N = int(line[-1])
            try:
                self.y = int(line[-3])
                self.x = int(line[-5])
            except:
                pass
            for i in range(self.N):
                self.particle_v[i] = []

            f.readline()

            for line in f:
                data = line.strip().split("(")
                timesteps = data[0].strip()
                self.time.append(int(timesteps))
                self.realTime.append(int(timesteps)*self.dt)

                for item_i in range(1, len(data)):
                    p_data = data[item_i].strip()
                    p_data = p_data[:-1]
                    p_data = p_data.split(',')
                    self.particle_v[item_i - 1].append((float(p_data[0]),float(p_data[1]),float(p_data[2])))
        f.close()

    def read_force(self, fname: str):
        """Read through a simulation dynamic force file and record all of
        the dynamic force magnitudes into a list"""

        with open(fname, "r") as f:
            self.dy_time = []
            self.force = []

            f.readline()
            f.readline()

            for line in f:
                data = line.strip().split()
                #print(data)
                time = data[0].strip()
                self.dy_time.append(float(time))
                dynamic_force = data[1].strip()
                self.force.append(float(dynamic_force))

    def read_kinetic(self, fname: str):
        """Read through a simulation kinetic energy file and record all of
        the enrgy magnitudes into a list"""

        run_num = len(self.kinetic_energy)
        with open(fname, "r") as f:
            self.kinetic_energy = []
            self.time = []
            self.realTime = []

            f.readline()

            for line in f:
                data = line.strip().split()
                #print(data)
                time = data[0].strip()
                self.time.append(float(time))
                self.realTime.append(float(time) * self.dt)
                kinetic_energy = data[1].strip()
                self.kinetic_energy.append(float(kinetic_energy))
        f.close()

        count = 0
        avg1 = 0
        avg2 = 0
        self.last_step = 0

        for timestep_i in range(1,len(self.time)):
            if self.time[timestep_i] % 1000 == 0:
                avg2 = avg(self.kinetic_energy, self.time, start=self.last_step, finish=timestep_i)
                #diff = abs(avg2 - avg1) / avg2
                diff = avg1 / avg2
                #print("--", diff, run_time[timestep_i], "avg1 =", avg1, "avg2=", avg2)
                if abs(1 - diff) < .105:
                    count += 1
                else:
                    count = 0
                if count == 5:
                    break
                avg1 = avg2
                self.last_step = timestep_i

        #print("stabalization cutoff =", run_time[self.last_step])

        average = avg(self.kinetic_energy, self.time, start=self.last_step, finish=-1)

        self.stabalization_step = self.last_step

        #print(average, run_time[(last_step-9)])
        #average = avg(run_kinetic_energy, run_time)

        #print("kinetic energy", run_kinetic_energy)
        self.average_kinetic.append(average)
        #print("Hoomd average =", average)

    def read_velocity_mult(self, fname: str):

        with open(fname, "r") as f:
            self.time = []
            line = f.readline().strip().split()
            self.dt = float(line[5])
            ### self.m = float(line[7])
            self.N = int(line[-1])
            for i in range(self.N):
                self.particle_v[i] = []

            f.readline()

            for line in f:
                data = line.strip().split("(")
                timesteps = data[0].strip()
                self.time.append(int(timesteps))

                for item_i in range(1, len(data)):
                    p_data = data[item_i].strip()
                    p_data = p_data[:-1]
                    p_data = p_data.split(',')
                    self.particle_v[item_i - 1].append((float(p_data[0]),float(p_data[1]),float(p_data[2])))

    def read_stiffnesses(self, fname: str):

        with open(fname, "r") as f:
            self.time = []
            line = f.readline().strip().split()
            self.dt = float(line[5])
	        #self.m = float(line[7])
            self.N = int(line[-1])
            for i in range(3):
                self.stiffness[i] = []

                f.readline()

                for line in f:
                    data = line.strip().split()
                    #print('data is ' , data)
                    timesteps = data[0].strip()
                    self.time.append(int(timesteps))

                    for item_i in range(1, len(data)):
                        p_data = data[item_i].strip()
	                    #print('p_data is ' , p_data)
	                    #p_data = p_data[:-1]
	                    #print('now p_data is ' , p_data)
	                    #p_data = p_data.split(',')
                        self.stiffness[item_i - 1].append(float(p_data))

    def kinetic_bin(self):
        self.bins = {}

        #create a bin for each collomn of particles (bin contians a tuple with list of particle indicies 0, and avg_k energy pos 2)
        for num in range(self.x):
            self.bins[num] = [[], 0]

        offset = 0
        for col_i in range(self.x):
            for p_i in range(2 * self.y):
                self.bins[col_i][0].append(p_i + offset)
            offset += (2 * self.y)

        for time_i in range(len(self.time[self.last_step:])):
            for bin_num in range(len(self.bins)):
                for p_i in self.bins[bin_num][0]:
                    velocity = self.particle_v[p_i][time_i + self.last_step]
                    V_mag_squared = (velocity[0]**2) + (velocity[1]**2) + (velocity[2]**2)
                    E_k = .5 * self.m * V_mag_squared
                    self.bins[bin_num][1] += E_k

        for bin_num in range(len(self.bins)):
            self.bins[bin_num][1] = self.bins[bin_num][1] / (self.time[-1] - self.time[self.last_step])

        bin_index = []
        bin_k = []

        avgk_total = 0
        for bin_num in range(len(self.bins)):
            avgk_total += self.bins[bin_num][1]

        print("my bin average =", avgk_total)
        #print(self.bins)

        for bin_num in range(len(self.bins)):
            bin_index.append(bin_num)
            bin_k.append(self.bins[bin_num][1])

        """
        plot_title = f"Kinetic energy of triangular lattice at different bins, force frequency of {100},\n Dynamic springs"

        g1 = plt.figure(1)
        plt.xlabel('bin')
        plt.ylabel('kinetic energy')
        plt.scatter(bin_index, bin_k)
        plt.title(plot_title, fontsize=10)
        plt.show()
        """

    def my_kinetic(self):
        #now matches the Hoomd. Therefore now need to change bins to match

        k_total = 0
        run_kinetic_energy = []
        for time_i in range(len(self.time)):
            run_kinetic_energy.append(0)
            for p_i in range(self.N):
                velocity = self.particle_v[p_i][time_i]
                V_mag_squared = (velocity[0]**2) + (velocity[1]**2) + (velocity[2]**2)
                E_k = .5 * self.m * V_mag_squared
                run_kinetic_energy[time_i] += E_k

        count = 0
        avg1 = 0
        avg2 = 0
        self.last_step = 0

        for timestep_i in range(1,len(self.time)):
            if self.time[timestep_i] % 1000 == 0:
                avg2 = avg(run_kinetic_energy, self.time, start=self.last_step, finish=timestep_i)
                diff = avg1 / avg2
                #print("--", diff, run_time[timestep_i])
                if abs(1 - diff) < .105:
                    count += 1
                else:
                    count = 0
                if count == 10:
                    break
                avg1 = avg2
                self.last_step = timestep_i

        print(self.last_step)

        #print("kinetic energy", run_kinetic_energy)
        #print(self.time)


        average = avg(run_kinetic_energy, self.time, start=self.last_step, finish=-1)

        #average = 0
        """
        for energy in run_kinetic_energy[self.last_step:]:
            average += energy
        average = average / (self.time[-1] - self.time[self.last_step])
        """

        #print(average, run_time[(last_step-9)])
        #average = avg(run_kinetic_energy, run_time)
        self.average_kinetic.append(average)
        print("average =", average)
        """
        k_total = 0
        for p_i in range(self.N):
            #need to pull out the velocities only after self.last_step
            for velocity in (self.particle_v[p_i][self.last_step:]):
                V_mag_squared = (velocity[0]**2) + (velocity[1]**2) + (velocity[2]**2)
                E_k = .5 * self.m * V_mag_squared
                k_total += E_k
        avgk_total = k_total / ((self.time[-1] - self.time[self.last_step]))

        print("my total average =", avgk_total)
        """

#Graphing Functions

    def heat_map(self, plot_title, display=True): #need to figure out how to pin edges of the lattice
        x_pos_list = []
        y_pos_list = []
        E_k_list = []

        for p_i in range(self.N):
            current_x = []
            current_y = []
            current_E = []


            for time_i in range(self.stabalization_step, len(self.time)):
                current_x.append(self.particle_p[p_i][time_i][0])
                current_y.append(self.particle_p[p_i][time_i][1])

                velocity = self.particle_v[p_i][time_i]
                V_mag_squared = (velocity[0]**2) + (velocity[1]**2) + (velocity[2]**2)
                current_E.append(.5 * self.m * V_mag_squared)

            avg_x_pos = stats.mean(current_x)
            avg_y_pos = stats.mean(current_y)
            avg_E = stats.mean(current_E)

            x_pos_list.append(avg_x_pos)
            y_pos_list.append(avg_y_pos)
            E_k_list.append(avg_E)


        H, xe, ye = np.histogram2d(x_pos_list, y_pos_list, [75,75], weights=E_k_list)
        H = H.T
        ax1 = plt.gca()
        plot = ax1.imshow(H, interpolation='nearest', origin='low', extent=[xe[0]*2.5, xe[-1]*2.5, ye[0], ye[-1]], cmap='Blues')
        ax1.set_xticks([])
        ax1.set_yticks([])

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        ax1.figure.colorbar(plot, cax=cax)
        plt.suptitle(plot_title)

        if display:
            plt.show()

        return((x_pos_list, y_pos_list, E_k_list))

    def graph_spring_line_particle_kinetic(self, plot_title, store_loc):
        E_k_list = []

        for p_i in range(self.N-1):

            current_E = []

            #for time_i in range(self.stabalization_step, len(self.time)):
            for time_i in range(int(1500000 / 500), len(self.time)):

                velocity = self.particle_v[p_i][time_i]
                V_mag_squared = (velocity[0]**2) + (velocity[1]**2) + (velocity[2]**2)
                current_E.append(.5 * self.m * V_mag_squared)

            avg_E = stats.mean(current_E)

            #E_k_list.append(math.log(avg_E))
            E_k_list.append(avg_E)

        plt.plot(E_k_list)
        plt.xlabel('particle')
        plt.ylabel('kinetic energy')
        plt.title(plot_title)

        string = plot_title.split(' ')
        force_freq = string[8].split('.')
        force_freq = force_freq[0] + "pt" + force_freq[1]
        # if string[8] == "0.1":
        #     pass
        if string[4] == "0,":
            save_title = f"1Dline_static_bonds_{force_freq}_forcefreq"
        else:
            save_title = f"1Dline_{string[4][:-1]}_bonds_{force_freq}_forcefreq"

        plt.savefig(save_title)
        shutil.copy(save_title + ".png", store_loc)

        plt.show()
        plt.clf()


        #fourier transform
        N = self.N - 100
        t = np.arange(100, self.N)
        fft = np.fft.fft(E_k_list[100:])
        T = t[1] - t[0] # sampling interval
        f = np.linspace(0, 1 / T, N)
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency [Hz]")
        #plt.plot(f[:N // 2], np.abs(fft)[:N // 2] * 1 / N)

        for value_i in range(len(np.abs(fft)[:N // 2] * 1 / N)):
            (np.abs(fft)[:N // 2] * 1 / N)[value_i] = math.log((np.abs(fft)[:N // 2] * 1 / N)[value_i])
        plt.plot(f[:N // 2], (np.abs(fft)[:N // 2] * 1 / N))

        #find peak of graph
        decreasing = True
        peak = (1, 100)
        for value_i in range(len(np.abs(fft)[:N // 2] * 1 / N)):
            amp = (np.abs(fft)[:N // 2] * 1 / N)[value_i]
            freq = f[:N // 2][value_i]

            if decreasing == True and peak[1] < amp:
                decreasing = False
            elif decreasing == True:
                peak = (freq, amp)
            elif decreasing == False and amp > peak[1]:
                peak = (freq, amp)

        peaks, _ = scipy.signal.find_peaks(np.abs(fft)[:N // 2] * 1 / N)
        plt.plot(f[:N // 2][peaks], (np.abs(fft)[:N // 2] * 1 / N)[peaks], "x")


        if string[4] == "0,":
            save_title = f"FourierPlot_static_bonds_{force_freq}_forcefreq"
        else:
            save_title = f"FourierPlot_{string[4][:-1]}_bonds_{force_freq}_forcefreq"
        plt.title(save_title + f"\nPeak at {peak[0]} Hz")
        plt.savefig(save_title)
        shutil.copy(save_title + ".png", store_loc)
        #plt.show()

        return peak

    def standing_wave_1D(self, p: int, A_p: float, B_p: float, save_loc: str, num_bonds: int):
        """Analysis function for 1D_line standing wave simulations. For the
        given p value the function creates position plots for each particle
        within the system (saving only the position plot for particle_i == p).
        From the position plot the function then preforms a Fourier transform to
        pull out the main ocilation frequency of the particles for the inputed
        p value.

        Returns - Tuple(avg_expected_freq, avg_simulation_freq)"""

        #calculate parameters depending on number of bond case
        w_p = 2 * math.sqrt(1 / self.m) * math.sin((p * math.pi) / (2 * self.N))
        if num_bonds == 1:
            w = .3044566161
            dk = .5
            k_0 = 1
            b = math.pi * p / 36
            gamma = dk / (2 * k_0)
            a = 4 * (b ** 2) * k_0 / (w ** 2)
            q = gamma * a

        #get current location for storing files
        current_loc = os.getcwd()

        #convert timesteps to realtime
        realtime = []
        for time in self.time:
            realtime.append(time * self.dt)

        expected_freq = []
        experimental_freq = []
        particle_list = []

        #run analysis on each particle for the selected p value
        for particle_i in range(1, self.N - 1):

            #for single dynamic bonds
            period = math.pi

            theoretical = []
            strob_time = []
            strob_pos = []
            z_time = []

            for time in realtime:
                #static theroy
                if num_bonds == 0:
                    theoretical.append((A_p * math.sin((p * math.pi * particle_i) / (self.N - 1)) * math.cos((w_p * time) + B_p)) - (float(self.N) / 2) + particle_i)
                elif num_bonds == 1:
                    z = w / 2 * time
                    z_time.append(z)
                    #theoretical.append(scipy.special.mathieu_cem(1, q, z)[0] - (float(self.N) / 2) + particle_i)
                    theoretical.append(0 - (float(self.N) / 2) + particle_i)

                #get strobed data
                if round(time % period, 1) == 0:
                    strob_time.append(time)
                    index = realtime.index(time)
                    strob_pos.append(self.particle_p[particle_i][index][0])

            #record x positions into lists
            p_i = []
            for position in self.particle_p[particle_i]:
                p_i.append(position[0])
                #p_i.append(math.log10(abs(position[0] - (float(self.N) / 2) + particle_i)))

            #create position plot
            if p == particle_i:
                g1 = plt.figure(1)
                #plt.xlabel('time')
                if num_bonds == 0:
                    plt.plot(realtime, theoretical)
                    plt.plot(realtime, p_i, "b")
                    plt.xlabel('time')
                elif num_bonds == 1:
                    plt.plot(z_time, p_i, "b")
                    plt.xlabel('z')
                plt.ylabel('position')
                #plt.semilogy(z_time, p_i)
                # plt.scatter(strob_time, strob_pos, color = 'r')
                plt.title(f"p: {p}, particle: {particle_i}")
                save_title = f"Position_plot_p_{int(p)}_Particle_{particle_i}"
                plt.savefig(save_title)
                shutil.move(current_loc + "/" + save_title + ".png", save_loc)
                # if p == 32:
                #     plt.show()
                plt.clf()

            #--------------------------------------

            #fourier transform
            N = len(p_i)
            fft = np.fft.fft(p_i)
            T = realtime[1] - realtime[0] # sampling interval
            f = np.linspace(0, 1 / T, N)
            plt.ylabel("Amplitude")
            plt.xlabel("Frequency [Hz]")

            for value_i in range(len(np.abs(fft)[:N // 2] * 1 / N)):
                (np.abs(fft)[:N // 2] * 1 / N)[value_i] = math.log((np.abs(fft)[:N // 2] * 1 / N)[value_i])
            plt.plot(f[:N // 2][1:], (np.abs(fft)[:N // 2] * 1 / N)[1:])
            #plt.xlim(.01, 2)

            #find peak of graph
            decreasing = True
            peak = (1, 100)
            for value_i in range(len(np.abs(fft)[:N // 2] * 1 / N)):
                amp = (np.abs(fft)[:N // 2] * 1 / N)[value_i]
                freq = f[:N // 2][value_i]

                if decreasing == True and peak[1] < amp:
                    decreasing = False
                elif decreasing == True:
                    peak = (freq, amp)
                elif decreasing == False and amp > peak[1]:
                    peak = (freq, amp)

            plt.title(f"\nPeak at {peak[0]} Hz")
            #plt.show()

            if p == particle_i:
                save_title = f"Fourier_plot_p_{int(p)}_Particle_{particle_i}"
                plt.savefig(save_title)
                shutil.move(current_loc + "/" + save_title + ".png", save_loc)
            plt.clf()

            expected_freq.append(w_p)
            experimental_freq.append(peak[0] * 2 * math.pi)
            particle_list.append(particle_i)

        return (stats.mean(expected_freq), stats.mean(experimental_freq))

    def standing_wave_1D_1spring(self, p, A_p, B_p, save_loc):

        #need to update the w value depending on the run paramaters
        w = .3044566161

        b = math.pi * p / 36
        gamma = dk / (2 * k_0)
        a = 4 * (b ** 2) * k_0 / (w ** 2)
        q = gamma * a

        current_loc = os.getcwd()

        #convert timesteps to realtime
        realtime = []
        for time in self.time:
            realtime.append(time * self.dt)

        expected_freq = []
        experimental_freq = []
        particle_list = []

        #for particle_i in specific_list:
        for particle_i in range(1, self.N - 1):

            w_p = 2 * math.sqrt(1 / self.m) * math.sin((p * math.pi) / (2 * self.N))

            #for single dynamic bonds
            period = math.pi

            theoretical = []
            strob_time = []
            strob_pos = []
            for time in realtime:
                #static theroy
                theoretical.append((A_p * math.sin((p * math.pi * particle_i) / (self.N - 1)) * math.cos((w_p * time) + B_p)) - (float(self.N) / 2) + particle_i)

                #1 spring theory
                z = w / (2 * time)

                #even (a solutions)
                theoretical.append(scipy.special.mathieu_cem())

                #get strobed data
                if round(time % period, 1) == 0:
                    strob_time.append(time)
                    index = realtime.index(time)
                    strob_pos.append(self.particle_p[particle_i][index][0])

            #record x positions into lists
            p_i = []
            for position in self.particle_p[particle_i]:
                p_i.append(position[0])

            g1 = plt.figure(1)
            plt.xlabel('time step')
            plt.ylabel('position')
            plt.plot(realtime, p_i, "b")
            plt.plot(realtime, theoretical)
            #print(len(strob_time), len(strob_pos))
            plt.scatter(strob_time, strob_pos, color = 'r')
            plt.title(f"p: {p}, particle: {particle_i}")
            #plt.show()
            #if int(p) in [2, 34]:
            if p == particle_i:
                save_title = f"Position_plot_p_{int(p)}_Particle_{particle_i}"
                plt.savefig(save_title)
                shutil.move(current_loc + "/" + save_title + ".png", save_loc)
            plt.clf()

            #fourier transform
            N = len(p_i)
            fft = np.fft.fft(p_i)
            T = realtime[1] - realtime[0] # sampling interval
            f = np.linspace(0, 1 / T, N)
            plt.ylabel("Amplitude")
            plt.xlabel("Frequency [Hz]")

            for value_i in range(len(np.abs(fft)[:N // 2] * 1 / N)):
                (np.abs(fft)[:N // 2] * 1 / N)[value_i] = math.log((np.abs(fft)[:N // 2] * 1 / N)[value_i])
            plt.plot(f[:N // 2][1:], (np.abs(fft)[:N // 2] * 1 / N)[1:])
            #plt.xlim(.01, 2)

            #find peak of graph
            decreasing = True
            peak = (1, 100)
            for value_i in range(len(np.abs(fft)[:N // 2] * 1 / N)):
                amp = (np.abs(fft)[:N // 2] * 1 / N)[value_i]
                freq = f[:N // 2][value_i]

                if decreasing == True and peak[1] < amp:
                    decreasing = False
                elif decreasing == True:
                    peak = (freq, amp)
                elif decreasing == False and amp > peak[1]:
                    peak = (freq, amp)

            plt.title(f"\nPeak at {peak[0]} Hz")
            #plt.show()
            #if p in [2.0, 34.0]:
            if p == particle_i:
                save_title = f"Fourier_plot_p_{int(p)}_Particle_{particle_i}"
                plt.savefig(save_title)
                shutil.move(current_loc + "/" + save_title + ".png", save_loc)
            plt.clf()

            expected_freq.append(w_p)
            experimental_freq.append(peak[0] * 2 * math.pi)
            particle_list.append(particle_i)

        return (stats.mean(expected_freq), stats.mean(experimental_freq))

    def graph_k_energy(self, run_num = 0, plot_title = "Kinetic_Engery_plot", store_loc = None):
        """Creates a graph of the kinetic energy of the system over the timesteps"""

        g1 = plt.figure(1)
        plt.xlabel('time')
        plt.ylabel('kinetic energy')
        plt.plot(self.realTime, self.kinetic_energy, "b")
        plt.title(plot_title, fontsize = 8)

        if store_loc != None:
            plt.savefig("Kinetic_Energy_plot")
            shutil.move("Kinetic_Energy_plot" + ".png", store_loc)
        else:
            plt.show()

        plt.clf()

    def graph_avg_k_energy(self):
        """Creates a graph that plots the average kinetic energy for each recorded simulation run
        over the angular frequency of the dynamic force in the simulation"""

        colors = []
        for w in self.dynamic_frequencies:
            if (self.spring_frequency % w == 0) or (w % self.spring_frequency == 0):
                colors.append('r')
            else:
                colors.append('b')

        g1 = plt.figure(1)
        plt.xlabel('angular frequency')
        plt.ylabel('average kinetic energy')
        #print(self.average_kinetic)
        plt.scatter(self.dynamic_frequencies, self.average_kinetic, c=colors)
        plt.show()

    def graph_pos_v_force(self, p_index):
        """creates graphs of the dynamic force magnitude on the system and the position and velocity of
        a select particle in the system"""

        g1 = plt.figure(1)
        plt.xlabel('time')
        plt.ylabel('force')
        plt.plot(self.dy_time, self.force, "b")

        p_i_a = []
        for velocity in self.particle_v[p_index]:
            p_i_a.append(velocity[0])
        g2 = plt.figure(2)
        plt.xlabel('time step')
        plt.ylabel('velocity')
        plt.plot(self.time, p_i_a, "b")

        p_i_b = []
        for position in self.particle_p[p_index]:
            p_i_b.append(position[0])
        g3 = plt.figure(3)
        plt.xlabel('time step')
        plt.ylabel('position')
        plt.plot(self.time, p_i_b, "b")

        plt.show()

    def graph_dynamic_force(self):
        """Creates and displays a graph of the magnitude of the dynamic force on the system over time"""

        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('force magnitude')
        plt.plot(self.dy_time, self.force, "b")
        plt.show()

    def graph_p_v(self, p_index: int):
        """Creates and displays a graph a particle (p_index) velocity over the timesteps of the simulation"""

        p_i = []
        for velocity in self.particle_v[p_index]:
            p_i.append(velocity[0])
        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('velocity')
        plt.plot(self.time, p_i, "b")
        plt.show()

    def graph_p_pos(self, p_index: int, store_loc = None):
        """Creates and displays a graph a particle p_index position over the timesteps of the simulation"""
        equilibriam_pos = (-self.N / 2) + p_index
        p_i = []
        for position in self.particle_p[p_index]:
            p_i.append(position[0] - equilibriam_pos)

        plot_title = f"Pos_plot_particle_{p_index}"
        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('position')
        plt.plot(self.time, p_i, "b")
        plt.title(plot_title)

        if store_loc != None:
            plt.savefig(plot_title)
            shutil.move(plot_title + ".png", store_loc)
        else:
            plt.show()

        plt.clf()

    def graph_PE(self):
        """To graph the potential energy the energy must be logged first"""

        data = numpy.genfromtxt(fname='PE_and_Temp.log', skip_header=True)
        plt.figure(figsize=(4, 2.2), dpi=140)
        plt.plot(data[10:, 0], data[10:, 1])
        plt.xlabel('time step')
        plt.ylabel('potential_energy')
        plt.show()

    def graph_Temp(self):
        """To graph the temperature must be logged first, done through logging potential energy"""

        data = numpy.genfromtxt(fname='PE_and_Temp.log', skip_header=True)
        plt.figure(figsize=(4, 2.2), dpi=140)
        plt.plot(data[10:, 0], data[10:, 2])
        plt.xlabel('time step')
        plt.ylabel('Thermal temperature')
        plt.show()

    def graph_mult_pos(self, center: int, spacing=1,dynamic=False):
        '''creates a stacked plot of the displacement vs time step of 4 particles to either side of the center particle
           separated by the specified spacing'''
        ax_list = []
        figure = plt.figure()
        height = 0.1

        #create axes for all 9 subplots
        for i in range (9):
            if i == 0:
                ax_list.append(figure.add_axes([0.1,0.05+(i*height),0.85,height],ylim=(-1.2,1.2),yticklabels=[]))
            else:
                ax_list.append(figure.add_axes([0.1,0.05+(i*height),0.85,height],xticklabels=[],ylim=(-1.2,1.2),yticklabels=[]))


            #populate list with particle displacement data
            p_i = []
            for position in self.particle_p[center+spacing*(i-4)]:
                p_i.append(position[0]-spacing*(i-4))

            #plot that particle's displacement vs time step
            ax_list[i].plot(self.time, p_i, ".")
            ax_list[i].set_ylabel(center+spacing*(i-4))

        if dynamic == True:
            plt.title('displacement from equilibrium vs time step (dynamic springs)')
        else:
            plt.title('displacement from equilibrium vs time step (static springs)')
        plt.show()

    def graph_mult_pos_3D(self, p_list: list, dynamic=False):
        '''creates a stacked plot of the displacement vs time step of the particles in p_list
           each plot displays lines of the particle's x, y, and z displacement from their
           initial positions.

           x displacement is the blue line
           y displacement is the red line
           z displacement is the green line '''

        initial_pos_x = {}
        initial_pos_y = {}
        initial_pos_z = {}
        ax_list = []
        figure = plt.figure()
        height = 0.12

        for p_i in p_list:
            position_x = self.particle_p[p_i][0][0]
            position_y = self.particle_p[p_i][0][1]
            position_z = self.particle_p[p_i][0][2]
            initial_pos_x[p_i] = position_x
            initial_pos_y[p_i] = position_y
            initial_pos_z[p_i] = position_z



        #create axes for all subplots
        counter = 0
        for p_i in p_list:
            if p_i == 0:
                ax_list.append(figure.add_axes([0.1,0.05+(counter*height),0.85,height],ylim=(-.2,.2),yticklabels=[]))
            else:
                ax_list.append(figure.add_axes([0.1,0.05+(counter*height),0.85,height],xticklabels=[],ylim=(-.2,.2),yticklabels=[]))

            #populate list with particle displacement data
            p_i_x = []
            p_i_y = []
            p_i_z = []
            for position in self.particle_p[p_i]:
                p_i_x.append(position[0]-initial_pos_x[p_i])
                p_i_y.append(position[1]-initial_pos_y[p_i])
                p_i_z.append(position[2]-initial_pos_z[p_i])


            #plot that particle's displacement vs time step
            #ax_list[counter].plot(self.time, p_i_x, "r", self.time, p_i_y, "b")
            ax_list[counter].plot(self.time, p_i_x, "b", self.time, p_i_y, "r", self.time, p_i_z, "g")
            ax_list[counter].set_ylabel(p_i)

            counter += 1

        if dynamic == True:
            plt.title('displacement from equilibrium vs time step (dynamic springs)')
        else:
            plt.title('displacement from equilibrium vs time step (static springs)')
        plt.show()

    def wave_packet(self, store_loc = None, plot_title = "Gausian plot", num_samples = 8):

        # times = [0]
        sample_period = (self.time[-1] + 500) / num_samples

        times = []
        for time_step_i in range(len(self.time)):
            time_step = self.time[time_step_i]
            if time_step % sample_period == 0:
                times.append(time_step_i)

        g1 = plt.figure(1)
        plt.xlabel('particle position')
        plt.ylabel('time')
        plt.yticks([])
        plt.title(plot_title, fontsize = 7)

        shift = 0

        #create and save waterfall plot of wave packet
        for i in times:
            time_step = self.time[i]
            p_list = []
            offsets = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_p[p_index][i][0] - equilibriam_pos + shift

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_p[p_index][i][0] > 0:
                        position = -self.particle_p[p_index][i][0] - equilibriam_pos + shift

                p_list.append(p_index)
                offsets.append(position)

            plt.plot(p_list, offsets)

            #for 12 particle
            shift -= .5
            # #for 7 particle
            # shift -= 2

        if store_loc != None:
            plt.savefig("Gausian_plot")
            shutil.move("Gausian_plot.png", store_loc)
        else:
            plt.show()
        plt.clf()

        #get fourier plot and data for either sample
        fourier_data_list = []
        for i in times:
            time_step = self.time[i]
            p_list = []
            offsets = []
            gausian = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_p[p_index][i][0] - equilibriam_pos

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_p[p_index][i][0] > 0:
                        position = -self.particle_p[p_index][i][0] - equilibriam_pos

                p_list.append(p_index)
                offsets.append(position)
                gausian.append(abs(position))

            #plt.plot(p_list, gausian)
            fourier_data = self.fourier_plot(p_list, offsets, store_loc = store_loc, plot_title = f"Fourier_plot_time={int(time_step * self.dt)}")
            fourier_data_list.append(fourier_data)

        #-----------------------------------------------------------------------
        mean_squared_error_list = []
        initial_data = fourier_data_list[0]
        for fourier_data in fourier_data_list:
            diff_sqrd_sum = 0
            for i in range(len(fourier_data[0])):
                diff_sqrd_sum += ((initial_data[1][i] - fourier_data[1][i]) ** 2)
            mean_squared_error_list.append(diff_sqrd_sum / len(fourier_data[0]))

        error_plot_title = "Mean Squared Error - (over course of Simulation)"
        plt.plot(mean_squared_error_list)
        plt.title(error_plot_title)
        plt.savefig(error_plot_title)
        if store_loc != None:
            shutil.move(f"{error_plot_title}.png", store_loc)
        plt.clf()


    def fourier_plot(self, x_data: list, y_data: list, store_loc = None, plot_title = "Fourier_plot"):
        """
        Function for ploting and returning the fourier data of a given inputed data set
        """

        x_data = np.array(x_data)
        y_data = np.array(y_data)

        fft = np.fft.fft(y_data)
        fft[0] = 0

        N = len(x_data)
        T = x_data[1] - x_data[0] # sampling interval
        freq = np.fft.fftfreq(N, d=T)

        plt.plot(abs(freq), abs(fft.real))
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency [1 / wave-length]")
        #plt.xlim(left=.025)
        #plt.ylim(-.5,.5)

        #find peak of graph
        decreasing = True
        peak = None
        for value_i in range(len(x_data)):
            amp = fft[value_i]
            frequency = freq[value_i]

            if peak == None:
                peak = (frequency, amp)
            elif decreasing == True and peak[1] < amp:
                decreasing = False
            elif decreasing == True:
                peak = (frequency, amp)
            elif decreasing == False and amp > peak[1]:
                peak = (frequency, amp)

        plt.title(f"{plot_title}\nPeak at {peak[0]} Hz")

        plt.savefig(plot_title)
        if store_loc != None:
            shutil.move(f"{plot_title}.png", store_loc)
        plt.clf()

        return (freq, fft.real)

    def graph_stiffness(self, spr_index: int):
        """Creates and displays a graph of spring spr_index stiffness over the timesteps of the simulation"""

        k_0 = []
        k_1 = []
        k_2 = []
        for stiffness in self.stiffness[0]:
            k_0.append(stiffness)
        for stiffness in self.stiffness[1]:
            k_1.append(stiffness)
        for stiffness in self.stiffness[2]:
            k_2.append(stiffness)
        #print('k_0 has ' , len(k_0), ' elements')
        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('stiffness')
        plt.plot(self.time, k_0, ".",self.time,k_1,'.',self.time,k_2,'.')
        plt.show()

    def contraction_expansion_tri(self):
        left_tags = []
        right_tags = []
        top_tags = []
        bottom_tags = []

        widths = []
        hights = []

        #left
        for i in range(1, (2*self.y), 2):
            left_tags.append(i)

        #right
        for i in range((2*self.x*self.y-(2*self.y)),(2*self.x*self.y), 2):
            right_tags.append(i)

        #top
        for i in range((2*self.y-2), (2*self.x*self.y), (2*self.y)):
            top_tags.append(i)

        #bottom
        for i in range(1, (2*self.x*self.y-(2*self.y)+3), (2*self.y)):
            bottom_tags.append(i)

        for i in range(len(self.time)):

            total = 0.0
            for p_i in left_tags:
                pos = self.particle_p[p_i][i][0]
                total += pos
            avg_left = total / len(left_tags)

            total = 0.0
            for p_i in right_tags:
                pos = self.particle_p[p_i][i][0]
                total += pos
            avg_right = total / len(right_tags)

            total = 0.0
            for p_i in top_tags:
                pos = self.particle_p[p_i][i][1]
                total += pos
            avg_top = total / len(top_tags)

            total = 0.0
            for p_i in bottom_tags:
                pos = self.particle_p[p_i][i][1]
                total += pos
            avg_bottom = total / len(bottom_tags)

            width = avg_right - avg_left
            hight = avg_top - avg_bottom

            widths.append(width)
            hights.append(hight)

        ratios = []
        for i in range(len(self.time)):
            delta_x = (widths[i] - widths[0])
            delta_y = (hights[i] - hights[0])
            if delta_x == 0:
                ratios.append(0)
            else:
                ratio = delta_y / delta_x
                ratios.append(ratio)


        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('width to hight ratio of lattice')
        #plt.plot(self.time, widths, "b", self.time, hights, "r")
        plt.plot(self.time, ratios)
        plt.show()

    def width_of_tri(self, force_freq: float, run_type: str, gamma: float, w: float):
        left_tags = []
        right_tags = []

        #left
        for i in range(1, (2*self.y), 2):
            left_tags.append(i)

        #right
        for i in range((2*self.x*self.y-(2*self.y)),(2*self.x*self.y), 2):
            right_tags.append(i)

        widths = []
        for i in range(len(self.time)):

            total = 0.0
            for p_i in left_tags:
                pos = self.particle_p[p_i][i][0]
                total += pos
            avg_left = total / len(left_tags)

            total = 0.0
            for p_i in right_tags:
                pos = self.particle_p[p_i][i][0]
                total += pos
            avg_right = total / len(right_tags)

            width = avg_right - avg_left

            widths.append(width)

        width_ratios = []
        width_initial = widths[0]
        average_width = 0
        #print(width_initial)
        for width in widths:
            ratio = abs((width - width_initial)) / width_initial * 100
            #print(f"width_initial is {width_initial}, width is {width}, ratio is {ratio}")
            width_ratios.append(ratio)
            #average_width += width

        real_time = []
        for time in self.time:
            real_time.append(time*self.dt)

        #print(len(widths))
        #print("calc_avg_width", self.last_step)
        average_width = avg(widths, real_time, start=self.last_step, finish=-1)
        #print(average_width)
        average_width = round(average_width, 8)
        #average_width = average_width / len(widths)


        average_k = round(self.average_kinetic[-1], 8)
        #print(width_ratios)
        g1 = plt.figure(1)
        plot_title = f"Fractional width change graph ({run_type} springs),\n spring frequency = {w}, spring mag = 1 ± .5, static force magnitude = .1, \n gamma ={gamma}, dt = .0001, average width = {average_width}, avg kinetic = {average_k}"

        plt.xlabel('time')
        plt.ylabel('fractional change of width of lattice')
        plt.title(plot_title, fontsize=10)
        #plt.plot(self.time, width_ratios)
        plt.plot(real_time, width_ratios)
        plt.savefig(f"fractional_width_change_graph_of_force_freq_{force_freq}_with_{run_type}_springs")
        # if save_location != " ":
        #     shutil.copy(f"binplot_force_freq{force_frequency}.pdf", save_location)
        #plt.show()
        plt.clf()



# SED Functions
    def SED(self, dof=1, Nk=15, Nw=100, cb=100, spring_const=4000):
        ## Nathan Villiger
        ## July 2019
        ## This code will take in velocity data from molecular dynamics simulations
        ## and compute the spectral energy density as a function of frequency and
        ## wavenumber.
        ##
        ## We will implement the Spectral Energy Density Method as outlined in
        ## Vila et al (2017) to generate dispersion diagrams from MD simulation data.
        ##
        ##
        ## This code will rely on the velocity data for each particle in the mass-spring
        ## chain throughout the simulation as well as the total time simulated and the
        ## time difference between successive data points, the mass of each particle,
        ## the total number of particles, the equilibrium spacing between adjacent
        ## masses, and the number of degrees of freedom within the unit cell.


        ## number of discrete points in frequency and wavenumber space that will be
        ## considered

        ## highest frequency to be considered
        w_max = np.sqrt(4*spring_const/self.m) * 1.1

        ## define vectors of points in frequency and wavenumber where SED will be
        ## calculated
        self.kspace_vec = np.linspace(-np.pi,np.pi,Nk)
        kspacelength = self.kspace_vec.size
        self.wspace_vec = np.linspace(0,w_max,Nw)
        wspacelength = self.wspace_vec.size

        ## number of degrees of freedom (e.g. different types of springs)
        dof = 1
        ''' would prefer to not hard program this here '''

        ## equilibrium spacing between adjacent masses
        a = 1
        ''' would prefer to get this from the simulation, not hard program '''

        ## number of unit cells
        #N = 25
        N = self.N / dof

        ## determine number of time steps and real time simulated
        tspacelength = self.time[-1]
        Time = tspacelength * self.dt

        ## define 2D array that will be filled with values of energy density
        SED = np.zeros((kspacelength,wspacelength), dtype=complex)

        ## define 4D array to be filled with fourier series representation of velocities
        fourier_rep = np.zeros((self.N,kspacelength,wspacelength,int(self.time[-1]/cb)), dtype=complex)

        ## also preallocate the preliminary matrices
        temp1 = np.zeros((self.N,kspacelength,wspacelength), dtype=complex)
        temp2 = np.zeros((dof,kspacelength,wspacelength), dtype=complex)

        print('cb is ' , cb)

        t0 = time.time()

        ## now go through and compute the energy density at each point in the SED matrix
        for q in range (kspacelength):
            #print(f"q{q}")
            for w in range (wspacelength):
                #print(f"w{w}")
                for num in range (self.N):
                    for t_index in range ((len(self.time)) - 1):
                        #v_mag = math.sqrt((self.particle_v[num][t_index][0])**2)
                        #+ (self.particle_v[num][t_index][1]**2)
                        #+ (self.particle_v[num][t_index][2]**2)

                        fourier_rep[num,q,w,int(self.time[t_index]/cb)] = self.particle_v[num][(t_index)][0] * np.exp(1j*(self.kspace_vec[q]*num*a - self.wspace_vec[w]*self.time[(t_index)]*self.dt))

        t1 = time.time()
        print(fourier_rep.nbytes/1e9 , 'gigabytes used by fourier_rep')
        print('fourier_rep has dimensions ' , fourier_rep.shape)
        print(t1-t0, ' seconds to compute fourier_rep')


        # perform the trapezoidal integral over time
        #temp1 = self.dt * np.sum(fourier_rep,axis=3)
        temp1 = np.trapz(fourier_rep,dx=cb,axis=3)
        t2 = time.time()
        print(t2-t1, " seconds to compute temp1")
        print(temp1.nbytes/1e9 , 'gigabytes used by temp1')

        for q in range (kspacelength):
            for w in range (wspacelength):

                # perform the sum over all unit cells, keeping track of which entries
                # correspond to which position in the unit cell
                for UCpos in range (dof):
                    temp2[UCpos,q,w] = np.sum(temp1[UCpos::dof,q,w],axis=0)

                # calculate spectral energy density
                SED[q,w] = (self.m/(2*(Time*N)**2)) * np.sum(np.abs(temp2[:,q,w])**2,axis=0)
        t3 = time.time()
        print(t3-t2, ' seconds to compute the SED')
        return SED

    def SED_2mass(self, dof=2, p_mass1=1, p_mass2=2, Nk=15, Nw=100, cb=100, k1=4000, k2=8000):
        ## Nathan Villiger
        ## July 2019
        ## This code will take in velocity data from molecular dynamics simulations
        ## and compute the spectral energy density as a function of frequency and
        ## wavenumber.
        ##
        ## We will implement the Spectral Energy Density Method as outlined in
        ## Vila et al (2017) to generate dispersion diagrams from MD simulation data.
        ##
        ##
        ## This code will rely on the velocity data for each particle in the mass-spring
        ## chain throughout the simulation as well as the total time simulated and the
        ## time difference between successive data points, the mass of each particle,
        ## the total number of particles, the equilibrium spacing between adjacent
        ## masses, and the number of degrees of freedom within the unit cell.


        mass_list = [p_mass1, p_mass2]

        ## highest frequency to be considered if different masses
        if p_mass1 != p_mass2:
            w_max = np.sqrt((k1/(p_mass1*p_mass2))*(p_mass1+p_mass2+np.sqrt(p_mass1**2+p_mass2**2+2*p_mass1*p_mass2))) * 1.1
        elif k1 != k2: ## highest frequency to be considered if different springs
            w_max = np.sqrt((k1+k2)/p_mass1+(1/p_mass1)*np.sqrt(k1**2+k2**2+2*k1*k2)) * 1.1

        ## number of degrees of freedom (e.g. different types of springs)
        dof = 2
        ''' would prefer to not hard program this here '''

        ## equilibrium spacing between adjacent masses
        a = 1
        ''' would prefer to get this from the simulation, not hard program '''

        ## define vectors of points in frequency and wavenumber where SED will be
        ## calculated
        self.kspace_vec = np.linspace(-np.pi/(a*dof),np.pi/(a*dof),Nk)
        kspacelength = self.kspace_vec.size
        self.wspace_vec = np.linspace(0,w_max,Nw)
        wspacelength = self.wspace_vec.size

        ## number of unit cells
        #N = 25
        N = self.N / dof

        ## determine number of time steps and real time simulated
        tspacelength = self.time[-1]
        Time = tspacelength * self.dt

        ## define 2D array that will be filled with values of energy density
        SED = np.zeros((kspacelength,wspacelength), dtype=complex)

        ## define 4D array to be filled with fourier series representation of velocities
        fourier_rep = np.zeros((self.N,kspacelength,wspacelength,int(self.time[-1]/cb)), dtype=complex)

        ## also preallocate the preliminary matrices
        temp1 = np.zeros((self.N,kspacelength,wspacelength), dtype=complex)
        temp2 = np.zeros((dof,kspacelength,wspacelength), dtype=complex)

        print('cb is ' , cb)

        t0 = time.time()

        ## now go through and compute the energy density at each point in the SED matrix
        for q in range (kspacelength):
            #print(f"q{q}")
            for w in range (wspacelength):
                #print(f"w{w}")
                for num in range (self.N):
                    for t_index in range ((len(self.time)) - 1):
                        #v_mag = math.sqrt((self.particle_v[num][t_index][0])**2)
                        #+ (self.particle_v[num][t_index][1]**2)
                        #+ (self.particle_v[num][t_index][2]**2)

                        fourier_rep[num,q,w,int(self.time[t_index]/cb)] = self.particle_v[num][(t_index)][0] * np.exp(1j*(self.kspace_vec[q]*num*a - self.wspace_vec[w]*self.time[(t_index)]*self.dt))

        t1 = time.time()
        print(fourier_rep.nbytes/1e9 , 'gigabytes used by fourier_rep')
        print('fourier_rep has dimensions ' , fourier_rep.shape)
        print(t1-t0, ' seconds to compute fourier_rep')


        # perform the trapezoidal integral over time
        #temp1 = self.dt * np.sum(fourier_rep,axis=3)
        temp1 = np.trapz(fourier_rep,dx=cb,axis=3)
        t2 = time.time()
        print(t2-t1, " seconds to compute temp1")
        print(temp1.nbytes/1e9 , 'gigabytes used by temp1')

        for q in range (kspacelength):
            for w in range (wspacelength):

                # perform the sum over all unit cells, keeping track of which entries
                # correspond to which position in the unit cell
                for UCpos in range (dof):
                    temp2[UCpos,q,w] = mass_list[UCpos] * np.sum(temp1[UCpos::dof,q,w],axis=0)

                # calculate spectral energy density
                SED[q,w] = (1/(2*(Time*N)**2)) * np.sum(np.abs(temp2[:,q,w])**2,axis=0)
        t3 = time.time()
        print(t3-t2, ' seconds to compute the SED')
        return SED

    def SED_3mass(self, dof=3, p_mass1=1, p_mass2=2, p_mass3=3, Nk=15, Nw=100, cb=100, k1=4000, k2=4000, k3=4000):
        ## Nathan Villiger
        ## July 2019
        ## This code will take in velocity data from molecular dynamics simulations
        ## and compute the spectral energy density as a function of frequency and
        ## wavenumber.
        ##
        ## We will implement the Spectral Energy Density Method as outlined in
        ## Vila et al (2017) to generate dispersion diagrams from MD simulation data.
        ##
        ##
        ## This code will rely on the velocity data for each particle in the mass-spring
        ## chain throughout the simulation as well as the total time simulated and the
        ## time difference between successive data points, the mass of each particle,
        ## the total number of particles, the equilibrium spacing between adjacent
        ## masses, and the number of degrees of freedom within the unit cell.


        mass_list = [p_mass1, p_mass2, p_mass3]

        ## highest frequency to be considered
        w_max = np.sqrt((k1+k2)/p_mass1+(1/p_mass1)*np.sqrt(k1**2+k2**2+2*k1*k2)) * 1.1

        ## number of degrees of freedom (e.g. different types of springs)
        dof = 3
        ''' would prefer to not hard program this here '''

        ## equilibrium spacing between adjacent masses
        a = 1
        ''' would prefer to get this from the simulation, not hard program '''

        ## define vectors of points in frequency and wavenumber where SED will be
        ## calculated
        self.kspace_vec = np.linspace(-np.pi/(a*dof),np.pi/(a*dof),Nk)
        kspacelength = self.kspace_vec.size
        self.wspace_vec = np.linspace(0,w_max,Nw)
        wspacelength = self.wspace_vec.size

        ## number of unit cells
        #N = 25
        N = self.N / dof

        ## determine number of time steps and real time simulated
        tspacelength = self.time[-1]
        Time = tspacelength * self.dt

        ## define 2D array that will be filled with values of energy density
        SED = np.zeros((kspacelength,wspacelength), dtype=complex)

        ## define 4D array to be filled with fourier series representation of velocities
        fourier_rep = np.zeros((self.N,kspacelength,wspacelength,int(self.time[-1]/cb)), dtype=complex)

        ## also preallocate the preliminary matrices
        temp1 = np.zeros((self.N,kspacelength,wspacelength), dtype=complex)
        temp2 = np.zeros((dof,kspacelength,wspacelength), dtype=complex)

        print('cb is ' , cb)

        t0 = time.time()

        ## now go through and compute the energy density at each point in the SED matrix
        for q in range (kspacelength):
            #print(f"q{q}")
            for w in range (wspacelength):
                #print(f"w{w}")
                for num in range (self.N):
                    for t_index in range (2000,(len(self.time)) - 1):
                        #v_mag = math.sqrt((self.particle_v[num][t_index][0])**2)
                        #+ (self.particle_v[num][t_index][1]**2)
                        #+ (self.particle_v[num][t_index][2]**2)

                        fourier_rep[num,q,w,int(self.time[t_index]/cb)] = self.particle_v[num][(t_index)][0] * np.exp(1j*(self.kspace_vec[q]*num*a - self.wspace_vec[w]*self.time[(t_index)]*self.dt))

        t1 = time.time()
        print(fourier_rep.nbytes/1e9 , 'gigabytes used by fourier_rep')
        print('fourier_rep has dimensions ' , fourier_rep.shape)
        print(t1-t0, ' seconds to compute fourier_rep')


        # perform the trapezoidal integral over time
        #temp1 = self.dt * np.sum(fourier_rep,axis=3)
        temp1 = np.trapz(fourier_rep,dx=cb,axis=3)
        t2 = time.time()
        print(t2-t1, " seconds to compute temp1")
        print(temp1.nbytes/1e9 , 'gigabytes used by temp1')

        for q in range (kspacelength):
            for w in range (wspacelength):

                # perform the sum over all unit cells, keeping track of which entries
                # correspond to which position in the unit cell
                for UCpos in range (dof):
                    temp2[UCpos,q,w] = mass_list[UCpos] * np.sum(temp1[UCpos::dof,q,w],axis=0)

                # calculate spectral energy density
                SED[q,w] = (1/(2*(Time*N)**2)) * np.sum(np.abs(temp2[:,q,w])**2,axis=0)
        t3 = time.time()
        print(t3-t2, ' seconds to compute the SED')
        return SED

    def Plot_SED1(self, SED, spring_const=4000):

        SED_mean = np.mean(SED,axis=2)

        plt.contourf(self.kspace_vec,self.wspace_vec,np.transpose(SED_mean))
        plt.plot(self.kspace_vec,np.abs(np.sqrt(4*spring_const/self.m) * np.sin(self.kspace_vec/2)),'--')
        plt.colorbar()
        plt.title('SED')
        #plt.show()

        X,Y = np.meshgrid(self.kspace_vec,self.wspace_vec)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,np.transpose(SED_mean))


        plt.show()

    def Plot_SED2(self, SED, p_mass1=1, p_mass2=2, k1=4000, k2=4000):

        SED_mean = np.mean(SED,axis=2)

        a=1

        plt.contourf(self.kspace_vec,self.wspace_vec,np.transpose(np.log10(SED_mean)))
        if p_mass1 != p_mass2: # use these if different masses
            plt.plot(self.kspace_vec,np.sqrt((spring_const/(p_mass1*p_mass2))*(p_mass1+p_mass2-np.sqrt(p_mass1**2+p_mass2**2+2*p_mass1*p_mass2*np.cos(self.kspace_vec*(2*a))))),'--')
            plt.plot(self.kspace_vec,np.sqrt((spring_const/(p_mass1*p_mass2))*(p_mass1+p_mass2+np.sqrt(p_mass1**2+p_mass2**2+2*p_mass1*p_mass2*np.cos(self.kspace_vec*(2*a))))),'--')
        elif k1 != k2: # use these if different springs
            plt.plot(self.kspace_vec,np.sqrt((k1+k2)/p_mass1-(1/p_mass1)*np.sqrt(k1**2+k2**2+2*k1*k2*np.cos(self.kspace_vec*2*a))),'--')
            plt.plot(self.kspace_vec,np.sqrt((k1+k2)/p_mass1+(1/p_mass1)*np.sqrt(k1**2+k2**2+2*k1*k2*np.cos(self.kspace_vec*2*a))),'--')

        plt.colorbar()
        plt.title('log10 SED')
        #plt.show()

        X,Y = np.meshgrid(self.kspace_vec,self.wspace_vec)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,np.transpose(SED_mean))


        plt.show()

    def Plot_SED3(self, SED, p_mass1=1, p_mass2=2, p_mass3=3, k1=4000, k2=4000, k3=4000):

        SED_mean = np.mean(SED,axis=2)

        plt.contourf(self.kspace_vec,self.wspace_vec,np.transpose(np.log10(SED_mean)))
        #plt.plot(self.kspace_vec,np.sqrt((spring_const/(p_mass1*p_mass2))*(p_mass1+p_mass2-np.sqrt(p_mass1**2+p_mass2**2+2*p_mass1*p_mass2*np.cos(self.kspace_vec)))),'--')
        #plt.plot(self.kspace_vec,np.sqrt((spring_const/(p_mass1*p_mass2))*(p_mass1+p_mass2+np.sqrt(p_mass1**2+p_mass2**2+2*p_mass1*p_mass2*np.cos(self.kspace_vec)))),'--')
        plt.colorbar()
        plt.title('log10 SED')
        #plt.show()

        X,Y = np.meshgrid(self.kspace_vec,self.wspace_vec)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,np.transpose(SED_mean))


        plt.show()

    def SED_3D(self, dof=3, p_mass1=1, p_mass2=1, p_mass3=1, Nk=15, Nw=100, cb=100, k1=4000, k2=4000, k3=4000):
        ## Nathan Villiger
        ## July 2019
        ## This code will take in velocity data from molecular dynamics simulations
        ## and compute the spectral energy density as a function of frequency and
        ## wavenumber.
        ##
        ## We will implement the Spectral Energy Density Method as outlined in
        ## Vila et al (2017) to generate dispersion diagrams from MD simulation data.
        ##
        ##
        ## This code will rely on the velocity data for each particle in the mass-spring
        ## chain throughout the simulation as well as the total time simulated and the
        ## time difference between successive data points, the mass of each particle,
        ## the total number of particles, the equilibrium spacing between adjacent
        ## masses, and the number of degrees of freedom within the unit cell.


        ## number of discrete points in frequency and wavenumber space that will be
        ## considered

        ## highest frequency to be considered
        w_max = np.sqrt(4*k1/p_mass1) * 1.1

        ## equilibrium spacing between adjacent masses
        a = 1
        ''' would prefer to get this from the simulation, not hard program '''

        ## number of degrees of freedom (e.g. different types of springs)
        #dof = 1
        ''' would prefer to not hard program this here '''

        mass_list = [p_mass1, p_mass2, p_mass3]

        ## define vectors of points in frequency and wavenumber where SED will be
        ## calculated
        self.kspace_vec = np.linspace(-np.pi/(dof * a),np.pi/(dof * a),Nk)
        kspacelength = self.kspace_vec.size
        self.wspace_vec = np.linspace(0,w_max,Nw)
        wspacelength = self.wspace_vec.size


        ## number of unit cells
        #N = 25
        N = self.N / dof

        ## determine number of time steps and real time simulated
        tspacelength = self.time[-1]
        Time = tspacelength * self.dt

        ## define 2D array that will be filled with values of energy density
        SED = np.zeros((kspacelength,wspacelength), dtype=complex)

        ## define 4D array to be filled with fourier series representation of velocities
        fourier_rep = np.zeros((self.N,kspacelength,wspacelength,int(self.time[-1]/cb),3), dtype=complex)

        ## also preallocate the preliminary matrices
        temp1 = np.zeros((self.N,kspacelength,wspacelength,3), dtype=complex)
        temp2 = np.zeros((dof,kspacelength,wspacelength,3), dtype=complex)

        print('cb is ' , cb)

        t0 = time.time()

        ## now go through and compute the energy density at each point in the SED matrix
        for q in range (kspacelength):
            #print(f"q{q}")
            for w in range (wspacelength):
                #print(f"w{w}")
                for num in range (self.N):
                    for t_index in range ((len(self.time)) - 1):
                        for dim in range (3):

                            fourier_rep[num,q,w,int(self.time[t_index]/cb),dim] = self.particle_v[num][(t_index)][dim] * np.exp(1j*(self.kspace_vec[q]*num*a - self.wspace_vec[w]*self.time[(t_index)]*self.dt))


        t1 = time.time()
        print(fourier_rep.nbytes/1e9 , 'gigabytes used by fourier_rep')
        print('fourier_rep has dimensions ' , fourier_rep.shape)
        print(t1-t0, ' seconds to compute fourier_rep')


        # perform the trapezoidal integral over time
        #temp1 = self.dt * np.sum(fourier_rep,axis=3)
        temp1 = np.trapz(fourier_rep,dx=cb,axis=3)
        t2 = time.time()
        print(t2-t1, " seconds to compute temp1")
        print(temp1.nbytes/1e9 , 'gigabytes used by temp1')


        for q in range (kspacelength):
            for w in range (wspacelength):
                for dim in range (3):

                    # perform the sum over all unit cells, keeping track of which entries
                    # correspond to which position in the unit cell
                    for UCpos in range (dof):
                        temp2[UCpos,q,w,dim] = mass_list[UCpos] * np.sum(temp1[UCpos::dof,q,w,dim],axis=0)


                # calculate spectral energy density
                SED[q,w] = (1/(2*(Time*N)**2)) * np.sum(np.abs(temp2[:,q,w,:])**2)
        t3 = time.time()
        print(t3-t2, ' seconds to compute the SED')
        return SED
