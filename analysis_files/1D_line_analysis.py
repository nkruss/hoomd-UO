"""
Author - Noah Kruss

File that contains the analysis class with the functions for analysising the
1D line test system
"""
#---------------IMPORT STATEMENTS------------------------------
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import os
import shutil
import statistics as stats
import scipy.signal
import scipy.special

#---------------Analysis Class-------------------------------
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
        self.kinetic_energy = []
        self.average_kinetic = []

    # Read Functions
    def read_pos(self, fname: str):
        """
        Read through a simulation position file and record all of the particle
        positions into a dictionary of tuples and store the simulation paramaters
        as analysis object porperties
        """

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
        """
        Read through a simulation velocity file and record all of the particle
        multi-dimensional velociy components into a dictionary of tuples and
        store the simulation paramaters as analysis object porperties
        """

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

    def read_kinetic(self, fname: str):
        """
        Read through a simulation kinetic energy file and record all of
        the enrgy magnitudes into a list
        """

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
                diff = avg1 / avg2
                if abs(1 - diff) < .105:
                    count += 1
                else:
                    count = 0
                if count == 5:
                    break
                avg1 = avg2
                self.last_step = timestep_i

        average = avg(self.kinetic_energy, self.time, start=self.last_step, finish=-1)

        self.stabalization_step = self.last_step
        self.average_kinetic.append(average)

    # Analysis Functions
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

    def graph_spring_line_particle_kinetic(self, plot_title, store_loc):
        """
        Function for plotting the kinetic energy within the spring mass chain by
        individual particle, and plotting the fourier transform of the system
        along with recording the peak location of the foureir transform
        """
        E_k_list = []

        for p_i in range(self.N-1):

            current_E = []
            for time_i in range(int(1500000 / 500), len(self.time)):

                velocity = self.particle_v[p_i][time_i]
                V_mag_squared = (velocity[0]**2) + (velocity[1]**2) + (velocity[2]**2)
                current_E.append(.5 * self.m * V_mag_squared)

            avg_E = stats.mean(current_E)

            E_k_list.append(avg_E)

        plt.plot(E_k_list)
        plt.xlabel('particle')
        plt.ylabel('kinetic energy')
        plt.title(plot_title)

        string = plot_title.split(' ')
        force_freq = string[8].split('.')
        force_freq = force_freq[0] + "pt" + force_freq[1]
        if string[4] == "0,":
            save_title = f"1Dline_static_bonds_{force_freq}_forcefreq"
        else:
            save_title = f"1Dline_{string[4][:-1]}_bonds_{force_freq}_forcefreq"

        plt.savefig(save_title)
        shutil.copy(save_title + ".png", store_loc)

        plt.show()
        plt.clf()

        #---------fourier transform----------------------
        N = self.N - 100
        t = np.arange(100, self.N)
        fft = np.fft.fft(E_k_list[100:])
        T = t[1] - t[0] # sampling interval
        f = np.linspace(0, 1 / T, N)
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency [Hz]")

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

            #create position plot
            if p == particle_i:
                g1 = plt.figure(1)

                if num_bonds == 0:
                    plt.plot(realtime, theoretical)
                    plt.plot(realtime, p_i, "b")
                    plt.xlabel('time')
                elif num_bonds == 1:
                    plt.plot(z_time, p_i, "b")
                    plt.xlabel('z')
                plt.ylabel('position')

                plt.title(f"p: {p}, particle: {particle_i}")
                save_title = f"Position_plot_p_{int(p)}_Particle_{particle_i}"
                plt.savefig(save_title)
                shutil.move(current_loc + "/" + save_title + ".png", save_loc)
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
