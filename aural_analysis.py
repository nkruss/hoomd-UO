import numpy as np
import matplotlib.pyplot as plt
import math
import time
import os
import shutil
import statistics as stats
import scipy.signal
import scipy.special
from scipy.stats import norm
import pandas as pd
import pywt

import gsd.pygsd as GSD_pygsd
import gsd.hoomd as GSD_hoomd

from scipy.optimize import curve_fit
from astropy.modeling import models, fitting

#fit functions
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

class Aural_Analysis():

    def __init__(self):

        self.particle_data = None

        self.dt = None
        self.m = None
        self.N = None

    def read_data(self, fname: str):

        #open data file
        f = GSD_pygsd.GSDFile(open(fname, 'rb'))
        t = GSD_hoomd.HOOMDTrajectory(f)

        self.particle_data = t
        self.N = len(t[0].particles.position)

    def wave_packet(self, dt, store_loc = None, plot_title = "Gausian plot", num_samples = 8, target_times = None):

        target_index = []
        if target_times == None:
            #select the time to plot the data at
            sample_period = (100000000) / num_samples

            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                if time_step % sample_period == 0:
                    target_index.append(time_step_i)

            target_index.append(len(self.particle_data) - 1)
        else:
            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                time = time_step * dt
                if (time in target_times):
                    target_index.append(time_step_i)

        #-----------------------------------------------------
        #create and save waterfall plot of wave packet
        g1 = plt.figure(1)
        plt.xlabel('particle position')
        plt.ylabel('time')
        plt.yticks([])
        plt.title(plot_title, fontsize = 7)

        shift = 0
        packet_amplitude_list = []
        fit_amp_list = []
        fit_std_list = []
        for i in target_index:
            time_step = i * 500
            #print(f"---{time_step}---")
            p_list = []
            offsets = []
            abs_offsets = []
            amplitude = 0
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos
                abs_position = abs(position)

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos
                        abs_position = abs(position)
                position += shift
                abs_position += shift

                #update amplitude
                if abs(position - shift) > amplitude:
                    amplitude = abs(position - shift)

                p_list.append(p_index)
                offsets.append(position)
                abs_offsets.append(abs_position)
            packet_amplitude_list.append(amplitude)

            plt.plot(p_list, offsets, color = "b")
            # plt.plot(p_list, abs_offsets, color = "g")

            #shift -= .5
            shift -= .175

        if store_loc != None:
            plt.savefig("Gausian_plot")
            shutil.move("Gausian_plot.png", store_loc)
        else:
            plt.show()
        plt.clf()

        df = pd.DataFrame(packet_amplitude_list)
        df.to_excel("amplitudes.xlsx")
        shutil.move("amplitudes.xlsx", store_loc)

        #get fourier plot and data for either sample
        fourier_data_list = []
        for i in target_index:
            time_step = i * 500
            p_list = []
            offsets = []
            gausian = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                p_list.append(p_index)
                offsets.append(position)
                gausian.append(abs(position))

            #plt.plot(p_list, gausian)
            fourier_data = self.fourier_plot(p_list, offsets, store_loc = store_loc, plot_title = f"Fourier_plot_time={int(time_step * .0001)}")
            fourier_data_list.append(fourier_data)


    def gaussian_fitting(self, dt, store_loc = None, plot_title = "Gausian plot", num_samples = 8, target_times = None):

        dimensionless_time = []
        w = 0.7853981633974483

        target_index = []
        if target_times == None:
            #select the time to plot the data at
            sample_period = (100000000) / num_samples

            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                if (time_step % sample_period == 0) and (time_step * dt * w > 200):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * dt * w)
            target_index.append(len(self.particle_data) - 1)
            dimensionless_time.append((len(self.particle_data) - 1) * 500 * w * dt)
        else:
            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                time = time_step * dt
                if (time in target_times):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w * dt)

        #-----------------------------------------------------
        g1 = plt.figure(1)

        shift = 0
        packet_amplitude_list = []
        fit_amp_list = []
        fit_std_list = []
        fit_cent_list = []
        unwrap_counter = 0
        for i in target_index:
            time_step = i * 500
            p_list = []
            offsets = []
            abs_offsets = []
            amplitude = 0
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos
                abs_position = abs(position)

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos
                        abs_position = abs(position)
                position += shift
                abs_position += shift

                #update amplitude
                if abs(position - shift) > amplitude:
                    amplitude = abs(position - shift)

                p_list.append(p_index)
                offsets.append(position)
                abs_offsets.append(abs_position)
            packet_amplitude_list.append(amplitude)


            #setting up gausian fit
            gausian_left_index = 0
            gausian_right_index = 0
            if abs_offsets[0] == shift:
                for i in range(0, len(abs_offsets)):
                    if abs_offsets[i] - shift != 0:
                        gausian_left_index = i
                        break
                for i in range(len(abs_offsets) - 1, 0, -1):
                    if abs_offsets[i] - shift != 0:
                        gausian_right_index = i
                        break
            else:
                count = 0
                for i in range(0, len(abs_offsets)):
                    if abs_offsets[i] - shift == 0:
                        count += 1
                    else:
                        gausian_right_index += 1
                        count = 0
                    if count >= 20:
                        break
                count = 0
                gausian_left_index = len(abs_offsets) - 1
                for i in range(len(abs_offsets) - 1, 0, -1):
                    if abs_offsets[i] - shift == 0:
                        count += 1
                    else:
                        gausian_left_index -= 1
                        count = 0
                    if count >= 20:
                        break

            #create a list of the amplitudes of the gausian without the periodic boundary conditions
            gausian_removed_periodic = []
            if gausian_left_index < gausian_right_index:
                for i in range(gausian_left_index, gausian_right_index, 1):
                    gausian_removed_periodic.append(abs_offsets[i] - shift)
                sigma1 = len(gausian_removed_periodic) / 4
            else:
                for i in range(gausian_left_index, len(abs_offsets), 1):
                    gausian_removed_periodic.append(abs_offsets[i] - shift)
                for i in range(gausian_right_index):
                    gausian_removed_periodic.append(abs_offsets[i] - shift)
                sigma1 = len(gausian_removed_periodic) / 4

            while(len(gausian_removed_periodic) != len(abs_offsets)):
                gausian_removed_periodic.append(0)

            amp1 = max(gausian_removed_periodic)
            cen1 = gausian_removed_periodic.index(amp1)

            popt_gauss, pcov_gauss = scipy.optimize.curve_fit(gauss_function, p_list, gausian_removed_periodic, p0=[amp1, cen1, sigma1])
            fit_amp_list.append(popt_gauss[0])
            fit_std_list.append(popt_gauss[2])

            fit_center = abs_offsets.index(amp1 + shift) + (unwrap_counter * self.N)
            if len(fit_cent_list) > 4:
                if(fit_cent_list[-4] > fit_center and
                   fit_cent_list[-3] > fit_center and
                   fit_cent_list[-2] > fit_center and
                   fit_cent_list[-1] > fit_center):
                    #print(f"wrapping on {abs_offsets.index(amp1 + shift) + (unwrap_counter * self.N)}, prev = {fit_cent_list[-2]}, unwrap_counter = {unwrap_counter}")
                    unwrap_counter += 1;
            fit_cent_list.append(abs_offsets.index(amp1 + shift) + (unwrap_counter * self.N))

            # print(f"{time_step} - ({gausian_left_index}, {gausian_right_index}), {popt_gauss[0]}, {popt_gauss[2]}")
            # plt.plot(p_list, gauss_function(p_list, *popt_gauss) + shift, color = "r")
            # plt.plot(p_list, abs_offsets, color = "g")

            # shift -= .5
            shift -= .175

        if store_loc != None:
            plt.savefig("Gaussian_Fits")
            shutil.move("Gaussian_Fits.png", store_loc)

            df = pd.DataFrame({"Dimensionless Time": dimensionless_time, "Amplitude": fit_amp_list, "STD": fit_std_list, "Center": fit_cent_list})
            df.to_excel("gaussian_fit_parameters.xlsx")
            shutil.move("gaussian_fit_parameters.xlsx", store_loc)
        else:
            plt.show()
        plt.clf()

        # df = pd.DataFrame(fit_amp_list)
        # df_2 = pd.DataFrame(fit_std_list)
        # df_3 = pd.DataFrame(fit_cent_list)
        # with pd.ExcelWriter("gausian_fit_parameters.xlsx") as writer:
        #     df.to_excel(writer, sheet_name = "Amp")
        #     df_2.to_excel(writer, sheet_name = "STD")
        #     df_3.to_excel(writer, sheet_name = "Center")
        # shutil.move("gausian_fit_parameters.xlsx", store_loc)

        #get fit parameter errors
        amp_perc_error_list = []
        std_perc_error_list = []
        for i in range(len(fit_amp_list)):
            amp_perc_error_list.append((fit_amp_list[i] - fit_amp_list[0]) / fit_amp_list[0])
            std_perc_error_list.append((fit_std_list[i] - fit_std_list[0]) / fit_std_list[0] * 100)

        #create error plots of fit parameters
        g1 = plt.figure(1)
        plt.xlabel('Time')
        plt.ylabel('Fit Amplitude Factor')
        plt.title("Gaussian Fit Amplitude Factor")
        plt.plot(dimensionless_time, amp_perc_error_list)
        if store_loc != None:
            plt.savefig("Gaussian_Fit_Amp")
            shutil.move("Gaussian_Fit_Amp.png", store_loc)
        else:
            plt.show()
        plt.clf()

        g1 = plt.figure(1)
        plt.xlabel('Time')
        plt.ylabel('Fit STD Percent Error')
        plt.title("Gaussian Fit STD Error")
        plt.plot(dimensionless_time, std_perc_error_list)
        if store_loc != None:
            plt.savefig("Gaussian_Fit_STD")
            shutil.move("Gaussian_Fit_STD.png", store_loc)
        else:
            plt.show()
        plt.clf()

        g1 = plt.figure(1)
        plt.xlabel('Time')
        plt.ylabel('Fit Center Position')
        plt.title("Gaussian Fit Center Position")
        plt.plot(dimensionless_time, fit_cent_list)
        if store_loc != None:
            plt.savefig("Gaussian_Fit_Center")
            shutil.move("Gaussian_Fit_Center.png", store_loc)
        else:
            plt.show()
        plt.clf()

        m, b = np.polyfit(dimensionless_time, fit_cent_list, 1)
        return(dimensionless_time, fit_cent_list, m, b, amp_perc_error_list)


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

        return (abs(freq), abs(fft.real))

    def RMS_error(self, dt, store_loc = None, plot_title = "Root Mean Squared Error", num_samples = 1000, target_times = None):

        dimensionless_time = []
        timesteps = []
        realtimes = []
        w = 0.7853981633974483

        target_index = []
        if target_times == None:
            #select the time to plot the data at
            sample_period = (100000000) / num_samples

            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                if time_step % sample_period == 0:
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w)
                    timesteps.append(time_step)
                    realtimes.append(time_step * .0001)
            target_index.append(len(self.particle_data) - 1)
            dimensionless_time.append((len(self.particle_data) - 1) * 500 * w)
            timesteps.append((len(self.particle_data) - 1) * 500)
        else:
            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                time = time_step * dt
                if (time in target_times):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w)
                    timesteps.append(time_step)
                    realtimes.append(time_step * dt)

        #get fourier plot and data for either sample
        fourier_data_list = []
        for i in target_index:
            time_step = i * 500
            p_list = []
            offsets = []
            gausian = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                p_list.append(p_index)
                offsets.append(position)
                gausian.append(abs(position))

            #plt.plot(p_list, gausian)
            fourier_data = self.fourier_plot(p_list, offsets)
            fourier_data_list.append(fourier_data)

        #-----------------------------------------------------------------------
        mean_squared_error_list = []
        initial_data = fourier_data_list[0]
        for fourier_data in fourier_data_list:
            diff_sqrd_sum = 0
            for i in range(len(fourier_data[0])):
                diff_sqrd_sum += ((initial_data[1][i] - fourier_data[1][i]) ** 2)
            mean_squared_error_list.append(math.sqrt(diff_sqrd_sum / len(fourier_data[0])))

        error_plot_title = "RMS Error - (over course of Simulation)"
        plt.plot(dimensionless_time, mean_squared_error_list)
        plt.title(plot_title, fontsize = 7)
        plt.savefig(error_plot_title)
        if store_loc != None:
            shutil.move(f"{error_plot_title}.png", store_loc)
        plt.clf()

        df = pd.DataFrame(mean_squared_error_list)
        df.to_excel("RMSE.xlsx")
        shutil.move("RMSE.xlsx", store_loc)

    def peak_error(self, dt, store_loc = None, plot_title = "Mean Peak Error", num_samples = 10, target_times = None):

        dimensionless_time = []
        w = 0.7853981633974483

        target_index = []
        if target_times == None:
            #select the time to plot the data at
            sample_period = (100000000) / num_samples

            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                if time_step % sample_period == 0:
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * dt * w)
            target_index.append(len(self.particle_data) - 1)
            dimensionless_time.append((len(self.particle_data) - 1) * 500 * w * dt)
        else:
            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                time = time_step * dt
                if (time in target_times):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w * dt)

        #get fourier plot and data for either sample
        fourier_data_list = []
        for i in target_index:
            time_step = i * 500
            p_list = []
            offsets = []
            gausian = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                p_list.append(p_index)
                offsets.append(position)
                gausian.append(abs(position))

            #plt.plot(p_list, gausian)
            fourier_data = self.fourier_plot(p_list, offsets)
            fourier_data_list.append(fourier_data)

        #-----------------------------------------------------------------------
        mean_peak_error_list = []
        peaks_initial = scipy.signal.find_peaks(fourier_data_list[0][1], height=.01)
        peak_pos_initial = fourier_data_list[0][0][peaks_initial[0]]
        shift = 15
        for fourier_data in fourier_data_list:
            diff_sum = 0

            for i in range(len(peak_pos_initial)):
                target_peak_index = peaks_initial[0][i]
                target_zone = fourier_data[1][target_peak_index - shift: target_peak_index + shift]
                #print("target zone = ", target_zone)
                peak = np.amax(target_zone)
                #print("peak = ", peak)
                peak_index = np.where(target_zone == peak)
                target_freq_zone = fourier_data[0][target_peak_index - shift:]
                peak_pos = target_freq_zone[peak_index[0][0]]
                #print(peak_pos_initial[i], peak_pos)

                diff_sum += (peak_pos_initial[i] - peak_pos) ** 2
            mean_peak_error_list.append(math.sqrt(diff_sum) / len(peaks_initial))
            #print()

        error_plot_title = "Fourier Peaks Mean Error - (over course of Simulation)"
        plt.plot(dimensionless_time, mean_peak_error_list)
        plt.title(error_plot_title)
        plt.xlabel("Dimensionless Time")
        plt.ylabel("Error")
        plt.savefig(error_plot_title)
        if store_loc != None:
            shutil.move(f"{error_plot_title}.png", store_loc)
        plt.clf()

        # #code for printing the peak points
        # plt.plot(fourier_data_list[0][0], fourier_data_list[0][1])
        # height = peaks_initial[1]['peak_heights']
        # peak_pos = fourier_data_list[0][0][peaks_initial[0]]
        # plt.scatter(peak_pos, height)

        df = pd.DataFrame(mean_peak_error_list)
        df.to_excel("Peak_pos.xlsx")
        shutil.move("Peak_pos.xlsx", store_loc)

    def normalized_error(self, dt, store_loc = None, plot_title = "Normalized Amplitude Error", num_samples = 1000, target_times = None):

        dimensionless_time = []
        w = 0.7853981633974483

        target_index = []
        if target_times == None:
            #select the time to plot the data at
            sample_period = (100000000) / num_samples

            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                #if time_step % sample_period == 0:
                if (time_step % sample_period == 0) and (time_step * dt * w > 200):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * dt * w)
            target_index.append(len(self.particle_data) - 1)
            dimensionless_time.append((len(self.particle_data) - 1) * 500 * w * dt)
        else:
            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                time = time_step * dt
                if (time in target_times):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w * dt)

        #get fourier plot and data for either sample
        fourier_data_list = []
        for i in target_index:
            time_step = i * 500
            p_list = []
            offsets = []
            gausian = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                p_list.append(p_index)
                offsets.append(position)
                gausian.append(abs(position))

            #plt.plot(p_list, gausian)
            fourier_data = self.fourier_plot(p_list, offsets)
            fourier_data_list.append(fourier_data)

        #-----------------------------------------------------------------------
        fractional_error_list = []
        initial_data = fourier_data_list[0]
        #initial_data = fourier_data_list[16]
        max_amp = max(initial_data[1])
        for fourier_data_i in range(len(fourier_data_list)):
            fourier_data = fourier_data_list[fourier_data_i]
            diff_sum = 0
            timestep = fourier_data_i * 500
            for i in range(len(fourier_data[0])):
                diff_sum += abs(initial_data[1][i] - fourier_data[1][i]) / max_amp

                # amp = fourier_data[1][i] / (math.e ** (0.0025 * timestep))
                # diff_sum += abs(initial_data[1][i] - amp) / max_amp

            fractional_error_list.append(diff_sum / len(fourier_data[0]))

        m, b = np.polyfit(dimensionless_time, fractional_error_list, 1)
        fit = []
        for time in dimensionless_time:
            fit.append(m*time + b)
        print(m)
        plt.plot(dimensionless_time, fit)
        plt.legend([f"m = {m}, b = {b}"])

        error_plot_title = "Foureir Amplitude Normalized Difference - (over course of Simulation)"
        plt.plot(dimensionless_time, fractional_error_list)
        plt.title(plot_title, fontsize = 7)
        plt.savefig(error_plot_title)
        if store_loc != None:
            shutil.move(f"{error_plot_title}.png", store_loc)

            df = pd.DataFrame(fractional_error_list)
            df.to_excel("normalized_RMSE.xlsx")
            shutil.move("normalized_RMSE.xlsx", store_loc)
        plt.clf()


        return (fractional_error_list, dimensionless_time)

    def integrety_test(self, dt, num_samples = 10, target_times = None):

        dimensionless_time = []
        w = 0.7853981633974483

        #pull out the indexes of the target times
        target_index = []
        if target_times == None:
            #select the time to plot the data at
            sample_period = (100000000) / num_samples

            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                if time_step % sample_period == 0:
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w)
            target_index.append(len(self.particle_data) - 1)
            dimensionless_time.append((len(self.particle_data) - 1) * 500 * w)
        else:
            for time_step_i in range(len(self.particle_data)):
                time_step = time_step_i * 500
                time = time_step * dt
                if (time in target_times):
                    target_index.append(time_step_i)
                    dimensionless_time.append(time_step * w)

        #get fourier plot and data for either sample
        fourier_data_list = []
        for i in target_index:
            time_step = i * 500
            p_list = []
            offsets = []
            gausian = []
            for p_index in range(self.N):
                equilibriam_pos = (-self.N / 2) + p_index
                position = self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                #check for particle 0 looping around due to boundery conditions
                if p_index == 0:
                    if self.particle_data[i].particles.position[p_index][0] > 0:
                        position = -self.particle_data[i].particles.position[p_index][0] - equilibriam_pos

                p_list.append(p_index)
                offsets.append(position)
                gausian.append(abs(position))

            fourier_data = self.fourier_plot(p_list, offsets)
            fourier_data_list.append(fourier_data)

        #-----------------------------------------------------------------------
        #set up for peak analysis
        mean_peak_error_list = []
        peaks_initial = scipy.signal.find_peaks(fourier_data_list[0][1], height=.01)
        peak_pos_initial = fourier_data_list[0][0][peaks_initial[0]]
        shift = 15

        #set up for amplitude analysis
        mean_applitude_error_list = []
        initial_data = fourier_data_list[0]

        for fourier_data in fourier_data_list:
            peak_diff_sum_sqrd = 0
            amp_diff_sum_sqrd = 0

            #sum sqrd errors on peak positions for current timeframe
            for i in range(len(peak_pos_initial)):
                target_peak_index = peaks_initial[0][i]
                target_zone = fourier_data[1][target_peak_index - shift: target_peak_index + shift]
                peak = np.amax(target_zone)
                peak_index = np.where(target_zone == peak)
                target_freq_zone = fourier_data[0][target_peak_index - shift:]
                peak_pos = target_freq_zone[peak_index[0][0]]

                peak_diff_sum_sqrd += (peak_pos_initial[i] - peak_pos) ** 2

            #sum sqrd errors on fourier amplitude for current timeframe
            for i in range(len(fourier_data[0])):
                amp_diff_sum_sqrd += ((initial_data[1][i] - fourier_data[1][i]) ** 2)

            #add difference values to apropriate lists
            mean_applitude_error_list.append(math.sqrt(amp_diff_sum_sqrd / len(fourier_data[0])))
            mean_peak_error_list.append(math.sqrt(peak_diff_sum_sqrd) / len(peaks_initial))

        return (mean_applitude_error_list, mean_peak_error_list, dimensionless_time)
