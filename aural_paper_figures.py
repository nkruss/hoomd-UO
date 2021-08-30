"""
Author - Noah Kruss

File that contains functions for creating the final figures for the Aural
metamaterial paper
"""

#---------------IMPORT STATEMENTS------------------------------
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
from matplotlib import ticker
import numpy as np
import os
import shutil
import sys

from analysis_files.aural_analysis import *

#------------------Set Figure Parameters---------------------------
fig_width = 3.46  # width in inches : 85 mm one-column format
fig_height = fig_width*3./4.      # height in inches
fig_size =  (fig_width,fig_height)
params = {
          'axes.labelsize': 9,
          'axes.labelpad':1,
          'axes.linewidth':.1,
          'font.size': 11,
          'legend.fontsize': 9,
          'xtick.labelsize': 7,
          'ytick.labelsize': 7,

          'figure.subplot.left' : 0.16,
          'figure.subplot.right' : 0.97,
          'figure.subplot.bottom' : 0.13,
          'figure.subplot.top' : 0.9,
          #~ 'figure.subplot.wspace' : 0.3,

          'legend.labelspacing' : 0.2,
          'legend.handletextpad' : 0.4,
          'legend.borderpad' : 0.,
          'legend.frameon': False,
          'legend.numpoints' : 1,
          'legend.columnspacing' : 1,

          'lines.linewidth' : 1.4,
          'lines.markersize' : 4,
          'lines.markeredgewidth' : .4,

          'xtick.major.pad' : 2,
          'ytick.major.pad' : 2,

          'figure.figsize': fig_size,

          #'text.usetex' : True,
          'font.family' : 'sans-serif',
          'font.serif' : 'DejaVu Sans',
          'mathtext.fontset' : 'dejavusans',
          # 'mathtext.it' : 'serif:italic',
          # 'mathtext.default' : 'it',

          'xtick.major.size' : 2,
          'ytick.major.size' : 2
          }
plt.rcParams.update(params)

#--------------------Helper Functions-----------------------------------
def normalize(z):
    z = z.copy()
    z -= z.min()
    z /= z.max()
    return z

#-------------------Figure Creation Functions--------------------------------
def gaussian_fiting_amp_analysis(dt = 0.0001):
    """
    Function for generating a plot of the fit amplitude error
    from gaussian fitting the system wave packet of various simulation trials
    """

    try:
        os.remove(os.getcwd() + "/aural_metamaterial/Gaussian_Fit_Amp_All_Data.pdf")
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

        with open(cond_file, "r") as f:
            time_cond = f.readline().strip()
            spring_cond = f.readline().strip()
            damp_cond = f.readline().strip()
            f.close()

        data = Aural_Analysis()
        data.read_data(gsd_file)

        gaussian_data = data.gaussian_fitting(dt, num_samples = 1000)
        gaussian_results.append(gaussian_data)

        xs.append(gaussian_data[0])
        ys.append(gaussian_data[4])
        yints.append(gaussian_data[4][500])

    ##set up figure axies and size
    plt.clf()
    plt.figure(figsize=fig_size)
    plt.xlabel('Time')
    plt.ylabel('Fit Amplitude Factor')
    #plt.title("Gaussian Fitting Amplitude Factor")

    ##plot data
    plt.plot(xs[0], ys[0], label=label_list[0], color='k')

    cmap = plt.get_cmap('winter')
    yints = np.asarray(yints)
    for x, y, color, label, dash in zip(xs[1:], ys[1:], normalize(yints[1:]), label_list[1:], dash_list[1:]):
        plt.plot(x, y, label=label, color=cmap(color), dashes=dash)
    plt.legend(loc = 'upper left', fontsize = 9)

    ##save figure
    plt.savefig('Gaussian_Fit_Amp_All_Data.pdf')
    shutil.move("Gaussian_Fit_Amp_All_Data.pdf", os.getcwd() + f"/aural_metamaterial")
