import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
import time

def avg(quantity: list, time: list, start=0, finish=-1):
    '''
    Calculates the mean value of the numbers in a list (quantity) over a some number of timesteps,
    by using the sum function to get a sum of the list and then
    deviding that by the number of timesteps

    '''
    if quantity == []: #checks for an empty list
        return 0.
    #print(start, finish, quantity)
    #print(quantity[start:finish])
    quantitySum = math.fsum(quantity[start:finish])
    timeSum = time[finish]
    #print(timeSum)
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
        ratio = (dynamic_avg - static_avg) / dynamic_avg
        ratio_list.append(ratio)

    colors = []

    g1 = plt.figure(1)
    plt.xlabel('angular frequency')
    plt.ylabel('average kinetic energy comparison')
    print(ratio_list)
    plt.scatter(static_data.dynamic_frequencies, ratio_list, c='b')
    plt.show()


class Analysis():

    def __init__(self):
        self.time = []
        self.particle_v = {}
        self.particle_p = {}
        self.stiffness = {}
        self.dt = None
        self.m = None
        self.N = None
        self.kspace_vec = None
        self.wspace_vec = None
        self.kinetic_energy = {}
        self.dynamic_frequencies = []
        self.average_kinetic = []
        self.spring_frequency = None

# Read Functions
    def read_pos(self, fname: str):
        """Read through a simulation position file and record all of the particle
        positions into a dictionary of tuples and store the simulation paramaters as analysis object porperties"""

        #open the position file
        with open(fname, "r") as f:
            #recorded the simulation general infomation
            self.time = []
            line = f.readline().strip().split()
            self.N = int(line[-1])
            #if the x, y replication for a 2 dimensional lattice was recorded store the infomation
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

                #for each particle record it's position as a tuple to the list stored in the dictionary
                #of particles
                for item_i in range(1, len(data)):
                    p_data = data[item_i].strip()
                    p_data = p_data[:-1]
                    p_data = p_data.split(',')
                    self.particle_p[item_i - 1].append((float(p_data[0]),float(p_data[1]),float(p_data[2])))

    def read_velocity(self, fname: str):
        """Read through a simulation velocity file and record all of the particle
        multi-dimensional velociy components into a dictionary of tuples and
        store the simulation paramaters as analysis object porperties"""

        with open(fname, "r") as f:
            self.time = []
            line = f.readline().strip().split()
            self.dt = float(line[5])
            self.m = float(line[7])
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
            run_kinetic_energy = []
            run_time= []

            f.readline()

            for line in f:
                data = line.strip().split()
                #print(data)
                time = data[0].strip()
                run_time.append(float(time))
                kinetic_energy = data[1].strip()
                run_kinetic_energy.append(float(kinetic_energy))
        self.kinetic_energy[run_num] = (run_time,run_kinetic_energy)

        count = 0
        avg1 = 0
        avg2 = 0
        last_step = 0

        for timestep_i in range(1,len(run_time)):
            if run_time[timestep_i] % 1000 == 0:
                avg2 = avg(run_kinetic_energy, run_time, start=last_step, finish=timestep_i)
                diff = abs(avg2 - avg1)
                #print(diff, run_time[timestep_i])
                if diff < .1:
                    count += 1
                else:
                    count = 0

                if count == 10:
                    break
                avg1 = avg2
                last_step = timestep_i

        average = avg(run_kinetic_energy, run_time, start=last_step, finish=-1)
        #print(average, run_time[(last_step-9)])
        #average = avg(run_kinetic_energy, run_time)
        self.average_kinetic.append(average)

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


#Graphing Functions

    def graph_k_energy(self, run_num: int):
        """Creates a graph of the kinetic energy of the system over the timesteps"""

        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('kinetic energy')
        plt.plot(self.kinetic_energy[run_num][0], self.kinetic_energy[run_num][1], "b")
        plt.show()

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

    def graph_p_pos(self, p_index: int):
        """Creates and displays a graph a particle p_index position over the timesteps of the simulation"""
        p_i = []
        for position in self.particle_p[p_index]:
            p_i.append(position[0])
        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('position')
        plt.plot(self.time, p_i, "b")
        plt.show()

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

#?
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
