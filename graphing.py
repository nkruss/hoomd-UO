import matplotlib.pyplot as plt
from matplotlib import interactive
import numpy as np

def create_spring_pos_graph(fname1: str, fname2: str):
    timestep = 0
    time = []
    p1_1 = []
    p2_1 = []
    p3_1 = []
    counter = 1
    with open(fname1) as f:
        for line in f:
            if counter == 3:
                data = line.strip().split()
                p1_1.append(float(data[0]))
            if counter == 4:
                data = line.strip().split()
                p2_1.append(float(data[0]))
            if counter == 5:
                data = line.strip().split()
                p3_1.append(float(data[0]))
                time.append(timestep)
                timestep += 50
                counter = 0
            counter += 1
    p1_2 = []
    p2_2 = []
    p3_2 = []
    counter = 1
    with open(fname2) as f:
        for line in f:
            if counter == 3:
                data = line.strip().split()
                p1_2.append(float(data[0]))
            if counter == 4:
                data = line.strip().split()
                p2_2.append(float(data[0]))
            if counter == 5:
                data = line.strip().split()
                p3_2.append(float(data[0]))
                counter = 0
            counter += 1

    dif_1 = []
    dif_2 = []
    dif_3 = []

    for i in range(len(p1_1)):
        dif = p1_1[i] - p1_2[i]
        dif_1.append(dif)
        dif = p2_1[i] - p2_2[i]
        dif_2.append(dif)
        dif = p3_1[i] - p3_2[i]
        dif_3.append(dif)

    plt.figure()
    plt.xlabel('time step')
    plt.ylabel("difference")
    plt.axis([0,50000,-.5,.5])
    plt.plot(time, dif_1)
    plt.show()

    """if show1:
        g1 = plt.figure(1)
        plt.xlabel('time step')
        plt.ylabel('position')
        plt.plot(time, p1, "b")

    if show2:
        g2 = plt.figure(2)
        plt.xlabel('time step')
        plt.ylabel('position')
        plt.plot(time, p2, "r")

    if show3:
        g3 = plt.figure(3)
        plt.xlabel('time step')
        plt.ylabel('position')
        plt.plot(time, p3, "g")

    plt.show()"""


def main():
    create_spring_pos_graph("100_rounded", "100_exact")
    create_spring_pos_graph("1000_rounded", "1000_exact")

main()
