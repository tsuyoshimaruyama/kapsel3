#!/usr/bin/python3.6

import h5py
import re
import os
import numpy as np
import pylab as py


def input_h5(filename):
    f = h5py.File(filename, 'r')
    a_group_key = list(f.keys())[0]
    data = np.array(list(f[a_group_key]))
    return (data)


def output_xyz(filename, param, data_x, data_y, data_z):
    f = open(filename, param)
    f.write(str(data_x) + '\t' + str(data_y) + '\t' + str(data_z) + '\n')
    f.close()


if __name__ == "__main__":
    #---------------------------------------
    # set boxsize
    lx = 32
    ly = 32
    lz = 32

    # number of data
    dnum = 11

    # time step interval
    interval = 1

    # average composition
    bar_psi = 0.0
    #---------------------------------------

    dir_name = "particle_fraction"
    os.makedirs(dir_name, exist_ok=True)

    for count in range(0, dnum, interval):

        psi = input_h5("orderparam_" + str(count) + ".h5")
        phi = input_h5("particle_" + str(count) + ".h5")
        phipsi = phi * psi

        n = np.sum(phi > 0.0)
        na = np.sum(phipsi > bar_psi)
        nb = np.sum(phipsi < bar_psi)
        print("%d: (a, b) = (%f, %f)" % (count, na / n, nb / n))
        # output
        output_xyz(dir_name + "/particle_fraction.txt",
                  "a", count, float(na)/float(n), float(nb)/float(n))
