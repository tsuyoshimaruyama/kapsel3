#!/usr/bin/python3.6

import h5py
import re
import os
import numpy as np
from scipy.fftpack import fftn, ifftn, fftshift

# spherical average
def spherical_average(field, rmax, dr):
    z, y, x = np.indices(field.shape)
    r = np.sqrt(x*x + y*y + z*z)
    imax = int(rmax / dr)

    data_x = np.arange(0)
    data_y = np.arange(0)

    for i in range(0, imax):
        ind = np.where((dr*float(i-0.5) <= r.flat) &
                       (r.flat < dr*(float(i)+0.5)))
        val_ext = field.flat[ind]
        n_ext = val_ext.size
        if (n_ext > 0):
            data_x = np.append(data_x, dr * float(i))
            data_y = np.append(data_y, np.sum(val_ext) / float(n_ext))
    return(data_x, data_y)

def calc_characteristic_q(data_x, data_y):
    qSq = data_x * data_y
    sum_qSq = np.sum(qSq)
    sum_Sq = np.sum(data_y)
    cq = sum_qSq / sum_Sq
    return (cq)

# input data
def input_h5(filename):
    f = h5py.File(filename, 'r')
    # List all groups
    a_group_key = list(f.keys())[0]
    data = np.array(list(f[a_group_key]))
    return (data)

def output_xy(filename, param,data_x, data_y):
    f = open(filename, param)
    for i in range(0, data_x.size):
        f.write(str(data_x[i]) + '\t' + str(data_y[i]) + '\n')
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
    #---------------------------------------

    N = float(lx * ly * lz)

    dir_name = "sq_fluid"
    os.makedirs(dir_name, exist_ok=True)

    for count in range(0, dnum, interval):

        data = input_h5("orderparam_" + str(count) + ".h5")
        # calculate power spectrum and pair correlation function
        F = fftn(data)
        # DC offset
        F[0][0][0] = 0

        AMP = np.abs(F)
        PS = AMP ** 2 / N

        data_x, data_y = spherical_average(PS, lx / 2.0, 1.0)

        data_x = data_x * (2. * np.pi / float(lx))

        cq = calc_characteristic_q(data_x, data_y)
        data_y = data_y / N

        print("%d: cq = %f" % (count, cq))
        #output
        output_xy(dir_name + "/cq.txt", "a", np.array([count]), np.array([cq]))
        output_xy(dir_name + "/sq_" + str(count) + ".txt", "w", data_x, data_y)
