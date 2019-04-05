#!/usr/bin/python3.6

import h5py
import re
import os
import numpy as np
from scipy.fftpack import fftn, ifftn, fftshift
import pylab as py

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
    a_group_key = list(f.keys())[0]
    data = np.array(list(f[a_group_key]))
    return (data)

def output_xy(filename, param,data_x, data_y):
    f = open(filename, param)
    for i in range(0, data_x.size):
        f.write(str(data_x[i]) + '\t' + str(data_y[i]) + '\n')
    f.close()

# simulation parameters
NX = 32
NY = NX
NZ = NX
N = float(NX * NY * NZ)

ix = 0
iy = 0
iz = 0

dir_name = "sq_data"
os.mkdir(dir_name)

for ts in range(0, 11):
    data = input_h5("orderparam_" + str(ts) + ".h5")
    F = fftn(data)
    F[0][0][0] = 0 # DC offset
    AMP = np.abs(F)
    PS = AMP ** 2 / N

    data_x, data_y = spherical_average(PS, NX / 2.0, 1.0)

    data_x = data_x * (2. * np.pi / float(NX))

    cq = calc_characteristic_q(data_x, data_y)
    data_x = data_x / cq
    data_y = data_y * (cq ** 3)

    print("ts: {0} cq: {1}".format(ts, cq))

    # plot
    py.figure(0)
    py.clf()
    py.loglog(data_x, data_y, '.-')
    py.title('$S(q)$ vs $q$ ' + str(ts) + 'steps')
    py.xlabel('$q$')
    py.ylabel('$S(q)$')
    py.show()

    #output
    output_xy(dir_name + "/cq.txt", "a", np.array([ts]), np.array([cq]))
    output_xy(dir_name + "/sq_" + str(ts) + ".txt", "w", data_x, data_y)
