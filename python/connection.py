#!/usr/bin/python3.6
import h5py
import itertools
import math
import numpy as np
import re
import os
from multiprocessing import Pool
import collections as cl
import matplotlib.pyplot as plt
import time

def output_xyz(filename, param, data_x, data_y, data_z):
    f = open(filename, param)
    for i in range(0, data_x.size):
        f.write(str(data_x[i]) + '\t' + str(data_y[i]) +
                '\t' + str(data_z[i]) + '\n')
    f.close()

def get_dist(i, j):
    dx = rx[i] - rx[j]
    if dx > 0.5 * lx:
        dx -= lx
    elif dx < - 0.5 * lx:
        dx += lx
    
    dy = ry[i] - ry[j]
    if dy > 0.5 * ly:
        dy -= ly
    elif dy < - 0.5 * ly:
        dy += ly   

    dz = rz[i] - rz[j]
    if dz > 0.5 * lz:
        dz -= lz
    elif dz < - 0.5 * lz :
        dz += lz

    return math.sqrt(dx*dx + dy*dy + dz*dz)

def get_dist_pbc(i, j, gt):
    # Caution: this function is called only when i and j would be connected across pbc y-plane. 
    dx = rx[i] - rx[j]
    dy = ry[i] - ry[j]

    if dy < 0:
        dx += gt * ly
    else:
        dx -= gt * ly

    if dx > 0.5 * lx:
        dx -= lx
    elif dx < - 0.5 * lx:
        dx += lx

    if dy > 0.5 * ly:
        dy -= ly
    elif dy < - 0.5 * ly:
        dy += ly   

    dz = rz[i] - rz[j]
    if dz > 0.5 * lz:
        dz -= lz
    elif dz < - 0.5 * lz :
        dz += lz

    return math.sqrt(dx*dx + dy*dy + dz*dz)
   
def check_pbc(i, j):
    dx = math.fabs(rx[i] - rx[j])
    dy = math.fabs(ry[i] - ry[j])
    dz = math.fabs(rz[i] - rz[j])

    if (dx < lx / 2.0) and (dy < ly / 2.0) and (dz < lz / 2.0) :
        retval = 0
    else:
        retval = 1 
    return retval

def adj(x, i, n):
    c = x + i
    if c < 0:
        ret = c + n
    elif c > n-1:
        ret = c - n
    else:
        ret = c
    return ret

def adj_pbc(x, i, n, pbc):

    if pbc == 0:
        ret = adj(x, i, n)
    else:
        dmyx = x
        if pbc < 0:
            dmyx += int(gt*ly*irmax)
        elif pbc > 0:
            dmyx -= int(gt*ly*irmax)

        if dmyx < 0:
            dmyx += nx
        elif dmyx > nx-1:
            dmyx -= nx
            
        ret = adj(dmyx, i, n)

    return ret

def get_particle_pair(i, adj_list):
    llstsub = []
    plstsub = []

    for m, pbc in adj_list:
        if pbc == 0:
            if get_dist(i,m) < rmax:
                if i > m:
                    llstsub.append([i,m])
                plstsub.append(m)
        else:
            if get_dist_pbc(i,m, gt) < rmax:
                if i > m:
                    llstsub.append([i,m])
                plstsub.append(m)

    return (llstsub, plstsub)          

def wrapper_get_particle_pair(args):
    return get_particle_pair(*args)

def cluster_labeling(vertex_flag, adj_particle_list):
    for i in pn:
        vertex_flag[i] = i
    
    while True:
        c = 0
        for i in pn:
            if adj_particle_list[i]: # is not empty?
                assoc_min = min([vertex_flag[j] for j in adj_particle_list[i]])
                if i < assoc_min:
                    if vertex_flag[i] != i:
                        vertex_flag[i] = i
                        c = 1
                else:
                    if vertex_flag[i] != assoc_min:
                        vertex_flag[i] = assoc_min
                        c = 1
        if c == 0:
            break
    
    
    cl_cluster = cl.Counter(vertex_flag)
    mc_cluster = cl_cluster.most_common()
    values, counts = zip(*cl_cluster.most_common())

    pnum = 0.0 # particle number
    cnum = len(counts) # cluster number
    wnum = 0.0
    for i in range(len(counts)):
        pnum = pnum + counts[i]
        wnum = wnum + counts[i] * counts[i]
    
    c_ave_n = float(pnum) / float(cnum)  # number average cluster weight
    c_ave_w = float(wnum) / float(pnum)  # weight average cluster weight
    print('cluster average: (number average, weight average) = (%f, %f) ' % (c_ave_n, c_ave_w))
    return (c_ave_n, c_ave_w) # weight_average

def make_xdmf(filename):
    f = open(filename, "w")
    f.write(r'<?xml version="1.0" ?>' + '\n')
    f.write(r'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' + '\n')
    f.write(r'<Xdmf Version="2.0">' + '\n')
    f.write(r'<Domain>' + '\n')
    f.write(r'<Grid Name="CellTime" GridType="Collection" CollectionType="Temporal">' + '\n')

    for count in range (0, dnum, interval):
        f.write('\n')
        f.write(r'<Grid Name="Structured grid" GridType="Uniform">' + '\n')
        f.write(r'<Time Value="'+ time_list[count] + r'" />' + '\n')
        if(nbond_list[int(count/interval)] > 0):
            f.write(r'<Topology Name="tp" TopologyType="Polyline" NumberOfElements="' + str(nbond_list[int(count/interval)]) + r'">' + '\n')
            f.write(r'<DataItem Format="HDF" Dimensions="' + str(nbond_list[int(count/interval)]) + r' 2">' + '\n')
            f.write(r'link_data_' + str(count) + r'.h5:/BOND' + '\n')
            f.write(r'</DataItem>' + '\n')

            f.write(r'</Topology>' + '\n')
        else:
            f.write(r'<Topology Name="tp" TopologyType="Polyline" NumberOfElements="' + str(1) + r'">' + '\n')
            f.write(r'<DataItem Dimensions="' + str(1) + r' 2">' + '\n')
            f.write('0 0\n')
            f.write(r'</DataItem>' + '\n')
            f.write(r'</Topology>' + '\n')            
        f.write('\n')
        f.write(r'<Geometry GeometryType="XYZ">' + '\n')
        f.write(r'<DataItem Dimensions="' + str(N) + r' 3" NumberType="Float" Precision="8" Format="HDF">' + '\n')
        f.write(r'link_data_' + str(count) + r'.h5:/POSITION' + '\n')
        f.write(r'</DataItem>' + '\n')
        f.write(r'</Geometry>' + '\n')

        #----
        f.write('\n')
        f.write(r'<Attribute Name = "bond flag" AttributeType = "Scalar">' + '\n')
        f.write(r'<DataItem Dimensions="' + str(N) + r' 1" NumberType="Float" Precision="8" Format="HDF">' + '\n')
        f.write(r'link_data_' + str(count) + r'.h5:/BOND_FLAG' + '\n')
        f.write(r'</DataItem>' + '\n')
        f.write(r'</Attribute>' + '\n')
        f.write('\n')
        f.write(r'<Attribute Name = "vertex flag" AttributeType = "Scalar">' + '\n')
        f.write(r'<DataItem Dimensions="' + str(N) + r' 1" NumberType="Float" Precision="8" Format="HDF">' + '\n')
        f.write(r'link_data_' + str(count) + r'.h5:/VERTEX_FLAG' + '\n')
        f.write(r'</DataItem>' + '\n')
        f.write(r'</Attribute>' + '\n')
        #----
        f.write(r'</Grid>' + '\n')

    f.write('\n')
    f.write(r'</Grid>' + '\n')
    f.write(r'</Domain>' + '\n')
    f.write(r'</Xdmf>' + '\n')

def make_xdmf_particle(filename):
    f = open(filename, "w")
    f.write(r'<?xml version="1.0" ?>' + '\n')
    f.write(r'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>' + '\n')
    f.write(r'<Xdmf Version="2.0">' + '\n')
    f.write(r'<Domain>' + '\n')
    f.write(r'<Grid Name="CellTime" GridType="Collection" CollectionType="Temporal">' + '\n')

    for count in range (0, dnum, interval):
        f.write('\n')
        f.write(r'<Grid Name="Structured grid" GridType="Uniform">' + '\n')
        f.write(r'<Time Value="'+ time_list[count] + r'" />' + '\n')
        f.write(r'<Topology Name="tp" TopologyType="Polyvertex" NumberOfElements="' + str(N) + r'"/>' + '\n')
        f.write(r'<Geometry GeometryType="XYZ">' + '\n')
        f.write(r'<DataItem Dimensions="' + str(N) + r' 3" NumberType="Float" Precision="8" Format="HDF">' + '\n')
        f.write(r'link_data_' + str(count) + r'.h5:/POSITION' + '\n')
        f.write(r'</DataItem>' + '\n')
        f.write(r'</Geometry>' + '\n')
        f.write('\n')
        f.write(r'<Attribute Name = "particle flag" AttributeType = "Scalar">' + '\n')
        f.write(r'<DataItem Dimensions="' + str(N) + r' 1" NumberType="Float" Precision="8" Format="HDF">' + '\n')
        f.write(r'link_data_' + str(count) + r'.h5:/VERTEX_FLAG' + '\n')
        f.write(r'</DataItem>' + '\n')
        f.write(r'</Attribute>' + '\n')
        f.write(r'</Grid>' + '\n')

    f.write('\n')
    f.write(r'</Grid>' + '\n')
    f.write(r'</Domain>' + '\n')
    f.write(r'</Xdmf>' + '\n')

def get_time(path):
    with open(path) as f:
        lines = f.readlines()
    
    lines_strip = [line.strip() for line in lines]
    time_ext = [line for line in lines_strip if 'Time Value' in line]
    time_ext_j = ' '.join(time_ext) 

    time_list = re.findall('\"(.+?)\"' , time_ext_j)
    return(time_list)


if __name__ == "__main__":
    #---------------------------------------
    # set boxsize
    lx = 32
    ly = 32
    lz = 32

    # set cut_off distance
    rmax = 2.0 ** (1.0 / 6.0) * 4 * 1.0

    # number of data
    dnum = 11

    # time step interval
    interval = 1

    # number of threads
    threadnum = 8
    #---------------------------------------
    
    dir_name = "connection"
    os.makedirs(dir_name, exist_ok=True)

    ilx = 1.0 / lx
    ily = 1.0 / ly
    ilz = 1.0 / lz

    irmax = 1.0 / (rmax)

    nx = int(lx * irmax)
    ny = int(ly * irmax)
    nz = int(lz * irmax)

    inx = 1.0 / float(nx)
    iny = 1.0 / float(ny)
    inz = 1.0 / float(nz)



    nbond_list = []
    nbond_list_minpath = []
    nbond_list_pos_minpath = []

    # loop--------------------------
    for count in range(0, dnum, interval):

        sw_start = time.time()

        # input hdf5
        fnamein = 'particle_data_' + str(count) + '.h5'
        fin = h5py.File(fnamein, 'r')

        # output hdf5
        fnameout = dir_name + '/link_data_' + str(count) + '.h5'
        fout = h5py.File(fnameout, 'w')

        rx = fin['RX']
        ry = fin['RY']
        rz = fin['RZ']
        pos = []
        pos.append(rx)
        pos.append(ry)
        pos.append(rz)
        tp = list(zip(*pos))

        # get number of particles
        N = rx.shape[0]

        # degree_oblique
        degree_oblique = fin['DEGREE_OBLIQUE']
        gt = degree_oblique[0]

        particle_position = [0] * N
        vertex_flag = [0] * N
        bond_flag = [0] * N
        cell_list = [[[[] for k in range(nz)] for j in range(ny)] for i in range(nx)]
        pn = list(range(N))

        # place particles in sub cells
        for i in pn:
            x = int(rx[i] * ilx * nx)
            y = int(ry[i] * ily * ny)
            z = int(rz[i] * ilz * nz)
            cell_list[x][y][z].append(i)
            particle_position[i] = (x,y,z)

        # set particle pair candidates in args
        args = []
        for i in pn:
            x = particle_position[i][0]
            y = particle_position[i][1]
            z = particle_position[i][2]
            
            adj_list = []
            for j, k, l in itertools.product(range(-1,2), range(-1,2), range(-1,2)):
                if(y == 0 and k == -1):
                    pbc = -1
                elif (y == ny-1 and k == 1):
                    pbc = 1
                else:
                    pbc = 0                
                
                for m in cell_list[adj_pbc(x,j,nx,pbc)][adj(y,k,ny)][adj(z,l,nz)]:
                    adj_list.append((m, pbc))
       
            args.append((i, adj_list))

        # get particle pair (multi processing)
        p = Pool(threadnum)
        res_pp = p.map(wrapper_get_particle_pair, args)
        p.close()

        # format results
        bond_list = []
        adj_particle_list = []
        for i in pn:
            bond_list.extend(res_pp[i][0])
            adj_particle_list.append(res_pp[i][1])

        # cluster labeling function
        c_ave_n, c_ave_w = cluster_labeling(vertex_flag, adj_particle_list)
        output_xyz(dir_name + '/cluster_size.txt', 'a',
                  np.array([count]), np.array([c_ave_n]), np.array([c_ave_w]))
        
        # deep copy vertex_flag -> bond_flag
        bond_flag = vertex_flag[:]


        # delete bonds which cross over pbc for visualization
        bond_list_pbc_del = []
        for i, j in bond_list:
            if (check_pbc(i, j) == 0): 
                bond_list_pbc_del.append((i,j))

        # write data in hdf5 list
        fout.create_dataset('BOND', data=bond_list_pbc_del)
        fout.create_dataset('POSITION', data=tp)
        fout.create_dataset('VERTEX_FLAG', data=vertex_flag)
        fout.create_dataset('BOND_FLAG', data=bond_flag)

        # make list of nbond
        nbond_list.append(len(bond_list_pbc_del))

        elapsed_time = time.time() - sw_start
        print('%d steps: timer %f sec' % (count, elapsed_time))

        print('-------------------')
    # loop--------------------------

    time_list = []
    time_list = get_time('particle_phase.xmf')
    make_xdmf(dir_name + '/link.xmf')
    make_xdmf_particle(dir_name + '/link_p.xmf')
