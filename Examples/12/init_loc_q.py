#Generating position coordinates of particles (quincke roller):

import random as rand
import numpy as np

radius = 2.0
z_max = 6.55
z_min = 6.45
Lx = 64 # system size (x)
Ly = 64
N = 65 # number of particles

e_dir = 2
# (X, Y, Z) = (0, 1, 2)

uvec_e = np.zeros(3)
uvec_e[e_dir] = 1.0


SIGMA = 2.0 * radius
overlap_length = SIGMA * 1.12
noise_x = Lx/2
noise_y = Ly/2

xp = np.zeros(N)
yp = np.zeros(N)
zp = np.zeros(N)

qs = np.zeros(N)
qx = np.zeros(N)
qy = np.zeros(N)
qz = np.zeros(N)

# location
for ip in range(N):
	while True:
		flag = 0
		xp[ip] = Lx / 2 + noise_x*rand.uniform(-1, 1)
		yp[ip] = Ly / 2 + noise_y*rand.uniform(-1, 1)
		zp[ip] = rand.uniform(z_min, z_max)
		for jp in range(0, ip, 1):
			sx = 1
			sy = 1
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = 1
			sy = 0
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = 1
			sy = -1
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = 0
			sy = 1
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = 0
			sy = 0
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = 0
			sy = -1
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = -1
			sy = 1
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = -1
			sy = 0
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
			sx = -1
			sy = -1
			u = np.array([xp[ip] + sx*Lx, yp[ip] + sy*Ly, zp[ip]]) - np.array([xp[jp], yp[jp], zp[jp]])
			d = np.linalg.norm(u)
			if d <= overlap_length: # if overlap
				flag = 1
				break
		if flag == 0:
			break;

#Quaternion
for ip in range(N):
	x = rand.uniform(0, 2.0 * np.pi)
	qs[ip] = np.cos(x/2.0)
	qx[ip] = uvec_e[0] * np.sin(x/2.0)
	qy[ip] = uvec_e[1] * np.sin(x/2.0)
	qz[ip] = uvec_e[2] * np.sin(x/2.0)

# Output
for ip in range(N):
	print('{{%.2f,%.2f,%.2f}{0.0,0.0,0.0}{%.6f,%.6f,%.6f,%.6f}{0.0,0.0,0.0}}' % (xp[ip],yp[ip],zp[ip], qs[ip],qx[ip],qy[ip],qz[ip]))

