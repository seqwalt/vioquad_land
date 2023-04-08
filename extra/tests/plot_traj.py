#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Load trajectory and keyframes
traj = np.loadtxt("test_data/test_traj.csv", delimiter=",", dtype=float)
keys = np.loadtxt("test_data/pos_keyframes.csv", delimiter=",", dtype=float)

# Split up trajectory data
t   = traj[:,0]
x   = traj[:,1]
y   = traj[:,2]
z   = traj[:,3]
yaw = traj[:,4]

# 2D Plots

# Position and yaw
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
axs[0, 0].plot(t, x)
axs[0, 0].set_title('x (m)')
axs[0, 1].plot(t, y, 'tab:orange')
axs[0, 1].set_title('y (m)')
axs[1, 0].plot(t, z, 'tab:green')
axs[1, 0].set_title('z (m)')
axs[1, 1].plot(t, yaw, 'tab:red')
axs[1, 1].set_title('yaw (rad)')

axs[1, 0].set(xlabel='time (s)')
axs[1, 1].set(xlabel='time (s)')
plt.subplots_adjust(bottom=0.1, left=0.05, right=0.95, top=0.9) # adjust plot margins

# x-z plot
fig, axz = plt.subplots()
axz.plot(x,z,'.', label='FOV constrained traj');
axz.plot(x[0], z[0], c='black', marker='o', label='initial position')
axz.plot(0, 0, c='black', marker='^', label='landmark')
axz.set_xlabel('X')
axz.set_ylabel('Z')
axz.axis('equal')
axz.legend()

# 3D Plot

ax3d = plt.figure().add_subplot(projection='3d')
ax3d.plot(x, y, z, label='reference')
ax3d.scatter(x[0], y[0], z[0], c='black', marker='o', label='initial position')
#ax3d.scatter(keys[0], keys[1], keys[2], c='green', marker='^', label='keyframe')
ax3d.set_title('Minimum Snap Trajectory')
ax3d.legend()
ax3d.set_xlabel('X')
ax3d.set_ylabel('Y')
ax3d.set_zlabel('Z')
ax3d.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z))) # set equal axis ratio

plt.show()
