#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Load trajectory
traj = np.loadtxt("test_traj.csv", delimiter=",", dtype=float)

# Split up data
x   = traj[:,0]
y   = traj[:,1]
z   = traj[:,2]
yaw = traj[:,3]

# Evaluation times
ti = 0
tf = 5
t = np.linspace(ti, tf, len(x))

# 2D Plots

# Position
#fig, axs = plt.subplots(1, 3, figsize=(12, 5))
#axs[0].plot(t, x)
#axs[0].set_title('x (m)')
#axs[1].plot(t, y, 'tab:orange')
#axs[1].set_title('y (m)')
#axs[2].plot(t, z, 'tab:green')
#axs[2].set_title('z (m)')

# Position and yaw
fig, axs = plt.subplots(1, 4, figsize=(16, 5))
axs[0].plot(t, x)
axs[0].set_title('x (m)')
axs[1].plot(t, y, 'tab:orange')
axs[1].set_title('y (m)')
axs[2].plot(t, z, 'tab:green')
axs[2].set_title('z (m)')
axs[3].plot(t, yaw, 'tab:red')
axs[3].set_title('yaw (rad)')

for ax in axs:
    ax.set(xlabel='time (s)')

plt.subplots_adjust(bottom=0.15, left=0.05, right=0.95, top=0.85)

# 3D Plot

ax3d = plt.figure().add_subplot(projection='3d')
ax3d.plot(x, y, z, label='reference')
ax3d.scatter(x[0], y[0], z[0], c='black', marker='o', label='initial position')
ax3d.set_title('Minimum Snap Trajectory')
ax3d.legend()

ax3d.set_xlabel('X')
ax3d.set_ylabel('Y')
ax3d.set_zlabel('Z')

plt.show()
