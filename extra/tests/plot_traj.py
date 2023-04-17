#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Load trajectory and keyframes
traj = np.loadtxt("test_data/test_traj.csv", delimiter=",", dtype=float)
keys = np.loadtxt("test_data/pos_keyframes.csv", delimiter=",", dtype=float)

# Split up trajectory data
t   = traj[:,0]
x = traj[:,1];    y = traj[:,2];    z = traj[:,3];    yaw = traj[:,4]      # position
#vx = traj[:,5];   vy = traj[:,6];   vz = traj[:,7];   vyaw = traj[:,8]      # velocity
#ax = traj[:,9];   ay = traj[:,10];  az = traj[:,11];  ayaw = traj[:,12]     # acceleration
#jx = traj[:,13];  jy = traj[:,14];  jz = traj[:,15];  jyaw = traj[:,16]     # jerk

# 2D Plots

# Position and yaw, velocity, acceleration, snap
titles = ['Position', 'Velocity', 'Acceleration', 'Jerk']
tex_dot = ['', '\dot', '\ddot', '\dddot']
pos_units = ['(m)', '(m/s)', '(m/s^2)', '(m/s^3)']
ang_units = ['(rad)', '(rad/s)', '(rad/s^2)', '(rad/s^3)']
plot_these = [0] # 0 = position, 1 = velocity, 2 = acceleration, 3 = jerk

for i in plot_these:
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    axs[0, 0].plot(t, traj[:, 1+4*i], 'tab:blue')
    axs[0, 0].set_title('$'+ tex_dot[i] +'{x}$ ' + pos_units[i])
    axs[0, 1].plot(t, traj[:, 2+4*i], 'tab:orange')
    axs[0, 1].set_title('$'+ tex_dot[i] +'{y}$ ' + pos_units[i])
    axs[1, 0].plot(t, traj[:, 3+4*i], 'tab:green')
    axs[1, 0].set_title('$'+ tex_dot[i] +'{z}$ ' + pos_units[i])
    axs[1, 1].plot(t, traj[:, 4+4*i], 'tab:red')
    axs[1, 1].set_title('$'+ tex_dot[i] +'{\psi}$ ' + ang_units[i])

    axs[1, 0].set(xlabel='time (s)')
    axs[1, 1].set(xlabel='time (s)')
    plt.subplots_adjust(bottom=0.1, left=0.05, right=0.95, top=0.9) # adjust plot margins
    fig.suptitle(titles[i],fontsize=16)
    plt.gcf().canvas.manager.set_window_title(titles[i])
    
# x-z plot
fig, axz = plt.subplots()
axz.plot(x,z,'.', label='FOV constrained traj');
axz.plot(x[0], z[0], c='black', marker='o', label='initial position')
axz.plot(0, 0, c='black', marker='^', label='landmark')
axz.set_xlabel('X')
axz.set_ylabel('Z')
axz.axis('equal')
axz.legend()
plt.gcf().canvas.manager.set_window_title('X-Z Plot')

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
plt.gcf().canvas.manager.set_window_title('3D Position')

plt.show()
