#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description='Plot a trajectory.')
parser.add_argument('trajectory_file', nargs=1, help='path to flatoutput trajectory csv')
parser.add_argument('-k', '--keys', metavar='keys_file', nargs=1, help='path to keyframes csv')
parser.add_argument('-l', '--landmark', metavar=('x','y','z'), type=float, nargs=3, help='3D position of landmark')
args = parser.parse_args()

# Load trajectory
traj_file = args.trajectory_file[0]
traj = np.loadtxt(traj_file, delimiter=",", dtype=float)
if traj.shape[1] != 17:
    print("Error: num columns != 17 (trajectory_file needs to contain flatoutput data.)")
    sys.exit()

# Load keyframes
if args.keys:
    keys_file = args.keys[0]
    keys = np.loadtxt(keys_file, delimiter=",", dtype=float)
    x_keys = keys[:,0]; y_keys = keys[:,1]; z_keys = keys[:,2];  # load in keyframes

# Load landmark
if args.landmark:
    landmark = args.landmark

if args.landmark or args.keys:
    x = traj[:,1];  y = traj[:,2];  z = traj[:,3]

# Trajectory times
t   = traj[:,0]

# 2D Plots

# Position and yaw, velocity, acceleration, snap
titles = ['Position', 'Velocity', 'Acceleration', 'Jerk']
tex_dot = ['', '\dot', '\ddot', '\dddot']
pos_units = ['(m)', '(m/s)', '(m/s^2)', '(m/s^3)']
ang_units = ['(rad)', '(rad/s)', '(rad/s^2)', '(rad/s^3)']
plot_these = [0, 1] # 0 = position, 1 = velocity, 2 = acceleration, 3 = jerk

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

if args.landmark:
    # x-z plot
    fig, axz = plt.subplots()
    axz.plot(x,z,'.', label='FOV constrained traj');
    axz.plot(x[0], z[0], c='black', marker='o', label='initial position')
    axz.plot(landmark[0], landmark[2], c='black', marker='^', label='landmark')
    axz.set_xlabel('X')
    axz.set_ylabel('Z')
    axz.axis('equal')
    axz.legend()
    plt.gcf().canvas.manager.set_window_title('X-Z Plot')

    # 3D Plot
if args.keys:
    ax3d = plt.figure().add_subplot(projection='3d')
    ax3d.plot(x, y, z, label='reference')
    ax3d.scatter(x[0], y[0], z[0], c='black', marker='o', label='initial position')
    ax3d.scatter(x_keys, y_keys, z_keys, c='green', marker='^', label='keyframes')
    ax3d.set_title('Minimum Snap Trajectory')
    ax3d.legend()
    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlabel('Z')
    ax3d.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z))) # set equal axis ratio
    plt.gcf().canvas.manager.set_window_title('3D Position')

plt.show()
