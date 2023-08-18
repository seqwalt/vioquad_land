/**
 * @file trajectory_gen_node.cpp
 * @brief Trajectory generation node for hardware flight using mavros and PX4
 */

#include "trajectory_gen.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "trajectory_gen_node");
    ros::NodeHandle nh;

    TrajectoryGen mavTraj;

    // Get trajectory
    string trajFile;
    nh.getParam("/trajectory_gen_node/trajectory_file", trajFile);
    mavTraj.csvToEigen(trajFile); // store trajectory in an Eigen matrix

    // Trajectory publisher
    mavTraj.traj_pub = nh.advertise<vioquad_land::FlatOutputs>("reference/flatoutputs", 1);

    // Timer that publishes setpoints at 100 Hz
    double freq = 100; // Hz
    bool autostart = false;
    mavTraj.traj_timer = nh.createTimer(ros::Duration(1./freq), &TrajectoryGen::trajCallback, &mavTraj, false, autostart);

    // Service server for sending the starting position and heading of the trajectory
    ros::ServiceServer init_traj_server = nh.advertiseService("initial_reference", &TrajectoryGen::initRefCallback, &mavTraj);

    // Service server for starting the traj_timer (i.e. to start publishing the trajectory)
    ros::ServiceServer send_traj_server = nh.advertiseService("stream_trigger", &TrajectoryGen::sendTrajCallback, &mavTraj);

    ros::spin();
    return 0;
}
