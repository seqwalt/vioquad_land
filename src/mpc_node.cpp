/**
 * @file mpc_node.cpp
 * @brief Model predictive control node for hardware flight using mavros and PX4
 */

#include "mpc.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "mpc_node");
    ros::NodeHandle nh;

    MPC mpc_ctrl;

    // Publisher of control input
    //mpc_ctrl.mpc_pub = nh.advertise<quad_control::FlatOutputs>("reference/flatoutputs", 1);
    mpc_ctrl.mpc_pub = nh.advertise<mavros_msgs::AttitudeTarget>
            ("mavros/setpoint_raw/attitude", 1);
    
    // Subscribers
    mpc_ctrl.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &MPC::mavPoseCallback, &mpc_ctrl, ros::TransportHints().tcpNoDelay());
    mpc_ctrl.vel_sub = nh.subscribe
            ("mavros/local_position/velocity", 1, &MPC::mavVelCallback, &mpc_ctrl, ros::TransportHints().tcpNoDelay());
    
    // Timer that publishes control inputs at 100 Hz
    double freq = 50; // Hz
    bool autostart = false;
    mpc_ctrl.mpc_timer = nh.createTimer(ros::Duration(1./freq), &MPC::mpcCallback, &mpc_ctrl, false, autostart);

    // Service server for sending the starting position and heading of the trajectory
    ros::ServiceServer init_mpc_server = nh.advertiseService("initial_reference", &MPC::initRefCallback, &mpc_ctrl);

    // Service server for starting the traj_timer (i.e. to start publishing the control inputs)
    ros::ServiceServer send_mpc_server = nh.advertiseService("stream_trigger", &MPC::streamMpcCallback, &mpc_ctrl);

    ros::spin();
    return 0;
}

