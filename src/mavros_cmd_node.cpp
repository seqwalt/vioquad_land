/**
 * @file mavros_cmd_node.cpp
 * @brief Offboard control node for hardware flight using mavros and PX4
 */

#include "mavros_cmd_node.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "mavros_cmd_node");
    ros::NodeHandle nh;

    MavrosCmd mavCmd;
    
    // Subscribers
    mavCmd.mode_sub = nh.subscribe
            ("mavros/state", 1, &MavrosCmd::modeCallback, &mavCmd); // http://wiki.ros.org/roscpp_tutorials/Tutorials/UsingClassMethodsAsCallbacks
    mavCmd.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &MavrosCmd::mavPoseCallback, &mavCmd, ros::TransportHints().tcpNoDelay());
    mavCmd.vel_sub = nh.subscribe
            ("mavros/local_position/velocity", 1, &MavrosCmd::mavVelCallback, &mavCmd, ros::TransportHints().tcpNoDelay());
    mavCmd.traj_sub = nh.subscribe
            ("reference/flatoutputs", 1, &MavrosCmd::mavRefCallback, &mavCmd, ros::TransportHints().tcpNoDelay());

    // Publishers
    mavCmd.pos_track_pub = nh.advertise<mavros_msgs::PositionTarget> // For control during trajectory tracking
            ("mavros/setpoint_raw/local", 1);
    mavCmd.pos_pub = nh.advertise<geometry_msgs::PoseStamped> // For takeoff, and going to start of trajectory
            ("mavros/setpoint_position/local", 1);
    mavCmd.att_pub = nh.advertise<mavros_msgs::AttitudeTarget>
            ("mavros/setpoint_raw/attitude", 1); // NOTE: setpoint_raw expects NED (while setpoint_attitude expects ENU)
//     mavCmd.att_pub = nh.advertise<geometry_msgs::PoseStamped>
//             ("mavros/setpoint_attitude/attitude", 1);
//     mavCmd.angVel_pub = nh.advertise<geometry_msgs::TwistStamped>
//             ("mavros/setpoint_attitude/cmd_vel", 1);
//     mavCmd.thrust_pub = nh.advertise<mavros_msgs::Thrust>
//             ("mavros/setpoint_attitude/thrust", 1);

    // Service clients for changing modes    
    mavCmd.arming_client = nh.serviceClient<mavros_msgs::CommandBool>("mavros/cmd/arming");
    mavCmd.set_mode_client = nh.serviceClient<mavros_msgs::SetMode>("mavros/set_mode");
    // Service client for requesting the starting position and heading of the trajectory
    mavCmd.init_traj_client = nh.serviceClient<quad_control::InitTraj>("initial_reference");
    // Service client for requesting that the trajectory be published by trajectory_gen_node
    mavCmd.send_traj_client = nh.serviceClient<std_srvs::Trigger>("send_trajectory");

    // Timers
    bool autostart = false;
    mavCmd.cmdloop_timer = nh.createTimer(ros::Duration(0.01), &MavrosCmd::cmdLoopCallback, &mavCmd, false, autostart); // send commands at 100 HZ
    mavCmd.setup_timer = nh.createTimer(ros::Duration(1), &MavrosCmd::setupCallback, &mavCmd, false, autostart); // send commands at 1 HZ

    // Parameters
    nh.param<string>("/mavros_cmd_node/ctrl_mode", mavCmd.ctrl_mode, "geometric");
    nh.param<bool>("/mavros_cmd_node/enable_sim", mavCmd.sim_enable, true); // used in setupCallback
    
    nh.param<double>("/mavros_cmd_node/max_err_acc", mavCmd.ctrl.max_err_acc, 20.0);     // largest magnitude (K_pos*err_pos + K_vel*err_vel) can have
    nh.param<double>("/mavros_cmd_node/thrust_const", mavCmd.ctrl.thrust_const, 0.05);   // normalized thrust = 
    nh.param<double>("/mavros_cmd_node/thrust_offset", mavCmd.ctrl.thrust_offset, 0.1);  //          thrust_const * Acc_des + thrust_offset
                                                                                         // normalized thrust means thrust that has magnitude in rng [0,1]
    nh.param<double>("/mavros_cmd_node/Kpos_x", mavCmd.ctrl.Kpos_x, 12.0);
    nh.param<double>("/mavros_cmd_node/Kpos_y", mavCmd.ctrl.Kpos_y, 12.0);
    nh.param<double>("/mavros_cmd_node/Kpos_z", mavCmd.ctrl.Kpos_z, 10.0);
    nh.param<double>("/mavros_cmd_node/Kvel_x", mavCmd.ctrl.Kvel_x, 3.0);
    nh.param<double>("/mavros_cmd_node/Kvel_y", mavCmd.ctrl.Kvel_y, 3.0);
    nh.param<double>("/mavros_cmd_node/Kvel_z", mavCmd.ctrl.Kvel_z, 3.3);
    nh.param<double>("/mavros_cmd_node/Katt_x", mavCmd.ctrl.Katt_x, 20.0);   // ang_rate = K_att * err_att
    nh.param<double>("/mavros_cmd_node/Katt_y", mavCmd.ctrl.Katt_y, 20.0);
    nh.param<double>("/mavros_cmd_node/Katt_z", mavCmd.ctrl.Katt_z, 20.0);
    
    // wait for FCU connection
    ros::Rate FCU_connect_rate(5.0);
    while(ros::ok() && !mavCmd.currentModes.connected){
        ros::spinOnce();
        FCU_connect_rate.sleep();
    }
    
    // Start command loop and setupInit
    mavCmd.cmdloop_timer.start();
    mavCmd.setup_timer.start();

    ros::spin();
    return 0;
}
