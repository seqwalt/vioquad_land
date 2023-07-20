/**
 * @file mavros_cmd_node.cpp
 * @brief Offboard control node for hardware flight using mavros and PX4
 */

#include "mavros_cmd.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "mavros_cmd_node");
    ros::NodeHandle nh;

    MavrosCmd mavCmd;

    // Parameters
    nh.param<bool>("/mavros_cmd_node/enable_sim", mavCmd.sim_enable, false); // simulation or not
    
    string ctrl_mode_str;
    nh.param<string>("/mavros_cmd_node/ctrl_mode", ctrl_mode_str, "position"); // default is position controller if not set from launch file
    if (ctrl_mode_str == "geometric"){
        if (mavCmd.sim_enable){
            mavCmd.ctrl_mode = MavrosCmd::GEOMETRIC;
        } else {
            ROS_ERROR_STREAM("Geometric control not yet supported in hardware. Shutting down mavros_cmd_node ...");
            ros::shutdown();
        }
    } else if (ctrl_mode_str == "position"){
        mavCmd.ctrl_mode = MavrosCmd::POSITION;
    } else if (ctrl_mode_str == "mpc") {
        if (mavCmd.sim_enable){
            mavCmd.ctrl_mode = MavrosCmd::MPC;
        } else {
            ROS_ERROR_STREAM("MPC control not yet supported in hardware. Shutting down mavros_cmd_node ...");
            ros::shutdown();
        }
    } else {
        ROS_ERROR_STREAM("mavros_cmd_node.cpp: No controller type named \'" << ctrl_mode_str << "\'. Shutting down mavros_cmd_node ...");
        ros::shutdown();
    }
    
    nh.param<double>("/mavros_cmd_node/max_err_acc", mavCmd.ctrl.max_err_acc, 20.0);     // largest magnitude (K_pos*err_pos + K_vel*err_vel) can have
    nh.param<double>("/mavros_cmd_node/Kpos_x", mavCmd.ctrl.K_pos[0], 12.0);
    nh.param<double>("/mavros_cmd_node/Kpos_y", mavCmd.ctrl.K_pos[1], 12.0);
    nh.param<double>("/mavros_cmd_node/Kpos_z", mavCmd.ctrl.K_pos[2], 10.0);
    nh.param<double>("/mavros_cmd_node/Kvel_x", mavCmd.ctrl.K_vel[0], 3.0);
    nh.param<double>("/mavros_cmd_node/Kvel_y", mavCmd.ctrl.K_vel[1], 3.0);
    nh.param<double>("/mavros_cmd_node/Kvel_z", mavCmd.ctrl.K_vel[2], 3.3);
    nh.param<double>("/mavros_cmd_node/Katt_x", mavCmd.ctrl.K_att[0], 20.0);   // ang_rate = K_att * err_att
    nh.param<double>("/mavros_cmd_node/Katt_y", mavCmd.ctrl.K_att[1], 20.0);
    nh.param<double>("/mavros_cmd_node/Katt_z", mavCmd.ctrl.K_att[2], 20.0);
    
    // Subscribers
    mavCmd.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &MavrosCmd::mavPoseCallback, &mavCmd, ros::TransportHints().tcpNoDelay());
    mavCmd.vel_sub = nh.subscribe
            ("mavros/local_position/velocity_local", 1, &MavrosCmd::mavVelCallback, &mavCmd, ros::TransportHints().tcpNoDelay());
    mavCmd.mode_sub = nh.subscribe
            ("mavros/state", 1, &MavrosCmd::modeCallback, &mavCmd); // http://wiki.ros.org/roscpp_tutorials/Tutorials/UsingClassMethodsAsCallbacks
    mavCmd.traj_sub = nh.subscribe
            ("reference/flatoutputs", 1, &MavrosCmd::mavRefCallback, &mavCmd, ros::TransportHints().tcpNoDelay());

    // Publishers
    mavCmd.pos_pub = nh.advertise<geometry_msgs::PoseStamped>        // For takeoff, and going to start of trajectory
            ("mavros/setpoint_position/local", 1);
    mavCmd.pos_track_pub = nh.advertise<mavros_msgs::PositionTarget> // For position controller during trajectory tracking
            ("mavros/setpoint_raw/local", 1);
    mavCmd.att_pub = nh.advertise<mavros_msgs::AttitudeTarget>       // For geometric controller during trajectory tracking
            ("mavros/setpoint_raw/attitude", 1);

    // Service clients for changing modes    
    mavCmd.arming_client = nh.serviceClient<mavros_msgs::CommandBool>("mavros/cmd/arming");
    mavCmd.set_mode_client = nh.serviceClient<mavros_msgs::SetMode>("mavros/set_mode");
    // Service client for requesting the starting position and heading of the trajectory
    mavCmd.init_setpnt_client = nh.serviceClient<quad_control::InitSetpoint>("initial_reference");
    // Service client for requesting that the trajectory be published by trajectory_gen_node
    //      OR that the mpc control inputs be published by mpc_node
    mavCmd.streamimg_client = nh.serviceClient<std_srvs::Trigger>("stream_trigger");

    // Timers
    bool autostart = false;
    // High-level logic loop (takoff, landing etc)
    mavCmd.cmdloop_timer = nh.createTimer(ros::Duration(0.1), &MavrosCmd::cmdLoopCallback, &mavCmd, false, autostart); // loop at 10 Hz
    // Setup loop (arm, disarm, offboard mode)
    mavCmd.setup_timer = nh.createTimer(ros::Duration(1), &MavrosCmd::setupCallback, &mavCmd, false, autostart);       // loop at  1 HZ
    
    // wait for FCU connection
    ros::Rate FCU_connect_rate(5.0); // 5 Hz
    while(ros::ok() && !mavCmd.currentModes.connected){
        if (ros::Time::now() - mavCmd.print_request > ros::Duration(3.0)){
            ROS_INFO("Waiting for FCU connection");
            mavCmd.print_request = ros::Time::now();
        }
        ros::spinOnce();
        FCU_connect_rate.sleep();
    }
    
    // Start command and setup loops
    mavCmd.cmdloop_timer.start();
    mavCmd.setup_timer.start();

    ros::spin();
    return 0;
}
