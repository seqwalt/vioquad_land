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
    mavCmd.pos_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mavros/setpoint_position/local", 10);
    mavCmd.att_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mavros/setpoint_attitude/attitude", 1);
    mavCmd.angVel_pub = nh.advertise<geometry_msgs::TwistStamped>
            ("mavros/setpoint_attitude/cmd_vel", 1);
    mavCmd.thrust_pub = nh.advertise<mavros_msgs::Thrust>
            ("mavros/setpoint_attitude/thrust", 1);

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
    // TODO: Implement the sim_enable parameter along with launch file
    //nh_private_.param<bool>("enable_sim", mavCmd.sim_enable, true);
    mavCmd.sim_enable = true; // used in setupCallback
    
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
