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

    // Parameters
    string searchTrajFile;
    nh.getParam("/mpc_node/search_traj_file", searchTrajFile);
    
    // Create MPC object
    MPC mpc_ctrl(searchTrajFile);
    nh.getParam("/mpc_node/enable_sim", mpc_ctrl.sim_enable);
    
    // Publishers
    mpc_ctrl.mpc_pub = nh.advertise<mavros_msgs::AttitudeTarget>
            ("mavros/setpoint_raw/attitude", 1); // publish control inputs
    mpc_ctrl.ref_total_pub = nh.advertise<nav_msgs::Path>
            ("mpc_total_reference", 1);    // total reference trajectory for mpc to track
    mpc_ctrl.ref_curr_pub = nh.advertise<nav_msgs::Path>
            ("mpc_curr_reference", 1);    // reference trajectory for current mpc iteration
    mpc_ctrl.pred_pub = nh.advertise<nav_msgs::Path>
            ("mpc_prediction", 1);    // reference trajectory for current mpc iteration     
    mpc_ctrl.est_pub = nh.advertise<nav_msgs::Path>
            ("state_estimate", 1); // estimated pose path
    mpc_ctrl.tag_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mpc_april_tag", 1); // estimated apriltag pose in map/world frame

    // Subscribers
    mpc_ctrl.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &MPC::mavPoseCallback, &mpc_ctrl, ros::TransportHints().tcpNoDelay());
    mpc_ctrl.vel_sub = nh.subscribe
            ("mavros/local_position/velocity_local", 1, &MPC::mavVelCallback, &mpc_ctrl, ros::TransportHints().tcpNoDelay());
    mpc_ctrl.imu_sub = nh.subscribe
            ("mavros/imu/data_raw", 1, &MPC::mavIMUCallback, &mpc_ctrl, ros::TransportHints().tcpNoDelay());
    mpc_ctrl.apriltag_sub = nh.subscribe
            ("tag_detections", 1, &MPC::mavAprilTagCallback, &mpc_ctrl, ros::TransportHints().tcpNoDelay());
            
    // Timer for publishing control inputs
    double mpc_freq = 50.0; // publish frequency in Hz
    bool autostart = false;
    mpc_ctrl.mpc_timer = nh.createTimer(ros::Duration(1.0/mpc_freq), &MPC::mpcCallback, &mpc_ctrl, false, autostart);
    
    // Timer for publishing various paths (reference, groundtruth)
    double path_freq = 10.0; // publish frequency in Hz
    mpc_ctrl.path_timer = nh.createTimer(ros::Duration(1.0/path_freq), &MPC::viewPathCallback, &mpc_ctrl, false, true);

    // Service server for sending the starting position and heading of the trajectory
    ros::ServiceServer init_mpc_server = nh.advertiseService("initial_reference", &MPC::initRefCallback, &mpc_ctrl);

    // Service server for starting the traj_timer (i.e. to start publishing the control inputs)
    ros::ServiceServer send_mpc_server = nh.advertiseService("stream_trigger", &MPC::streamMpcCallback, &mpc_ctrl);

    cout << fixed;
    cout << setprecision(4);
    
    ros::spin();
    return 0;
}

