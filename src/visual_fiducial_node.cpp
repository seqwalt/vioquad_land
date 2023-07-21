/**
 * @file visual_fiducial_node.cpp
 * @brief Node for estimating the postion of the landing pad fiducial in world space
 */

#include "visual_fiducial.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "visual_fiducial_node");
    ros::NodeHandle nh;

    // Parameters
    
    // Create VisFid object
    VisFid vis_fid;
    nh.param<double>("/visual_fiducial_node/tag_smoothing_factor", vis_fid.tag_smoothing_factor, 0.5);
    nh.param<double>("/mpc_node/tran_BC_x", vis_fid.tran_BC_x, 0.108);
    nh.param<double>("/mpc_node/tran_BC_y", vis_fid.tran_BC_y, 0.0);
    nh.param<double>("/mpc_node/tran_BC_z", vis_fid.tran_BC_z, 0.0);
    
    // Publishers
    vis_fid.tag_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("fiducial_pose_est", 1); // estimated apriltag pose in map/world frame
    vis_fid.est_pub = nh.advertise<nav_msgs::Path>
            ("state_estimate", 1);  // estimated pose path (collects mavros/local_position/pose into a path for rviz)

    // Subscribers
    vis_fid.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &VisFid::mavPoseCallback, &vis_fid, ros::TransportHints().tcpNoDelay());
    vis_fid.apriltag_sub = nh.subscribe
            ("tag_detections", 1, &VisFid::mavAprilTagCallback, &vis_fid, ros::TransportHints().tcpNoDelay());
    
    // Timer for publishing various paths (reference, groundtruth)
    double path_freq = 10.0; // publish frequency in Hz
    vis_fid.path_timer = nh.createTimer(ros::Duration(1.0/path_freq), &VisFid::viewPathCallback, &vis_fid, false, true);
            
    cout << fixed;
    cout << setprecision(4);
    
    ros::spin();
    return 0;
}


