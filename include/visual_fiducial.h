#pragma once

#include <ros/ros.h>
#include <mavros/mavros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <apriltag_ros/AprilTagDetectionArray.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <random>

using namespace std;

class VisFid {
    public:
        ros::Publisher est_pub;
        ros::Publisher tag_pub;
        ros::Subscriber pose_sub;
        ros::Subscriber apriltag_sub;
        ros::Timer path_timer;
        
        // VisFid class constructor
        VisFid(){
            tag_visible = false;
        }

        // Publish tragectories and poses for rviz visualization
        void viewPathCallback(const ros::TimerEvent &event){
            if (tag_visible){
                // Publish estimated AprilTag pose in map/world frame
                tag_pub.publish(tagPose);
            }
            // Publish estimated path for visualization (vio or ekf2, depending on direct_from_vio param)
            addPoseToPath(curPose, pathEst);
            est_pub.publish(pathEst);
        }

        // Pass along current mav pose data, wrt world frame.
        void mavPoseCallback(const geometry_msgs::PoseStamped &msg) {
            curPose = msg.pose;
            Eigen::Quaterniond quadQuat = msg2quaternion(curPose.orientation);  // convert quadQuat from message to Eigen::Quateniond
            curYaw = quaternion2yaw(quadQuat);
        }
        
        // Compute AprilTag pose in World frame (compute transformation T_WA)
        void mavAprilTagCallback(const apriltag_ros::AprilTagDetectionArray &msg) {
            // Check if tag is visible
            tag_visible = !msg.detections.empty();
            if (tag_visible){
                // Form T_WB (load from state estimation: mocap, VIO, etc.)
                Eigen::Isometry3d T_WB = Eigen::Isometry3d::Identity();
                Eigen::Vector3d tran_WB(curPose.position.x, curPose.position.y, curPose.position.z); // translation part
                Eigen::Quaterniond quat_WB = msg2quaternion(curPose.orientation); // [w, x, y, z] rotation part
                T_WB.translate(tran_WB);
                T_WB.rotate(quat_WB);
                
                // Form T_BC (transformation from quadcopter body frame to down-facing camera frame)
                // TODO: load T_BC info in from launch file
                Eigen::Isometry3d T_BC = Eigen::Isometry3d::Identity();
                Eigen::Vector3d tran_BC(0.108, 0.0, 0.0); // translation part
                Eigen::Matrix3d rot_BC;
                // TODO: fix sdf to apriltag stuff
                //vector<double> sdf_rpy{0.0, 1.5708, 0.0}; // roll pitch yaw (XYZ form), given by iris_downward_camera.sdf
                //vector<double> apriltag_rpy{sdf_rpy[2], -sdf_rpy[0], sdf_rpy[1]}; // rotation from sdf to apriltag representaiton (https://github.com/AprilRobotics/apriltag/wiki/AprilTag-User-Guide#coordinate-system).
                vector<double> apriltag_rpy{3.1415, 0.0, -1.5708};
                rot_BC = Eigen::AngleAxisd(apriltag_rpy[2], Eigen::Vector3d::UnitZ())
                       * Eigen::AngleAxisd(apriltag_rpy[1], Eigen::Vector3d::UnitY())
                       * Eigen::AngleAxisd(apriltag_rpy[0], Eigen::Vector3d::UnitX()); // XYZ rotation http://sdformat.org/tutorials?tut=specify_pose
                T_BC.translate(tran_BC);
                T_BC.rotate(rot_BC);
                
                // Compute T_CA (transformation from down-facing camera frame to AprilTag frame)
                geometry_msgs::Pose AprilTagPose_Cam = msg.detections[0].pose.pose.pose; // pose of AprilTag in camera frame
                Eigen::Vector3d tran_CA(AprilTagPose_Cam.position.x, AprilTagPose_Cam.position.y, AprilTagPose_Cam.position.z); // translation part
                Eigen::Quaterniond quat_CA = msg2quaternion(AprilTagPose_Cam.orientation); // rotation part
                Eigen::Isometry3d T_CA = Eigen::Isometry3d::Identity();
                T_CA.translate(tran_CA);
                T_CA.rotate(quat_CA); // transformation matrix representing the AprilTag in the camera frame
                
                // Compute T_WA
                Eigen::Isometry3d T_WA = T_WB * T_BC * T_CA;
                tagPos = T_WA.translation();
                Eigen::Matrix3d tagRot = T_WA.rotation();
                Eigen::Quaterniond tagQuat = Eigen::Quaterniond(tagRot);
                Eigen::Quaterniond quatDiff = quat_WB.inverse() * tagQuat;
                yawDiff = quaternion2yaw(quatDiff); // yaw difference with correct sign
                 
                // Fill tagPose
                geometry_msgs::PoseStamped tagPoseTemp;
                tagPoseTemp.pose.position.x = tagPos(0);
                tagPoseTemp.pose.position.y = tagPos(1);
                tagPoseTemp.pose.position.z = tagPos(2);
                tagPoseTemp.pose.orientation.w = tagQuat.coeffs().w();
                tagPoseTemp.pose.orientation.x = tagQuat.coeffs().x();
                tagPoseTemp.pose.orientation.y = tagQuat.coeffs().y();
                tagPoseTemp.pose.orientation.z = tagQuat.coeffs().z();
                tagPoseTemp.header.frame_id = "map";
                tagPoseTemp.header.stamp = ros::Time::now();
                tagPose = tagPoseTemp;
                
                // TODO for Apriltag as state est
                // Use AprilTag as position estimation of quadcopter
//                 if (first_april_tag){
//                     first_april_tag = false;
//                     T_WA_initial = T_WA;
//                 }
//                 Eigen::Isometry3d T_WB_tag = T_WA_initial * T_CA.inverse() * T_BC.inverse(); // T_WB estimate, using apriltag as reference
//                 OdomData odom_data;
//                 Eigen::Vector3d pos_WB_tag = T_WB_tag.translation();
//                 odom_data.pose.position.x = pos_WB_tag(0);
//                 odom_data.pose.position.x = pos_WB_tag(1);
//                 odom_data.pose.position.x = pos_WB_tag(2);
//                 odom_data.source = APRIL_TAG;
//                 odomMsgSelector(odom_data);
//                 ROS_INFO_STREAM_ONCE("Using AprilTag as position estimate.");
            } else {
                // AprilTag is not detected, reset
                // first_april_tag = true; // TODO: for Apriltag as state est
            }
        }

    private:
        // Pose data
        double curYaw;
        double yawDiff;
        geometry_msgs::Pose curPose;
        geometry_msgs::PoseStamped tagPose;
        Eigen::Vector3d tagPos;
        nav_msgs::Path pathEst;       // estimate (vio of ekf2 depending on direct_from_vio param)
        
        // Landing info
        // Eigen::Isometry3d T_WA_initial;
        //bool first_april_tag; // TODO: for Apriltag as state est
        bool tag_visible;
        
        // Method to convert Eigen matrix data into a ROS nav_msgs/Path message
        void eigen2PathMsg(const Eigen::Matrix<double, Eigen::Dynamic, 11>& data, nav_msgs::Path& path) {
            int num_times = data.rows();
            path.header.stamp = ros::Time::now(); // Set the timestamp of the Path message
            path.header.frame_id = "map";
            
            // Clear the poses array in the Path message
            path.poses.clear();
            // Create a PoseStamped message
            geometry_msgs::PoseStamped pose_stamped;
            
            // Loop through each row of data
            for (int i = 0; i < num_times; i++) {
                pose_stamped.header.stamp = ros::Time(data(i, 0));
                pose_stamped.pose.position.x = data(i, 1);
                pose_stamped.pose.position.y = data(i, 2);
                pose_stamped.pose.position.z = data(i, 3);
                pose_stamped.pose.orientation.w = data(i, 7);
                pose_stamped.pose.orientation.x = data(i, 8);
                pose_stamped.pose.orientation.y = data(i, 9);
                pose_stamped.pose.orientation.z = data(i, 10);

                // Add the PoseStamped message to the poses array in the Path message
                path.poses.push_back(pose_stamped);
            }
        }
        
        // Method to append pose data to existing mav_msgs/Path message
        void addPoseToPath(const geometry_msgs::Pose& pose, nav_msgs::Path& path) {
            // Create pose stamped message
            geometry_msgs::PoseStamped pose_stamped;
            pose_stamped.header.stamp = ros::Time::now();
            pose_stamped.pose = pose;
            
            // Append the pose to the path
            path.header.frame_id = "map";
            path.poses.push_back(pose_stamped);
        }
        
        // Convert quaternion to yaw, for ZYX euler order
        double quaternion2yaw(const Eigen::Quaterniond& quat) {
            vector<double> q{quat.coeffs().w(), quat.coeffs().x(), quat.coeffs().y(), quat.coeffs().z()};
            double yaw = atan2(2*(q[0]*q[3] + q[1]*q[2]), pow(q[0],2) + pow(q[1],2) - pow(q[2],2) - pow(q[3],2));
            return yaw;
        }
        
        Eigen::Quaterniond msg2quaternion(const geometry_msgs::Quaternion& q) {
            return Eigen::Quaterniond{q.w, q.x, q.y, q.z};
        }
};

