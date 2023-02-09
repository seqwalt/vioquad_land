#ifndef OFFB_NODE_H
#define OFFB_NODE_H

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <Eigen/Core>
#include <cmath>
#include <iostream>

using namespace std;

enum FlightState { WAITING_FOR_HOME_POSE, MISSION_EXECUTION, LANDING, LANDED } flight_state;

mavros_msgs::State current_state;
mavros_msgs::CommandBool arm_cmd;
mavros_msgs::SetMode offb_set_mode;
geometry_msgs::Pose home_pose;
geometry_msgs::PoseStamped hover_msg;
geometry_msgs::Pose hover_pose;
geometry_msgs::Pose hover_offset;

ros::Time last_request;

ros::Publisher pos_pub;

ros::Timer cmdloop_timer;
ros::Timer simInit_timer;

ros::ServiceClient arming_client;
ros::ServiceClient set_mode_client;

bool sim_enable;
bool received_home_pose;

template <class T>
void waitForPose(const T *first_pose, const std::string &msg, double hz = 2.0) {
    ros::Rate wait_rate(hz);
    ROS_INFO_STREAM(msg);
    while (ros::ok() && !(*first_pose)) {
        ros::spinOnce();
        wait_rate.sleep();
    }
};

#endif
