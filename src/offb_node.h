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

class OffbCtrl {
    public:
        OffbCtrl(){
            flight_state = WAITING_FOR_HOME_POSE;
        }
        
        bool sim_enable;
        geometry_msgs::Pose ref_offset;
        mavros_msgs::State current_state;

        ros::Subscriber state_sub;
        ros::Subscriber pose_sub;
        ros::Publisher pos_pub;
        ros::ServiceClient arming_client;
        ros::ServiceClient set_mode_client;
        ros::Timer cmdloop_timer;
        ros::Timer sim_timer;
        
        void stateCallback(const mavros_msgs::State::ConstPtr& msg);    // Pass along mav state info
        void cmdloopCallback(const ros::TimerEvent& event);             // Publish flight commands according to current flight state
        void mavposeCallback(const geometry_msgs::PoseStamped &msg);    // Pass along mav pose data. Sets reference pose upon first call.
        void simCallback(const ros::TimerEvent& event);                 // Enable OFFBOARD mode and ARM, used for simulation

    private:
        enum FlightState { WAITING_FOR_HOME_POSE, TAKEOFF, MISSION, RETURNING, LANDING } flight_state;
        
        mavros_msgs::CommandBool arm_cmd;
        mavros_msgs::SetMode offb_set_mode;

        geometry_msgs::Pose home_pose;
        geometry_msgs::PoseStamped ref_msg;
        geometry_msgs::Pose ref_pose;
        geometry_msgs::Pose curr_pose;

        ros::Time last_request;
        ros::Time land_time;
        ros::Time print_request;
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
};



#endif
