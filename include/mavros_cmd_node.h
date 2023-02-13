#ifndef MAVROS_CMD_NODE_H
#define MAVROS_CMD_NODE_H

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <Eigen/Core>
#include <cmath>
#include <iostream>

#include "controller.h"

using namespace std;

class MavrosCmd {
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

        // Private Function
        template <class T>
        void waitForPose(const T *first_pose, const std::string &msg, double hz = 2.0) {
            ros::Rate wait_rate(hz);
            ROS_INFO_STREAM(msg);
            while (ros::ok() && !(*first_pose)) {
                ros::spinOnce();
                wait_rate.sleep();
            }
        };

    public:
        MavrosCmd(){
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
        ros::Timer setup_timer;
        
        // List of Public functions:
        /*
        void stateCallback(const mavros_msgs::State::ConstPtr& msg);    // Pass along mav state info
        void cmdloopCallback(const ros::TimerEvent& event);             // Publish flight commands according to current flight state
        void mavposeCallback(const geometry_msgs::PoseStamped &msg);    // Pass along mav pose data. Sets reference pose upon first call.
        void setupCallback(const ros::TimerEvent& event);               // Enable OFFBOARD mode and ARM, used for simulation
        */
        
        // ---------- Public Functions ----------- //

        // Pass along mav state info
        void stateCallback(const mavros_msgs::State::ConstPtr& msg){
            current_state = *msg;
        }

        // Publish flight commands according to current flight state
        void cmdloopCallback(const ros::TimerEvent& event) {
        switch (flight_state) {
            case WAITING_FOR_HOME_POSE: {
            waitForPose(&received_home_pose, "Waiting for home pose...");
            ROS_INFO("Got pose! Drone Ready to be armed.");
            flight_state = TAKEOFF;
            break;
            }
            case TAKEOFF: {
            // hover over home position
            geometry_msgs::PoseStamped takeoff_msg;
            takeoff_msg.header.stamp = ros::Time::now();
            takeoff_msg.pose = home_pose;
            takeoff_msg.pose.position.z += 0.5;
            double x_diff = takeoff_msg.pose.position.x - curr_pose.position.x;
            double y_diff = takeoff_msg.pose.position.y - curr_pose.position.y;
            double z_diff = takeoff_msg.pose.position.z - curr_pose.position.z;
            double pos_error = sqrt(pow(x_diff,2.) + pow(y_diff,2.) + pow(z_diff,2.)); // in meters
            if (pos_error > 0.1){
                pos_pub.publish(takeoff_msg);
            } else {
                flight_state = MISSION;
                ros::spinOnce();
            }
            break;
            }
            case MISSION: {
            double x_diff = ref_pose.position.x - curr_pose.position.x;
            double y_diff = ref_pose.position.y - curr_pose.position.y;
            double z_diff = ref_pose.position.z - curr_pose.position.z;
            double pos_error = sqrt(pow(x_diff,2.) + pow(y_diff,2.) + pow(z_diff,2.)); // in meters
            if (ros::Time::now() - print_request > ros::Duration(2.0)){
                cout << "x diff: " << x_diff << endl;
                cout << "y diff: " << y_diff << endl;
                cout << "z diff: " << z_diff << endl;
                cout << "total error: " << pos_error << endl << endl;
                print_request = ros::Time::now();
            }
            if (pos_error > 0.08) {
                ref_msg.header.stamp = ros::Time::now();
                ref_msg.pose = ref_pose;
                pos_pub.publish(ref_msg);
            } else {
                flight_state = RETURNING;
                pos_pub.publish(ref_msg);
            }
            break;
            }
            case RETURNING: {
            // hover over home position
            geometry_msgs::PoseStamped return_msg;
            return_msg.header.stamp = ros::Time::now();
            return_msg.pose = home_pose;
            return_msg.pose.position.z = home_pose.position.z + 0.5;
            double x_diff = return_msg.pose.position.x - curr_pose.position.x;
            double y_diff = return_msg.pose.position.y - curr_pose.position.y;
            double z_diff = return_msg.pose.position.z - curr_pose.position.z;
            double pos_error = sqrt(pow(x_diff,2.) + pow(y_diff,2.) + pow(z_diff,2.)); // in meters
            if (pos_error > 0.1){
                pos_pub.publish(return_msg);
            } else {
                flight_state = LANDING;
                ros::spinOnce();
            }
            break;
            }
            case LANDING: {
            // auto land
            offb_set_mode.request.custom_mode = "AUTO.LAND";
            if (current_state.mode != "AUTO.LAND" && (ros::Time::now() - last_request > ros::Duration(2.0))) {
                if (set_mode_client.call(offb_set_mode) && offb_set_mode.response.mode_sent) {
                ROS_INFO("AUTO.LAND enabled");
                ROS_INFO("Once landed, switch to position mode and disarm.");
                cmdloop_timer.stop();
                }
                last_request = ros::Time::now();
            }
            break;
            }
            default: {
            ROS_INFO("Error: Flight State Not Defined");
            }
        }
        }

        // Pass along mav pose data. Sets reference pose upon first call.
        void mavposeCallback(const geometry_msgs::PoseStamped &msg) {
        if (!received_home_pose) {
            received_home_pose = true;
            home_pose = msg.pose;
            ROS_INFO_STREAM("Home pose initialized to: " << home_pose);
            
            geometry_msgs::Pose temp = home_pose;
            temp.position.x += ref_offset.position.x;
            temp.position.y += ref_offset.position.y;
            temp.position.z += ref_offset.position.z;
            ref_pose = temp;
            ROS_INFO_STREAM("Reference offset initialized to: " << ref_offset);
        }
        curr_pose = msg.pose;
        }

        // Enable OFFBOARD mode and ARM, used for simulation
        void setupCallback(const ros::TimerEvent &event) {
            if (sim_enable) {
                switch(flight_state) {
                    case LANDING: {
                        // try to disarm
                        arm_cmd.request.value = false;
                        if (current_state.armed && (ros::Time::now() - last_request > ros::Duration(2.0))) {
                            if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                                ROS_INFO("Vehicle disarmed");
                                setup_timer.stop();
                            }
                            last_request = ros::Time::now();
                        }
                    }
                    default: {
                        arm_cmd.request.value = true;
                        offb_set_mode.request.custom_mode = "OFFBOARD";
                        if (current_state.mode != "OFFBOARD" && (ros::Time::now() - last_request > ros::Duration(2.0))) {
                            if (set_mode_client.call(offb_set_mode) && offb_set_mode.response.mode_sent) {
                                ROS_INFO("Offboard enabled");
                            }
                            last_request = ros::Time::now();
                        } else {
                            if (!current_state.armed && (ros::Time::now() - last_request > ros::Duration(2.0))) {
                                if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                                    ROS_INFO("Vehicle armed");
                                }
                                last_request = ros::Time::now();
                            }
                        }
                    }
                }
            } else { // not using sim
                setup_timer.stop();
            }
        }

};



#endif
