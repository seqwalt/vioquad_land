#pragma once

#ifndef MAVROS_CMD_NODE_H
#define MAVROS_CMD_NODE_H

#include <ros/ros.h>
#include <mavros/mavros.h>
#include <geometry_msgs/PoseStamped.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/Thrust.h>
#include <std_srvs/Trigger.h>

#include <iostream>

#include "quad_control/InitTraj.h"     // custom service 
#include "quad_control/FlatOutputs.h"  // custom message
#include "controller.h"

using namespace std;

class MavrosCmd {
    public:
        MavrosCmd(){
            flight_state = INITIALIZATION;
            received_home_pose = false;
            trajDone = false;
            
            mode_request = ros::Time::now();
            land_request = ros::Time::now();
            traj_request = ros::Time::now();
            print_request = ros::Time::now();
        }
        
        bool sim_enable;
        mavros_msgs::State currentModes;

        ros::Subscriber mode_sub;
        ros::Subscriber pose_sub;
        ros::Subscriber vel_sub;
        ros::Subscriber traj_sub;
        
        ros::Publisher pos_pub;
        ros::Publisher att_pub;
        ros::Publisher angVel_pub;
        ros::Publisher thrust_pub;
        
        ros::ServiceClient arming_client;
        ros::ServiceClient set_mode_client;
        ros::ServiceClient init_traj_client;
        ros::ServiceClient send_traj_client;
        
        ros::Timer cmdloop_timer;
        ros::Timer setup_timer;
        
        // List of Public functions:
        /*
        void modeCallback(const mavros_msgs::State::ConstPtr& msg);     // Pass along mav mode info
        void cmdloopCallback(const ros::TimerEvent& event);             // Publish flight commands according to current flight state
        void mavposeCallback(const geometry_msgs::PoseStamped &msg);    // Pass along mav pose data. Sets reference pose upon first call.
        void setupCallback(const ros::TimerEvent& event);               // Enable OFFBOARD mode and ARM, used for simulation
        */
        
        // ---------- Public Functions ----------- //

        // Pass along mav state info
        void modeCallback(const mavros_msgs::State::ConstPtr& msg){
            currentModes = *msg;
        }

        // Publish flight commands according to current flight state
        void cmdLoopCallback(const ros::TimerEvent& event) {
            switch (flight_state) {
                case INITIALIZATION: {
                    if (!received_home_pose) {
                        waitForData(&received_home_pose, "Waiting for home pose...");
                        ROS_INFO("Got home pose.");
                    } else if(ros::Time::now() - traj_request > ros::Duration(2.0)) {
                        traj_request = ros::Time::now();
                        ROS_INFO("Waiting for initial trajectory pose");
                        if (init_traj_client.call(init_traj) && init_traj.response.success) {
                            // ROS_INFO("Received initial trajectory pose. Commencing takeoff...");
                            ROS_INFO("Received.");
                            ROS_INFO_STREAM("Initial traj position: " << init_traj.response.position);
                            ROS_INFO_STREAM("Initial traj yaw     : " << init_traj.response.yaw);
                            flight_state = TAKEOFF;
                        }
                    }
                    break;
                }
                case TAKEOFF: {
                    // hover at trajectory start pose
                    geometry_msgs::Point takeoff_pos = init_traj.response.position;
                    double takeoff_yaw = init_traj.response.yaw;
                    Eigen::AngleAxisd rotation_vector(takeoff_yaw, Eigen::Vector3d(0,0,1));  // rotate about z-axis
                    Eigen::Quaterniond takeoff_quat(rotation_vector); // compute yaw rotation as a quaternion
                    
                    geometry_msgs::PoseStamped takeoff_msg;
                    takeoff_msg.header.stamp = ros::Time::now();
                    takeoff_msg.pose.orientation = tf2::toMsg(takeoff_quat);
                    takeoff_msg.pose.position = takeoff_pos;
                    takeoff_msg.pose.position.z = takeoff_pos.z + home_pose.position.z;
                    
                    // Error
                    double x_diff = takeoff_msg.pose.position.x - curPose.position.x;
                    double y_diff = takeoff_msg.pose.position.y - curPose.position.y;
                    double z_diff = takeoff_msg.pose.position.z - curPose.position.z;
                    double pos_error = sqrt(pow(x_diff,2.) + pow(y_diff,2.) + pow(z_diff,2.)); // in meters
                    
                    Eigen::Quaterniond curQuat = ctrl.quatToEigen(curPose.orientation);
                    double curYaw = mavros::ftf::quaternion_get_yaw(curQuat);
                    double yaw_error = abs(takeoff_yaw - curYaw);
                    cout << abs(z_diff) << endl;
                    /*
                    if (abs(z_diff) > 0.1){
                        dispOnce("ZDIFFFFFFF", sw1);
                        geometry_msgs::PoseStamped vert_takeoff_msg;
                        vert_takeoff_msg.header.stamp = ros::Time::now();
                        vert_takeoff_msg.pose = home_pose;
                        vert_takeoff_msg.pose.position.z += 1;
                        pos_pub.publish(vert_takeoff_msg);
                    } else
                    */
                    if (pos_error > 0.25 || yaw_error > 0.2){
                        dispOnce("pos_eROROR", sw2);
                        pos_pub.publish(takeoff_msg);
                    } else {
                        pos_pub.publish(takeoff_msg);
                        flight_state = MISSION;
                        ROS_INFO("Starting mission.");
                        ros::spinOnce();
                    }
                    break;
                }
                case MISSION: {
                    if (!trajDone) {
                        // request the trajectory to be sent
                        if(!send_traj.response.success && (ros::Time::now() - traj_request > ros::Duration(2.0))){
                            send_traj_client.call(send_traj);
                            traj_request = ros::Time::now();
                        }
                        // TODO: cout the absolute error btw curr pose and desired pose
                    } else {
                        flight_state = LANDING;
                        pos_pub.publish(curPose);
                    }
                    break;
                }
                case LANDING: {
                    // auto land
                    mavros_set_mode.request.custom_mode = "AUTO.LAND";
                    if (currentModes.mode != "AUTO.LAND" && (ros::Time::now() - land_request > ros::Duration(2.0))) {
                        if (set_mode_client.call(mavros_set_mode) && mavros_set_mode.response.mode_sent) {
                        ROS_INFO("AUTO.LAND enabled");
                        ROS_INFO("Once landed, switch to position mode and disarm.");
                        cmdloop_timer.stop();
                        }
                        land_request = ros::Time::now();
                    }
                    break;
                }
                default: {
                    ROS_INFO("Error: Flight State Not Defined");
                }
            }
        }

        // Pass along current mav pose data. Sets reference pose upon first call.
        void mavPoseCallback(const geometry_msgs::PoseStamped &msg) {
            if (!received_home_pose) {
                received_home_pose = true;
                home_pose = msg.pose;
                ROS_INFO_STREAM("Home pose initialized to: " << home_pose);
            }
            curPose = msg.pose;
        }
        
        // Pass along current mav velocity data.
        void mavVelCallback(const geometry_msgs::TwistStamped &msg) {
            curVel = msg.twist.linear;
        }
        
        // Control Inputs are sent, after recieving reference setpoint
        void mavRefCallback(const quad_control::FlatOutputs &msg) {
            // Reference info
            flatRef.position = msg.position;
            flatRef.velocity = msg.velocity;
            flatRef.acceleration = msg.acceleration;
            flatRef.yaw = msg.yaw;
            trajDone = msg.trajectoryDone;
            
            // Current State info
            curState.pose = curPose;
            curState.velocity = curVel;
            
            // Compute and send control commands
            ctrl.Geometric(ctrlInputs, flatRef, curState);
            att_pub.publish(ctrlInputs.attitude);
            angVel_pub.publish(ctrlInputs.cmd_vel);
            thrust_pub.publish(ctrlInputs.norm_thrust);
        }

        // Enable OFFBOARD mode and ARM, used for simulation
        void setupCallback(const ros::TimerEvent &event) {
            if (sim_enable) {
                switch(flight_state) {
                    case LANDING: {
                        // try to disarm
                        arm_cmd.request.value = false;
                        if (currentModes.armed && (ros::Time::now() - mode_request > ros::Duration(2.0))) {
                            if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                                ROS_INFO("Vehicle disarmed");
                                setup_timer.stop();
                            }
                            mode_request = ros::Time::now();
                        }
                    }
                    default: {
                        arm_cmd.request.value = true;
                        mavros_set_mode.request.custom_mode = "OFFBOARD";
                        if (currentModes.mode != "OFFBOARD" && (ros::Time::now() - mode_request > ros::Duration(2.0))) {
                            if (set_mode_client.call(mavros_set_mode) && mavros_set_mode.response.mode_sent) {
                                ROS_INFO("Offboard enabled");
                            }
                            mode_request = ros::Time::now();
                        } else {
                            if (!currentModes.armed && (ros::Time::now() - mode_request > ros::Duration(2.0))) {
                                if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                                    ROS_INFO("Vehicle armed");
                                }
                                mode_request = ros::Time::now();
                            }
                        }
                    }
                }
            } else { // not using sim
                setup_timer.stop();
            }
        }

    private:
        enum FlightState { INITIALIZATION, TAKEOFF, MISSION, LANDING } flight_state;
        
        mavros_msgs::CommandBool arm_cmd;
        mavros_msgs::SetMode mavros_set_mode;
        quad_control::InitTraj init_traj;
        std_srvs::Trigger send_traj;
        
        geometry_msgs::Pose home_pose;
        geometry_msgs::Pose curPose;
        geometry_msgs::Vector3 curVel;
        
        Controller ctrl;
        Controller::FlatReference flatRef;
        Controller::ControlInputs ctrlInputs;
        Controller::State curState;
        bool trajDone;
        
        ros::Time mode_request;
        ros::Time land_request;
        ros::Time traj_request;
        ros::Time print_request;
        bool received_home_pose;

        // Private Function
        template <class T>
        void waitForData(const T *data, const std::string &msg, double hz = 100.0) {
            ros::Rate wait_rate(hz);
            ROS_INFO_STREAM(msg);
            while (ros::ok() && !(*data)) {
                ros::spinOnce();
                wait_rate.sleep();
            }
        };
        
        bool sw1, sw2, sw3 = true;
        void dispOnce(const string &str, bool &sw){
            if (sw){
                cout << str << endl;
                sw = false;
            }
        }
};



#endif
