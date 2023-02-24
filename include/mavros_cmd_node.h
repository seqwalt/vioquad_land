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

        Controller ctrl;
        
        bool sim_enable;
        string ctrl_mode;
        mavros_msgs::State currentModes;

        ros::Subscriber mode_sub;
        ros::Subscriber pose_sub;
        ros::Subscriber vel_sub;
        ros::Subscriber traj_sub;

        ros::Publisher pos_pub;
        ros::Publisher att_pub;

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
                    // Vertical takeoff
                    geometry_msgs::PoseStamped takeoff_msg;
                    takeoff_msg.header.stamp = ros::Time::now();
                    takeoff_msg.pose = home_pose;
                    takeoff_msg.pose.position.z += 1.0;

                    // Error
                    double z_diff = takeoff_msg.pose.position.z - curPose.position.z;
                    if (ros::Time::now() - print_request > ros::Duration(3.0)){
                        ROS_INFO_STREAM("z_diff = " << z_diff);
                        print_request = ros::Time::now();
                    }
                    if (abs(z_diff) > 0.1){
                        pos_pub.publish(takeoff_msg);
                        dispOnce("Trying vertical takeoff. ", switch1);
                    } else {
                        pos_pub.publish(takeoff_msg);
                        flight_state = TO_START_POSE;
                        dispOnce("Done vertical takeoff. ", switch4);
                        ros::spinOnce();
                    }
                    break;
                }
                case TO_START_POSE: {
                    //dispOnce("Start pose checkpoint 1", switch6);
                    // hover at trajectory start pose
                    geometry_msgs::Point firstPosition = init_traj.response.position;
                    double takeoff_yaw = init_traj.response.yaw;
                    Eigen::AngleAxisd rotation_vector(takeoff_yaw, Eigen::Vector3d(0,0,1));  // rotate about z-axis
                    Eigen::Quaterniond firstQuat(rotation_vector); // compute yaw rotation as a quaternion

                    firstPose_msg.header.stamp = ros::Time::now();
                    firstPose_msg.pose.orientation = tf2::toMsg(firstQuat);
                    firstPose_msg.pose.position = firstPosition;

                    // Error
                    double x_diff = firstPose_msg.pose.position.x - curPose.position.x;
                    double y_diff = firstPose_msg.pose.position.y - curPose.position.y;
                    double z_diff = firstPose_msg.pose.position.z - curPose.position.z;
                    double pos_error = sqrt(pow(x_diff,2.) + pow(y_diff,2.) + pow(z_diff,2.)); // in meters

                    Eigen::Quaterniond curQuat = ctrl.quatToEigen(curPose.orientation);
                    double curYaw = mavros::ftf::quaternion_get_yaw(curQuat);
                    double yaw_error = abs(takeoff_yaw - curYaw);
                    if (ros::Time::now() - print_request > ros::Duration(1.0)){
                        ROS_INFO_STREAM("pos_err = " << pos_error);
                        ROS_INFO_STREAM("yaw_err = " << yaw_error);
                        print_request = ros::Time::now();
                    }
                    if (pos_error > 0.1 || yaw_error > 0.05){
                        // Go to first trajectory pose
                        //dispOnce("Going to first trajectory pose.", switch2);
                        pos_pub.publish(firstPose_msg);
                    } else {
                        pos_pub.publish(firstPose_msg);
                        flight_state = MISSION;
                        ros::spinOnce();
                    }
                    break;
                }
                case MISSION: {
                    if (!trajDone) {
                        // request the trajectory to be sent
                        if(!send_traj.response.success){
                            pos_pub.publish(firstPose_msg); // fixed hover while waiting for trajectory data
                            if(ros::Time::now() - traj_request > ros::Duration(2.0)){
                                ROS_INFO("Requesting Trajectory.");
                                send_traj_client.call(send_traj);
                                traj_request = ros::Time::now();
                            }
                        } else {
                            pos_pub.publish(firstPose_msg);
                            ros::spinOnce();
                        }
                        // TODO: cout the absolute error btw curr pose and desired pose
                    } else {
                        geometry_msgs::PoseStamped hover_msg;
                        hover_msg.header.stamp = ros::Time::now();
                        hover_msg.pose = curPose;
                        pos_pub.publish(hover_msg);
                        flight_state = LANDING;
                        ros::spinOnce();
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
            } // end switch statement
        } // end cmdLoopCallback function

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
        void mavRefCallback(const quad_control::FlatOutputs &flatRefMsg) {
            // Reference info
            trajDone = flatRefMsg.trajectoryDone;

            // Current State info
            curState.pose = curPose;
            curState.velocity = curVel;

            // Compute and send control commands

            // Geometric controller TODO: Choose mode based on ctrl_mode string
            //ctrl.Geometric(attInputs, flatRefMsg, curState);
            //attInputs.header.stamp = ros::Time::now();
            //attInputs.type_mask = 128;  // ignore attitude
            //att_pub.publish(attInputs); // set attitude, body rate and thrust to mavros

            // Position/Yaw controller
            ctrl.PosYaw(posYawInputs, flatRefMsg);
            posYawInputs.header.stamp = ros::Time::now();
            pos_pub.publish(posYawInputs);
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
        enum FlightState { INITIALIZATION, TAKEOFF, TO_START_POSE, MISSION, LANDING } flight_state;

        mavros_msgs::CommandBool arm_cmd;
        mavros_msgs::SetMode mavros_set_mode;
        quad_control::InitTraj init_traj;
        std_srvs::Trigger send_traj;

        geometry_msgs::Pose home_pose;
        geometry_msgs::PoseStamped firstPose_msg;
        geometry_msgs::Pose curPose;
        geometry_msgs::Vector3 curVel;

        Controller::State curState;
        mavros_msgs::AttitudeTarget attInputs;
        mavros_msgs::PositionTarget posYawInputs;
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

        bool switch1, switch2, switch3, switch4, switch5, switch6, switch7, switch8 = true;
        void dispOnce(const string &str, bool &sw){
            if (sw){
                ROS_INFO_STREAM(str);
                sw = false;
            }
        }
};



#endif
