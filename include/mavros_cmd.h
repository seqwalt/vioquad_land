#ifndef MAVROS_CMD_H_INCLUDED
#define MAVROS_CMD_H_INCLUDED

#include <ros/ros.h>
#include <mavros/mavros.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/Thrust.h>
#include <std_srvs/Trigger.h>

#include <iostream>

#include "vioquad_land/InitSetpoint.h"  // custom service
#include "vioquad_land/FlatOutputs.h"   // custom message
#include "tracking_controller.h"

using namespace std;

class MavrosCmd {
    public:
        MavrosCmd(){
            flight_state = INITIALIZATION;
            received_home_pose = false;
            mission_done = false;

            mode_request = ros::Time::now();
            land_request = ros::Time::now();
            setpnt_request = ros::Time::now();
            print_request = ros::Time::now();
        }

        Controller ctrl;
        
        enum ControllerType { POSITION, GEOMETRIC, MPC } ctrl_mode; // chosen in launch file
        bool sim_enable;
        mavros_msgs::State currentModes;

        ros::Subscriber mode_sub;
        ros::Subscriber pose_sub;
        ros::Subscriber vel_sub;
        ros::Subscriber traj_sub;

        ros::Publisher pos_pub;
        ros::Publisher pos_track_pub;
        ros::Publisher att_pub;

        ros::ServiceClient arming_client;
        ros::ServiceClient set_mode_client;
        ros::ServiceClient init_setpnt_client;
        ros::ServiceClient streamimg_client;

        ros::Timer cmdloop_timer;
        ros::Timer modes_timer;
        
        ros::Time print_request;

        // Pass along mav state info
        void modeCallback(const mavros_msgs::State::ConstPtr& msg) {
            currentModes = *msg;
        }

        // Switch between flight states (e.g. takeoff, landing) and
        // perform tasks for each flight state (e.g. request trajectory, tell quadcopter to takeoff)
        void cmdLoopCallback(const ros::TimerEvent& event) {
            switch (flight_state) {
                case INITIALIZATION: {
                    if (!received_home_pose) {
                        waitForData(&received_home_pose, "Waiting for home pose.");
                        ROS_INFO("Got home pose.");
                    } else if(ros::Time::now() - setpnt_request > ros::Duration(2.0)) {
                        setpnt_request = ros::Time::now();
                        //ROS_INFO("Waiting for initial trajectory pose.");
                        if (init_setpnt_client.call(init_setpnt) && init_setpnt.response.success) {
                            ROS_INFO_STREAM("Initial setpoint: " << endl << init_setpnt.response);
                            flight_state = TAKEOFF;
                        }
                    }
                    break;
                }
                case TAKEOFF: {
                    // Vertical takeoff
                    
                    // Only start offboard mode after publishing to pos_pub for a second
                    static ros::Time takeoff_time = ros::Time::now();
                    if (ros::Time::now() - takeoff_time > ros::Duration(1.0)) modes_timer.start();
                    
                    geometry_msgs::PoseStamped takeoff_msg;
                    takeoff_msg.header.stamp = ros::Time::now();
                    takeoff_msg.pose = home_pose;
                    double z_ref = 1.0;
                    takeoff_msg.pose.position.z += z_ref;

                    // Error
                    double z_diff = takeoff_msg.pose.position.z - curPose.position.z;
                    //if (ros::Time::now() - print_request > ros::Duration(3.0)){
                    //    ROS_INFO_STREAM("z_diff = " << z_diff);
                    //    print_request = ros::Time::now();
                    //}
                    if (abs(z_diff) > 0.8){
                        pos_pub.publish(takeoff_msg);
                    } else {
                        pos_pub.publish(takeoff_msg);
                        flight_state = TO_START_POSE;
                        ros::spinOnce();
                    }
                    break;
                }
                case TO_START_POSE: {
                    // Go to start of trajectory then hover
                    geometry_msgs::Point firstPosition = init_setpnt.response.position;
                    double takeoff_yaw = init_setpnt.response.yaw;
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
                    /*if (ros::Time::now() - print_request > ros::Duration(1.0)){
                        ROS_INFO_STREAM("pos_err = " << pos_error);
                        ROS_INFO_STREAM("yaw_err = " << yaw_error);
                        print_request = ros::Time::now();
                    }*/
                    if (pos_error > 0.08 || yaw_error > 0.05){
                        // Go to first trajectory pose
                        ROS_INFO_ONCE("Going to start of trajectory.");
                        pos_pub.publish(firstPose_msg);
                    } else {
                        pos_pub.publish(firstPose_msg);
                        flight_state = MISSION;
                        ros::spinOnce();
                    }
                    break;
                }
                case MISSION: {
                    // Request trajectory data. Each time trajectory data is sent, 
                    // mavRefCallback() will send a control input
                    if (!mission_done) {
                        // request the trajectory to be sent
                        if(!start_stream.response.success){
                            pos_pub.publish(firstPose_msg); // fixed hover while waiting for trajectory data
                            if(ros::Time::now() - setpnt_request > ros::Duration(2.0)){
                                //ROS_INFO("Requesting trajectory.");
                                streamimg_client.call(start_stream);
                                setpnt_request = ros::Time::now();
                            }
                        } else {
                            ROS_INFO_ONCE("Executing mission.");
                        }
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
                    if (currentModes.mode != "AUTO.LAND" && (ros::Time::now() - land_request > ros::Duration(0.5))) {
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
                    ROS_ERROR("Flight state not defined, landing");
                    flight_state = LANDING;
                    break;
                }
            } // end switch statement
        } // end of cmdLoopCallback()

        // Pass along current mav pose data. Sets reference pose upon first call.
        // Required for feedback control of quadcopter
        void mavPoseCallback(const geometry_msgs::PoseStamped &msg) {
            if (!received_home_pose) {
                received_home_pose = true;
                home_pose = msg.pose;
                ROS_INFO_STREAM("Home pose initialized to: " << endl << home_pose);
            }
            curPose = msg.pose;
        }

        // Pass along current mav velocity data.
        // Required for feedback control of quadcopter
        void mavVelCallback(const geometry_msgs::TwistStamped &msg) {
            curVel = msg.twist.linear;
        }

        // Control Inputs are sent, after recieving reference setpoint
        // ctrl_mode is a parameter set during node start-up
        void mavRefCallback(const vioquad_land::FlatOutputs &flatRefMsg) {
            // Reference info
            mission_done = flatRefMsg.reached_goal;

            // Current State info
            curState.pose = curPose;
            curState.velocity = curVel;

            // Compute and send control commands
            switch(ctrl_mode) {
                case GEOMETRIC: { // Geometric controller
                    ctrl.Geometric(attInputs, flatRefMsg, curState);
                    attInputs.header.stamp = ros::Time::now();
                    att_pub.publish(attInputs); // set attitude, body rate and thrust to mavros
                    break;
                }
                case POSITION: { // Position/Yaw controller
                    ctrl.PosYaw(posYawInputs, flatRefMsg);
                    posYawInputs.header.stamp = ros::Time::now();
                    posYawInputs.coordinate_frame = 1; // corresponds to MAV_FRAME_LOCAL_NED
                    pos_track_pub.publish(posYawInputs);
                    break;
                }
                default: {
                    ROS_ERROR("Controller type not defined, landing");
                    flight_state = LANDING;
                    break;
                }
            }
        } // end of mavRefCallback()

        // Enable OFFBOARD mode. Arm/Disarm if in simulation
        void modesCallback(const ros::TimerEvent &event) {
            switch(flight_state) {
                case LANDING: {
                    // try to disarm if in simulation
                    if (sim_enable && currentModes.armed && (ros::Time::now() - mode_request > ros::Duration(2.0))) {
                        arm_cmd.request.value = false;
                        if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                            ROS_INFO("Vehicle disarmed");
                            modes_timer.stop();
                        }
                        mode_request = ros::Time::now();
                    }
                    break;
                }
                default: {
                    // Enter offboard mode only if trying to take off 
                    if (flight_state == TAKEOFF && currentModes.mode != "OFFBOARD" && (ros::Time::now() - mode_request > ros::Duration(0.5))) {
                        mavros_set_mode.request.custom_mode = "OFFBOARD";
                        if (set_mode_client.call(mavros_set_mode) && mavros_set_mode.response.mode_sent) {
                            ROS_INFO("Offboard enabled");
                        }
                        mode_request = ros::Time::now();
                    }
                    // try to arm if in simulation
                    if (try_to_arm && sim_enable && !currentModes.armed && (ros::Time::now() - mode_request > ros::Duration(0.5))) {
                        arm_cmd.request.value = true;
                        if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                            ROS_INFO("Vehicle armed for takeoff");
                            try_to_arm = false;
                        }
                        mode_request = ros::Time::now();
                    }
                    break;
                }
            }
        } // end of setupCallback()

    private:
        enum FlightState { INITIALIZATION, TAKEOFF, TO_START_POSE, MISSION, LANDING } flight_state;

        mavros_msgs::CommandBool arm_cmd;
        mavros_msgs::SetMode mavros_set_mode;
        vioquad_land::InitSetpoint init_setpnt;
        std_srvs::Trigger start_stream;

        geometry_msgs::Pose home_pose;
        geometry_msgs::PoseStamped firstPose_msg;
        geometry_msgs::Pose curPose;
        geometry_msgs::Vector3 curVel;

        Controller::State curState;
        mavros_msgs::AttitudeTarget attInputs;
        mavros_msgs::PositionTarget posYawInputs;
        bool mission_done;
        bool try_to_arm = true;

        ros::Time mode_request;
        ros::Time land_request;
        ros::Time setpnt_request;
        bool received_home_pose;

        template <class T>
        void waitForData(const T *data, const std::string &msg, double hz = 100.0) {
            ros::Rate wait_rate(hz);
            ROS_INFO_STREAM(msg);
            while (ros::ok() && !(*data)) {
                ros::spinOnce();
                wait_rate.sleep();
            }
        };
};

#endif // MAVROS_CMD_H_INCLUDED
