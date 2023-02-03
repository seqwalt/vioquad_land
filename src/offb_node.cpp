/**
 * @file offb_node.cpp
 * @brief Offboard control example node, written with MAVROS version 0.19.x, PX4 Pro Flight
 * Stack and tested in Gazebo SITL
 */

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>

mavros_msgs::State current_state;
void state_cb(const mavros_msgs::State::ConstPtr& msg){
    current_state = *msg;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "offb_node");
    ros::NodeHandle nh;

    ros::Subscriber state_sub = nh.subscribe<mavros_msgs::State>
            ("mavros/state", 10, state_cb);
    // ros::Subscriber pose_sub = nh.subscribe<geometry_msgs::PoseStamped>
    //         ("mavros/state", 10, state_cb);
    ros::Publisher local_pos_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mavros/setpoint_position/local", 10);
    ros::ServiceClient arming_client = nh.serviceClient<mavros_msgs::CommandBool>
            ("mavros/cmd/arming");
    ros::ServiceClient set_mode_client = nh.serviceClient<mavros_msgs::SetMode>
            ("mavros/set_mode");

    //the setpoint publishing rate MUST be faster than 2Hz
    ros::Rate rate(20.0);

    // wait for FCU connection
    while(ros::ok() && !current_state.connected){
        ros::spinOnce();
        rate.sleep();
    }

    geometry_msgs::PoseStamped pose;
    pose.pose.position.x = 0;
    pose.pose.position.y = 0;
    pose.pose.position.z = 2;

    //send a few setpoints before starting
    for(int i = 100; ros::ok() && i > 0; --i){
        local_pos_pub.publish(pose);
        ros::spinOnce();
        rate.sleep();
    }

    // service request variable for offboard mode and landing mode
    mavros_msgs::SetMode offb_set_mode, land_set_mode;
    offb_set_mode.request.custom_mode = "OFFBOARD";
    land_set_mode.request.custom_mode = "AUTO.LAND";

    // service request variable for arming
    mavros_msgs::CommandBool arm_cmd;
    arm_cmd.request.value = true;

    ros::Time last_request = ros::Time::now();

    bool landing = false;

    // Arm
    while(ros::ok() && (!current_state.armed || current_state.mode != "OFFBOARD")){
      if( current_state.mode != "OFFBOARD" &&
            (ros::Time::now() - last_request > ros::Duration(5.0))){
            if( set_mode_client.call(offb_set_mode) &&
                offb_set_mode.response.mode_sent){
                ROS_INFO("Offboard enabled");
            }
            last_request = ros::Time::now();
      } else {
            if( !current_state.armed &&
                (ros::Time::now() - last_request > ros::Duration(5.0))){
                std::cout << "try arm" << std::endl;
                if( arming_client.call(arm_cmd) &&
                    arm_cmd.response.success){
                    ROS_INFO("Vehicle armed");
                }
                last_request = ros::Time::now();
            }
      }

      local_pos_pub.publish(pose);

      ros::spinOnce();
      rate.sleep();
    }

    ros::Time start_time = ros::Time::now();

    while(ros::ok()){

        // Send set point
        if( landing == false ){
          local_pos_pub.publish(pose);
        }

        // land after 15 seconds
        if( current_state.mode == "OFFBOARD" &&
            (ros::Time::now() - start_time > ros::Duration(15.0)) &&
            (ros::Time::now() - last_request > ros::Duration(5.0)) ){
            if( set_mode_client.call(land_set_mode) &&
                land_set_mode.response.mode_sent){
                ROS_INFO("Landing enabled");
                landing == true;
            }
            last_request = ros::Time::now();
        }

        ros::spinOnce();
        rate.sleep();
    }

    return 0;
}
