/**
 * @file offb_node.cpp
 * @brief Offboard control example node, written with MAVROS version 0.19.x, PX4 Pro Flight
 * Stack and tested in Gazebo SITL
 */

#include "offb_node.h"

void stateCallback(const mavros_msgs::State::ConstPtr& msg);
void cmdloopCallback(const ros::TimerEvent& event);
void mavposeCallback(const geometry_msgs::PoseStamped &msg);
void simInitCallback(const ros::TimerEvent& event);

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "offb_node");
    ros::NodeHandle nh;

    // Subscribers
    ros::Subscriber state_sub = nh.subscribe
            ("mavros/state", 1, stateCallback);
    ros::Subscriber pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, mavposeCallback, ros::TransportHints().tcpNoDelay());
    
    // Publishers
    pos_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mavros/setpoint_position/local", 10);
    
    // Services        
    arming_client = nh.serviceClient<mavros_msgs::CommandBool>
            ("mavros/cmd/arming");
    set_mode_client = nh.serviceClient<mavros_msgs::SetMode>
            ("mavros/set_mode");
    
    // Timers
    cmdloop_timer = nh.createTimer(ros::Duration(0.01), cmdloopCallback); // send commands at 100 HZ
    simInit_timer = nh.createTimer(ros::Duration(1), simInitCallback); // send commands at 1 HZ
            
    // Parameters
    // TODO: Implement the sim_enable parameter along with launch file
    //nh_private_.param<bool>("enable_sim", sim_enable, true);
    sim_enable = true;
    
    // wait for FCU connection
    ros::Rate FCU_connect_rate(5.0);
    while(ros::ok() && !current_state.connected){
        ros::spinOnce();
        FCU_connect_rate.sleep();
    }

    // goal pose wrt home pose
    ref_offset.position.x = 0;
    ref_offset.position.y = 0;
    ref_offset.position.z = 1.5;

    ros::spin();
    return 0;
}

// ---------- Functions ----------- //

// Pass along mav state info
void stateCallback(const mavros_msgs::State::ConstPtr& msg){
    current_state = *msg;
}

void cmdloopCallback(const ros::TimerEvent& event) {
  switch (flight_state) {
    case WAITING_FOR_HOME_POSE: {
      waitForPose(&received_home_pose, "Waiting for home pose...");
      ROS_INFO("Got pose! Drone Ready to be armed.");
      flight_state = MISSION_EXECUTION;
      break;
    }

    case MISSION_EXECUTION: {
      double x_diff = ref_pose.position.x - curr_pose.position.x;
      double y_diff = ref_pose.position.y - curr_pose.position.y;
      double z_diff = ref_pose.position.z - curr_pose.position.z;
      double pos_error = sqrt(pow(x_diff,2.) + pow(y_diff,2.) + pow(z_diff,2.)); // in meters
        
      if (pos_error > 0.01) {
        ref_msg.header.stamp = ros::Time::now();
        ref_msg.pose = ref_pose;
        pos_pub.publish(ref_msg);
      } else {
        flight_state = LANDING;
      }
      break;
    }

    case LANDING: {
      geometry_msgs::PoseStamped landingmsg;
      landingmsg.header.stamp = ros::Time::now();
      landingmsg.pose = home_pose;
      landingmsg.pose.position.z = landingmsg.pose.position.z + 1.0;
      pos_pub.publish(landingmsg);
      flight_state = LANDED;
      ros::spinOnce();
      break;
    }
    case LANDED: {
      ROS_INFO("Landed. Please set to position control and disarm.");
      cmdloop_timer.stop();
      break;
    }
  }
}

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
  }
  curr_pose = msg.pose;
}

void simInitCallback(const ros::TimerEvent &event) {
    if (sim_enable) {
        ROS_INFO("sim init");
        // Enable OFFBoard mode and arm automatically
        // ONLY run if the vehicle is simulated
        arm_cmd.request.value = true;
        offb_set_mode.request.custom_mode = "OFFBOARD";
        if (current_state.mode != "OFFBOARD" && (ros::Time::now() - last_request > ros::Duration(5.0))) {
            if (set_mode_client.call(offb_set_mode) && offb_set_mode.response.mode_sent) {
            ROS_INFO("Offboard enabled");
            }
            last_request = ros::Time::now();
        } else {
            if (!current_state.armed && (ros::Time::now() - last_request > ros::Duration(5.0))) {
            if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                ROS_INFO("Vehicle armed");
            }
            last_request = ros::Time::now();
            }
        }
    } else {
        simInit_timer.stop();
    }
}
