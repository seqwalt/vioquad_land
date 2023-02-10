/**
 * @file offb_node.cpp
 * @brief Offboard control node for hardware flight using mavros and PX4
 */

#include "offb_node.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "offb_node");
    ros::NodeHandle nh;

    OffbCtrl offbCtrl;
    
    // goal pose wrt home pose
    offbCtrl.ref_offset.position.x = 1;
    offbCtrl.ref_offset.position.y = 1;
    offbCtrl.ref_offset.position.z = 1.5;
    
    // Subscribers
    offbCtrl.state_sub = nh.subscribe
            ("mavros/state", 1, &OffbCtrl::stateCallback, &offbCtrl); // http://wiki.ros.org/roscpp_tutorials/Tutorials/UsingClassMethodsAsCallbacks
    offbCtrl.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &OffbCtrl::mavposeCallback, &offbCtrl, ros::TransportHints().tcpNoDelay());
    
    // Publishers
    offbCtrl.pos_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mavros/setpoint_position/local", 10);
    
    // Services        
    offbCtrl.arming_client = nh.serviceClient<mavros_msgs::CommandBool>
            ("mavros/cmd/arming");
    offbCtrl.set_mode_client = nh.serviceClient<mavros_msgs::SetMode>
            ("mavros/set_mode");
    
    // Timers
    bool oneshot = false;
    bool autostart = false;
    offbCtrl.cmdloop_timer = nh.createTimer(ros::Duration(0.01), &OffbCtrl::cmdloopCallback, &offbCtrl, oneshot, autostart); // send commands at 100 HZ
    offbCtrl.sim_timer = nh.createTimer(ros::Duration(1), &OffbCtrl::simCallback, &offbCtrl, oneshot, autostart); // send commands at 1 HZ
            
    // Parameters
    // TODO: Implement the sim_enable parameter along with launch file
    //nh_private_.param<bool>("enable_sim", offbCtrl.sim_enable, true);
    offbCtrl.sim_enable = true;
    
    // wait for FCU connection
    ros::Rate FCU_connect_rate(5.0);
    while(ros::ok() && !offbCtrl.current_state.connected){
        ros::spinOnce();
        FCU_connect_rate.sleep();
    }
    
    // Start command loop and simInit
    offbCtrl.cmdloop_timer.start();
    offbCtrl.sim_timer.start();

    ros::spin();
    return 0;
}

// ---------- Functions ----------- //

// Pass along mav state info
void OffbCtrl::stateCallback(const mavros_msgs::State::ConstPtr& msg){
    current_state = *msg;
}

// Publish flight commands according to current flight state
void OffbCtrl::cmdloopCallback(const ros::TimerEvent& event) {
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
      takeoff_msg.pose.position.z += 1.0;
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
      if (pos_error > 0.06) {
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
void OffbCtrl::mavposeCallback(const geometry_msgs::PoseStamped &msg) {
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
void OffbCtrl::simCallback(const ros::TimerEvent &event) {
    if (sim_enable) {
        switch(flight_state) {
            case LANDING: {
                // try to disarm
                arm_cmd.request.value = false;
                if (current_state.armed && (ros::Time::now() - last_request > ros::Duration(2.0))) {
                    if (arming_client.call(arm_cmd) && arm_cmd.response.success) {
                        ROS_INFO("Vehicle disarmed");
                        sim_timer.stop();
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
        sim_timer.stop();
    }
}
