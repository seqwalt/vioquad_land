/**
 * @file mavros_cmd_node.cpp
 * @brief Offboard control node for hardware flight using mavros and PX4
 */

#include "mavros_cmd_node.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "mavros_cmd_node");
    ros::NodeHandle nh;

    MavrosCmd mavCmd;
    
    // goal pose wrt home pose
    mavCmd.ref_offset.position.x = 1;
    mavCmd.ref_offset.position.y = 1;
    mavCmd.ref_offset.position.z = 1.5;
    
    // Subscribers
    mavCmd.state_sub = nh.subscribe
            ("mavros/state", 1, &MavrosCmd::stateCallback, &mavCmd); // http://wiki.ros.org/roscpp_tutorials/Tutorials/UsingClassMethodsAsCallbacks
    mavCmd.pose_sub = nh.subscribe
            ("mavros/local_position/pose", 1, &MavrosCmd::mavposeCallback, &mavCmd, ros::TransportHints().tcpNoDelay());
    
    // Publishers
    mavCmd.pos_pub = nh.advertise<geometry_msgs::PoseStamped>
            ("mavros/setpoint_position/local", 10);
    
    // Services        
    mavCmd.arming_client = nh.serviceClient<mavros_msgs::CommandBool>
            ("mavros/cmd/arming");
    mavCmd.set_mode_client = nh.serviceClient<mavros_msgs::SetMode>
            ("mavros/set_mode");
    
    // Timers
    bool autostart = false;
    mavCmd.cmdloop_timer = nh.createTimer(ros::Duration(0.01), &MavrosCmd::cmdloopCallback, &mavCmd, false, autostart); // send commands at 100 HZ
    mavCmd.setup_timer = nh.createTimer(ros::Duration(1), &MavrosCmd::setupCallback, &mavCmd, false, autostart); // send commands at 1 HZ
            
    // Parameters
    // TODO: Implement the sim_enable parameter along with launch file
    //nh_private_.param<bool>("enable_sim", mavCmd.sim_enable, true);
    mavCmd.sim_enable = true; // used in setupCallback
    
    // wait for FCU connection
    ros::Rate FCU_connect_rate(5.0);
    while(ros::ok() && !mavCmd.current_state.connected){
        ros::spinOnce();
        FCU_connect_rate.sleep();
    }
    
    // Start command loop and setupInit
    mavCmd.cmdloop_timer.start();
    mavCmd.setup_timer.start();

    ros::spin();
    return 0;
}
