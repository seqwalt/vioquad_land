#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <mavros_msgs/Thrust.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

using namespace std;

class Controller {
    public:
        struct Inputs {         // Control inputs to be sent to FCU through mavros
            geometry_msgs::PoseStamped attitude;
            geometry_msgs::TwistStamped cmd_vel;
            mavros_msgs::Thrust thrust;
        };
        
        struct FlatReference {  // Tracking reference (for using FlatOutputs message)
            geometry_msgs::Point position;
            geometry_msgs::Vector3 velocity;
            geometry_msgs::Vector3 acceleration;
            double yaw;
        };
        
        struct States {         // Quadcopter states
            geometry_msgs::Pose pose;
            geometry_msgs::Vector3 velocity;
        };
        
        geometry_msgs::Point err_pos;
        geometry_msgs::Vector3 err_vel;
        
        // Geometric tracking controller
        void Geometric(&Inputs inputs, const &FlatReference ref, const &States cur) {
            // Generate control inputs to use with mavros, using the paper:
            // Geometric Tracking Control of a Quadrotor UAV on SE(3) (Lee et al., 2010)
            // Implementation inspired by https://github.com/Jaeyoung-Lim/mavros_controllers
            
            // Position and velocity errors
            err_pos = cur.position - ref.pose.position;    // position error
            err_vel = cur.velocity - ref.velocity;         // velocity error
            
            // Convert current attitude to Eigen rotation matrix R_curr
            geometry_msgs::Quaternion q = cur.pose.orientation;
            Eigen::Matrix3d R_curr = Eigen::Quaterniond(q.x, q.y, q.z, q.w).toRotationMatrix();
            
            // Compute reference attitude R_ref (Note: need NEU )
            Eigen::Vector3d zb; // unit vector in direction of thrust
            
            err_att = 0.5 * ; // attitude error
        }
};

#endif















