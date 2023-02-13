#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <mavros_msgs/Thrust.h>
#include <tf/tf.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

using namespace std;

class Controller {
    private:
        const double g = 9.81;
        const Eigen::Vector3d e3(0,0,1);
        
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
        geometry_msgs::Vector3 err_att;
        
        // Geometric tracking controller
        void Geometric(Inputs& inputs, const FlatReference& ref, const States& cur) {
            // Generate control inputs to use with mavros, using the paper:
            // Geometric Tracking Control of a Quadrotor UAV on SE(3) (Lee et al., 2010)
            // Implementation inspired by https://github.com/Jaeyoung-Lim/mavros_controllers
            
            // Position and velocity errors
            err_pos = cur.pose.position - ref.position;    // position error
            err_vel = cur.velocity - ref.velocity;         // velocity error
            
            // Convert current attitude to Eigen rotation matrix R_curr
            geometry_msgs::Quaternion q = cur.pose.orientation;
            Eigen::Matrix3d R_curr = Eigen::Quaterniond(q.x, q.y, q.z, q.w).toRotationMatrix();
            
            // Compute reference attitude R_ref
            // Note: want ENU coord system, but the Lee paper uses NED,
            // so we make the changes: z <-> -z, and x <-> y
            geometry_msgs::Vector3 a = ref.acceleration;
            Eigen::Vector3d acc_ref(a.x, a.y, a.z);
            
            Eigen::Vector3d thrust_vect = - K_pos*err_pos - K_vel*err_vel - g*e3 + acc_ref; // vector in direction of desired thrust
            Eigen::Vector3d zb = thrust_vect / thrust_vect.norm();
            Eigen::Vector3d y_head(-sin(ref.yaw), cos(ref.yaw), 0); // heading vector
            Eigen::Vector3d xb = y_head.cross(zb) / (y_head.cross(zb)).norm();
            Eigen::Vector3d yb = zb.cross(xb);
            Eigen::Matrix3d R_ref;
            
            R_ref.col(0) = xb;
            R_ref.col(1) = yb;
            R_ref.col(2) = zb;
            
            // --- Attitude setpoint --- //
            Eigen::Quaterniond quat_ref(R_ref);
            tf::quaternionEigenToMsg(quat_ref, inputs.attitude.pose.orientation); // fill inputs.attitude.pose.orientation
            
            // Attitude error
            err_att = 0.5 * vee(R_ref.transpose()*R_curr - R_curr.transpose()*R_ref); // attitude error
            
            // --- Angular rate setpoint (K_att = 20.0) --- //
            Eigen::Vector3d ang_rate = K_att*err_att;
            tf::vectorEigenToMsg(ang_rate, inputs.cmd_vel.twist.angular); // fill inputs.cmd_vel.twist.angular

            // --- Thrust setpoint --- //
            double Thrust = thrust_vect.dot(R_curr*e3);
            inputs.thrust.thrust = Thrust;
        }
        
        // The "vee" operator converts a skew-symm matrix to a vector
        Eigen::Vector3d vee(const Eigen::Matrix3d& R){
            assert(R == -R.transpose()); // Check skew symmetry
            return Eigen::Vector3(R(2,1), R(0,2), R(1,0));
        }
};

#endif















