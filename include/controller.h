#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <mavros_msgs/Thrust.h>
#include <tf2_eigen/tf2_eigen.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

using namespace std;

class Controller {
    public:
        struct AttitudeInputs {      // Control inputs to be sent to FCU through mavros
            geometry_msgs::PoseStamped attitude;
            geometry_msgs::TwistStamped cmd_vel;
            mavros_msgs::Thrust norm_thrust;
        };
        
        struct PoseInputs {         // Control inputs to be sent to FCU through mavros
            geometry_msgs::PoseStamped pose;
        };
        
        struct FlatReference {      // Tracking reference (for using FlatOutputs message)
            geometry_msgs::Point position;
            geometry_msgs::Vector3 velocity;
            geometry_msgs::Vector3 acceleration;
            double yaw;
        };
        
        struct State {             // Quadcopter states
            geometry_msgs::Pose pose;
            geometry_msgs::Vector3 velocity;
        };
        
        Eigen::Vector3d err_pos;
        Eigen::Vector3d err_vel;
        Eigen::Vector3d err_att;
        
        Eigen::Vector3d K_pos = Eigen::Vector3d(8.0, 8.0, 10.0);
        Eigen::Vector3d K_vel = Eigen::Vector3d(1.5, 1.5, 3.0);
        Eigen::Vector3d K_att = Eigen::Vector3d(20.0, 20.0, 20.0);
        double thrust_const1 = 0.05;
        double thrust_const2 = 0.1;
        
        // Position + heading tracking controller
        // Directly feed through the setpoint positions and heading (yaw)
        void PosYaw(PoseInputs& inputs, const FlatReference& ref) {
            inputs.pose.pose.position = ref.position;
            
            Eigen::AngleAxisd temp_angAxis(ref.yaw, Eigen::Vector3d(0,0,1));  // rotate about z-axis
            Eigen::Quaterniond quat_ref(temp_angAxis);
            inputs.pose.pose.orientation = tf2::toMsg(quat_ref);
        }
        
        // Geometric tracking controller
        void Geometric(AttitudeInputs& inputs, const FlatReference& ref, const State& cur) {
            // Generate control inputs to use with mavros, using the paper:
            // Geometric Tracking Control of a Quadrotor UAV on SE(3) (Lee et al., 2010)
            // Implementation inspired by https://github.com/Jaeyoung-Lim/mavros_controllers
            
            // Position and velocity errors
            err_pos = vectToEigen(cur.pose.position) - vectToEigen(ref.position);    // position error
            err_vel = vectToEigen(cur.velocity) - vectToEigen(ref.velocity);         // velocity error
            
            // Convert current attitude to Eigen rotation matrix R_curr
            geometry_msgs::Quaternion q = cur.pose.orientation;
            Eigen::Matrix3d R_curr = Eigen::Quaterniond(q.x, q.y, q.z, q.w).toRotationMatrix();
            
            // Compute reference attitude R_ref
            // Note: want ENU coord system, but the Lee paper uses NED,
            // so we make the changes: z <-> -z, and x <-> y
            geometry_msgs::Vector3 a = ref.acceleration;
            Eigen::Vector3d acc_ref(a.x, a.y, a.z);
            
            Eigen::Vector3d thrust_vect = acc_ref - K_pos.asDiagonal()*err_pos - K_vel.asDiagonal()*err_vel - g*e3; // vector in direction of desired thrust
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
            inputs.attitude.pose.orientation = tf2::toMsg(quat_ref);
            
            // Attitude error
            err_att = 0.5 * vee(R_ref.transpose()*R_curr - R_curr.transpose()*R_ref); // attitude error
            
            // --- Angular rate setpoint (K_att = 20.0) --- //
            Eigen::Vector3d ang_rate = K_att.asDiagonal()*err_att;
            tf2::toMsg(ang_rate, inputs.cmd_vel.twist.angular); // fill inputs.cmd_vel.twist.angular

            // --- Thrust setpoint --- //
            double Thrust = thrust_vect.dot(R_curr*e3);
            inputs.norm_thrust.thrust = std::max(0.0, std::min(1.0, thrust_const1 * Thrust + thrust_const2));
        }
        
        template <class T>
        Eigen::Quaterniond quatToEigen(const T& msg){
            Eigen::Quaterniond q;
            q.w() = msg.w;
            q.x() = msg.x;
            q.y() = msg.y;
            q.z() = msg.z;
            return q;
        };
        
    private:
        const double g = 9.81;
        const Eigen::Vector3d e3 = Eigen::Vector3d(0,0,1);
        
        // The "vee" operator converts a skew-symm matrix to a vector
        Eigen::Vector3d vee(const Eigen::Matrix3d& R){
            assert(R == -R.transpose()); // Check skew symmetry
            return Eigen::Vector3d(R(2,1), R(0,2), R(1,0));
        }
        
        template <class T>
        Eigen::Vector3d vectToEigen(const T& msg){
            return Eigen::Vector3d(msg.x, msg.y, msg.z);
        };
};

#endif















