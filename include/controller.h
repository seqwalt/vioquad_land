#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <mavros_msgs/Thrust.h>
#include <mavros_msgs/PositionTarget.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <tf2_eigen/tf2_eigen.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

using namespace std;

class Controller {
    public:        
        struct State {             // Quadcopter states
            geometry_msgs::Pose pose;
            geometry_msgs::Vector3 velocity;
        };
        
        Eigen::Vector3d err_pos;
        Eigen::Vector3d err_vel;
        Eigen::Vector3d err_att;
        
        double Kpos_x, Kpos_y, Kpos_z;
        double Kvel_x, Kvel_y, Kvel_z;
        double Katt_x, Katt_y, Katt_z;
        double thrust_const, thrust_offset, max_err_acc;

        Eigen::Vector3d K_pos = Eigen::Vector3d(Kpos_x, Kpos_y, Kpos_z);
        Eigen::Vector3d K_vel = Eigen::Vector3d(Kvel_x, Kvel_y, Kvel_z);
        Eigen::Vector3d K_att = Eigen::Vector3d(Katt_x, Katt_y, Katt_z);
        
        // Position + heading tracking controller
        // Directly feed through the setpoint positions and heading (yaw)
//         void PosYaw(geometry_msgs::PoseStamped& inputs, const FlatReference& ref) {
//             inputs.pose.position = ref.position;
//             
//             Eigen::AngleAxisd temp_angAxis(ref.yaw, Eigen::Vector3d(0,0,1));  // rotate about z-axis
//             Eigen::Quaterniond quat_ref(temp_angAxis);
//             inputs.pose.orientation = tf2::toMsg(quat_ref);
//         }
        
        // Position + heading tracking controller
        // Give PX4 the setpoint positions/velocities/accels and heading (yaw)
        // PX4 uses the following control scheme:
        //      The velocity and the acceleration setpoints are used as feedforwards;
        //      the velocity setpoint is added to the output of the position controller
        //      and the result is used as the input to the velocity controller;
        //      the acceleration setpoint is added to the output of the velocity controller
        //      and the result used to compute the thrust vector
        void PosYaw(mavros_msgs::PositionTarget& inputs, const quad_control::FlatOutputs &ref) {
            inputs.position = ref.position;
            inputs.velocity = ref.velocity;
            inputs.acceleration_or_force = ref.acceleration;
            inputs.yaw = ref.yaw;
        }
        
        // Geometric tracking controller
        void Geometric(mavros_msgs::AttitudeTarget& inputs, const quad_control::FlatOutputs &ref, const State& cur) {
            // Generate control inputs to use with mavros, using the paper:
            // Geometric Tracking Control of a Quadrotor UAV on SE(3) (Lee et al., 2010)
            // Implementation inspired by https://github.com/Jaeyoung-Lim/mavros_controllers
            
            // Rotation matrix ENU <-> NED. Note symmetry
            Eigen::Matrix3d ENUtoNED;
            ENUtoNED << 0.0, 1.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0,-1.0;
            
            // Position and velocity errors
            err_pos = vectToEigen(cur.pose.position) - vectToEigen(ref.position);    // position error
            err_vel = vectToEigen(cur.velocity) - vectToEigen(ref.velocity);         // velocity error
            
            // Convert current attitude to Eigen rotation matrix R_curr
            geometry_msgs::Quaternion q = cur.pose.orientation;
            Eigen::Matrix3d R_curr = ENUtoNED * Eigen::Quaterniond(q.x, q.y, q.z, q.w).toRotationMatrix(); // convert to NED
            
            // Compute reference attitude R_ref
            // We use the ENU frame here, then convert to NED at the end
            geometry_msgs::Vector3 a = ref.acceleration;
            Eigen::Vector3d acc_ref(a.x, a.y, a.z); // convert from ENU to NED
            
            Eigen::Vector3d acc_err = K_pos.asDiagonal()*err_pos + K_vel.asDiagonal()*err_vel;
            Eigen::Vector3d acc_des = acc_ref + g*e3 - clip(acc_err, max_err_acc); // vector in direction of desired thrust
            Eigen::Vector3d zb = acc_des / acc_des.norm();
            Eigen::Vector3d y_heading = Eigen::Vector3d(-sin(ref.yaw), cos(ref.yaw), 0); // heading vector
            Eigen::Vector3d xb = y_heading.cross(zb) / (y_heading.cross(zb)).norm();
            Eigen::Vector3d yb = zb.cross(xb);
            Eigen::Matrix3d R_ref;
            
            R_ref.col(0) = xb;
            R_ref.col(1) = yb;
            R_ref.col(2) = zb;
            
            // --- Attitude setpoint --- //
            Eigen::Quaterniond quat_ref(ENUtoNED * R_ref); // convert to NED
            inputs.orientation = tf2::toMsg(quat_ref);
            
            // Attitude error
            err_att = 0.5 * vee(R_ref.transpose()*R_curr - R_curr.transpose()*R_ref); // attitude error
            
            // --- Angular rate setpoint (K_att = 20.0) --- //
            Eigen::Vector3d ang_rate = K_att.asDiagonal()*err_att;
            tf2::toMsg(ang_rate, inputs.body_rate); // fill inputs.body_rate

            // --- Thrust setpoint --- //
            double Acc_des = (ENUtoNED * acc_des).dot(R_curr*e3);
            inputs.thrust = (float)std::max(0.0, std::min(1.0, thrust_const * Acc_des + thrust_offset));
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
            return Eigen::Vector3d(msg.y, msg.x, -msg.z);
        };
        
        Eigen::Vector3d clip(const Eigen::Vector3d& vec, const double& mag) {
            // IF vec has magnitude greater than mag, then
            // scale it's magnitude to mag and return scaled vec.
            // ELSE, return unscaled vec.
            return (vec.norm() > mag) ? vec*(mag/vec.norm()) : vec;
        }
};

#endif















