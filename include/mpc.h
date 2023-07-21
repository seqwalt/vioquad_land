#pragma once

#include <ros/ros.h>
#include <mavros/mavros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <sensor_msgs/Imu.h>
#include <apriltag_ros/AprilTagDetectionArray.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/Thrust.h>
#include <mavros_msgs/AttitudeTarget.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>
#include <stdio.h>
#include <random>

#include <std_srvs/Trigger.h>
#include "quad_control/InitSetpoint.h" // custom service
#include "quad_control/FlatOutputs.h"  // custom message

// for trajectory generation
#include "min_snap_traj.h"
#include "filesystem.hpp"
#include <iomanip>
#include <chrono>
#include <ctime>

// for mpc
#include "acado_common.h"
#include "acado_auxiliary_functions.h"
#include "acado_messages.h" // custom header (does not come with ACADO)

using namespace std;

// Some convenient private ACADO definitions
#define NX          ACADO_NX  // Number of differential state variables
#define NXA         ACADO_NXA // Number of algebraic variables
#define NU          ACADO_NU  // Number of control inputs
#define NOD         ACADO_NOD // Number of online data values

#define NY          ACADO_NY  // Number of measurements/references on nodes 0..N - 1
#define NYN         ACADO_NYN // Number of measurements/references on node N

#define N           ACADO_N   // Number of intervals in the horizon
#define NUM_STEPS   3         // Number of real-time iterations
#define VERBOSE     1         // Show iterations: 1, silent: 0

// global variables used by the solver
ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

class MPC {
    public:
        ros::Publisher mpc_pub;
        ros::Publisher ref_total_pub;
        ros::Publisher ref_curr_pub;
        ros::Publisher pred_pub;
        ros::Publisher est_pub;
        ros::Publisher tag_pub;
        
        ros::Subscriber pose_sub;
        ros::Subscriber vel_sub;
        ros::Subscriber imu_sub;
        ros::Subscriber apriltag_sub;
        ros::Timer mpc_timer;
        ros::Timer path_timer;
        
        quad_control::FlatOutputs ref;
        
        // Parameters
        bool sim_enable;
        double tag_smoothing_factor;
        double tran_BC_x;
        double tran_BC_y;
        double tran_BC_z;
        double land_spd_factor;
        double land_height;
        bool Do_Fov;
        bool Use_Percep_Cost; // TODO make use_percep_cost a roslaunch param

        // MPC class constructor
        MPC(const std::string& searchTrajFileName){
            // Clear solver memory.
            memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
            memset(&acadoVariables, 0, sizeof( acadoVariables ));

            // Initialize the solver
            acado_initializeSolver();

            // Initialize acadoVariables for horizontal search path
            traj = loadSearchTraj(searchTrajFileName);
            bool use_percep_cost = false;
            initAcadoVariables(use_percep_cost);
            
            if( VERBOSE ) acado_printHeader();
            ROS_INFO_STREAM("MPC solver initialized.");

            first_mpc_call = true;
            near_landing_pad = false;
            has_landed = false;            
            traj_gen_time = chrono::steady_clock::now();
            prev_time2land = 0;
            first_april_tag = true; // TODO: for Apriltag as state est
            tag_visible = false;
            
            // Prepare first step
            acado_preparationStep();
        }

        // Send initial pose back to mavros_cmd_node, upon request
        bool initRefCallback(quad_control::InitSetpoint::Request &req, quad_control::InitSetpoint::Response &res){
            res.success = true;
            res.position.x = traj(0,1);
            res.position.y = traj(0,2);
            res.position.z = traj(0,3);
            res.yaw = quaternion2yaw(Eigen::Quaterniond(traj(0,7), traj(0,8), traj(0,9), traj(0,10)));
            ROS_INFO_STREAM("Initial MPC pose sent.");
            return true;
        }

        // Start ROS timer for MPC control (mpcCallback), upon request
        bool streamMpcCallback(std_srvs::Trigger::Request &req, std_srvs::Trigger::Response &res){
            // Align first trajectory quaternion with current quadcopter quaternion
            Eigen::Quaterniond first_traj_quat(traj(0,7), traj(0,8), traj(0,9), traj(0,10));
            Eigen::Quaterniond curQuat = msg2quaternion(curPose.orientation);
            if (quatsDiffSign(first_traj_quat, curQuat)) flipQuatSigns(traj);
            
            mpc_timer.start();  // start publishing the MPC control values
            if(mpc_timer.hasStarted()){
                res.success = true;
            }
            return true;
        }

        // Publish tragectories and poses for rviz visualization
        void viewPathCallback(const ros::TimerEvent &event){
            // Publish reference trajectory for visualization
            ref_total_pub.publish(mpcTotalRef);
            if (tag_visible){
                // Publish estimated AprilTag pose in map/world frame
                tag_pub.publish(tagPose);
            }
            
            if(!first_mpc_call){
                // Publish estimated path for visualization (vio or ekf2, depending on direct_from_vio param)
                addPoseToPath(curPose, mpcEst);
                est_pub.publish(mpcEst);
                
                // Publish reference for current mpc iteration
                acadoArr2PathMsg(acadoVariables.y, NY, mpcCurrRef);
                ref_curr_pub.publish(mpcCurrRef);
                
                // Publish state prediction for current mpc iteration
                acadoArr2PathMsg(acadoVariables.x, NX, mpcPred);
                pred_pub.publish(mpcPred);
            }
        }
        
        // Track the trajectory with an acado MPC by
        // publishing thrust and angular rates to mavros
        void mpcCallback(const ros::TimerEvent &event){
            if (!near_landing_pad){
                // Not close enough to landing pad to turn off motors
                //acado_tic( &t ); // for timing mpc

                // ------------------ Update Acado Variables ------------------ //
                
                // Current quadcopter state
                acadoVariables.x0[0] = curPose.position.x;
                acadoVariables.x0[1] = curPose.position.y;
                acadoVariables.x0[2] = curPose.position.z;
                acadoVariables.x0[3] = curPose.orientation.w;
                acadoVariables.x0[4] = curPose.orientation.x;
                acadoVariables.x0[5] = curPose.orientation.y;
                acadoVariables.x0[6] = curPose.orientation.z;
                acadoVariables.x0[7] = curVel.x;
                acadoVariables.x0[8] = curVel.y;
                acadoVariables.x0[9] = curVel.z;
                
                double horizon_start; // start time of mpc time horizon
                int num_steps;
                if (first_mpc_call){
                    first_mpc_call = false;
                    num_steps = 3;
                    mpc_start_time = ros::Time::now().toSec();
                    horizon_start = 0.0;
                } else {
                    num_steps = NUM_STEPS;
                    double max_start = traj(traj.rows()-1,0) - traj(0,0);
                    horizon_start = min(ros::Time::now().toSec() - mpc_start_time, max_start);
                }
                
                updateAcadoReferences(horizon_start);
                
                // ------------------ MPC optimization ------------------ //
                int iter, status;
                // Real-time iteration (RTI) loop
                for(iter = 0; iter < num_steps; ++iter){
                    // Perform the feedback step
                    status = acado_feedbackStep();

                    if (status != 0){
                        cout << "Iteration:" << iter << ", QP problem! QP status: " << status << endl;
                        ACADO_MSG::showStatus(status);
                        ros::shutdown();
                        break;
                    }

                    // Prepare for the next step
                    acado_preparationStep();
                }

                // Extract control input
                T  = acadoVariables.u[0]; // mass-normalized thrust
                wx = acadoVariables.u[1];
                wy = acadoVariables.u[2];
                wz = acadoVariables.u[3];
                
            } else {
                // Near the landing pad, so set control inputs to zero
                T = 0.0;
                wx = 0.0; wy = 0.0; wz = 0.0;
                
                has_landed = true;
                ROS_INFO_STREAM_ONCE(endl << " ---- The quadcopter has landed! ---- " << endl);
            }
            
            // Apply control input
            mpcInputs.header.stamp = ros::Time::now();
            mpcInputs.type_mask = 128; // Ignore attitude messages
            mpcInputs.body_rate.x = wx;
            mpcInputs.body_rate.y = wy;
            mpcInputs.body_rate.z = wz;
            mpcInputs.thrust = thrustMap(T);
            mpc_pub.publish(mpcInputs); // set attitude, body rate and thrust to mavros
            
            //acado_printDifferentialVariables();
            //acado_printControlVariables();        // debugging info
            //float te = acado_toc( &t );           // for timing mpc
            //cout << "mpcCallback time: " << te << "sec" << endl;
        }

        // Pass along current mav pose data, wrt world frame.
        void mavPoseCallback(const geometry_msgs::PoseStamped &msg) {
            OdomData data;
            data.source = POSE_BY_EKF2;
            data.pose = msg.pose;
            odomMsgSelector(data);
        }

        // Pass along current mav velocity data, wrt world frame.
        void mavVelCallback(const geometry_msgs::TwistStamped &msg) {
            OdomData data;
            data.source = VEL_BY_EKF2;
            data.velocity = msg.twist.linear;
            odomMsgSelector(data);
        }
        
        // Pass along current raw mav imu data, wrt body frame.
        void mavIMUCallback(const sensor_msgs::Imu &msg) {
            curLinAcc = msg.linear_acceleration; // linear acceleration in imu body frame
            curYawVel = msg.angular_velocity.z;
        }
        
        // Compute AprilTag pose in World frame (compute transformation T_WA)
        void mavAprilTagCallback(const apriltag_ros::AprilTagDetectionArray &msg) {
            // Check if tag is visible
            tag_visible = !msg.detections.empty();
            if (tag_visible && !has_landed){
                // Form T_WB (load from state estimation: mocap, VIO, etc.)
                Eigen::Isometry3d T_WB = Eigen::Isometry3d::Identity();
                Eigen::Vector3d tran_WB(curPose.position.x, curPose.position.y, curPose.position.z); // translation part
                Eigen::Quaterniond quat_WB = msg2quaternion(curPose.orientation); // [w, x, y, z] rotation part
                T_WB.translate(tran_WB);
                T_WB.rotate(quat_WB);
                
                // Form T_BC (transformation from quadcopter body frame to down-facing camera frame)
                // TODO: load T_BC info in from launch file
                Eigen::Isometry3d T_BC = Eigen::Isometry3d::Identity();
                Eigen::Vector3d tran_BC(tran_BC_x, tran_BC_y, tran_BC_z); // translation part
                Eigen::Matrix3d rot_BC;
                // TODO: fix sdf to apriltag stuff
                //vector<double> sdf_rpy{0.0, 1.5708, 0.0}; // roll pitch yaw (XYZ form), given by iris_downward_camera.sdf
                //vector<double> apriltag_rpy{sdf_rpy[2], -sdf_rpy[0], sdf_rpy[1]}; // rotation from sdf to apriltag representaiton (https://github.com/AprilRobotics/apriltag/wiki/AprilTag-User-Guide#coordinate-system).
                vector<double> apriltag_rpy{3.1415, 0.0, -1.5708};
                rot_BC = Eigen::AngleAxisd(apriltag_rpy[2], Eigen::Vector3d::UnitZ())
                       * Eigen::AngleAxisd(apriltag_rpy[1], Eigen::Vector3d::UnitY())
                       * Eigen::AngleAxisd(apriltag_rpy[0], Eigen::Vector3d::UnitX()); // XYZ rotation http://sdformat.org/tutorials?tut=specify_pose
                T_BC.translate(tran_BC);
                T_BC.rotate(rot_BC);
                
                // Compute T_CA (transformation from down-facing camera frame to AprilTag frame)
                geometry_msgs::Pose AprilTagPose_Cam = msg.detections[0].pose.pose.pose; // pose of AprilTag in camera frame
                Eigen::Vector3d tran_CA(AprilTagPose_Cam.position.x, AprilTagPose_Cam.position.y, AprilTagPose_Cam.position.z); // translation part
                Eigen::Quaterniond quat_CA = msg2quaternion(AprilTagPose_Cam.orientation); // rotation part
                Eigen::Isometry3d T_CA = Eigen::Isometry3d::Identity();
                T_CA.translate(tran_CA);
                T_CA.rotate(quat_CA); // transformation matrix representing the AprilTag in the camera frame
                near_landing_pad = (abs(tran_CA[2]) < land_height); // check if near landing pad
                
                // Compute T_WA
                Eigen::Isometry3d T_WA = T_WB * T_BC * T_CA;
                Eigen::Vector3d tag_pos_raw = T_WA.translation();
                Eigen::Matrix3d tag_rot_raw = T_WA.rotation();
                Eigen::Quaterniond tag_quat_raw = Eigen::Quaterniond(tag_rot_raw);
                static Eigen::Quaterniond tagQuat = tag_quat_raw;  // initialization for filter
                
                // Exponential smoothing of AprilTag pose
                if(first_april_tag){
                    tagPos = T_WA.translation(); // initialization for filter
                    first_april_tag = false;
                }
                Eigen::Vector3d z_tag = tag_quat_raw*Eigen::Vector3d::UnitZ();
                Eigen::Vector3d z_world = Eigen::Vector3d::UnitZ();
                bool tag_outlier = z_world.dot(z_tag) < cos(M_PI/8);  // assume tag is facing close to up
                if (!tag_outlier) {
                    double a = tag_smoothing_factor; // (0,1] smoothing factor
                    assert(a >= 0.0 && a < 1.0);
                    tagPos = a*tag_pos_raw + (1 - a)*tagPos; // s_t = a*x_t + (1-a)*s_{t-1}
                    if (tagQuat.coeffs().dot(tag_quat_raw.coeffs()) < 0) tag_quat_raw.coeffs() *= -1; // correct for "double cover"
                    tagQuat.coeffs() = a*tag_quat_raw.coeffs() + (1 - a)*tagQuat.coeffs();
                } else {
                    return;
                }
                
                Eigen::Quaterniond quatDiff = quat_WB.inverse() * tagQuat;
                yawDiff = quaternion2yaw(quatDiff); // yaw difference with correct sign
                 
                // Fill tagPose
                geometry_msgs::PoseStamped tagPoseTemp;
                tagPoseTemp.pose.position.x = tagPos(0);
                tagPoseTemp.pose.position.y = tagPos(1);
                tagPoseTemp.pose.position.z = tagPos(2);
                tagPoseTemp.pose.orientation.w = tagQuat.coeffs().w();
                tagPoseTemp.pose.orientation.x = tagQuat.coeffs().x();
                tagPoseTemp.pose.orientation.y = tagQuat.coeffs().y();
                tagPoseTemp.pose.orientation.z = tagQuat.coeffs().z();
                tagPoseTemp.header.frame_id = "map";
                tagPoseTemp.header.stamp = ros::Time::now();
                tagPose = tagPoseTemp;
                
                // TODO for Apriltag as state est
                // Use AprilTag as position estimation of quadcopter
//                 if (first_april_tag){
//                     first_april_tag = false;
//                     T_WA_initial = T_WA;
//                 }
//                 Eigen::Isometry3d T_WB_tag = T_WA_initial * T_CA.inverse() * T_BC.inverse(); // T_WB estimate, using apriltag as reference
//                 OdomData odom_data;
//                 Eigen::Vector3d pos_WB_tag = T_WB_tag.translation();
//                 odom_data.pose.position.x = pos_WB_tag(0);
//                 odom_data.pose.position.x = pos_WB_tag(1);
//                 odom_data.pose.position.x = pos_WB_tag(2);
//                 odom_data.source = APRIL_TAG;
//                 odomMsgSelector(odom_data);
//                 ROS_INFO_STREAM_ONCE("Using AprilTag as position estimate.");
                
                
                // Gerenerate landing trajectory
                double h_land = tran_CA[2]; // height above landing pad
                double vel_base = 0.43;     // average speed if curVel.z = 0
                double scale_val = land_spd_factor;    // (0.45 for gt, 0.35 for VIO) increasing scale_val causes shorter land times when curVel.z < 0, and longer times when curVel.z > 0
                double avg_spd = abs(vel_base*exp(-scale_val*min(max(curVel.z, -1.0), 1.0))); // average spd (m/s)
                double time2land = h_land/avg_spd; // time to land
                chrono::steady_clock::time_point curr_time = chrono::steady_clock::now();
                chrono::duration<double> time_since_last = chrono::duration_cast<chrono::duration<double>>(curr_time - traj_gen_time);
                updateAcadoOnlineData(Use_Percep_Cost); // needs to be before the following if statement, for current mpc
                if ((double)time_since_last.count() > 0.6*prev_time2land
                    && time2land > 0.6
                    && checkFOVFeasible(Do_Fov))
                {
                    prev_time2land = time2land;
                    traj_gen_time = chrono::steady_clock::now();
                    traj.resize(0,0);
                    traj = genLandingTraj(time2land, Do_Fov); // generate trajectory. each row has t, x, y, z, vx, vy, vz, qw, qx, qy, qz
                    initAcadoVariables(Use_Percep_Cost);
                    first_mpc_call = true; // reset for new trajectory
                    
                    eigen2PathMsg(traj, mpcTotalRef); // update ROS Trajectory message
                    //saveTraj("landing_traj", traj);   // save trajectory to file TODO: make bool ros param for whether or not to save
                }
            } else {
                // AprilTag is not detected, reset
                first_april_tag = true; // TODO: for Apriltag as state est
            }
        }

    private:
        // Pose data
        double curYaw;
        double yawDiff;
        double curYawVel;
        geometry_msgs::Pose curPose;
        geometry_msgs::Vector3 curVel;
        geometry_msgs::Vector3 curLinAcc;
        geometry_msgs::PoseStamped tagPose;
        Eigen::Vector3d tagPos;
        enum OdomSource {
            APRIL_TAG,
            POSE_BY_EKF2,
            VEL_BY_EKF2
        };
        struct OdomData {
            OdomSource source;
            geometry_msgs::Pose pose;
            geometry_msgs::Vector3 velocity;
        };
        
        // MPC data and info
        mavros_msgs::AttitudeTarget mpcInputs;
        nav_msgs::Path mpcTotalRef;  // total reference
        nav_msgs::Path mpcCurrRef;   // current iteration reference
        nav_msgs::Path mpcPred;      // mpc prediction
        nav_msgs::Path mpcEst;       // estimate from ekf2 as a path
        double mpc_start_time;
        bool first_mpc_call;
        const double mpc_time_horizon = 2.0; // From ../acado/PAMPC/quadrotor_pampc.cpp
        double T, wx, wy, wz;
        acado_timer t;
        
        // Trajectory data, info, data struct
        Eigen::MatrixXd traj;
        chrono::steady_clock::time_point traj_gen_time;
        struct MinSnapData{
            MinSnapTraj::vectOfMatrix24d bounds;
            vector<MinSnapTraj::keyframe> keyframes;
            MinSnapTraj::FOVdata fov_data;
        };
        
        // Landing info
        bool near_landing_pad;
        double prev_time2land;
        bool has_landed;
        Eigen::Isometry3d T_WA_initial;
        bool first_april_tag; // TODO: for Apriltag as state est
        bool tag_visible;
        
        // Initialize the acado variables
        void initAcadoVariables(bool use_percep_cost){
            // Check if trajectory is too short
            double traj_duration = traj(traj.rows()-1,0) - traj(0,0);
            if (traj_duration < mpc_time_horizon){
                ROS_ERROR_STREAM("Need traj_duration > mpc_time_horizon. Shutting down node ...");
                ros::shutdown(); 
            }
            
            // Initialize the references for start of trajectory (time 0.0)
            updateAcadoReferences(0.0);
            
            // Initialize the controls
            initAcadoControls();

            // Initialize the states and state-feedback vector
            initAcadoStates();

            // Initialize the weight acado MPC weight matrix
            initAcadoWeights(use_percep_cost);
            
            // NOTE acado online data is updated separately in mavAprilTagCallback(),
            //   since it updates the PAMPC point of interest based on the AprilTag location
        }

        // Initialize the acado state variables and state-feedback vector
        void initAcadoStates() {
            int row1, row2, row_try;
            double s;
            double dt_mpc = mpc_time_horizon/N; // horizon/(num intervals)
            double dt_traj = traj(1,0) - traj(0,0);

            // Update the states
            for (int i = 0; i < N + 1; ++i){
                row_try = (int)floor(i*dt_mpc/dt_traj);
                row1 = min((int)traj.rows()-1, row_try);
                row2 = min((int)traj.rows()-1, row1 + 1);
                s = (i*dt_mpc/dt_traj) - row_try; // interpolation parameter
                acadoVariables.x[i * NX + 0] = lerp(traj(row1, 1), traj(row2, 1), s); // x
                acadoVariables.x[i * NX + 1] = lerp(traj(row1, 2), traj(row2, 2), s); // y
                acadoVariables.x[i * NX + 2] = lerp(traj(row1, 3), traj(row2, 3), s); // z
                acadoVariables.x[i * NX + 3] = lerp(traj(row1, 7), traj(row2, 7), s); // qw
                acadoVariables.x[i * NX + 4] = lerp(traj(row1, 8), traj(row2, 8), s); // qx
                acadoVariables.x[i * NX + 5] = lerp(traj(row1, 9), traj(row2, 9), s); // qy
                acadoVariables.x[i * NX + 6] = lerp(traj(row1,10), traj(row2,10), s); // qz
                acadoVariables.x[i * NX + 7] = lerp(traj(row1, 4), traj(row2, 4), s); // vx
                acadoVariables.x[i * NX + 8] = lerp(traj(row1, 5), traj(row2, 5), s); // vy
                acadoVariables.x[i * NX + 9] = lerp(traj(row1, 6), traj(row2, 6), s); // vz
            }
            
            // Initialize state-feedback vector
            #if ACADO_INITIAL_STATE_FIXED
                for(int i = 0; i < NX; ++i)
                    acadoVariables.x0[i] = acadoVariables.x[i];
            #endif
        }
        
        // Update Acado references, given a trajectory start time (horizon_start) 
        void updateAcadoReferences(double horizon_start) {
            int row1, row2, row_try;
            double s;
            double dt_mpc = mpc_time_horizon/N; // horizon/(num intervals)
            double dt_traj = traj(1,0) - traj(0,0);
            
            // Update the references
            for (int i = 0; i < N; ++i){ // iterate over time steps within horizon
                row_try = (int)floor((horizon_start + i*dt_mpc)/dt_traj);
                row1 = min((int)traj.rows()-1, row_try);
                row2 = min((int)traj.rows()-1, row1 + 1);
                s = ((horizon_start + i*dt_mpc)/dt_traj) - row_try; // interpolation parameter
                acadoVariables.y[i * NY + 0] = lerp(traj(row1, 1), traj(row2, 1), s); // x
                acadoVariables.y[i * NY + 1] = lerp(traj(row1, 2), traj(row2, 2), s); // y
                acadoVariables.y[i * NY + 2] = lerp(traj(row1, 3), traj(row2, 3), s); // z
                acadoVariables.y[i * NY + 3] = lerp(traj(row1, 7), traj(row2, 7), s); // qw
                acadoVariables.y[i * NY + 4] = lerp(traj(row1, 8), traj(row2, 8), s); // qx
                acadoVariables.y[i * NY + 5] = lerp(traj(row1, 9), traj(row2, 9), s); // qy
                acadoVariables.y[i * NY + 6] = lerp(traj(row1,10), traj(row2,10), s); // qz 
                acadoVariables.y[i * NY + 7] = lerp(traj(row1, 4), traj(row2, 4), s); // vx
                acadoVariables.y[i * NY + 8] = lerp(traj(row1, 5), traj(row2, 5), s); // vy
                acadoVariables.y[i * NY + 9] = lerp(traj(row1, 6), traj(row2, 6), s); // vz
                acadoVariables.y[i * NY + 10] = 0; // horiz landmark projection
                acadoVariables.y[i * NY + 11] = 0; // vert landmark projection
                acadoVariables.y[i * NY + 12] = 9.8; // Thrust
                acadoVariables.y[i * NY + 13] = 0; // wx
                acadoVariables.y[i * NY + 14] = 0; // wy
                acadoVariables.y[i * NY + 15] = 0; // wz
            }
            row_try = (int)floor((horizon_start + N*dt_mpc)/dt_traj);
            row1 = min((int)traj.rows()-1, row_try);
            row2 = min((int)traj.rows()-1, row1 + 1);
            s = ((horizon_start + N*dt_mpc)/dt_traj) - row_try; // interpolation parameter
            acadoVariables.yN[0] = lerp(traj(row1, 1), traj(row2, 1), s); // x
            acadoVariables.yN[1] = lerp(traj(row1, 2), traj(row2, 2), s); // y
            acadoVariables.yN[2] = lerp(traj(row1, 3), traj(row2, 3), s); // z
            acadoVariables.yN[3] = lerp(traj(row1, 7), traj(row2, 7), s); // qw
            acadoVariables.yN[4] = lerp(traj(row1, 8), traj(row2, 8), s); // qx
            acadoVariables.yN[5] = lerp(traj(row1, 9), traj(row2, 9), s); // qy
            acadoVariables.yN[6] = lerp(traj(row1,10), traj(row2,10), s); // qz
            acadoVariables.yN[7] = lerp(traj(row1, 4), traj(row2, 4), s); // vx
            acadoVariables.yN[8] = lerp(traj(row1, 5), traj(row2, 5), s); // vy
            acadoVariables.yN[9] = lerp(traj(row1, 6), traj(row2, 6), s); // vz
            acadoVariables.yN[10] = 0; // horiz landmark projection
            acadoVariables.yN[11] = 0; // vert landmark projection
        }

        // Initialize acado weight matrices
        void initAcadoWeights(bool use_percep_cost) {
            float p_cost = 0.0f;
            if (use_percep_cost) p_cost = 100.0f;
            
            int blk_size = NY*NY; // block size
            float xCostExp = 1.0f;  // state cost scaling
            float uCostExp = 1.0f;  // input cost scaling
            float xExp, uExp;
            for(int i = 0; i < N; ++i){
                xExp = exp(-((float)i/(float)N) * xCostExp);
                uExp = exp(-((float)i/(float)N) * uCostExp);
                acadoVariables.W[i*blk_size] = 500.0f * xExp; // x gain
                acadoVariables.W[i*blk_size + NY + 1] = 500.0f * xExp; // y gain
                acadoVariables.W[i*blk_size + NY*2 + 2] = 300.0f * xExp; // z gain
                acadoVariables.W[i*blk_size + NY*3 + 3] = 100.0f * xExp; // qw gain
                acadoVariables.W[i*blk_size + NY*4 + 4] = 100.0f * xExp; // qx gain
                acadoVariables.W[i*blk_size + NY*5 + 5] = 100.0f * xExp; // qy gain
                acadoVariables.W[i*blk_size + NY*6 + 6] = 100.0f * xExp; // qz gain
                acadoVariables.W[i*blk_size + NY*7 + 7] = 5.0f * xExp; // v_x gain
                acadoVariables.W[i*blk_size + NY*8 + 8] = 5.0f * xExp; // v_y gain
                acadoVariables.W[i*blk_size + NY*9 + 9] = 5.0f * xExp; // v_z gain
                acadoVariables.W[i*blk_size + NY*10 + 10] = p_cost * xExp; // horiz perception cost
                acadoVariables.W[i*blk_size + NY*11 + 11] = p_cost * xExp; // vert perception cost
                acadoVariables.W[i*blk_size + NY*12 + 12] = 1.0f * uExp; // T gain
                acadoVariables.W[i*blk_size + NY*13 + 13] = 15.0f * uExp; // w_x gain
                acadoVariables.W[i*blk_size + NY*14 + 14] = 15.0f * uExp; // w_y gain
                acadoVariables.W[i*blk_size + NY*15 + 15] = 1.0f * uExp; // w_z gain
            }
            
            xExp = exp(-xCostExp);
            //xExp = 1.0f;
            acadoVariables.WN[0] = 500.0f * xExp; // x gain
            acadoVariables.WN[NYN+1] = 500.0f * xExp; // y gain
            acadoVariables.WN[NYN*2 + 2] = 300.0f; // z gain (not xExp bc want landing height to be more exact)
            acadoVariables.WN[NYN*3 + 3] = 100.0f * xExp; // qw gain
            acadoVariables.WN[NYN*4 + 4] = 100.0f * xExp; // qx gain
            acadoVariables.WN[NYN*5 + 5] = 100.0f * xExp; // qy gain
            acadoVariables.WN[NYN*6 + 6] = 100.0f * xExp; // qz gain
            acadoVariables.WN[NYN*7 + 7] = 5.0f * xExp; // vz gain
            acadoVariables.WN[NYN*8 + 8] = 5.0f * xExp; // vz gain
            acadoVariables.WN[NYN*9 + 9] = 5.0f * xExp; // vz gain
            acadoVariables.WN[NYN*10 + 10] = p_cost * xExp; // horiz perception cost
            acadoVariables.WN[NYN*11 + 11] = p_cost * xExp; // vert perception cost
        }
        
        // Update acado online data
        void updateAcadoOnlineData(bool use_percep_cost) {
            if (use_percep_cost){
                double tagYaw = curYaw + yawDiff;
                for (int i = 0; i < N + 1; ++i){
                    acadoVariables.od[i * NOD + 0] = tagPos(0) - 0.108*cos(tagYaw); // landmark x position
                    acadoVariables.od[i * NOD + 1] = tagPos(1) - 0.108*sin(tagYaw); // landmark y position
                    acadoVariables.od[i * NOD + 2] = tagPos(2); // landmark z position
                    acadoVariables.od[i * NOD + 3] = 0.0; // NOT USED translation body to cam x   
                    acadoVariables.od[i * NOD + 4] = 0.0; // NOT USED translation body to cam y
                    acadoVariables.od[i * NOD + 5] = 0.0; // NOT USED translation body to cam z
                    acadoVariables.od[i * NOD + 6] = 0.0; // quat w body to cam
                    acadoVariables.od[i * NOD + 7] = 0.0; // quat x body to cam
                    acadoVariables.od[i * NOD + 8] = 1.0; // quat y body to cam
                    acadoVariables.od[i * NOD + 9] = 0.0; // quat z body to cam
                }
            }
        }
        
        // Initialize acado control inputs
        void initAcadoControls() {
            for (int i = 0; i < N; ++i){ // iterate over time steps within horizon
                acadoVariables.u[i * NU + 0] = 9.8; // Thrust (mass-normalized)
                acadoVariables.u[i * NU + 1] = 0.0; // w_x
                acadoVariables.u[i * NU + 2] = 0.0; // w_y
                acadoVariables.u[i * NU + 3] = 0.0; // w_z
            }
        }
        
        // Load-in a pre-made search trajectory
        Eigen::MatrixXd loadSearchTraj(const std::string& filename) {
            
            // find file
            std::ifstream file(filename);
            if (!file){
                ROS_ERROR_STREAM("The file named \'" << filename << "\' could not be opened. Shutting down node ...");
                ros::shutdown();
            }

            // store file data in a vector of vectors
            std::vector<std::vector<double>> data;
            std::string line;
            while (std::getline(file, line)) {
                std::vector<double> row;
                std::stringstream ss(line);
                std::string token;
                while (std::getline(ss, token, ',')) {
                    row.push_back(std::stod(token));
                }
                data.push_back(row);
            }

            // transfer the data to an Eigen matrix
            const int num_rows = data.size();
            const int num_cols = data[0].size();
            Eigen::MatrixXd eval = Eigen::MatrixXd(num_rows, num_cols);
            for (int i = 0; i < num_rows; ++i) {
                for (int j = 0; j < num_cols; ++j) {
                    eval(i, j) = data[i][j];
                }
            }
                        
            // ------------------ Load ROS Trajectory message ------------------ //
            eigen2PathMsg(eval, mpcTotalRef);
            
            return eval;
        }
        
        // Generate a landing trajectory once the landing pad (AprilTag) is spotted
        Eigen::MatrixXd genLandingTraj(double totalTime, bool doFOV){
            // ------------------ Parameters ------------------ //
            int order = 8;          // order of piecewise polynomials (must be >= 4 for min snap) (works well when order >= numFOVtimes)
            int num_intervals = 1;  // number of time intervals (must be >= 1) (setting to 1 is fine if not using keyframes, and only using FOV constraints)
            int num_seg_times = num_intervals + 1;
            double ros_start = ros::Time::now().toSec();
            vector<double> seg_times = MinSnapTraj::linspace(ros_start, totalTime + ros_start, num_seg_times); // time points for the polynomial segments
            
            // ------------------ Prepare the solver (boundary conditions etc.) ------------------ //
            MinSnapData data = prepLandingSolver(doFOV);

            // ------------------ Solve the trajectory ------------------ //
            MinSnapTraj prob(order, seg_times, data.bounds, data.keyframes, data.fov_data);  // create object

            chrono::steady_clock::time_point tic, toc;                // time the solver
            tic = chrono::steady_clock::now();

            MinSnapTraj::TrajSolution sol = prob.solveTrajectory();   // solve

            toc = chrono::steady_clock::now();
            showDuration("Init solve time: ", toc - tic);

            // ------------------ Save trajectory in Eigen matrix ------------------ //
            int eval_freq = 10; // evaluation data points per second (approximate, due to linspace spacing)
            int num_eval_times = eval_freq*totalTime + 1;
            vector<double> eval_times = MinSnapTraj::linspace(ros_start, totalTime + ros_start, num_eval_times); // time points to evaluate the trajectory
            Eigen::MatrixXd eval = prob.QuaternionTraj(sol.coeffs, eval_times);
            
            // Align first trajectory quaternion with current quadcopter quaternion
            Eigen::Quaterniond first_quat(eval(0,7), eval(0,8), eval(0,9), eval(0,10));
            Eigen::Quaterniond curQuat = msg2quaternion(curPose.orientation);
            if (quatsDiffSign(first_quat, curQuat)) flipQuatSigns(eval);
            
            // Concatonate extra rows if totalTime <= mpc_time_horizon,
            //   allowing acado mpc to still function
            if (totalTime <= mpc_time_horizon){
                double dt_eval = eval(1,0) - eval(0,0);
                double start_time = eval_times.back() + dt_eval; // start one sample past last time
                int num_extra_times = (int)ceil((mpc_time_horizon - totalTime)/dt_eval); // end one sample past mpc_time_horizon
                double final_time = start_time + dt_eval*(num_extra_times - 1);
                
                vector<double> tspan_extra = MinSnapTraj::linspace(start_time, final_time, num_extra_times);
                Eigen::MatrixXd eval_extra(tspan_extra.size(), eval.cols());  // matrix for loading extra rows
                for (size_t i = 0; i < tspan_extra.size(); i++){
                    eval_extra.row(i) = eval.row(eval.rows()-1);  // copy the last row of eval
                    eval_extra(i, 0) = tspan_extra[i];            // update the time
                }
                Eigen::MatrixXd eval_total(eval.rows() + eval_extra.rows(), eval.cols());
                eval_total << eval, eval_extra;  // concatenate both matrices
                eval = eval_total;
            }
            
            return eval;
        }
        
        // Prepare boundary conditions and other data for trajectory generation
        MinSnapData prepLandingSolver(bool doFOV) {
            // ------------------ Boundary conditions ------------------ //
            MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
            geometry_msgs::Point curPos = curPose.position;
            double tagYaw = curYaw + yawDiff;
            // pos/yaw:  x          y          z          yaw
            p0_bounds << curPos.x,  curPos.y,  curPos.z,  curYaw, // initial
                         tagPos(0) - 0.108*cos(tagYaw), tagPos(1) - 0.108*sin(tagYaw), tagPos(2), tagYaw; // final

            MinSnapTraj::Matrix24d p1_bounds; // p=1 means 1st derivative
            // velocity: vx   vy   vz   vyaw
            p1_bounds << curVel.x, curVel.y, curVel.z, curYawVel, // initial
                        0.0, 0.0, 0.0, 0.0; // final

            MinSnapTraj::Matrix24d p2_bounds; // p=2 means 2nd derivative
            // accel:    ax         ay      az      ayaw
            p2_bounds << curLinAcc.x, curLinAcc.y, curLinAcc.z-9.81, 0.0, // initial
                        0.0, 0.0, 1.0, 0.0; // final

            MinSnapData data;
            data.bounds.push_back(p0_bounds);
            data.bounds.push_back(p1_bounds);
            data.bounds.push_back(p2_bounds);

            // ------------------ Keyframe/Waypoints ------------------ //
            data.keyframes = {}; // no keyframes

            // ------------------ FOV data ------------------ //
            data.fov_data.do_fov = doFOV;
            data.fov_data.l = vector<double> {tagPos(0) - 0.108*cos(tagYaw), tagPos(1) - 0.108*sin(tagYaw), tagPos(2)};  // 3D landmark to keep in FOV (i.e. AprilTag center)
            data.fov_data.alpha_x = 0.5;    // half of vertical fov (radians), orig: 0.59
            data.fov_data.alpha_y = 0.7;    // orig: 0.79
            
            return data;
        }
        
        // Check if a FOV-constrained trajectory if feasible
        bool checkFOVFeasible(bool doFOV) {
            if (doFOV) {
                MinSnapData data = prepLandingSolver(doFOV);
                vector<double> l = data.fov_data.l;
                double alpha_x = data.fov_data.alpha_x;
                double alpha_y = data.fov_data.alpha_y; 
                
                double x0 = data.bounds[0](0,0); double xm = data.bounds[0](1,0);
                double y0 = data.bounds[0](0,1); double ym = data.bounds[0](1,1);
                double z0 = data.bounds[0](0,2); double zm = data.bounds[0](1,2);           
                double theta_x0 = atan2(l[0] - x0, z0 - l[2]); // initial theta for ax az constraint
                double theta_xm = atan2(l[0] - xm, zm - l[2]); // final theta for ax az constraint
                double theta_y0 = atan2(l[1] - y0, z0 - l[2]); // initial theta for ay az constraint
                double theta_ym = atan2(l[1] - ym, zm - l[2]); // final theta for ay az constraint
                
                // Check FOV constraint satisfied at initial and final times
                bool theta_x_BC_FOV = ( (alpha_x - abs(theta_x0) > 1e-3) && (alpha_x - abs(theta_xm) > 1e-3) );
                bool theta_y_BC_FOV = ( (alpha_y - abs(theta_y0) > 1e-3) && (alpha_y - abs(theta_ym) > 1e-3) );
                
                return (theta_x_BC_FOV && theta_y_BC_FOV);
            } else {
                return true;
            }
        }
        
        void odomMsgSelector(const OdomData& data) {
            assert(data.source == APRIL_TAG || data.source == POSE_BY_EKF2 || data.source == VEL_BY_EKF2); // make sure data.source is filled
            // TODO for Apriltag as state est
//             if (tag_visible) { // Use AprilTag for position estimate if available
//                 if(data.source == APRIL_TAG){
//                     curPose.position = data.pose.position; // position of quad, using apriltag estimate
//                 } else if (data_has_pose){
//                     curPose.orientation = data.pose.orientation;
//                     geometry_msgs::Quaternion curQuat = curPose.orientation;
//                     Eigen::Quaterniond quadQuat(curQuat.w, curQuat.x, curQuat.y, curQuat.z);  // convert quadQuat from message to Eigen::Quateniond
//                     curYaw = quaternion2yaw(quadQuat);
//                 }
//             } else if (data_has_pose) {
//                 curPose = data.pose;
//                 geometry_msgs::Quaternion curQuat = curPose.orientation;
//                 Eigen::Quaterniond quadQuat(curQuat.w, curQuat.x, curQuat.y, curQuat.z);  // convert quadQuat from message to Eigen::Quateniond
//                 curYaw = quaternion2yaw(quadQuat);
//             }
//             if (data_has_vel){
//                 curVel = data.velocity;
//             }
            if (data.source == POSE_BY_EKF2) {
                curPose = data.pose;
                //negativeQuat(curPose.orientation);   // to avoid acado MPC yaw issues
                Eigen::Quaterniond quadQuat = msg2quaternion(curPose.orientation);  // convert quadQuat from message to Eigen::Quateniond
                curYaw = quaternion2yaw(quadQuat);
            }
            if (data.source == VEL_BY_EKF2){
                curVel = data.velocity;
            }
        }
        
        // Convert mass-normalized thrust to PX4 normalized thrust
        float thrustMap(double th){
            // input: mass-normalized thrust (m/s^2)
            // output: PX4-normalized thrust ([0,1])
            double a, b, c, d;
            
            if (sim_enable){
                a = 0.00013850414341400538;
                b = -0.005408755617324549;
                c = 0.11336157680888627;
                d = -0.0022807568577082674;
            } else {
                ROS_ERROR_STREAM("No thrust map for real quadcopter.");
                ros::shutdown();
            }
            
            double norm_th = a*th*th*th + b*th*th + c*th + d;
            if (th == 0.0) {
                return 0.0;
            } else {
                return (float)std::max(0.0, std::min(1.0, norm_th ));
            }
        }
        
        // Method to convert Eigen matrix data into a ROS nav_msgs/Path message
        void eigen2PathMsg(const Eigen::Matrix<double, Eigen::Dynamic, 11>& data, nav_msgs::Path& path) {
            int num_times = data.rows();
            path.header.stamp = ros::Time::now(); // Set the timestamp of the Path message
            path.header.frame_id = "map";
            
            // Clear the poses array in the Path message
            path.poses.clear();
            // Create a PoseStamped message
            geometry_msgs::PoseStamped pose_stamped;
            
            // Loop through each row of data
            for (int i = 0; i < num_times; i++) {
                pose_stamped.header.stamp = ros::Time(data(i, 0));
                pose_stamped.pose.position.x = data(i, 1);
                pose_stamped.pose.position.y = data(i, 2);
                pose_stamped.pose.position.z = data(i, 3);
                pose_stamped.pose.orientation.w = data(i, 7);
                pose_stamped.pose.orientation.x = data(i, 8);
                pose_stamped.pose.orientation.y = data(i, 9);
                pose_stamped.pose.orientation.z = data(i, 10);

                // Add the PoseStamped message to the poses array in the Path message
                path.poses.push_back(pose_stamped);
            }
        }
        
        // Method to convert ACADO array (x or y) into a ROS nav_msgs/Path message
        void acadoArr2PathMsg(const float* A, int stride, nav_msgs::Path& path) {
            path.header.stamp = ros::Time::now(); // Set the timestamp of the Path message
            path.header.frame_id = "map";
            path.poses.clear();
            
            // Create a PoseStamped message
            geometry_msgs::PoseStamped pose_stamped;
            
            for (int i = 0; i < N; ++i){ // iterate over time steps within horizon
                pose_stamped.pose.position.x = A[i * stride + 0]; // x
                pose_stamped.pose.position.y = A[i * stride + 1]; // y
                pose_stamped.pose.position.z = A[i * stride + 2]; // z
                pose_stamped.pose.orientation.w = A[i * stride + 3]; // qw
                pose_stamped.pose.orientation.x = A[i * stride + 4]; // qx
                pose_stamped.pose.orientation.y = A[i * stride + 5]; // qy
                pose_stamped.pose.orientation.z = A[i * stride + 6]; // qz
                
                // Add the PoseStamped message to the poses array in the Path message
                path.poses.push_back(pose_stamped);
            }
        }
        
        // Method to append pose data to existing mav_msgs/Path message
        void addPoseToPath(const geometry_msgs::Pose& pose, nav_msgs::Path& path) {
            // Create pose stamped message
            geometry_msgs::PoseStamped pose_stamped;
            pose_stamped.header.stamp = ros::Time::now();
            pose_stamped.pose = pose;
            
            // Append the pose to the path
            path.header.frame_id = "map";
            path.poses.push_back(pose_stamped);
        }
        
        // Convert quaternion to yaw, for ZYX euler order
        double quaternion2yaw(const Eigen::Quaterniond& quat) {
            vector<double> q{quat.coeffs().w(), quat.coeffs().x(), quat.coeffs().y(), quat.coeffs().z()};
            double yaw = atan2(2*(q[0]*q[3] + q[1]*q[2]), pow(q[0],2) + pow(q[1],2) - pow(q[2],2) - pow(q[3],2));
            return yaw;
        }
        
        Eigen::Quaterniond msg2quaternion(const geometry_msgs::Quaternion& q) {
            return Eigen::Quaterniond{q.w, q.x, q.y, q.z};
        }
        
        // save trajectory data for debugging
        void saveTraj(string name, Eigen::MatrixXd eval){
            string file_path = __FILE__;
            string dir_path = file_path.substr(0, file_path.rfind("mpc.h")); // get directory containing source file
            dir_path = dir_path + "../extra/traj/";
            string data_dir = "gen_landing_traj/";
            ghc::filesystem::create_directories(dir_path + data_dir); // create new directory for data
            string stamp = getTimestamp();
            string filename = name + stamp + ".csv";
            saveMatrix(dir_path + data_dir + filename, eval);
            cout << "Trajectory data saved at " << endl << dir_path + data_dir << endl << endl;
        }
        
        // store trajectory in file
        void saveMatrix(string fileName, const Eigen::MatrixXd& matrix) {
            //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
            const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");
        
            ofstream file(fileName);
            if (file.is_open()){
                file << matrix.format(CSVFormat);
                file.close();
            }
        }
        
        // return string time stamp in format mmdd_HHMMSS.mmm
        std::string getTimestamp() {
            auto now = std::chrono::system_clock::now();
            auto time = std::chrono::system_clock::to_time_t(now);
            auto subSec = now.time_since_epoch() % std::chrono::seconds(1);
            std::stringstream ss;
            ss << std::put_time(std::localtime(&time), "%m-%d_%H-%M-%S");
            ss << '.' << std::setfill('0') << std::setw(3) << std::chrono::duration_cast<std::chrono::milliseconds>(subSec).count();
            return ss.str();
        }
        
        // Linear interpolation between a and b
        float lerp(double a, double b, double t){
            assert(t >= 0.0);
            assert(t <= 1.0);
            return (float)(a + t * (b - a));
        }
        
        // Flip the signs of all quaternions in a trajectory
        void flipQuatSigns(Eigen::MatrixXd& M) {
            M.block(0,7,M.rows(),4) *= -1; // start at (0,7), block of size (M.rows(),4)
        }
        
        // Check if two nearby quaternions (e.g. same rotations) use opposite signage
        // E.g. quatsDiffSign(q1, q2) is false if
        //      q1 = (0, -0.707, 0, 0.707), q2 = (0, 0.707, 0, -0.707)
        bool quatsDiffSign(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) {
            auto coeffs = q1.coeffs() + q2.coeffs();
            double abs_sum = abs(coeffs.w()) + abs(coeffs.x()) + abs(coeffs.y()) + abs(coeffs.z());
            double eps = 0.5;
            return abs_sum < eps;
        }
        
        // Print message and chrono duration in ms 
        void showDuration(string message, chrono::duration<double> t_diff){
            chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t_diff);
            cout << message << time_used.count()*1000. << " ms" << endl;
        }
};
