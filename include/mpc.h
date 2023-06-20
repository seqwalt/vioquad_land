#ifndef MPC_H_INCLUDED
#define MPC_H_INCLUDED

#include <ros/ros.h>
#include <mavros/mavros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/Imu.h>
#include <apriltag_ros/AprilTagDetectionArray.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/Thrust.h>
#include <mavros_msgs/AttitudeTarget.h>

#include <nav_msgs/Path.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <random>

#include <std_srvs/Trigger.h>
#include "quad_control/InitSetpoint.h" // custom service
#include "quad_control/FlatOutputs.h"  // custorm message

// for min snap traj gen
#include "min_snap_traj.h"
#include <iomanip>
#include <chrono>

// for loading a trajectory in

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
        // Some convenient public ACADO definitions
        ros::Publisher mpc_pub;
        ros::Publisher ref_total_pub;
        ros::Publisher ref_curr_pub;
        ros::Publisher pred_pub;
        ros::Publisher gt_pub;
        ros::Publisher tag_pub;
        
        ros::Subscriber pose_sub;
        ros::Subscriber vel_sub;
        ros::Subscriber imu_sub;
        ros::Subscriber apriltag_sub;
        ros::Timer mpc_timer;
        ros::Timer path_timer;
        
        quad_control::FlatOutputs ref;
        unsigned int iter = 0;

        //MPC(const std::string file){ // constructor
        MPC(const std::string& searchTrajFileName){ // constructor
            // Clear solver memory.
            memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
            memset(&acadoVariables, 0, sizeof( acadoVariables ));

            // Initialize the solver
            acado_initializeSolver();

            // Initialize acadoVariables for horizontal search path
            traj = loadSearchTraj(searchTrajFileName);
            init_acadoVariables();
            
            if( VERBOSE ) acado_printHeader();
            ROS_INFO_STREAM("MPC solver initialized.");

            first_mpc_call = true;
            firstTime_tagVisible = true;
            near_landing_pad = false;
            second_traj_solve = true;
            
            // Prepare first step
            acado_preparationStep();
        }

        bool initRefCallback(quad_control::InitSetpoint::Request &req, quad_control::InitSetpoint::Response &res){
            res.success = true;
            res.position.x = traj(0,1);
            res.position.y = traj(0,2);
            res.position.z = traj(0,3);
            res.yaw = quaternion2yaw(Eigen::Quaterniond(traj(0,7), traj(0,8), traj(0,9), traj(0,10)));
            ROS_INFO_STREAM("Initial MPC pose sent.");
            return true;
        }

        bool streamMpcCallback(std_srvs::Trigger::Request &req, std_srvs::Trigger::Response &res){
            mpc_timer.start();  // start publishing the MPC control values
            if(mpc_timer.hasStarted()){
                res.success = true;
                //ROS_INFO_STREAM("Publishing MPC control inputs.");
            }
            return true;
        }

        void viewPathCallback(const ros::TimerEvent &event){
            // Publish reference trajectory for visualization
            ref_total_pub.publish(mpcTotalRef);
            // Publish estimated AprilTag pose in map/world frame
            tag_pub.publish(tagPose);
            
            if(!first_mpc_call){
                // Publish ground truth path for visualization
                addPoseToPath(curPose, mpcGT);
                gt_pub.publish(mpcGT);
                
                // Publish reference for current mpc iteration
                acadoArr2PathMsg(acadoVariables.y, NY, mpcCurrRef);
                ref_curr_pub.publish(mpcCurrRef);
                
                // Publish state prediction for current mpc iteration
                acadoArr2PathMsg(acadoVariables.x, NX, mpcPred);
                pred_pub.publish(mpcPred);
            }
        }
        
        void mpcCallback(const ros::TimerEvent &event){
            if (!near_landing_pad){
                // Not close enough to landing pad to turn off motors
                //acado_tic( &t );
                // ref.reached_goal = false;

                // Current quadcopter state
                static int callback_iter = 0;

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

                if (first_mpc_call){
                    first_mpc_call = false;
                    start_time = ros::Time::now().toSec();
                    duration = 0.0;
                } else {
                    duration = min(ros::Time::now().toSec() - start_time, traj(traj_num_times-1,0));
                }
                
                int row1, row2;
                int start_row = (int)floor(duration/traj_time_step);
                double s = (duration - traj(start_row,0))/traj_time_step; // interpolation parameter
                for (int i = 0; i < N; ++i){ // iterate over time steps within horizon
                    // Update the references
                    row1 = min(traj_num_times-1, start_row + i*numIntersampleTimes);
                    row2 = min(traj_num_times-1, start_row + i*numIntersampleTimes + 1);
                    acadoVariables.y[i * NY + 0] = (float)lerp(traj(row1, 1), traj(row2, 1), s); // x
                    acadoVariables.y[i * NY + 1] = (float)lerp(traj(row1, 2), traj(row2, 2), s); // y
                    acadoVariables.y[i * NY + 2] = (float)lerp(traj(row1, 3), traj(row2, 3), s); // z
                    acadoVariables.y[i * NY + 3] = (float)lerp(traj(row1, 7), traj(row2, 7), s); // qw
                    acadoVariables.y[i * NY + 4] = (float)lerp(traj(row1, 8), traj(row2, 8), s); // qx
                    acadoVariables.y[i * NY + 5] = (float)lerp(traj(row1, 9), traj(row2, 9), s); // qy
                    acadoVariables.y[i * NY + 6] = (float)lerp(traj(row1,10), traj(row2,10), s); // qz 
                    acadoVariables.y[i * NY + 7] = (float)lerp(traj(row1, 4), traj(row2, 4), s); // vx
                    acadoVariables.y[i * NY + 8] = (float)lerp(traj(row1, 5), traj(row2, 5), s); // vy
                    acadoVariables.y[i * NY + 9] = (float)lerp(traj(row1, 6), traj(row2, 6), s); // vz
                    acadoVariables.y[i * NY + 10] = 0; // horiz landmark projection
                    acadoVariables.y[i * NY + 11] = 0; // vert landmark projection
                    acadoVariables.y[i * NY + 12] = 9.8; // Thrust
                    acadoVariables.y[i * NY + 13] = 0; // wx
                    acadoVariables.y[i * NY + 14] = 0; // wy
                    acadoVariables.y[i * NY + 15] = 0; // wz
                }
                row1 = min(traj_num_times-1, start_row + N*numIntersampleTimes);
                row2 = min(traj_num_times-1, start_row + N*numIntersampleTimes + 1);
                acadoVariables.yN[0] = (float)lerp(traj(row1, 1), traj(row2, 1), s); // x
                acadoVariables.yN[1] = (float)lerp(traj(row1, 2), traj(row2, 2), s); // y
                acadoVariables.yN[2] = (float)lerp(traj(row1, 3), traj(row2, 3), s); // z
                acadoVariables.yN[3] = (float)lerp(traj(row1, 7), traj(row2, 7), s); // qw
                acadoVariables.yN[4] = (float)lerp(traj(row1, 8), traj(row2, 8), s); // qx
                acadoVariables.yN[5] = (float)lerp(traj(row1, 9), traj(row2, 9), s); // qy
                acadoVariables.yN[6] = (float)lerp(traj(row1,10), traj(row2,10), s); // qz
                acadoVariables.yN[7] = (float)lerp(traj(row1, 4), traj(row2, 4), s); // vx
                acadoVariables.yN[8] = (float)lerp(traj(row1, 5), traj(row2, 5), s); // vy
                acadoVariables.yN[9] = (float)lerp(traj(row1, 6), traj(row2, 6), s); // vz
                acadoVariables.yN[10] = 0; // horiz landmark projection
                acadoVariables.yN[11] = 0; // vert landmark projection

                callback_iter += 1;
                //std::cout << "callback iteration: " << callback_iter << std::endl;

                int iter, status;
                // MPC optimization (real-time iteration (RTI) loop)
                for(iter = 0; iter < NUM_STEPS; ++iter){
                    // Perform the feedback step
                    status = acado_feedbackStep();

                    if (status != 0){
                        std::cout << "Iteration:" << iter << ", QP problem! QP status: " << status << std::endl;
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
            }
            
            // Apply control input
            mpcInputs.header.stamp = ros::Time::now();
            mpcInputs.type_mask = 128; // Ignore attitude messages
            mpcInputs.body_rate.x = wx;
            mpcInputs.body_rate.y = wy;
            mpcInputs.body_rate.z = wz;
            mpcInputs.thrust = thrust_map(T);
            mpc_pub.publish(mpcInputs); // set attitude, body rate and thrust to mavros
            
            //acado_printDifferentialVariables();
            //acado_printControlVariables();
            //float te = acado_toc( &t );
            //cout << "mpcCallback time: " << te << "sec" << endl;
        }

        // Pass along current mav pose data, wrt world frame.
        void mavPoseCallback(const geometry_msgs::PoseStamped &msg) {
            curPose = msg.pose;
            
            geometry_msgs::Quaternion curQuat = curPose.orientation;
            Eigen::Quaterniond quadQuat(curQuat.w, curQuat.x, curQuat.y, curQuat.z);  // convert quadQuat from message to Eigen::Quateniond
            curYaw = quaternion2yaw(quadQuat);
        }

        // Pass along current mav velocity data, wrt world frame.
        void mavVelCallback(const geometry_msgs::TwistStamped &msg) {
            curVel = msg.twist.linear;
        }
        
        // Pass along current raw mav imu data,wrt body frame.
        void mavIMUCallback(const sensor_msgs::Imu &msg) {
            curImuAcc = msg.linear_acceleration; // linear acceleration in imu body frame
        }
        
        // Compute AprilTag pose in World frame (compute transformation T_WA)
        void mavAprilTagCallback(const apriltag_ros::AprilTagDetectionArray &msg) {
            // Check if tag is visible
            tagVisible = !msg.detections.empty();            
            if (tagVisible){                
                // Form T_WB (load from state estimation: mocap, VIO, etc.)
                Eigen::Isometry3d T_WB = Eigen::Isometry3d::Identity();
                Eigen::Vector3d tran_WB(curPose.position.x, curPose.position.y, curPose.position.z); // translation part
                Eigen::Quaterniond quat_WB(curPose.orientation.w, curPose.orientation.x, curPose.orientation.y, curPose.orientation.z); // [w, x, y, z] rotation part
                T_WB.translate(tran_WB);
                T_WB.rotate(quat_WB);
                
                // Form T_BC (transformation from quadcopter body frame to down-facing camera frame)
                // TODO: load T_BC info in from launch file
                Eigen::Isometry3d T_BC = Eigen::Isometry3d::Identity();
                Eigen::Vector3d tran_BC(0.108, 0.0, 0.0); // translation part
                Eigen::Matrix3d rot_BC;
                vector<double> sdf_rpy{0.0, 1.5708, 3.1415}; // roll pitch yaw (XYZ form), given by iris_downward_camera.sdf
                vector<double> apriltag_rpy{sdf_rpy[2], -sdf_rpy[0], sdf_rpy[1]}; // rotation from sdf to apriltag representaiton (https://github.com/AprilRobotics/apriltag/wiki/AprilTag-User-Guide#coordinate-system).
                rot_BC = Eigen::AngleAxisd(apriltag_rpy[2], Eigen::Vector3d::UnitZ())
                       * Eigen::AngleAxisd(apriltag_rpy[1], Eigen::Vector3d::UnitY())
                       * Eigen::AngleAxisd(apriltag_rpy[0], Eigen::Vector3d::UnitX()); // XYZ rotation http://sdformat.org/tutorials?tut=specify_pose
                T_BC.translate(tran_BC);
                T_BC.rotate(rot_BC);
                
                // Compute T_CA (transformation from down-facing camera frame to AprilTag frame)
                geometry_msgs::Pose AprilTagPose_Cam = msg.detections[0].pose.pose.pose; // pose of AprilTag in camera frame
                Eigen::Vector3d tran_CA(AprilTagPose_Cam.position.x, AprilTagPose_Cam.position.y, AprilTagPose_Cam.position.z); // translation part
                Eigen::Quaterniond quat_CA(AprilTagPose_Cam.orientation.w, AprilTagPose_Cam.orientation.x, AprilTagPose_Cam.orientation.y, AprilTagPose_Cam.orientation.z); // rotation part
                Eigen::Isometry3d T_CA = Eigen::Isometry3d::Identity();
                T_CA.translate(tran_CA);
                T_CA.rotate(quat_CA); // transformation matrix representing the AprilTag in the camera frame
                near_landing_pad = (abs(tran_CA[2]) < land_height); // check if near landing pad
                
                // Compute T_WA
                Eigen::Isometry3d T_WA = T_WB * T_BC * T_CA;
                tagPos = T_WA.translation();
                Eigen::Matrix3d tagRot = T_WA.rotation();
                Eigen::Quaterniond tagQuat = Eigen::Quaterniond(tagRot);
                tagYaw = quaternion2yaw(tagQuat);
                
                geometry_msgs::PoseStamped tagPoseTemp;
                tagPoseTemp.pose.position.x = tagPos(0);
                tagPoseTemp.pose.position.y = tagPos(1);
                tagPoseTemp.pose.position.z = tagPos(2);
                tagPoseTemp.pose.orientation.w = tagQuat.coeffs()(3);
                tagPoseTemp.pose.orientation.x = tagQuat.coeffs()(0);
                tagPoseTemp.pose.orientation.y = tagQuat.coeffs()(1);
                tagPoseTemp.pose.orientation.z = tagQuat.coeffs()(2);
                tagPoseTemp.header.frame_id = "map";
                tagPoseTemp.header.stamp = ros::Time::now();
                tagPose = tagPoseTemp;
                
                // Gerenerate landing trajectory for first time
                // TODO create landing traj after continuously seeing AprilTag for some amount of time
                if (firstTime_tagVisible){
                    // start a timer
                    land_start_time = chrono::steady_clock::now();
                    firstTime_tagVisible = false;
                    numIntersampleTimes = 10;
                    traj.resize(0,0);
                    traj = genLandingTraj(numIntersampleTimes, 3.5, false); // generate trajectory. each row has t, x, y, z, vx, vy, vz, qw, qx, qy, qz
                    init_acadoVariables();
                    first_mpc_call = true; // reset for new trajectory
                } else {
                    // TODO fix this shifting
                    //shiftTraj2AprilTag();
                    chrono::steady_clock::time_point curr_time = chrono::steady_clock::now();
                    chrono::duration<double> dur = chrono::duration_cast<chrono::duration<double>>(curr_time - land_start_time);
                    if ((double)dur.count() > 1.4 && second_traj_solve){
                        second_traj_solve = false;
                        traj = genLandingTraj(numIntersampleTimes, 3.5 - (double)dur.count(), false); // generate trajectory. each row has t, x, y, z, vx, vy, vz, qw, qx, qy, qz
                        init_acadoVariables();
                        first_mpc_call = true;
                    }
                }
            } else {
                // No AprilTags detected
            }
        }

    private:
        geometry_msgs::Pose curPose;
        double curYaw;
        geometry_msgs::Vector3 curVel;
        geometry_msgs::Vector3 curImuAcc;
        Eigen::Vector3d tagPos;
        double tagYaw;
        geometry_msgs::PoseStamped tagPose;
        mavros_msgs::AttitudeTarget mpcInputs;
        nav_msgs::Path mpcTotalRef;  // total reference
        nav_msgs::Path mpcCurrRef;  // current iteration reference
        nav_msgs::Path mpcPred;     // mpc prediction
        nav_msgs::Path mpcGT;   // ground truth
        int numIntersampleTimes;
        Eigen::MatrixXd traj;
        bool first_mpc_call;
        bool tagVisible;  // visibility of AprilTag
        bool firstTime_tagVisible;
        bool near_landing_pad;
        bool second_traj_solve;
        double start_time, duration; // in seconds
        double traj_time_step;
        int traj_num_times;
        
        const double land_height = 0.1; // stop motors when this high (m) above landing pad
        chrono::steady_clock::time_point land_start_time;

        double T, wx, wy, wz;
        acado_timer t;

        void init_acadoVariables(){
            // generate and store trajectory
            traj_time_step = traj(1,0) - traj(0,0);
            traj_num_times = traj.rows();

            // Initialize the states TODO: fix the dependency on numIntersampleTimes
            int row;
            for (int i = 0; i < N + 1; ++i){
                row = i*numIntersampleTimes;
                acadoVariables.x[i * NX + 0] = (float)traj(row, 1); // x
                acadoVariables.x[i * NX + 1] = (float)traj(row, 2); // y
                acadoVariables.x[i * NX + 2] = (float)traj(row, 3); // z
                acadoVariables.x[i * NX + 3] = (float)traj(row, 7); // qw
                acadoVariables.x[i * NX + 4] = (float)traj(row, 8); // qx
                acadoVariables.x[i * NX + 5] = (float)traj(row, 9); // qy
                acadoVariables.x[i * NX + 6] = (float)traj(row,10); // qz
                acadoVariables.x[i * NX + 7] = (float)traj(row, 4); // vx
                acadoVariables.x[i * NX + 8] = (float)traj(row, 5); // vy
                acadoVariables.x[i * NX + 9] = (float)traj(row, 6); // vz
            }

            // Initialize the controls and references
            for (int i = 0; i < N; ++i){ // iterate over time steps within horizon
                // Initize the controls
                acadoVariables.u[i * NU + 0] = 9.8; // Thrust (mass-normalized)
                acadoVariables.u[i * NU + 1] = 0.0; // w_x
                acadoVariables.u[i * NU + 2] = 0.0; // w_y
                acadoVariables.u[i * NU + 3] = 0.0; // w_z

                // Initialize the references
                row = i*numIntersampleTimes;
                //cout << "traj time: " << traj(row, 0) << endl;
                acadoVariables.y[i * NY + 0] = (float)traj(row, 1); // x
                acadoVariables.y[i * NY + 1] = (float)traj(row, 2); // y
                acadoVariables.y[i * NY + 2] = (float)traj(row, 3); // z
                acadoVariables.y[i * NY + 3] = (float)traj(row, 7); // qw
                acadoVariables.y[i * NY + 4] = (float)traj(row, 8); // qx
                acadoVariables.y[i * NY + 5] = (float)traj(row, 9); // qy
                acadoVariables.y[i * NY + 6] = (float)traj(row,10); // qz
                acadoVariables.y[i * NY + 7] = (float)traj(row, 4); // vx
                acadoVariables.y[i * NY + 8] = (float)traj(row, 5); // vy
                acadoVariables.y[i * NY + 9] = (float)traj(row, 6); // vz
                acadoVariables.y[i * NY + 10] = 0; // horiz landmark projection
                acadoVariables.y[i * NY + 11] = 0; // vert landmark projection
                acadoVariables.y[i * NY + 12] = 9.8; // Thrust
                acadoVariables.y[i * NY + 13] = 0; // wx
                acadoVariables.y[i * NY + 14] = 0; // wy
                acadoVariables.y[i * NY + 15] = 0; // wz
            }
            // Initialize the final reference
            row = N*numIntersampleTimes;
            
            acadoVariables.yN[0] = (float)traj(row, 1); // x
            acadoVariables.yN[1] = (float)traj(row, 2); // y
            acadoVariables.yN[2] = (float)traj(row, 3); // z
            acadoVariables.yN[3] = (float)traj(row, 7); // qw
            acadoVariables.yN[4] = (float)traj(row, 8); // qx
            acadoVariables.yN[5] = (float)traj(row, 9); // qy
            acadoVariables.yN[6] = (float)traj(row,10); // qz
            acadoVariables.yN[7] = (float)traj(row, 4); // vx
            acadoVariables.yN[8] = (float)traj(row, 5); // vy
            acadoVariables.yN[9] = (float)traj(row, 6); // vz
            acadoVariables.yN[10] = 0; // horiz landmark projection
            acadoVariables.yN[11] = 0; // vert landmark projection

            /* MPC: initialize the current state feedback. */
            #if ACADO_INITIAL_STATE_FIXED
                acadoVariables.x0[0] = acadoVariables.x[0];
                acadoVariables.x0[1] = acadoVariables.x[1];
                acadoVariables.x0[2] = acadoVariables.x[2];
                acadoVariables.x0[3] = acadoVariables.x[3];
                acadoVariables.x0[4] = acadoVariables.x[4];
                acadoVariables.x0[5] = acadoVariables.x[5];
                acadoVariables.x0[6] = acadoVariables.x[6];
                acadoVariables.x0[7] = acadoVariables.x[7];
                acadoVariables.x0[8] = acadoVariables.x[8];
                acadoVariables.x0[9] = acadoVariables.x[9];
            #endif

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
                acadoVariables.W[i*blk_size + NY*10 + 10] = 0.0f * xExp; // 5 horiz perception cost
                acadoVariables.W[i*blk_size + NY*11 + 11] = 0.0f * xExp; // 5 vert perception cost
                acadoVariables.W[i*blk_size + NY*12 + 12] = 1.0f * uExp; // T gain
                acadoVariables.W[i*blk_size + NY*13 + 13] = 4.0f * uExp; // w_x gain
                acadoVariables.W[i*blk_size + NY*14 + 14] = 4.0f * uExp; // w_y gain
                acadoVariables.W[i*blk_size + NY*15 + 15] = 4.0f * uExp; // w_z gain
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
            acadoVariables.WN[NYN*10 + 10] = 0.0f * xExp; // 5 horiz perception cost
            acadoVariables.WN[NYN*11 + 11] = 0.0f * xExp; // 5 vert perception cost
            
            // Initialize online data
            for (int i = 0; i < N + 1; ++i){
                acadoVariables.od[i * NOD + 0] = 0.0; // landmark x position
                acadoVariables.od[i * NOD + 1] = 0.0; // landmark y position
                acadoVariables.od[i * NOD + 2] = -0.3; // landmark y position
                acadoVariables.od[i * NOD + 3] = 0.0; // translation body to cam x   
                acadoVariables.od[i * NOD + 4] = 0.0; // translation body to cam y
                acadoVariables.od[i * NOD + 5] = 0.0; // translation body to cam z
                acadoVariables.od[i * NOD + 6] = 0.0; // quat w body to cam
                acadoVariables.od[i * NOD + 7] = 0.7071068; // quat x body to cam
                acadoVariables.od[i * NOD + 8] = 0.7071068; // quat y body to cam
                acadoVariables.od[i * NOD + 9] = 0.0; // quat z body to cam
            }
        }

        float thrust_map(double th){
            // determined for SIMULATION iris quadcopter only
            // input: mass-normalized thrust (m/s^2)
            // output: PX4-normalized thrust ([0,1])
            double a = 0.00013850414341400538;
            double b = -0.005408755617324549;
            double c = 0.11336157680888627;
            double d = -0.0022807568577082674;
            double norm_th = a*th*th*th + b*th*th + c*th + d;
            if (th == 0.0) {
                return 0.0;
            } else {
                return (float)std::max(0.0, std::min(1.0, norm_th ));
            }
        }

        double lerp(double a, double b, double t){
            // linear interpolation between a and b
            assert(t >= 0.0);
            assert(t <= 1.0);
            return a + t * (b - a);
        }

        float lerp(float a, float b, float c, float t){
            // linear interpolation between a, b and c
            assert(t >= 0.0f);
            assert(t <= 1.0f);
            if (t <= 0.5f) return a + t * (b - a);
            else return b + t * (c - b);
        }

        double rand_num(double min, double max) {
            // Making rng static ensures that it stays the same
            // between different invocations of the function
            static std::default_random_engine rng;
            std::uniform_real_distribution<double> dist(min, max);
            return dist(rng);
        }

        void showDuration(string message, chrono::duration<double> t_diff){
            chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t_diff);
            cout << message << time_used.count()*1000. << " ms" << endl;
        }

        // Load in a pre-made search trajectory
        Eigen::MatrixXd loadSearchTraj(const std::string& filename){
            
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
            numIntersampleTimes = 10;
            Eigen2PathMsg(eval, mpcTotalRef);
            
            return eval;
        }
        
        // Generate a landing trajectory once the landing pad (AprilTag) is spotted
        Eigen::MatrixXd genLandingTraj(int numIntersampleTimes, double totalTime, bool doFOV){
            // Parameters
            int order = 8;         // order of piecewise polynomials (must be >= 4 for min snap) (works well when order >= numFOVtimes)
            int numIntervals = 1;  // number of time intervals (must be >= 1) (setting to 1 is fine if not using keyframes, and only using FOV constraints)
            double T = totalTime;  // 3.5 duration of trajectory in seconds (must be > 0.0)
            vector<double> times = MinSnapTraj::linspace(0.0, T, numIntervals + 1); // times for the polynomial segments
            
            // ------------------ Get acceleration in world frame ------------------ //
            Eigen::Vector3d imuAcc(curImuAcc.x, curImuAcc.y, curImuAcc.z);  // convert acc from message to Eigen 3D vector
            geometry_msgs::Quaternion curQuat = curPose.orientation;
            Eigen::Quaterniond quadQuat(curQuat.w, curQuat.x, curQuat.y, curQuat.z);  // convert quadQuat from message to Eigen::Quateniond
            Eigen::Vector3d accW = quadQuat * imuAcc - Eigen::Vector3d(0.0, 0.0, 9.81);  // rotate imuAcc by quadQuat to get acc wrt world frame, using Eigen overloaded "*" syntax, then account for gravity term
            
            // ------------------ Boundary conditions ------------------ //
            MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
            geometry_msgs::Point curPos = curPose.position;
            // pos/yaw:  x          y          z          yaw
            p0_bounds << curPos.x,  curPos.y,  curPos.z,  curYaw, // initial
                         tagPos(0) - 0.108*cos(tagYaw), tagPos(1) - 0.108*sin(tagYaw), tagPos(2), tagYaw; // final

            MinSnapTraj::Matrix24d p1_bounds; // p=1 means 1st derivative
            // velocity: vx   vy   vz   vyaw
            p1_bounds << curVel.x, curVel.y, curVel.z, 0.0, // initial
                        0.0, 0.0, 0.0, 0.0; // final

            MinSnapTraj::Matrix24d p2_bounds; // p=2 means 2nd derivative
            // accel:    ax         ay      az      ayaw
            p2_bounds << accW(0), accW(1), accW(2), 0.0, // initial
                        0.0, 0.0, 0.0, 0.0; // final

            MinSnapTraj::vectOfMatrix24d BC;
            BC.push_back(p0_bounds);
            BC.push_back(p1_bounds);
            BC.push_back(p2_bounds);

            // ------------------ Keyframe/Waypoints ------------------ //
            vector<MinSnapTraj::keyframe> Keyframes {}; // no keyframes

            // ------------------ FOV data ------------------ //
            MinSnapTraj::FOVdata fov_data;
            fov_data.do_fov = doFOV;
            fov_data.l = vector<double> {tagPos(0) - 0.108*cos(tagYaw), tagPos(1) - 0.108*sin(tagYaw), tagPos(2)};  // 3D landmark to keep in FOV (i.e. AprilTag center)
            fov_data.alpha_x = 0.5925;                // half of vertical fov (radians)
            fov_data.alpha_y = 0.79;                  // half of horizontal fov (radians)

            // ------------------ Solve the trajectory ------------------ //
            MinSnapTraj prob(order, times, BC, Keyframes, fov_data);  // create object

            chrono::steady_clock::time_point tic, toc;                // time the solver
            tic = chrono::steady_clock::now();

            MinSnapTraj::TrajSolution sol = prob.solveTrajectory();   // solve

            toc = chrono::steady_clock::now();
            showDuration("Init solve time: ", toc - tic);

            // ------------------ Save trajectory in Eigen matrix ------------------ //
            double d_time = 0.1; // discretization time for the mpc
            double num_times = 1.0 + (double)numIntersampleTimes*(T/d_time);

            vector<double> tspan = MinSnapTraj::linspace(0.0, T, num_times); // time points to evaluate the trajectory
            Eigen::MatrixXd eval = prob.QuaternionTraj(sol.coeffs, tspan);
            
            // ------------------ Load ROS Trajectory message ------------------ //
            Eigen2PathMsg(eval, mpcTotalRef);
            
            return eval;
        }
        
        // Shift the end of the landing trajectory to the AprilTag marker, to adjust for inaccuracy/motion of AprilTag pose estimate
        void shiftTraj2AprilTag() {
            Eigen::MatrixXd eval = traj;
            
            int num_times = eval.rows();
            for (int i = 0; i < num_times; i++) {
                // Shift position
                eval(i, 1) += tagPos(0) - eval(num_times-1, 1); // x
                eval(i, 2) += tagPos(1) - eval(num_times-1, 2); // y
                eval(i, 3) += tagPos(2) - eval(num_times-1, 3); // z
                
                // TODO: Shift yaw
                
            }
            
            traj = eval;
            Eigen2PathMsg(eval, mpcTotalRef);
        }
        
        // Method to convert Eigen matrix data into a ROS nav_msgs/Path message
        void Eigen2PathMsg(const Eigen::Matrix<double, Eigen::Dynamic, 11>& data, nav_msgs::Path& path) {
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
        double quaternion2yaw(const Eigen::Quaterniond& quat){
            vector<double> q{quat.coeffs()(3), quat.coeffs()(0), quat.coeffs()(1), quat.coeffs()(2)};
            double yaw = atan2(2*(q[0]*q[3] + q[1]*q[2]), pow(q[0],2) + pow(q[1],2) - pow(q[2],2) - pow(q[3],2));
            return yaw;
        }
};

#endif // MPC_H_INCLUDED
