#ifndef MPC_H_INCLUDED
#define MPC_H_INCLUDED

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/Thrust.h>
#include <mavros_msgs/AttitudeTarget.h>

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <random>

#include <std_srvs/Trigger.h>
#include "quad_control/InitSetpoint.h" // custom service 
#include "quad_control/FlatOutputs.h"  // custorm message

// TODO configure CMakeLists.txt for acado
// TODO move acado/qpoases folders/files into filesystem
#include "acado_common.h"
#include "acado_auxiliary_functions.h"

using namespace std;

// Some convenient private ACADO definitions
#define NX          ACADO_NX  // Number of differential state variables
#define NXA         ACADO_NXA // Number of algebraic variables
#define NU          ACADO_NU  // Number of control inputs
//#define NOD         ACADO_NOD // Number of online data values

#define NY          ACADO_NY  // Number of measurements/references on nodes 0..N - 1
#define NYN         ACADO_NYN // Number of measurements/references on node N

#define N           ACADO_N   // Number of intervals in the horizon
#define NUM_STEPS   10        // Number of real-time iterations
#define VERBOSE     1         // Show iterations: 1, silent: 0

// global variables used by the solver
ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

class MPC {
    public:
        // Some convenient public ACADO definitions


        ros::Publisher mpc_pub;
        ros::Subscriber pose_sub;
        ros::Subscriber vel_sub;
        ros::Timer mpc_timer;
        
        chrono::steady_clock::time_point t1, t2;
        
        quad_control::FlatOutputs ref;
        Eigen::MatrixXd trajMatrix;
        unsigned int iter = 0;
        
        //MPC(const std::string file){ // constructor
        MPC(){ // constructor
            // Clear solver memory.
            memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
            memset(&acadoVariables, 0, sizeof( acadoVariables ));
            
            // load parameters
            //TODO (low priority) paramMPC par = loadParams(file);
            
            // set acado online data (cost matrix Q, camera orientation, landing pad position estimate, etc)
            // TODO (low priority) acadoVariables.od[0] = paramMPC.diagQ[0];
            // try acadoVariables.od[0:M] = paramMPC.diagQ;
            
            // Initialize the solver
            acado_initializeSolver();

            // Initialize acadoVariables
            init_acadoVariables();
            
            if( VERBOSE ) acado_printHeader();
            ROS_INFO_STREAM("MPC solver initialized.");
            
            // Prepare first step
            acado_preparationStep();
        }
        
        bool initRefCallback(quad_control::InitSetpoint::Request &req, quad_control::InitSetpoint::Response &res){
            res.success = true;
            res.position.x = 2.0;
            res.position.y = 2.0;
            res.position.z = 1.5;
            res.yaw        = 0.0;
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
        
        void mpcCallback(const ros::TimerEvent &event){
//             acado_tic( &t );
            
            // ref.reached_goal = false;
            
            // TODO use acado_shiftStates and acado_shiftControls instead of this:
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
            
            int iter;
            // MPC optimization (real-time iteration (RTI) loop)
            for(iter = 0; iter < NUM_STEPS; ++iter){
                // Perform the feedback step
                acado_feedbackStep();
                // Shift the initialization (look at acado_common.h)
//                 acado_shiftStates(2, 0, 0);
//                 acado_shiftControls(0);
                // Prepare for the next step
                acado_preparationStep();
            }
            
            // Extract control input
            T  = acadoVariables.u[0]; // mass-normalized thrust
            wx = acadoVariables.u[1];
            wy = acadoVariables.u[2];
            wz = acadoVariables.u[3];
            
            // Apply control input
            mpcInputs.header.stamp = ros::Time::now();
            mpcInputs.type_mask = 128; // Ignore attitude messages
            mpcInputs.body_rate.x = wx;
            mpcInputs.body_rate.y = wy;
            mpcInputs.body_rate.z = wz;
            mpcInputs.thrust = thrust_map(T);
            mpc_pub.publish(mpcInputs); // set attitude, body rate and thrust to mavros
            
            acado_printDifferentialVariables();
            acado_printControlVariables();
//             real_t te = acado_toc( &t );
//             cout << "mpcCallback time: " << te << "sec" << endl;
        }
        
        // Pass along current mav pose data.
        void mavPoseCallback(const geometry_msgs::PoseStamped &msg) {
            curPose = msg.pose;
        }

        // Pass along current mav velocity data.
        void mavVelCallback(const geometry_msgs::TwistStamped &msg) {
            curVel = msg.twist.linear;
        }

    private:
        geometry_msgs::Pose curPose;
        geometry_msgs::Vector3 curVel;
        mavros_msgs::AttitudeTarget mpcInputs;
        
        double T, wx, wy, wz;
        acado_timer t;
        
        // TODO (low priority) setup yaml file reading
//         // MPC parameters from yaml file
//         struct paramMPC {
//             std::vector<double> diagQ;  // diagonal elemetns of cost matrix
//         };
//         
//         paramMPC loadParams(const std::string file){
//             
//         }
        
        void init_acadoVariables(){
            int i;
            float s;
            // Initialize the states
            for (i = 0; i < N + 1; ++i){
                s = (float)i/(float)N; // linear interp parameter in range [0,1]
                acadoVariables.x[i * NX + 0] = lerp(2.0f, 0.0f, s); // x; 2 if s==0, 0 if s==1
                acadoVariables.x[i * NX + 1] = lerp(2.0f, 0.0f, s); // y
                acadoVariables.x[i * NX + 2] = lerp(1.5f, 1.0f, s); // z
                acadoVariables.x[i * NX + 3] = 1.0; // qw
                acadoVariables.x[i * NX + 4] = 0.0; // qx
                acadoVariables.x[i * NX + 5] = 0.0; // qy
                acadoVariables.x[i * NX + 6] = 0.0; // qz
                acadoVariables.x[i * NX + 7] = 0.0; // vx
                acadoVariables.x[i * NX + 8] = 0.0; // vy
                acadoVariables.x[i * NX + 9] = 0.0; // vz
            }
            // Initialize the controls and references
            for (i = 0; i < N; ++i){
                // Initize the controls
                acadoVariables.u[i * NU + 0] = 9.8; // Thrust (mass-normalized)
                acadoVariables.u[i * NU + 1] = 0.0; // w_x
                acadoVariables.u[i * NU + 2] = 0.0; // w_y
                acadoVariables.u[i * NU + 3] = 0.0; // w_z

                // Initialize the references
                acadoVariables.y[i * NY + 0] = 0.0; // x
                acadoVariables.y[i * NY + 1] = 0.0; // y
                acadoVariables.y[i * NY + 2] = 1.0; // z
                acadoVariables.y[i * NY + 3] = 1.0; // qw
                acadoVariables.y[i * NY + 4] = 0.0; // qx
                acadoVariables.y[i * NY + 5] = 0.0; // qy
                acadoVariables.y[i * NY + 6] = 0.0; // qz
                acadoVariables.y[i * NY + 7] = 0.0; // vx
                acadoVariables.y[i * NY + 8] = 0.0; // vy
                acadoVariables.y[i * NY + 9] = 0.0; // vz
            }
            // Initialize the final reference
            acadoVariables.yN[0] = 0.0; // x
            acadoVariables.yN[1] = 0.0; // y
            acadoVariables.yN[2] = 1.0; // z
            acadoVariables.yN[3] = 1.0; // qw
            acadoVariables.yN[4] = 0.0; // qx
            acadoVariables.yN[5] = 0.0; // qy
            acadoVariables.yN[6] = 0.0; // qz
            acadoVariables.yN[7] = 0.0; // vx
            acadoVariables.yN[8] = 0.0; // vy
            acadoVariables.yN[9] = 0.0; // vz
            
            /* MPC: initialize the current state feedback. */
        #if ACADO_INITIAL_STATE_FIXED
            for (i = 0; i < NX; ++i) acadoVariables.x0[i] = acadoVariables.x[i];
        #endif
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
            return (float)std::max(0.0, std::min(1.0, norm_th ));
        }
        
        double rand_num() {
            // Making rng static ensures that it stays the same
            // between different invocations of the function
            static std::default_random_engine rng;
            std::uniform_real_distribution<double> dist(-1.0, 1.0); 
            return dist(rng); 
        }
        
        float lerp(float a, float b, float t){
            assert(t <= 1.0f);
            assert(t >= 0.0f);
            return a + t * (b - a);
        }
};

#endif // MPC_H_INCLUDED















