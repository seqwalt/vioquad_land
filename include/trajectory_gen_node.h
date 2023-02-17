#ifndef TRAJ_GEN_NODE
#define TRAJ_GEN_NODE

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/CommandTOL.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/State.h>
#include <mavros_msgs/Thrust.h>

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <std_srvs/Trigger.h>
#include "quad_control/InitTraj.h"     // custom service 
#include "quad_control/FlatOutputs.h"  // custorm message

using namespace std;

class TrajectoryGen{
    public:
        ros::Publisher traj_pub;
        ros::Timer traj_timer;
        
        quad_control::FlatOutputs ref;
        Eigen::MatrixXd trajMatrix;
        double z_offset = 1; // start trajectory slightly off the ground
        unsigned int iter = 0;
        unsigned int numRows;
        
        void trajCallback(const ros::TimerEvent &event){
            if (iter < numRows){
                ref.trajectoryDone = false;
                ref.position.x = trajMatrix(iter,1);
                ref.position.y = trajMatrix(iter,2);
                ref.position.z = trajMatrix(iter,3);
                ref.velocity.x = trajMatrix(iter,5);
                ref.velocity.y = trajMatrix(iter,6);
                ref.velocity.z = trajMatrix(iter,7);
                ref.acceleration.x = trajMatrix(iter,9);
                ref.acceleration.y = trajMatrix(iter,10);
                ref.acceleration.z = trajMatrix(iter,11);
                ref.yaw = trajMatrix(iter,4);
                traj_pub.publish(ref);
                iter += 1;
            } else {
                ref.trajectoryDone = true;
                ref.position.x = trajMatrix(numRows-1,1);
                ref.position.y = trajMatrix(numRows-1,2);
                ref.position.z = trajMatrix(numRows-1,3);
                ref.velocity.x = 0.0;
                ref.velocity.y = 0.0;
                ref.velocity.z = 0.0;
                ref.acceleration.x = 0.0;
                ref.acceleration.y = 0.0;
                ref.acceleration.z = 0.0;
                ref.yaw = trajMatrix(numRows-1,4);
                traj_pub.publish(ref);
            }
        }
        
        bool sendTrajCallback(std_srvs::Trigger::Request &req, std_srvs::Trigger::Response &res){
            traj_timer.start();  // start publishing the trajectory
            if(traj_timer.hasStarted()){
                res.success = true;
                ROS_INFO_STREAM("Publishing trajectory.");
            }
            return true;
        }
        
        bool initRefCallback(quad_control::InitTraj::Request &req, quad_control::InitTraj::Response &res){
            res.success = true;
            res.position.x = trajMatrix(0,1);
            res.position.y = trajMatrix(0,2);
            res.position.z = trajMatrix(0,3);
            res.yaw        = trajMatrix(0,4);
            ROS_INFO_STREAM("Initial trajectory pose sent.");
            return true;
        }
        
        void csvToEigen(const std::string& filename)
        {
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
            const int num_cols = 18;
            trajMatrix = Eigen::MatrixXd(num_rows, num_cols);
            for (int i = 0; i < num_rows; ++i) {
                for (int j = 0; j < num_cols; ++j) {
                    trajMatrix(i, j) = data[i][j];
                }
            }
            
            // start trajectory slightly off the ground
            trajMatrix.col(3) = trajMatrix.col(3).array() + z_offset; // add an offset to the z column (column 3)
            numRows = static_cast<unsigned int>(num_rows);
            ROS_INFO_STREAM("Loaded trajectory file named \'" << filename << "\'." );
        }
};

#endif














