/**
 * @file compare_traj_test.cpp
 * @brief Sanity check comparing MinSnapTraj::FlatOutputTraj() and MinSnapTraj::QuaternionTraj() trajectories. Should contain same info.
 * @note executable located at .../vioquad_ws/devel/.private/quad_control/lib/quad_control
 */

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <iostream>
#include <string>
#include <vector>
#include "include/filesystem.hpp"

using std::string;
using std::cout;
using std::endl;

// Load csv matrix into Eigen::MatrixXd
void csvToEigen(const std::string& filename, Eigen::MatrixXd& trajMatrix) {
    // find file
    std::ifstream file(filename);
    if (!file) throw std::runtime_error{"The file could not be opened."};

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
    trajMatrix = Eigen::MatrixXd(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            trajMatrix(i, j) = data[i][j];
        }
    }
}

// Convert quaternion to yaw, for ZYX euler order
double quaternion2yaw(const Eigen::Quaterniond& quat) {
    std::vector<double> q{quat.coeffs()(3), quat.coeffs()(0), quat.coeffs()(1), quat.coeffs()(2)};
    double yaw = atan2(2*(q[0]*q[3] + q[1]*q[2]), pow(q[0],2) + pow(q[1],2) - pow(q[2],2) - pow(q[3],2));
    return yaw;
}

int main(){
    // File names to compare
    string quat_file = "quat_traj_20230711191509.csv";
    string flat_file = "flat_traj_20230711191509.csv";
    
    // Get full file path
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("compare_traj_test.cpp")); // get directory containing source file
    dir_path = dir_path + "../traj/gen_landing_traj/";
    string flat_path = dir_path + flat_file;
    string quat_path = dir_path + quat_file;
    
    // Load Eigen matrices
    Eigen::MatrixXd Mf, Mq;
    csvToEigen(quat_path, Mq); // time, x, y, z, vx, vy, vz, qw, qx, qy, qz
    csvToEigen(flat_path, Mf); // time, x, y, z, yaw, vx, vy, vz, vyaw, ax, ay, az, ayaw, jx, jy, jz, jyaw
    assert(Mf.rows() == Mq.rows() && Mf(0,0) == Mq(0,0) && Mf(Mf.rows()-1,0) == Mq(Mq.rows()-1,0));
    
    // Convert quaternion to yaw
    double max = -1.0;
    double yawDiffAbs;
    for (int i = 0; i < Mq.rows(); ++i) {
        yawDiffAbs = std::abs(quaternion2yaw(Eigen::Quaterniond(Mq(i,7), Mq(i,8), Mq(i,9), Mq(i,10))) - Mf(i,4));
        max = (yawDiffAbs > max ? yawDiffAbs : max);
    }
    cout << "Max abs yaw difference: " << max << endl;
}
