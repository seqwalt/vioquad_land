/**
 * @file min_snap_traj_test.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 * @note executable located at .../vioquad_ws/devel/.private/quad_control/lib/quad_control
 */

#include "min_snap_traj.h"

#include <iomanip>
#include <chrono>
#include "include/filesystem.hpp"

using namespace std;

// store trajectory in file
void saveMatrix(string fileName, Eigen::MatrixXd matrix){
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");
 
    ofstream file(fileName);
    if (file.is_open()){
        file << matrix.format(CSVFormat);
        file.close();
    }
}

void showDuration(string message, chrono::duration<double> t_diff){
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t_diff);
    cout << message << time_used.count()*1000. << " ms" << endl;
}

int main(int argc, char **argv)
{
    // Parameters
    int order = 4;         // order of piecewise polynomials (must be >= 4 for min snap)
    int numIntervals = 1;  // number of time intervals (must be >= 1)
    double T = 3.0;        // duration of trajectory in seconds (must be > 0.0)
    vector<double> times = MinSnapTraj::linspace(0.0, T, numIntervals + 1); // times for the polynomial segments
    
// ------------------ Position and velocity boundary conditions ------------------ //
    MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
    // pos/yaw:   x    y    z   yaw
    p0_bounds << -2.3, -2.3, 2.5, 0.0, // initial
                 0.0, 0.0, 0.1, 0.0; // final
                 
    MinSnapTraj::Matrix24d p1_bounds; // p=1 means 1st derivative
    // velocity: vx   vy   vz   vyaw
    p1_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final
    
    MinSnapTraj::vectOfMatrix24d BC;
    BC.push_back(p0_bounds);
    BC.push_back(p1_bounds);
    
// ------------------ Keyframe/Waypoints ------------------ //
    // - Fix x, y, z and/or yaw at a desired time. Can fix just x, or just z and yaw, for example.
    // - Keyframes lie at meeting point between two polynomial segments
    
    // Fix at position pos at time times[time_ind]
    int time_ind = 2;
    Eigen::Matrix<double, 1, 3> pos {2.0, 3.0, 1.0};
    MinSnapTraj::keyframe k1x {time_ind, 0, pos(0)}; // time_ind, flat_ind, value
    MinSnapTraj::keyframe k1y {time_ind, 1, pos(1)};
    MinSnapTraj::keyframe k1z {time_ind, 2, pos(2)};
    //vector<MinSnapTraj::keyframe> Keyframes {k1x, k1y, k1z}; // store the keyframe
    vector<MinSnapTraj::keyframe> Keyframes; // store the keyframe
    
// ------------------ Solve the trajectory ------------------ //
    MinSnapTraj prob(order, times, BC, Keyframes);          // create object
    
    chrono::steady_clock::time_point tic, toc;              // time the solver
    tic = chrono::steady_clock::now();
    
    MinSnapTraj::TrajSolution sol = prob.solveTrajectory(); // solve
    
    toc = chrono::steady_clock::now();
    showDuration("Init solve time: ", toc - tic);
    
// ------------------ Analyze/save trajectory ------------------ //
    vector<double> tspan = MinSnapTraj::linspace(0.0, T, 100); // time points to evaluate the trajectory
    Eigen::MatrixXd eval = prob.polyEval(sol.coeffs, tspan);
    
    // Print trajectory
    //   For more advanced printing: https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    //cout << "---- Trajectory: x, y, z, yaw ----" << endl;
    //cout << eval << endl;
    
    // Save trajectory in a file within dir_path directory (i.e. extra/tests folder)
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("min_snap_traj_test.cpp")); // get directory containing source file
    string data_dir = "test_data/";
    ghc::filesystem::create_directories(dir_path + data_dir); // create new directory for data
    
    string traj_file = "test_traj.csv";
    saveMatrix(dir_path + data_dir + traj_file, eval);
    
    string key_file = "pos_keyframes.csv";
    saveMatrix(dir_path + data_dir + key_file, pos);
    
    cout << "Trajectory data saved at " << endl << dir_path + data_dir << endl << endl;
    
    return 0;
}


