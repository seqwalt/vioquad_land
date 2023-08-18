/**
 * @file min_snap_traj_test.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 * @note executable located at .../vioquad_ws/devel/.private/vioquad_land/lib/vioquad_land
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
    
// ------------------ Position and velocity boundary conditions ------------------ //
    MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
    // pos/yaw:   x    y    z   yaw
    p0_bounds << -2.4, -2.4, 2.5, 0.0, // initial
                  0.0,  0.0, 0.1, 0.0; // final
                 
    MinSnapTraj::Matrix24d p1_bounds; // p=1 means 1st derivative
    // velocity: vx   vy   vz   vyaw
    p1_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final
                 
    MinSnapTraj::Matrix24d p2_bounds; // p=2 means 2nd derivative
    // accel:    ax   ay   az   ayaw
    p2_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final
    
    MinSnapTraj::vectOfMatrix24d BC;
    BC.push_back(p0_bounds);
    BC.push_back(p1_bounds);
    BC.push_back(p2_bounds);
    
// ------------------ Keyframe/Waypoints ------------------ //
    vector<MinSnapTraj::keyframe> empty_Keyframes {}; // no keyframes
    
// ------------------ FOV data ------------------ //
    MinSnapTraj::FOVdata fov_data;
    fov_data.do_fov = true;
    fov_data.l = vector<double> {0.0,0.0,0.0};  // 3D landmark to keep in FOV
    fov_data.alpha_x = M_PI/4;                  // half of horizontal fov (radians)
    fov_data.alpha_y = M_PI/4;                  // half of vertical fov (radians)
    
// ------------------ Initial solve of the FOV constrained trajectory ------------------ //
    
    // Parameters
    int order = 8;         // order of piecewise polynomials (must be >= 4 for min snap) (works well when order >= numFOVtimes)
    int numIntervals = 1;  // number of time intervals (must be >= 1) (setting to 1 is fine if not using keyframes, and only using FOV constraints)
    double T = 5.0;        // duration of trajectory in seconds (must be > 0.0)
    vector<double> times = MinSnapTraj::linspace(0.0, T, numIntervals + 1); // times for the polynomial segments
    
    chrono::steady_clock::time_point tic, toc;                // time the solver
    tic = chrono::steady_clock::now();
    
    MinSnapTraj fov_prob(order, times, BC, empty_Keyframes, fov_data);  // create object
    MinSnapTraj::TrajSolution sol = fov_prob.solveTrajectory();   // solve
    
    toc = chrono::steady_clock::now();
    showDuration("Init solve time: ", toc - tic);
    
// ------------------ Second solve of the FOV constrained trajectory ------------------ //
    
    MinSnapTraj::vectOfMatrix24d bounds_new = BC;
    bounds_new[0](0,0) = -1; // set initial x position to -1
    fov_prob.updateBounds(bounds_new);  // update bounds
    
    tic = chrono::steady_clock::now();
    
    sol = fov_prob.solveTrajectory();   // solve
    
    toc = chrono::steady_clock::now();
    showDuration("Second solve time: ", toc - tic);
    
    fov_prob.deleteQpData();
    
// ------------------ Analyze/save trajectory ------------------ //
    double time_step = 0.1;
    int numIntersampleTimes = 10;
    int num_times = (int)(1.0 + (double)numIntersampleTimes*(T/time_step));
    vector<double> tspan = MinSnapTraj::linspace(0.0, T, num_times); // time points to evaluate the trajectory
    //vector<double> tspan = MinSnapTraj::linspace(0.0, T, 50); // time points to evaluate the trajectory
    //Eigen::MatrixXd eval_flat = prob.FlatOutputTraj(sol.coeffs, tspan);
    Eigen::MatrixXd eval_quat = fov_prob.QuaternionTraj(sol.coeffs, tspan);
    
    // Print trajectory data
    //   For more advanced printing: https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    //cout << "---- Trajectory: t, x, y, z, yaw, vx, vy, vz, vyaw, ax, ay, az, ayaw, sx, sy, sz, syaw ----" << endl;
    //cout << eval << endl;
    
    // Save trajectory in a file within dir_path directory (i.e. extra/tests folder)
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("min_snap_traj_test.cpp")); // get directory containing source file
    string data_dir = "test_data/";
    ghc::filesystem::create_directories(dir_path + data_dir); // create new directory for data
    
    string traj_file = "test_traj.csv";
    //saveMatrix(dir_path + data_dir + traj_file, eval_flat);
    saveMatrix(dir_path + data_dir + traj_file, eval_quat);
    
    //string key_file = "pos_keyframes.csv";
    //saveMatrix(dir_path + data_dir + key_file, pos);
    
    cout << "Trajectory data saved at " << endl << dir_path + data_dir << endl << endl;

    return 0;
}


