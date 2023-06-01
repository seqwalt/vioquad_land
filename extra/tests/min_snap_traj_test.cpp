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
    
// ------------------ Initial solve of the search path trajectory ------------------ //
    /*
    // keyframes
    // - Fix x, y, z and/or yaw at a desired time. Can fix just x, or just z and yaw, for example.
    // - Keyframes lie at meeting point between two polynomial segments
    // Fix at position pos at time times[time_ind]
    
    // For horizontal sweep search path
    int numHlines = 5;          // number of staight horizontal lines
    int numKeys = numHlines*2 - 2;  // number of 3D keypoints minus initial and final positions
    double xWid = 4.0;
    double yWid = 4.0;
    vector<MinSnapTraj::keyframe> Keyframes {}; // store the keyframe
    //MinSnapTraj::keyframe k1x, k1y, k2x, k2y;
//     for (int i = 1; i <= numKeys/2; i++) {
//         // only need to specify x, y since z is const.
//         k1x.t_ind = 2*i-1; k1x.flat_ind = 0; k1x.value = xWid*(double)(i%2); // time_ind, flat_ind, value
//         k1y.t_ind = 2*i-1; k1y.flat_ind = 1; k1y.value = (double)(i-1)*yWid/((double)numHlines-1.0);
//         k2x.t_ind = 2*i;   k2x.flat_ind = 0; k2x.value = xWid*(double)(i%2);
//         k2y.t_ind = 2*i;   k2y.flat_ind = 1; k2y.value = (double)i*yWid/((double)numHlines-1.0);
//         Keyframes.push_back(k1x);
//         Keyframes.push_back(k1y);
//         Keyframes.push_back(k2x);
//         Keyframes.push_back(k2y);
//     }

    MinSnapTraj::keyframe k1x {1, 0, 2.0}; // time_ind, flat_ind, value
    MinSnapTraj::keyframe k1y {1, 1, 0.0};
    MinSnapTraj::keyframe k1z {1, 2, 1.0};
    Keyframes.push_back(k1x);
    Keyframes.push_back(k1y);
    Keyframes.push_back(k1z);
    
    MinSnapTraj::keyframe k2x {2, 0, 2.0}; // time_ind, flat_ind, value
    MinSnapTraj::keyframe k2y {2, 1, 1.0};
    MinSnapTraj::keyframe k2z {2, 2, 2.0};
    Keyframes.push_back(k2x);
    Keyframes.push_back(k2y);
    Keyframes.push_back(k2z);
    
    MinSnapTraj::keyframe k3x {3, 0, 3.0}; // time_ind, flat_ind, value
    MinSnapTraj::keyframe k3y {3, 1, 2.0};
    MinSnapTraj::keyframe k3z {3, 2, 2.0};
    Keyframes.push_back(k3x);
    Keyframes.push_back(k3y);
    Keyframes.push_back(k3z);
    
    MinSnapTraj::keyframe k4x {4, 0, 3.0}; // time_ind, flat_ind, value
    MinSnapTraj::keyframe k4y {4, 1, 2.0};
    MinSnapTraj::keyframe k4z {4, 2, 2.0};
    Keyframes.push_back(k4x);
    Keyframes.push_back(k4y);
    Keyframes.push_back(k4z);

//     cout << "num keyframes: " << Keyframes.size() << endl << endl;
//     for (size_t i = 0; i < Keyframes.size(); i++){
//         cout << "Keyframe " << i << ": " << endl;
//         cout << " time ind: " << Keyframes[i].t_ind << endl;
//         cout << " flat_ind: " << Keyframes[i].flat_ind << endl;
//         cout << " value: " << Keyframes[i].value << endl << endl;
//     }
    
    // boundary conditions
    // pos/yaw:   x     y    z   yaw
    p0_bounds <<  0.0,  0.0, 1.5, 0.0, // initial
                 xWid, yWid, 1.5, 0.0; // final

    // velocity: vx   vy   vz   vyaw
    p1_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final

    // accel:    ax   ay   az   ayaw
    p2_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final

    MinSnapTraj::vectOfMatrix24d search_BC;
    search_BC.push_back(p0_bounds);
    search_BC.push_back(p1_bounds);
    //search_BC.push_back(p2_bounds);
    
    // Parameters
    int search_order = 5;         // order of piecewise polynomials (must be >= 4 for min snap) (works well when order >= numFOVtimes)
    //int search_numIntervals = 2*numHlines-1;  // number of time intervals (must be >= 1) (setting to 1 is fine if not using keyframes, and only using FOV constraints)
    int search_numIntervals = 5;
    T = 25.0;        // duration of trajectory in seconds (must be > 0.0)
    vector<double> search_times = MinSnapTraj::linspace(0.0, T, search_numIntervals + 1); // times for the polynomial segments
    
    // No FOV constraints
    MinSnapTraj::FOVdata no_fov_data;
    no_fov_data.do_fov = false;
    
    MinSnapTraj search_prob(search_order, search_times, search_BC, Keyframes, no_fov_data);  // create object
    tic = chrono::steady_clock::now();
    
    MinSnapTraj::TrajSolution search_sol = search_prob.solveTrajectory();   // solve
    
    toc = chrono::steady_clock::now();
    showDuration("Search solve time: ", toc - tic);
    */
// ------------------ Analyze/save trajectory ------------------ //
    double time_step = 0.1;
    int numIntersampleTimes = 10;
    double num_times = 1 + (double)numIntersampleTimes*(T/time_step);
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


