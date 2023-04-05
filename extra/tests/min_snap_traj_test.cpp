/**
 * @file min_snap_traj_test.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 * @note executable located at .../vioquad_ws/devel/.private/quad_control/lib/quad_control
 */

#include "min_snap_traj.h"
#include <iomanip>

using namespace std;

// store trajectory in file
void saveTraj(string fileName, Eigen::MatrixXd matrix){
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");
 
    ofstream file(fileName);
    if (file.is_open()){
        file << matrix.format(CSVFormat);
        file.close();
    }
}

int main(int argc, char **argv)
{
    // Parameters
    int order = 4;         // order of piecewise polynomials (must be >= 4 for min snap)
    int numIntervals = 1;  // number of time intervals (must be >= 1)
    double T = 5.0;        // duration of trajectory in seconds (must be > 0.0)

// ------------------ Position and velocity boundary conditions ------------------ //
    MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
    // pos/yaw:   x    y    z   yaw
    p0_bounds << 1.0, 2.0, 0.5, 0.0, // initial
                 0.0, 0.0, 2.5, 0.0; // final
                 
    MinSnapTraj::Matrix24d p1_bounds; // p=1 means 1st derivative
    // velocity: vx   vy   vz   vyaw
    p1_bounds << 0.0, 0.0, 0.0, 0.5, // initial
                 0.0, 0.0, 0.5, 0.0; // final
    
    MinSnapTraj::vectOfMatrix24d BC;
    BC.push_back(p0_bounds);
    BC.push_back(p1_bounds);
    
// ------------------ Keyframe/Waypoints ------------------ //
    // - Fix x, y, z and/or yaw at a desired time. Can fix just x, or just z and yaw, for example.
    // - Keyframes define the meeting point between two polynomial segments
    struct keyframe{
        double time_ratio; // {keyframe time} = T*time_ratio + {initial time} (require 0.0 < time_ratio < 1.0)
        int flat_ind;      // flat_ind = 0 for x, 1 for y, 2 for z, 3 for yaw
        double value;      // The desired value of the state corresponding to flat_ind
    };
    
    // Fix position at time_ratio 0.5
    double t1_ratio = 0.5;
    keyframe k1x {t1_ratio, 0, 2.0}; // time_ratio, flat_ind, value
    keyframe k1y {t1_ratio, 1, 3.0};
    keyframe k1z {t1_ratio, 2, 1.0};
    vector<keyframe> Keyframes {k1x, k1y, k1z}; // store the keyframe
    
// ------------------ Solve the trajectory ------------------ //
    MinSnapTraj prob(order, numIntervals, T, BC, Keyframes); // create object
    MinSnapTraj::TrajSolution sol = prob.solveTrajectory();
    
    // Evaluate trajectory
    vector<double> tspan = prob.linspace(0.0, T, 100); // time points to evaluate the trajectory
    Eigen::MatrixXd eval = prob.polyEval(sol.coeffs, tspan);
    
    // Print trajectory
    //   For more advanced printing: https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    //cout << "---- Trajectory: x, y, z, yaw ----" << endl;
    //cout << eval << endl;
    
    // Save trajectory in a file within dir_path directory (i.e. extra/tests folder)
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("min_snap_traj_test.cpp")); // get directory containing source file
    string csv_name = "test_traj.csv";
    saveTraj(dir_path + csv_name, eval);
    cout << "Trajectory saved at " << endl << dir_path + csv_name << endl << endl;
    
    return 0;
}


