/**
 * @file min_snap_traj_test3.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 * @note executable located at .../vioquad_ws/devel/.private/vioquad_land/lib/vioquad_land
 */

#include "min_snap_traj.h"

#include <iomanip>
#include <chrono>
#include "include/filesystem.hpp"

using namespace std;

// store trajectory in file
void saveMatrix(string fileName, const Eigen::MatrixXd& matrix){
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");

    ofstream file(fileName);
    if (file.is_open()){
        file << matrix.format(CSVFormat);
        file.close();
    }
}

void saveKeyframes(string fileName, const vector<MinSnapTraj::keyframe>& Keyframes){
    // assumes 3D position data for keyframes
    Eigen::MatrixXd matrix(Keyframes.size()/3, 3);
    int temp_t_ind = Keyframes[0].t_ind;
    int row = 0;
    // store keyframe values in a matrix
    for (size_t i = 0; i < Keyframes.size(); i++){
        if (temp_t_ind != Keyframes[i].t_ind){
            // go to next row in the matrix
            row += 1;
        }
        temp_t_ind = Keyframes[i].t_ind;
        matrix(row, Keyframes[i].flat_ind) = Keyframes[i].value;
    }
    // save matrix as a csv file
    saveMatrix(fileName, matrix);
}

void showDuration(string message, chrono::duration<double> t_diff){
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t_diff);
    cout << message << time_used.count()*1000. << " ms" << endl;
}

int main(int argc, char **argv)
{

// ------------------ Initial solve of the search path trajectory ------------------ //

    // keyframes
    // - Fix x, y, z and/or yaw at a desired time. Can fix just x, or just z and yaw, for example.
    // - Keyframes lie at meeting point between two polynomial segments
    // Fix at position pos at time times[time_ind]

    vector<MinSnapTraj::keyframe> Keyframes {}; // store the keyframe
    MinSnapTraj::keyframe kx, ky, kz;
    vector<double> kxVals{0,1,1,1,0,-1,-1,-1,0,1,0};
    vector<double> kyVals{0,1,0,-1,-1,-1,0,1,1,1,0};
    vector<double> kzVals{0.5,1,1.7,1,1.7,1,1.7,1,1.7,1,0.5};
    for (int i = 1; i < kxVals.size()-1; ++i){
        kx.t_ind = i; ky.t_ind = i; kz.t_ind = i;
        kx.flat_ind = 0; ky.flat_ind = 1; kz.flat_ind = 2;
        kx.value = kxVals[i]; ky.value = kyVals[i]; kz.value = kzVals[i];
        Keyframes.push_back(kx);
        Keyframes.push_back(ky);
        Keyframes.push_back(kz);
    }


    // boundary conditions
    MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
    // pos/yaw:   x                 y             z         yaw
    p0_bounds << kxVals[0],     kyVals[0],     kzVals[0],     0, // initial
                 kxVals.back(), kyVals.back(), kzVals.back(), 0; // final

    MinSnapTraj::Matrix24d p1_bounds; // p=0 means 0th derivative
    // velocity: vx   vy   vz   vyaw
    p1_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final

    MinSnapTraj::Matrix24d p2_bounds; // p=0 means 0th derivative
    // accel:    ax   ay   az   ayaw
    p2_bounds << 0.0, 0.0, 0.0, 0.0, // initial
                 0.0, 0.0, 0.0, 0.0; // final

    MinSnapTraj::vectOfMatrix24d search_BC;
    search_BC.push_back(p0_bounds);
    search_BC.push_back(p1_bounds);
    search_BC.push_back(p2_bounds);

    // Parameters
    int search_order = 5;         // order of piecewise polynomials (must be >= 4 for min snap) (works well when order >= numFOVtimes)
    double T = 25.0;
    vector<double> search_times = MinSnapTraj::linspace(0.0, T, kxVals.size()); // time points to solve the trajectory
    
    // No FOV constraints
    MinSnapTraj::FOVdata no_fov_data;
    no_fov_data.do_fov = false;

    MinSnapTraj search_prob(search_order, search_times, search_BC, Keyframes, no_fov_data);  // create object

    chrono::steady_clock::time_point tic, toc;                // time the solver
    tic = chrono::steady_clock::now();

    MinSnapTraj::TrajSolution search_sol = search_prob.solveTrajectory();   // solve

    toc = chrono::steady_clock::now();
    showDuration("Search solve time: ", toc - tic);

// ------------------ Analyze/save trajectory ------------------ //
    double time_step = 0.1;
    int numIntersampleTimes = 10;
    int num_times = (int)(1.0 + (double)numIntersampleTimes*(T/time_step));
    vector<double> tspan = MinSnapTraj::linspace(0.0, T, num_times); // time points to evaluate the trajectory
    Eigen::MatrixXd eval_quat = search_prob.FlatOutputTraj(search_sol.coeffs, tspan);

    // Print trajectory data
    //   For more advanced printing: https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    //cout << "---- Trajectory: t, x, y, z, yaw, vx, vy, vz, vyaw, ax, ay, az, ayaw, sx, sy, sz, syaw ----" << endl;
    //cout << eval << endl;

    // Save trajectory in a file within dir_path directory (i.e. extra/tests folder)
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("min_snap_traj_test3.cpp")); // get directory containing source file
    string data_dir = "test_data/";
    ghc::filesystem::create_directories(dir_path + data_dir); // create new directory for data

    string traj_file = "traj_squiggle.csv";
    //saveMatrix(dir_path + data_dir + traj_file, eval_flat);
    saveMatrix(dir_path + data_dir + traj_file, eval_quat);

    string key_file = "keyframes_squiggle.csv";
    saveKeyframes(dir_path + data_dir + key_file, Keyframes);

    cout << "Trajectory data saved at " << endl << dir_path + data_dir << endl << endl;

    return 0;
}
