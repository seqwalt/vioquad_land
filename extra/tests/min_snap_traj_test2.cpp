/**
 * @file min_snap_traj_test2.cpp
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

    // For horizontal sweep search path
    int numHlines = 5;          // number of staight horizontal lines
    int numKeys = numHlines*2 - 2;  // number of 3D keypoints minus initial and final positions
    double xWid = 4.0;
    double yWid = 4.0;
    double zInit = 1.5;
    vector<MinSnapTraj::keyframe> Keyframes {}; // store the keyframe
    MinSnapTraj::keyframe k1x, k1y, k1z, k2x, k2y, k2z;
    for (int i = 1; i <= numKeys/2; i++) {
        // only need to specify x, y since z is const.
        k1x.t_ind = 2*i-1; k1x.flat_ind = 0; k1x.value = xWid*(double)(i%2); // time_ind, flat_ind, value
        k1y.t_ind = 2*i-1; k1y.flat_ind = 1; k1y.value = (double)(i-1)*yWid/((double)numHlines-1.0);
        k1z.t_ind = 2*i-1; k1z.flat_ind = 2; k1z.value = 1.5;
        k2x.t_ind = 2*i;   k2x.flat_ind = 0; k2x.value = xWid*(double)(i%2);
        k2y.t_ind = 2*i;   k2y.flat_ind = 1; k2y.value = (double)i*yWid/((double)numHlines-1.0);
        k2z.t_ind = 2*i;   k2z.flat_ind = 2; k2z.value = 1.5;
        Keyframes.push_back(k1x);
        Keyframes.push_back(k1y);
        Keyframes.push_back(k1z);
        Keyframes.push_back(k2x);
        Keyframes.push_back(k2y);
        Keyframes.push_back(k2z);
    }

    // boundary conditions
    MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
    // pos/yaw:   x     y    z   yaw
    p0_bounds <<  0.0,  0.0, zInit, 3.14, // initial
                 xWid, yWid, zInit, 3.14; // final

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

    double spd = 1; // rough speed estimate. units: m/s
    double dtHlines = xWid/spd;
    double dtVlines = 2.0*(yWid/(numHlines-1))/spd;

    double T = numHlines*dtHlines + (numHlines-1)*dtVlines;

    cout << "dtHlines: " << dtHlines << endl;
    cout << "dtVlines: " << dtVlines << endl;
    cout << "       T: " << T << endl;

    vector<double> search_times; // times for the polynomial segments
    search_times.push_back(0.0);
    double curr_time = 0.0;
    for (int i = 0; i < numKeys/2; i++){
        curr_time += dtHlines;
        search_times.push_back(curr_time);
        curr_time += dtVlines;
        search_times.push_back(curr_time);
    }
    search_times.push_back(T);

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
    //vector<double> tspan = MinSnapTraj::linspace(0.0, T, 50); // time points to evaluate the trajectory
    //Eigen::MatrixXd eval_flat = prob.FlatOutputTraj(sol.coeffs, tspan);
    Eigen::MatrixXd eval_quat = search_prob.QuaternionTraj(search_sol.coeffs, tspan);

    // Print trajectory data
    //   For more advanced printing: https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    //cout << "---- Trajectory: t, x, y, z, yaw, vx, vy, vz, vyaw, ax, ay, az, ayaw, sx, sy, sz, syaw ----" << endl;
    //cout << eval << endl;

    // Save trajectory in a file within dir_path directory (i.e. extra/tests folder)
    string file_path = __FILE__;
    string dir_path = file_path.substr(0, file_path.rfind("min_snap_traj_test2.cpp")); // get directory containing source file
    string data_dir = "test_data/";
    ghc::filesystem::create_directories(dir_path + data_dir); // create new directory for data

    string traj_file = "test_traj.csv";
    //saveMatrix(dir_path + data_dir + traj_file, eval_flat);
    saveMatrix(dir_path + data_dir + traj_file, eval_quat);

    string key_file = "keyframes.csv";
    saveKeyframes(dir_path + data_dir + key_file, Keyframes);

    cout << "Trajectory data saved at " << endl << dir_path + data_dir << endl << endl;

    return 0;
}
