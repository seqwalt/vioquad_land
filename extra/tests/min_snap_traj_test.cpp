/**
 * @file min_snap_traj_test.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 */

#include "min_snap_traj.h"

using namespace std;

void printVectorAsMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char **argv)
{
    // Params
    int order = 4;         // order of piecewise polynomials
    int numIntervals = 5; // number of time intervals
    double duration = 5.0; // duration of trajectory in seconds
    
    
    // Position and velocity boundary conditions
    MinSnapTraj::Matrix24d p0_bounds; // p=0 means 0th derivative
    p0_bounds << 1.0, 2.0, 4.0, 0.0, // initial x, y, z, yaw
                 0.0, 0.0, 1.0, 0.0; // final values
    MinSnapTraj::Matrix24d p1_bounds; // p=1 means 1st derivative
    p1_bounds << 0.0, 0.0, 0.0, 0.0, // initial velocities
                 0.0, 0.0, 0.0, 0.0; // final velocities
    MinSnapTraj::vectOfMatrix24d bounds;
    bounds.push_back(p0_bounds);
    bounds.push_back(p1_bounds);
    
    MinSnapTraj prob(order, numIntervals, duration, bounds); // create object

    MinSnapTraj::TrajSolution sol = prob.solveTrajectory();
    
    vector<double> tspan = prob.linspace(0.0, 5.0, 100);
    vector<vector<double>> eval = prob.polyEval(sol.coeffs, tspan);
    printVectorAsMatrix(eval);
    
    return 0;
}


