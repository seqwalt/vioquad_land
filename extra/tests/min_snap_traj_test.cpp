/**
 * @file min_snap_traj_test.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 */

#include "min_snap_traj.h"

using namespace std;

int main(int argc, char **argv)
{
    int numIntervals = 2;
    double duration = 5.0; // seconds
    MinSnapTraj obj(numIntervals, duration); // create object

    int stats = obj.getCoeffs();
    cout << "my status: " << stats << endl;
    
    return 0;
}


