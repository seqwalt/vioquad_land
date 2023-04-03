/**
 * @file min_snap_traj_test.cpp
 * @brief Test for minimum snap trajector generation with FOV constraints
 */

#include "min_snap_traj.h"

using namespace std;

int main(int argc, char **argv)
{

    float duration = 5.0f; // seconds
    MinSnapTraj obj(duration); // create object
    obj.getCoeffs();

    return 0;
}


