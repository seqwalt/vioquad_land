#ifndef MIN_SNAP_TRAJ_H_INCLUDED
#define MIN_SNAP_TRAJ_H_INCLUDED

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <random>

#include <QProblem.hpp>

using namespace std;
using namespace Eigen;

class MinSnapTraj {
    public:
        // Variables/Datatypes
        const static int order = 4; // order of piecewise polynomials
        const static int numFlatOut = 4; // x, y, z, yaw
          // Create datatype to hold coefficients for all time intervals
        typedef vector< Matrix<float, order, numFlatOut>, aligned_allocator<Matrix<float, order, numFlatOut>> > flatOutputCoeffs; 
        
        // Constructor
        MinSnapTraj(double duration){
            // Initialize
            finalTime = duration;
            //numCons = 4*numSteps + 5*(numSteps + 1); // 4 for FOV constraint, 5 for boundary conditions 
            numCons = 12;
            numSteps = 2;
            numVars = numSteps*numFlatOut*order;
            coeffs.reserve(numSteps);              // human-friendly coeff storage, type flatOutputCoeffs
            c.reserve(numVars);  // decision variable c, type vector<float>
        }
        
        void getCoeffs() {
            real_t H[numVars*numVars] = {0.0};  // size(H) = [numVars, numVars], init to zeros
            real_t A[numCons*numVars];          // size(A) = [numCons, numVars]
            real_t g[numVars] = {0.0};          // size(g) = [numVars, 1], set to zeros
            real_t lbA[numCons];
            real_t ubA[numCons];
            
            // Keyframes
            double keyframes[] = {
                1.0, 2.0, 4.0, 0.0,  // x, y, z, yaw at initial time
                1.0, 1.0, 3.0, 0.0,  // middle time
                0.0, 0.0, 1.0, 0.0}; // final time
            int array_len = (int)(sizeof(keyframes)/sizeof(keyframes[0]));
            int numKey = array_len / numFlatOut; // number of keyframes
            numSteps = numKey - 1;
            
            // Time points (one for each keyframe)
            vector<double> time = linspace(0.0, finalTime, numKey);
            
            
            // Load values into H
            int k_r = 4; // min snap for pva trajectory
            int k_yaw = 2;
            assert(k_r <= order && k_yaw <= 2);
            float w_r;   // r weight
            float w_yaw; // yaw weight
            
            for(int j = 0; j <= numSteps; j++){      // iterate over time points
                for(int i = k_r; i <= order; i++){   // order
                    // time j, order i
                    int elem_x   = (order*j + i + j)*numFlatOut;
                    int elem_y   = elem_x + 1;
                    int elem_z   = elem_x + 2;
                    
                    w_r = (float)(fact(i)/fact(i - k_r)) * (float)pow(time[j],i-k_r);
                    cout << "w_r^2: " << w_r*w_r << endl;
                    H[numVars*elem_x + elem_x] = w_r*w_r;
                    H[numVars*elem_y + elem_y] = w_r*w_r;
                    H[numVars*elem_z + elem_z] = w_r*w_r;
                }
                for(int i = k_yaw; i <= order; i++){  // order
                    // time j, order i
                    int elem_yaw = (order*j + i + j)*numFlatOut + 3;
                    
                    w_yaw = (float)(fact(i)/fact(i - k_yaw)) * (float)pow(time[j],i-k_yaw);
                    cout << "w_yaw^2: " << w_yaw*w_yaw << endl;
                    H[numVars*elem_yaw + elem_yaw] = w_yaw*w_yaw;
                }
            }
            
            // Position + yaw boundary conditions
            int currCon = 0; // current constraint
            for(int k = 0; k <= numFlatOut - 1; k++){
                for(int j = 0; j <= numSteps - 1; j++){  // iterate over time intervals
                    // First element of current constraint
                    int firstElement = currCon*numVars;
                    
                    for(int i = 0; i <= order; i++){
                        int elem  = (order*j + i + j)*numFlatOut + firstElement;
                        A[elem + k] = pow(time[j], i);
                    }
                    
                    // Set start of interval (time t_j)
                    lbA[currCon] = keyframes[k + j*numFlatOut];
                    ubA[currCon] = keyframes[k + j*numFlatOut];
                    currCon += 1;
                    
                    // Set final conditions
                    if(j == numSteps - 1){
                        // First element of current constraint
                        int firstElement = currCon*numVars;
                        
                        for(int i = 0; i <= order; i++){
                            int elem  = (order*j + i + j)*numFlatOut + firstElement;
                            A[elem + k] = pow(time[j + 1], i);
                        }
                        
                        // Set start of interval (time t_j)
                        lbA[currCon] = keyframes[k + numSteps*numFlatOut];
                        ubA[currCon] = keyframes[k + numSteps*numFlatOut];
                        currCon += 1;
                    }
                }
            }
            //cout << "Total pos + yaw constraints: " << currCon << endl;
            
            // Solve the QP
            QProblem testQP(numVars, numCons);
            testQP.setHessianType(HST_SEMIDEF);
            
            //my_printmatrix("H", H, numVars, numVars);
            
            //testQP.setPrintLevel(PL_HIGH);
            int nWSR = 10; //  maximum number of working set recalculations to be performed during the initial homotopy
            real_t* nullP = NULL;
            returnValue status = testQP.init(H, g, A, nullP, nullP, lbA, ubA, nWSR, 0, 0); // lb, ub set to null pointers
            
            printStatus(status);
            
//             float cOpt[numVars];
//             testQP.getPrimalSolution(cOpt);
//             for(int i = 0; i <= numVars - 1; i++){
//                 cout << cOpt[i] << endl;
//             }           

        }
        
        void my_printmatrix(char *name, float *A, int m, int n) {
            int i, j;

            printf("%s = [...\n", name);
            for (i = 0; i < m; i++) {
                for (j = 0; j < n; j++)
                    printf("  % 9.4f", A[i*n+j]);
                printf(",\n");
            }
            printf("];\n");
        }
        
    private:
        static bool initialized;
        flatOutputCoeffs coeffs;
        vector<float> c;
        int numVars; // number of desicion variables
        int numCons;  // number of constraints
        int numSteps;  // number of time steps
        double finalTime;  // number of time steps
        
        // Factorial function
        int fact(int i){
            assert(i >= 0);
            if (i == 0 || i == 1){
                return 1;
            } else{
                return i*fact(i-1);
            }
        }
        
        template<typename T>
        vector<double> linspace(T start_in, T end_in, int num){
            // https://stackoverflow.com/questions/27028226/python-linspace-in-c
            vector<double> linspaced;
            double start = static_cast<double>(start_in);
            double end = static_cast<double>(end_in);

            if (num == 0) return linspaced;
            if (num == 1){
                linspaced.push_back(start);
                return linspaced;
            }

            double delta = (end - start) / (double)(num - 1);

            for(int i = 0; i < num-1; ++i){
                linspaced.push_back(start + delta * i);
            }
            linspaced.push_back(end); // Ensure that start and end
                                      // are exactly the same as the input
            return linspaced;
        }
        
        template<typename T>
        void printArray(const T (&arr)[]){
            size_t len = sizeof(arr)/sizeof(arr[0]);
            for(size_t i = 0; i <= len - 1; i++){
                cout << arr[i] << endl;
            }
        }
        
        int printStatus(returnValue status){
            string enum_status;
            switch(status) {
                case SUCCESSFUL_RETURN:
                    enum_status = "SUCCESSFUL_RETURN"; break;
                case RET_INIT_FAILED:
                    enum_status = "RET_INIT_FAILED"; break;
                case RET_INIT_FAILED_CHOLESKY:
                    enum_status = "RET_INIT_FAILED_CHOLESKY"; break;
                case RET_INIT_FAILED_TQ:
                    enum_status = "RET_INIT_FAILED_TQ"; break;
                case RET_INIT_FAILED_HOTSTART:
                    enum_status = "RET_INIT_FAILED_HOTSTART"; break;
                case RET_INIT_FAILED_INFEASIBILITY:
                    enum_status = "RET_INIT_FAILED_INFEASIBILITY"; break;
                case RET_INIT_FAILED_UNBOUNDEDNESS:
                    enum_status = "RET_INIT_FAILED_UNBOUNDEDNESS"; break;
                case RET_MAX_NWSR_REACHED:
                    enum_status = "RET_MAX_NWSR_REACHED"; break;
                case RET_INVALID_ARGUMENTS:
                    enum_status = "RET_INVALID_ARGUMENTS"; break;
                case RET_INACCURATE_SOLUTION:
                    enum_status = "RET_INACCURATE_SOLUTION"; break;
                case RET_NO_SOLUTION:
                    enum_status = "RET_NO_SOLUTION"; break;
                default:
                    cout << "Unknown status value: " << status << endl;
                    return 1;
            }
            cout << "Status value: " << status << ". Meaning: " << enum_status << endl;
            return 0;    
        }
};

#endif // MIN_SNAP_TRAJ_H_INCLUDED
















