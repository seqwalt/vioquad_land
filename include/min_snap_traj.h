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
        const static int numSteps = 2; // num time intervals (time points - 1)
          // Create datatype to hold coefficients for all time intervals
        typedef vector< Matrix<float, order, numFlatOut>, aligned_allocator<Matrix<float, order, numFlatOut>> > flatOutputCoeffs; 
        
        // Constructor
        MinSnapTraj(float duration){
            // Initialize
            finalTime = duration;
            //numCons = 4*numSteps + 5*(numSteps + 1); // 4 for FOV constraint, 5 for boundary conditions 
            numCons = 12;
            numVars = (numSteps + 1)*numFlatOut*(order+1);
            coeffs.reserve(numSteps);              // human-friendly coeff storage, type flatOutputCoeffs
            c.reserve(numVars);  // decision variable c, type vector<float>
        }
        
        void getCoeffs() {
            //real_t H[numVars*numVars] = {0.0};  // size(H) = [numVars, numVars], init to zeros
            //real_t A[numCons*numVars];          // size(A) = [numCons, numVars]
            //real_t g[numVars] = {0.0};          // size(g) = [numVars, 1], set to zeros
            //real_t lbA[numCons];
            //real_t ubA[numCons];
            
            // Keyframes
            float keyframes[] = {
                1.0f, 2.0f, 4.0f, 0.0f,  // x, y, z, yaw at initial time
                1.0f, 1.0f, 3.0f, 0.0f,  // middle time
                0.0f, 0.0f, 1.0f, 0.0f}; // final time
            int array_len = (int)(sizeof(keyframes)/sizeof(keyframes[0]));
            int numKey = array_len / numFlatOut; // number of keyframes
            //numSteps = numKey - 1;
            
            // Time points (one for each keyframe)
            vector<float> time = linspace(0.0f, finalTime, numKey);
            
            
            // ---- Load values into H ---- //
            int k_r = 4; // min snap for pva trajectory
            int k_yaw = 2;
            assert(k_r <= order && k_yaw <= 2);
            Matrix<float, (numSteps+1)*numFlatOut*(order+1), (numSteps+1)*numFlatOut*(order+1)> H;
            H = Matrix<float, (numSteps+1)*numFlatOut*(order+1), (numSteps+1)*numFlatOut*(order+1)>::Zero();
            
            // Load x part
            float factorial_term;
            float temp;
            float integral_term;
            int H_ind;
            Matrix<float, order+1, order+1> Hx;
            for(int ti = 0; ti <= numSteps; ti++){      // iterate over time points
                // temp matrix to store current time point block diagram for x coeffs
                Hx = Matrix<float, order+1, order+1>::Zero();
                for(int i = k_r; i <= order; i++){      // iterate over orders
                    for(int j = k_r; j <= order; j++){
                        factorial_term = (float)fact(i)*fact(j)/(fact(i-k_r)*fact(j-k_r));
                        temp = (float)(i+j-2*k_r+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hx(i,j) = factorial_term*integral_term;
                        } else {
                            Hx(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(order+1);
                H.block<order+1, order+1>(H_ind, H_ind) = Hx;
                cout << "Hx: " << endl
                     << Hx << endl << endl;
            }
            
            // Load y part
            Matrix<float, order+1, order+1> Hy;
            for(int ti = 0; ti <= numSteps; ti++){      // iterate over time points
                // temp matrix to store current time point block diagram for x coeffs
                Hy = Matrix<float, order+1, order+1>::Zero();
                for(int i = k_r; i <= order; i++){      // iterate over orders
                    for(int j = k_r; j <= order; j++){
                        factorial_term = (float)fact(i)*fact(j)/(fact(i-k_r)*fact(j-k_r));
                        temp = (float)(i+j-2*k_r+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hy(i,j) = factorial_term*integral_term;
                        } else {
                            Hy(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(order+1) + (numSteps+1)*(order+1);
                H.block<order+1, order+1>(H_ind, H_ind) = Hy;
                cout << "Hy: " << endl
                     << Hy << endl << endl;
            }
            
            // Load z part
            Matrix<float, order+1, order+1> Hz;
            for(int ti = 0; ti <= numSteps; ti++){      // iterate over time points
                // temp matrix to store current time point block diagram for x coeffs
                Hz = Matrix<float, order+1, order+1>::Zero();
                for(int i = k_r; i <= order; i++){      // iterate over orders
                    for(int j = k_r; j <= order; j++){
                        factorial_term = (float)fact(i)*fact(j)/(fact(i-k_r)*fact(j-k_r));
                        temp = (float)(i+j-2*k_r+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hz(i,j) = factorial_term*integral_term;
                        } else {
                            Hz(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(order+1) + 2*(numSteps+1)*(order+1);
                H.block<order+1, order+1>(H_ind, H_ind) = Hz;
                cout << "Hz: " << endl
                     << Hz << endl << endl;
            }
            
            // Load yaw part
            Matrix<float, order+1, order+1> Hyaw;
            for(int ti = 0; ti <= numSteps; ti++){      // iterate over time points
                // temp matrix to store current time point block diagram for x coeffs
                Hyaw = Matrix<float, order+1, order+1>::Zero();
                for(int i = k_yaw; i <= order; i++){      // iterate over orders
                    for(int j = k_yaw; j <= order; j++){
                        factorial_term = (float)fact(i)*fact(j)/(fact(i-k_yaw)*fact(j-k_yaw));
                        temp = (float)(i+j-2*k_yaw+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hyaw(i,j) = factorial_term*integral_term;
                        } else {
                            Hyaw(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(order+1) + 3*(numSteps+1)*(order+1);
                H.block<order+1, order+1>(H_ind, H_ind) = Hyaw;
                cout << "Hyaw: " << endl
                     << Hyaw << endl << endl;
            }
            
            cout << "H: " << endl
                 << H << endl;
                 
            // TODO: convert eigen H matrix to 1-D real_t array
            // TODO: create A, ubA, lbA with new desicion var format
            // Solve the QP
            /*
            QProblem testQP(numVars, numCons);
            testQP.setHessianType(HST_SEMIDEF);
            
            //my_printmatrix("H", H, numVars, numVars);
            
            //testQP.setPrintLevel(PL_HIGH);
            int nWSR = 10; //  maximum number of working set recalculations to be performed during the initial homotopy
            real_t* nullP = NULL;
            returnValue status = testQP.init(H, g, A, nullP, nullP, lbA, ubA, nWSR, 0, 0); // lb, ub set to null pointers
            
            printStatus(status);
            */
            
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
        //int numSteps;  // number of time steps
        float finalTime;  // number of time steps
        
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
        vector<float> linspace(T start_in, T end_in, int num){
            // https://stackoverflow.com/questions/27028226/python-linspace-in-c
            vector<float> linspaced;
            float start = static_cast<float>(start_in);
            float end = static_cast<float>(end_in);

            if (num == 0) return linspaced;
            if (num == 1){
                linspaced.push_back(start);
                return linspaced;
            }

            float delta = (end - start) / (float)(num - 1);

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
















