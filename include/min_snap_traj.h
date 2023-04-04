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

//#include <QProblem.hpp>
#include <qpOASES.hpp>

using namespace std;
USING_NAMESPACE_QPOASES

class MinSnapTraj {
    public:
        
        
        // Variables/Datatypes
        const static int n = 4; // highest order of the piecewise polynomials
        const static int numFlatOut = 4; // x, y, z, yaw
        typedef vector< Eigen::Matrix<double, 2, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 4>> > vectOfMatrix24d;
        typedef Eigen::Matrix<double, 2, 4> Matrix24d;

        // Constructor
        MinSnapTraj(int numIntervals, double duration){
            // Initialize
            finalTime = duration;
            m = numIntervals;
            numVars = m*numFlatOut*(n+1);
        }
        
        int getCoeffs() {
            
            // Position and velocity boundary conditions (BC)
            Matrix24d p0_bounds; // p=0 means 0th derivative
            p0_bounds << 1.0f, 2.0f, 4.0f, 0.0f, // initial x, y, z, yaw
                         0.0f, 0.0f, 1.0f, 0.0f; // final values
            Matrix24d p1_bounds; // p=1 means 1st derivative
            p1_bounds << 0.0f, 0.0f, 0.0f, 0.0f, // initial velocities
                         0.0f, 0.0f, 0.0f, 0.0f; // final velocities
            vectOfMatrix24d bounds;
            bounds.push_back(p0_bounds);
            bounds.push_back(p1_bounds);
            int maxBoundOrder = bounds.size() - 1; // 0 if only position BC, 2 if up to acceleration BC
            cout << "maxBoundOrder: " << maxBoundOrder << endl << endl;
            
            // Time points (one for each keyframe)
            vector<double> time = linspace(0.0, finalTime, m+1); // number of time points = numIntervals + 1
            cout << "times: " << endl;
            for(size_t i = 0; i < time.size(); i++) cout << time[i] << " ";
            cout << endl << endl;
            
            // ---- Load values into H ---- //
            int k_r = 4; // min snap for pva trajectory
            int k_yaw = 2;
            vector<int> k = {k_r, k_r, k_r, k_yaw}; // derivative for objective function for x,y,z,yaw
            
            cout << "k: " << endl;
            for(size_t i = 0; i < k.size(); i++) cout << k[i] << " ";
            cout << endl << endl;
            
            assert(k_r <= n && k_yaw <= 2);
            Eigen::MatrixXd H_eig(numVars, numVars);
            H_eig = Eigen::MatrixXd::Zero(numVars, numVars);

            // Load x part
            double factorial_term;
            double temp;
            double integral_term;
            int H_ind;
            Eigen::Matrix<double, n+1, n+1> Hx;
            for(int ti = 0; ti <= m-1; ti++){      // iterate over time intervals
                // temp matrix to store current time point block diagram for x coeffs
                Hx = Eigen::Matrix<double, n+1, n+1>::Zero();
                for(int i = k_r; i <= n; i++){      // iterate over orders
                    for(int j = k_r; j <= n; j++){
                        factorial_term = fact(i)*fact(j)/(fact(i-k_r)*fact(j-k_r));
                        temp = (double)(i+j-2*k_r+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hx(i,j) = factorial_term*integral_term;
                        } else {
                            Hx(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(n+1);
                H_eig.block<n+1, n+1>(H_ind, H_ind) = Hx;
                //cout << "Hx: " << endl
                //    << Hx << endl << endl;
            }
            
            // Load y part
            Eigen::Matrix<double, n+1, n+1> Hy;
            for(int ti = 0; ti <= m-1; ti++){      // iterate over time intervals
                // temp matrix to store current time point block diagram for x coeffs
                Hy = Eigen::Matrix<double, n+1, n+1>::Zero();
                for(int i = k_r; i <= n; i++){      // iterate over orders
                    for(int j = k_r; j <= n; j++){
                        factorial_term = fact(i)*fact(j)/(fact(i-k_r)*fact(j-k_r));
                        temp = (double)(i+j-2*k_r+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hy(i,j) = factorial_term*integral_term;
                        } else {
                            Hy(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(n+1) + m*(n+1);
                H_eig.block<n+1, n+1>(H_ind, H_ind) = Hy;
                //cout << "Hy: " << endl
                //    << Hy << endl << endl;
            }
            
            // Load z part
            Eigen::Matrix<double, n+1, n+1> Hz;
            for(int ti = 0; ti <= m-1; ti++){      // iterate over time intervals
                // temp matrix to store current time point block diagram for x coeffs
                Hz = Eigen::Matrix<double, n+1, n+1>::Zero();
                for(int i = k_r; i <= n; i++){      // iterate over orders
                    for(int j = k_r; j <= n; j++){
                        factorial_term = fact(i)*fact(j)/(fact(i-k_r)*fact(j-k_r));
                        temp = (double)(i+j-2*k_r+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hz(i,j) = factorial_term*integral_term;
                        } else {
                            Hz(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(n+1) + 2*m*(n+1);
                H_eig.block<n+1, n+1>(H_ind, H_ind) = Hz;
                //cout << "Hz: " << endl
                //    << Hz << endl << endl;
            }
            
            // Load yaw part
            Eigen::Matrix<double, n+1, n+1> Hyaw;
            for(int ti = 0; ti <= m-1; ti++){      // iterate over time intervals
                // temp matrix to store current time point block diagram for x coeffs
                Hyaw = Eigen::Matrix<double, n+1, n+1>::Zero();
                for(int i = k_yaw; i <= n; i++){      // iterate over orders
                    for(int j = k_yaw; j <= n; j++){
                        factorial_term = fact(i)*fact(j)/(fact(i-k_yaw)*fact(j-k_yaw));
                        temp = (double)(i+j-2*k_yaw+1);
                        integral_term = (1/temp)*(pow(time[ti+1],temp) - pow(time[ti],temp));
                        if (i == j){
                            Hyaw(i,j) = factorial_term*integral_term;
                        } else {
                            Hyaw(i,j) = factorial_term*0.5*integral_term;
                        }
                        
                    }
                }
                H_ind = ti*(n+1) + 3*m*(n+1);
                H_eig.block<n+1, n+1>(H_ind, H_ind) = Hyaw;
                //cout << "Hyaw: " << endl
                //    << Hyaw << endl << endl;
            }
                 
            // Constraints
            numCons = 34;
            double* A = new double[numCons*numVars] ();
            double* lbA = new double[numCons] ();
            double* ubA = new double[numCons] ();
            
            int con = 0; // current constraint index
            int A_ind, A_ind_seg1, A_ind_seg2;   // 1-D array indices of constraint matrix A
            for(int ti = 0; ti <= m; ti++){      // iterate over time points
                if(ti == 0){ // initial time
                    // Fix initial positions/yaw and velocities
                    for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){ // iterate over x,y,z,yaw
                        for(int p = 0; p <= maxBoundOrder; p++){ // iterate over derivatives (p=0 means positions, p=1 means velocities)
                            // create the constraint
                            for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                                A_ind = i + flat_ind*(n+1)*m + con*numVars;
                                A[A_ind] = pow(time[ti],i-p)*fact(i)/fact(i-p);
                            }
                            lbA[con] = bounds[p](0,flat_ind); // initial conditions for p-th derivative
                            ubA[con] = lbA[con];
                            con += 1;
                        }
                    }
                } else if(ti == m){ // final time
                    // Fix final positions/yaw and velocities
                    for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){ // iterate over x,y,z,yaw
                        for(int p = 0; p <= maxBoundOrder; p++){ // iterate over derivatives (p=0 means positions, p=1 means velocities)
                            // create the constraint
                            for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                                A_ind = (ti-1)*(n+1) + i + flat_ind*(n+1)*m + con*numVars;
                                A[A_ind] = pow(time[ti],i-p)*fact(i)/fact(i-p);
                            }
                            lbA[con] = bounds[p](1,flat_ind); // final conditions for p-th derivative
                            ubA[con] = lbA[con];
                            con += 1;
                        }
                    }
                } else { // other times
                    // Enforce continuity of the first k_r derivates of position
                    //  and the first k_yaw derivates of yaw
                    for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){ // iterate over x,y,z,yaw
                        for(int p = 0; p <= k[flat_ind]; p++){ // p = 0 means continuity of 0-th derivative
                            // create the constraint
                            for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                                A_ind_seg1 = (ti-1)*(n+1) + i + flat_ind*(n+1)*m + con*numVars;
                                A_ind_seg2 =     ti*(n+1) + i + flat_ind*(n+1)*m + con*numVars;
                                A[A_ind_seg1] = pow(time[ti],i-p)*fact(i)/fact(i-p);
                                A[A_ind_seg2] = -A[A_ind_seg1];
                            }
                            lbA[con] = 0.0f;
                            ubA[con] = 0.0f;
                            con += 1;
                        }
                    }
                }
            }
            
            // QP problem setup
            double* H = matrixToArray(H_eig);   // convert H to array
            double* g = new double[numVars] (); // set g to zeros
            double* nullP = NULL;               // null pointer for lb and ub
            int nWSR = 100;                     //  maximum number of working set recalculations to be performed during the initial homotopy
            
            QProblem testQP(numVars, numCons, HST_SEMIDEF);     // specify QP params (H is semidefinite)
            //testQP.setPrintLevel(PL_HIGH);
            
            // Solve the QP
            returnValue status = testQP.init(H, g, A, nullP, nullP, lbA, ubA, nWSR); // lb, ub set to null pointers
            
            // Print matrices
            printStatus(status);
            
            cout << "number of constraints: " << con << endl << endl;
            cout << "number of decision vars: " << numVars << endl << endl;
            
            cout << "H: " << endl;
            printMatrix(H,numVars,numVars);
            
            cout << endl << "A: " << endl;
            printMatrix(A,numCons,numVars);
            
            cout << endl << "lbA: " << endl;
            printMatrix(lbA,1, numCons);
            
            cout << endl << "ubA: " << endl;
            printMatrix(ubA,1, numCons);
            
            double* cOpt = new double[numVars];
            testQP.getPrimalSolution(cOpt);
            cout << endl << "c: " << endl;
            printMatrix(cOpt, 1, numVars);
            
            // Manually clear memory
            delete[] H;
            delete[] A;
            delete[] lbA;
            delete[] ubA;
            delete[] g;
            delete[] cOpt;
            
            return status;
        }
        
        void printMatrix(double *A, int rows, int cols) {
            int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++)
                    printf("  % 9.2f", A[i*cols+j]);
                printf(",\n");
            }
        }
        
    private:
        int numVars; // number of desicion variables
        int numCons;  // number of constraints
        int m;            // number of time intervals
        double finalTime;  // number of time steps
        
        // Factorial function
        double fact(int i){
            assert(i >= 0);
            if (i == 0 || i == 1){
                return 1.0f;
            } else{
                return (double)i*fact(i-1);
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
        
        double* matrixToArray(const Eigen::MatrixXd& matrix) {
            // Get the number of rows and columns in the matrix
            int rows = matrix.rows();
            int cols = matrix.cols();

            // Allocate memory for the array
            double* array = new double[rows * cols];

            // Copy the matrix data into the array
            int index = 0;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    array[index++] = matrix(i, j);
                }
            }
            return array;
        }
        
        int printStatus(qpOASES::returnValue status){
            USING_NAMESPACE_QPOASES
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
                //case RET_INACCURATE_SOLUTION:
                //    enum_status = "RET_INACCURATE_SOLUTION"; break;
                //case RET_NO_SOLUTION:
                //    enum_status = "RET_NO_SOLUTION"; break;
                default:
                    cout << "Unknown status value: " << status << endl;
                    return 1;
            }
            cout << "Status value: " << status << ". Meaning: " << enum_status << endl;
            return 0;    
        }
};

#endif // MIN_SNAP_TRAJ_H_INCLUDED
















