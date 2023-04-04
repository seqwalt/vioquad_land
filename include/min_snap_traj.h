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
#include <algorithm>
#include <numeric>

#include <qpOASES.hpp>

using namespace std;
USING_NAMESPACE_QPOASES

class MinSnapTraj {
    public:
        
        // Datatypes
        struct TrajSolution {
            vector<vector<double>> coeffs;
            int status;
        };
        typedef Eigen::Matrix<double, 2, 4> Matrix24d;
        typedef vector< Eigen::Matrix<double, 2, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 4>> > vectOfMatrix24d;

        // Constructor
        MinSnapTraj(int order, int numIntervals, double duration, const vectOfMatrix24d& bound_conds){
            // Initialize
            finalTime = duration;
            m = numIntervals;
            n = order;
            numVars = m*numFlatOut*(n+1); // numbr of desicion variables
            bounds = bound_conds;
        }
        
        TrajSolution solveTrajectory() {
            
// -------- Position and velocity boundary conditions (BC) -------- //
            int maxBoundOrder = bounds.size() - 1; // 0 if only position BC, 2 if up to acceleration BC
            
            // Time points (one for each keyframe)
            time = linspace(0.0, finalTime, m+1); // number of time points = numIntervals + 1
            
// -------- Hessian: Build H matrix -------- //
            // Specify the trajectory is min snap (minimize 4-th derivative of r)
            int k_r = 4;
            int k_yaw = 2;
            assert(k_r <= n && k_yaw <= 2);
            vector<int> k = {k_r, k_r, k_r, k_yaw}; // store for later use
            
            // Init some variables
            double factorial_term;
            double index_term;
            double integral_term;
            int H_ind;
            
            // Init some matrices
            Eigen::MatrixXd H_eig(numVars, numVars);
            H_eig = Eigen::MatrixXd::Zero(numVars, numVars);
            Eigen::MatrixXd H_block(n+1, n+1);
            Eigen::MatrixXd zero_block = Eigen::MatrixXd::Zero(n+1, n+1);
            
            // Load data into H_eig
            for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){    // iterate over x,y,z,yaw
                for(int ti = 0; ti <= m-1; ti++){                           // iterate over time intervals
                    // Store current block matrix for given flat output (x, y ...)
                    H_block = zero_block;
                    for(int i = k[flat_ind]; i <= n; i++){                  // iterate over degrees
                        for(int j = k[flat_ind]; j <= n; j++){
                            factorial_term = fact(i)*fact(j)/(fact(i-k[flat_ind])*fact(j-k[flat_ind]));
                            index_term = (double)(i+j-2*k[flat_ind]+1);
                            integral_term = (1/index_term)*(pow(time[ti+1],index_term) - pow(time[ti],index_term));
                            if (i == j){
                                H_block(i,j) = factorial_term*integral_term;
                            } else {
                                H_block(i,j) = factorial_term*0.5*integral_term;
                            }
                            
                        }
                    }
                    H_ind = ti*(n+1) + flat_ind*m*(n+1);
                    H_eig.block(H_ind, H_ind, n+1, n+1) = H_block;
                    //cout << "H_block: " << endl
                    //    << H_block << endl << endl;
                }
            }
            
            // Convert H to array that can be used by qpOASES
            double* H = matrixToArray(H_eig);
                 
// -------- Constraints: Load A, lbA and ubA -------- //
            // Constraints: pos/vel boundary consitions, enforced continuity b/w time intervals
            numCons = numFlatOut*(maxBoundOrder+1)*2 + (numFlatOut-1)*(m-1)*(k_r+1) + (m-1)*(k_yaw+1);
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
                            lbA[con] = 0.0;
                            ubA[con] = 0.0;
                            con += 1;
                        }
                    }
                }
            }
              
// -------- Solve the QP -------- //
            // QP problem setup
            double* g = new double[numVars] (); // set g to zeros
            double* nullP = NULL;               // null pointer for lb and ub
            int nWSR = 300;                     //  maximum number of working set recalculations to be performed during the initial homotopy
            
            QProblem testQP(numVars, numCons, HST_SEMIDEF);     // specify QP params (H is semidefinite)
            testQP.setPrintLevel(PL_LOW);                       // Just print errors
            
            // Solve the QP and get the solution
            returnValue status = testQP.init(H, g, A, nullP, nullP, lbA, ubA, nWSR); // lb, ub set to null pointers
            double* cOpt = new double[numVars];
            testQP.getPrimalSolution(cOpt);
            
            // Print stuff
            assert(con == numCons);
            cout << endl << "Number of constraints: " << numCons << endl;
            cout << "Number of decision variables: " << numVars << endl;

            cout << "Time points (sec): ";
            for(size_t i = 0; i < time.size(); i++) cout << time[i] << " ";
            cout << endl << endl;
            
            //cout << "H: " << endl;
            //printMatrix(H,numVars,numVars);
            
            //cout << endl << "A: " << endl;
            //printMatrix(A,numCons,numVars);
            
            //cout << endl << "lbA: " << endl;
            //printMatrix(lbA,1, numCons);
            
            //cout << endl << "ubA: " << endl;
            //printMatrix(ubA,1, numCons);
            
            //cout << endl << "optimal solution: " << endl;
            //printMatrix(cOpt, 1, numVars);
            
            // Load return value
            TrajSolution sol;
            sol.status = status;
            // want sol.coeffs indexed by time, instead of by flat_ind, so change ordering:
            vector<double> ti_vect;
            for(int ti = 0; ti <= m; ti++){
                ti_vect.clear();
                for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){
                    for(int p = 0; p <= n; p++){
                        ti_vect.push_back(cOpt[flat_ind*(n+1)*m + ti*(n+1) + p]);
                    }
                }
                sol.coeffs.push_back(ti_vect);
            }
            
            // Manually clear memory
            delete[] H;
            delete[] A;
            delete[] lbA;
            delete[] ubA;
            delete[] g;
            delete[] cOpt;
            
            return sol;
        }
        
        // Given coeffs, and a vector of times, output the evaluated polynomial
        vector<vector<double>> polyEval(const vector<vector<double>>& coeffs, const vector<double>& tspan){
            assert(tspan[0] >= time[0] && *tspan.end() <= time[m]);
            vector<vector<double>> eval;
            
            // Collect all time power terms
            vector<vector<double>> t_pow_mat;
            vector<double> t_pow_vect;
            for(size_t i = 0; i <= tspan.size()-1; i++){
                t_pow_vect.clear();
                for(int p = 0; p <= n; p++){ t_pow_vect.push_back(pow(tspan[i],p));}
                t_pow_mat.push_back(t_pow_vect);
            }
            
            int j = 1;
            double sum_x, sum_y, sum_z, sum_yaw;
            vector<double> flatout;
            for(size_t i = 0; i <= tspan.size()-1; i++){
                while(true){
                    if (tspan[i] > time[j]){
                        j += 1;
                    } else {
                        // eval j-th segment at ti
                        sum_x   = std::inner_product(begin(t_pow_mat[i]), end(t_pow_mat[i]), begin(coeffs[j]), 0);
                        sum_y   = std::inner_product(begin(t_pow_mat[i]), end(t_pow_mat[i]), begin(coeffs[j]) + n+1, 0);
                        sum_z   = std::inner_product(begin(t_pow_mat[i]), end(t_pow_mat[i]), begin(coeffs[j]) + 2*(n+1), 0);
                        sum_yaw = std::inner_product(begin(t_pow_mat[i]), end(t_pow_mat[i]), begin(coeffs[j]) + 3*(n+1), 0);
                        flatout.clear();
                        flatout.push_back(sum_x); flatout.push_back(sum_y);
                        flatout.push_back(sum_z); flatout.push_back(sum_yaw);
                        eval.push_back(flatout);
                        break;
                    }
                }
            }
            return eval; // eval[i] gives the flatoutputs at time tspan[i]
        }
        
        // Mimic numpy linspace function
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
        
    private:
        int numVars;                // number of desicion variables
        int numCons;                // number of constraints
        int n;                      // order of the polynomials
        int m;                      // number of time intervals
        double finalTime;           // number of time steps
        vector<double> time;        // time points
        const int numFlatOut = 4;   // x, y, z, yaw
        vectOfMatrix24d bounds;     // boundary conditions
        
        // Factorial function
        double fact(int i){
            assert(i >= 0);
            if (i == 0 || i == 1){
                return 1.0f;
            } else{
                return (double)i*fact(i-1);
            }
        }
        
        // Print array as matrix
        void printMatrix(double *A, int rows, int cols) {
            int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++)
                    printf("  % 9.2f", A[i*cols+j]);
                printf(",\n");
            }
        }
        
        // Convert Eigen::MatrixXd to array
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
};

#endif // MIN_SNAP_TRAJ_H_INCLUDED
















