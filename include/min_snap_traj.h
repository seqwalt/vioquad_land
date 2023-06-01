#ifndef MIN_SNAP_TRAJ_H_INCLUDED
#define MIN_SNAP_TRAJ_H_INCLUDED

#include <Eigen/Core>
#include <Eigen/Geometry>

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

class MinSnapTraj {
    public:
        
        // Datatypes
        struct TrajSolution {
            vector<vector<double>> coeffs;
            int status;
        };
        struct keyframe{
            int t_ind;    // time index. t_ind \in {1, ..., m-1}
            int flat_ind; // flat_ind = 0 for x, 1 for y, 2 for z, 3 for yaw
            double value; // The desired value of the state corresponding to flat_ind
        };
        struct FOVdata{
            bool do_fov;        // true to do fov constraint, false otherwise
            vector<double> l;   // landmark 3D position
            double alpha_x;     // half of horizontal FOV (radians)
            double alpha_y;     // half of vertical FOV (radians)
        };
        typedef Eigen::Matrix<double, 2, 4> Matrix24d;
        typedef vector< Eigen::Matrix<double, 2, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 2, 4>> > vectOfMatrix24d;

        // Constructor
        MinSnapTraj(int order, vector<double> time_in, vectOfMatrix24d& bounds_in, vector<keyframe> keys_in, FOVdata fov){
            // Initialize
            assert(order >= 4);
            n = order;                    // polynomial order
            time = time_in;               // time points
            m = time.size() - 1;          // number of time segments
            
            do_fov_constraint = fov.do_fov;
            if (do_fov_constraint){
                numVars = m*numFlatOut*(n+1) + 2; // number of desicion variables (with 2 slack variables)
                l = fov.l;
                alpha_x = fov.alpha_x;
                alpha_y = fov.alpha_y;
            } else {
                numVars = m*numFlatOut*(n+1); // number of desicion variables (no slack variables)
            }
            
            H = buildHessian();
            g = buildGvector();
            
            bounds = bounds_in;           // boundary conditions
            keys = keys_in;               // keyframes
        }
        
        TrajSolution solveTrajectory() {
            
            buildConstraints();
              
// -------- Solve the QP -------- //
            // QP problem setup
            double* nullP = NULL;               // null pointer for lb and ub
            int nWSR = 10000;                     //  maximum number of working set recalculations to be performed during the initial homotopy
            
            qpOASES::QProblem testQP(numVars, numCons, qpOASES::HST_SEMIDEF);     // specify QP params (H is semidefinite)
            testQP.setPrintLevel(qpOASES::PL_LOW);                       // Just print errors
            
//             qpOASES::Options options;
//             options.setToReliable();
//             testQP.setOptions(options);
            
            // Solve the QP and get the solution
            qpOASES::returnValue status = testQP.init(H, g, A, nullP, nullP, lbA, ubA, nWSR); // lb, ub set to null pointers
            cOpt = new double[numVars];
            testQP.getPrimalSolution(cOpt);
            
            // Print stuff
            cout << endl << "Number of constraints: " << numCons << endl;
            cout << "Number of decision variables: " << numVars << endl;

            cout << "Segment times (seconds): ";
            cout << "[";
            for(size_t i = 0; i < time.size()-1; i++) cout << time[i] << ", ";
            cout << time[time.size()-1] << "]" << endl;
            
            //cout << "H: " << endl;
            //printArrayAsMatrix(H,numVars,numVars);
            
            //cout << endl << "A: " << endl;
            //printArrayAsMatrix(A,numCons,numVars);
            
            //cout << endl << "lbA: " << endl;
            //printArrayAsMatrix(lbA,1, numCons);
            
            //cout << endl << "ubA: " << endl;
            //printArrayAsMatrix(ubA,1, numCons);
            
            //cout << "SIZE of original optimal solution: " << numVars << endl;
            //cout << endl << "original optimal solution: " << endl;
            //printArrayAsMatrix(cOpt, 1, numVars);
            
// -------- Load return value -------- //
           
            // NOTE: The decision variable c_ (not including slack vars) used for the return of solveTrajectory has the form:
            // c_ = [x_{0,0}, ..., x_{n,0}, y_{0,0}, ..., y_{n,0}, . . ., yaw_{n,0}, x_{0,1}, ..., x_{n,1}, y_{0,1}, . . .],
            // where x_{i,j} is the x coefficient of order i, at time index j.
            // The array c_ is "ordered by time" meaning time index 0 comes first, then 1 and so on.
            
            TrajSolution sol;
            sol.status = status;
            // want sol.coeffs indexed by time, instead of by state (see NOTEs), so change ordering:
            vector<double> ti_vect;
            for(int ti = 0; ti <= m-1; ti++){
                ti_vect.clear();
                for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){
                    for(int p = 0; p <= n; p++){
                        ti_vect.push_back(cOpt[flat_ind*(n+1)*m + ti*(n+1) + p]); // no collection of slack variables
                    }
                }
                sol.coeffs.push_back(ti_vect);
            }
            
            //cout << "SIZE of re-ordered optimal solution: " << sol.coeffs.size()*sol.coeffs[0].size() << endl;
            cout << "Optimal coefficients (one row per time segment): " << endl;
            printVectorOfVectors(sol.coeffs, 2);
            
            if(do_fov_constraint){
                cout << "Slack variables: " << cOpt[numVars - 2] << " " << cOpt[numVars - 1] << endl;
            }
            
            return sol;
        }
        
        void updateBounds(vectOfMatrix24d& bounds_new) {
            bounds = bounds_new;
        }
        
        void deleteQpData() {
            delete[] H;
            delete[] A;
            delete[] g;
            delete[] lbA;
            delete[] ubA;
            delete[] cOpt;
        }
        
        //  Input: coeffs output by solveTrajectory (i.e. rows correspond to polynomial segments), and a vector of times tspan
        // Output: matrix with each row: time, x, y, z, yaw, vx, vy, vz, vyaw, ax, ay, az, ayaw, jx, jy, jz, jyaw
        Eigen::MatrixXd FlatOutputTraj(const vector<vector<double>>& coeffs, const vector<double>& tspan){
            assert(tspan[0] >= time[0] && tspan.back() <= time[m]);
            
            Eigen::MatrixXd eval = Eigen::MatrixXd::Zero(tspan.size(), 17); // row: t, x, y, z, yaw, vx, vy, vz, vyaw, ax, ay, az, ayaw, jx, jy, jz, jyaw
            
            vector<int> seg_ind = polySegmentInds(tspan); // get vector of polynomial segments, given vector of times 
            double time_term;
            for(size_t ti = 0; ti <= tspan.size()-1; ti++){ // iterate over tspan times
                eval(ti, 0) = tspan[ti];
                for(int p = 0; p <= 3; p++){      // iterate over orders (pos is p=0, vel is p=1, acc is p=2, jerk is p=3)
                    for(int i = p; i <= n; i++){  // sum over polynomial terms for given order
                        time_term = pow(tspan[ti], i-p) * (fact(i)/fact(i-p));
                        eval(ti, 1 + p*4) += time_term * coeffs[seg_ind[ti]][i];           // x, vx, ax etc.
                        eval(ti, 2 + p*4) += time_term * coeffs[seg_ind[ti]][i + n+1];     // y, vy, ay etc.
                        eval(ti, 3 + p*4) += time_term * coeffs[seg_ind[ti]][i + 2*(n+1)]; // z, vz, az etc.
                        eval(ti, 4 + p*4) += time_term * coeffs[seg_ind[ti]][i + 3*(n+1)]; // yaw, vyaw, ayaw etc.
                    }
                }
            }
            
            return eval;
        }
        
        //  Input: coeffs output by solveTrajectory (i.e. rows correspond to polynomial segments), and a vector of times tspan
        // Output: matrix with each row: time, x, y, z, vx, vy, vz, qw, qx, qy, qz
        Eigen::MatrixXd QuaternionTraj(const vector<vector<double>>& coeffs, const vector<double>& tspan){
            assert(tspan[0] >= time[0] && tspan.back() <= time[m]);
            vector<int> seg_ind = polySegmentInds(tspan); // get vector of polynomial segments, given vector of times 
            
            Eigen::MatrixXd eval = Eigen::MatrixXd::Zero(tspan.size(), 11); // time, x, y, z, vx, vy, vz, qw, qx, qy, qz
            double time_term, ax, ay, az, yaw;
            Eigen::Quaterniond quat;
            Eigen::Vector3d acc;

            for(size_t ti = 0; ti <= tspan.size()-1; ti++){ // iterate over tspan times
                eval(ti, 0) = tspan[ti];
                // Position and Velocity
                for(int p = 0; p <= 1; p++){      // iterate over orders (pos is p=0, vel is p=1)
                    for(int i = p; i <= n; i++){  // sum over polynomial terms for given order
                        time_term = pow(tspan[ti], i-p) * (fact(i)/fact(i-p));
                        eval(ti, 1 + p*3) += time_term * coeffs[seg_ind[ti]][i];           // x, vx
                        eval(ti, 2 + p*3) += time_term * coeffs[seg_ind[ti]][i + n+1];     // y, vy
                        eval(ti, 3 + p*3) += time_term * coeffs[seg_ind[ti]][i + 2*(n+1)]; // z, vz
                    }
                }
                // Acceleration
                int p = 2; // acc is p=2
                ax = 0; ay = 0; az = 0; // reset
                for(int i = p; i <= n; i++){  // sum over polynomial terms for acceleration
                    time_term = pow(tspan[ti], i-p) * (fact(i)/fact(i-p));
                    ax += time_term * coeffs[seg_ind[ti]][i];           // ax
                    ay += time_term * coeffs[seg_ind[ti]][i + n+1];     // ay
                    az += time_term * coeffs[seg_ind[ti]][i + 2*(n+1)]; // az
                }
                // Yaw
                p = 0; // pos/yaw is p=0
                yaw = 0; // reset
                for(int i = p; i <= n; i++){  // sum over polynomial terms for acceleration
                    time_term = pow(tspan[ti], i-p) * (fact(i)/fact(i-p));
                    yaw += time_term * coeffs[seg_ind[ti]][i + 3*(n+1)]; // yaw
                }
                // Compute Quaternion from accleration and yaw
                acc << ax, ay, az;
                quat = acc2quaternion(acc, yaw);
                eval(ti, 7) = -quat.coeffs().w();
                eval(ti, 8) = -quat.coeffs().x();
                eval(ti, 9) = -quat.coeffs().y();
                eval(ti,10) = -quat.coeffs().z();
            }
            
            return eval;
        }
        
        // Return vector of polynomial segment indices,
        //   corrseponding to given vector of times
        vector<int> polySegmentInds(const vector<double>& new_times){
            int seg_ind = 0;
            vector<int> seg_ind_vect;
            for(size_t i = 0; i <= new_times.size()-1; i++){
                while(true){                  
                    if (new_times[i] > time[seg_ind + 1]){
                        seg_ind += 1;
                    } else {
                        // store segment index
                        seg_ind_vect.push_back(seg_ind);
                        break;
                    }
                }
            }
            
            // check seg_ind_vect correctness:
            //cout << endl << "seg_ind check: " << endl;
            //for(size_t i = 0; i <= new_times.size()-1; i++){
            //    assert(time[seg_ind_vect[i]] <= new_times[i]);
            //    assert(new_times[i] <= time[seg_ind_vect[i] + 1]);
            //    
            //    // print check
            //    cout << time[seg_ind_vect[i]] << " <= " << new_times[i] << " <= " << time[seg_ind_vect[i] + 1] << endl;
            //}
            
            return seg_ind_vect;
        }
        
        // Mimic numpy linspace function
        template<typename T>
        static vector<double> linspace(T start_in, T end_in, int num){
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
        vector<double> time;        // time points
        const int numFlatOut = 4;   // x, y, z, yaw
        vectOfMatrix24d bounds;     // boundary conditions
        vector<keyframe> keys;      // keyframes
        vector<double> l;           // landmark 3D position
        double alpha_x;             // half of horizontal FOV (radians)
        double alpha_y;             // half of vertical FOV (radians)
        bool do_fov_constraint;     // true if doing fov constraint, false otherwise
        
        int k_r = 4;       // order of polynomial for position
        int k_yaw = 2;     // order of polynomial for yaw

        // Variables for use by qpOases
        double* H;
        double* A;
        double* lbA;
        double* ubA;
        double* g;
        double* cOpt;
        
        double* buildHessian() {
            // -------- Hessian: Build H matrix -------- //
            
            // NOTE: The decision variable c used to solve the QP has the following form:
            // c = [x_{0,0}, x_{1,0}, ..., x_{n,0}, x_{0,1}, x_{1,1}, ..., x_{n,1}, . . ., x_{n,m-1}, y_{0,0}, y_{1,0}, . . .],
            // where x_{i,j} is the x coefficient of order i, at time index j.
            // The array c is "ordered by state" meaning x comes first, then y, z, yaw.
            
            // Specify the trajectory is min snap (minimize 4-th (i.e. k_r-th) derivative of r)
            assert(k_r <= n && k_yaw <= 2);
            vector<int> k = {k_r, k_r, k_r, k_yaw}; // store for later use
            
            // Init some variables
            double factorial_term;
            double index_term;
            double integral_term;
            int H_ind;
            
            // Init some matrices
            Eigen::MatrixXd H_min_snap(numVars, numVars);
            Eigen::MatrixXd H_min_z(numVars, numVars);
            H_min_snap = Eigen::MatrixXd::Zero(numVars, numVars);
            H_min_z = Eigen::MatrixXd::Zero(numVars, numVars);
            Eigen::MatrixXd H_block(n+1, n+1);
            Eigen::MatrixXd zero_block = Eigen::MatrixXd::Zero(n+1, n+1);
            
            // Load data into H_min_snap and H_min_z
            for(int ti = 0; ti <= m-1; ti++){                           // iterate over time intervals
                // Minimize Snap
                for(int flat_ind = 0; flat_ind <= numFlatOut-1; flat_ind++){    // iterate over x,y,z,yaw
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
                    H_min_snap.block(H_ind, H_ind, n+1, n+1) = H_block;
                    //cout << "H_block: " << endl
                    //    << H_block << endl << endl;
                }
                
                // Minimize z (flat_ind = 2)
                if(do_fov_constraint){
                    H_block = zero_block;
                    for(int i = 0; i <= n; i++){ // iterate over degrees
                        for(int j = 0; j <= n; j++){
                            index_term = (double)(i+j+1);
                            integral_term = (1/index_term)*(pow(time[ti+1],index_term) - pow(time[ti],index_term));
                            if (i == j){
                                H_block(i,j) = integral_term;
                            } else {
                                H_block(i,j) = 0.5*integral_term;
                            }
                            
                        }
                    }
                    H_ind = ti*(n+1) + 2*m*(n+1);
                    H_min_z.block(H_ind, H_ind, n+1, n+1) = H_block;
                }
            }
            
            // Convert from Eigen matrix to array that can be used by qpOASES
            Eigen::MatrixXd H_eig(numVars, numVars);
            double min_z_weight = 8.0;
            H_eig = H_min_snap + min_z_weight*H_min_z;
            H = matrixToArray(H_eig);
            return H;
        }
        
        void buildConstraints() {
            // -------- Constraints: Load A, lbA and ubA -------- //
            
            // Constraint types:
            // - pos/vel boundary consitions
            // - enforced continuity b/w time intervals
            // - keyframes/waypoints
            // - camera FOV constraints
            
            int maxBoundOrder = bounds.size() - 1; // 0 if only position boundary condition (BC), 2 if up to acceleration BC, for example
            int numFOVtimes;
            if (do_fov_constraint){
                numFOVtimes = 8;  // (>= 2) number of time instants to enforce FOV constraint (use 4 (resp. 8) when theta values are fixed (resp. have +-0.25 constraint relaxation))
            } else {
                numFOVtimes = 2;  // setting to 2 is same as no fov constraints, becaues the initial and final conditions are already fixed
            }
            
            int posConSides = 8;    // number of sides of cone in position constraint
            
            int numBoundCons = numFlatOut*(maxBoundOrder+1)*2;                  // number of boundary condition constraints
            int numContCons  = (numFlatOut-1)*(m-1)*(k_r+1) + (m-1)*(k_yaw+1);  // number of "enforced continuity b/w time interval" constraints
            int numKeyCons   = keys.size();                                     // number of keyframe/waypoint constraints (not including boundary conditions)
            int numFovCons   = 4*(numFOVtimes - 2);                             // number of field-of-view (FOV) constraints
            int numPosCons   = (posConSides + 2)*(numFOVtimes - 2);             // number of position constraints
            numCons = numBoundCons + numContCons + numKeyCons + numFovCons + numPosCons;
            
            A = new double[numCons*numVars] ();
            lbA = new double[numCons] ();
            ubA = new double[numCons] ();
            vector<int> k = {k_r, k_r, k_r, k_yaw};
            
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
            
            // Waypoint/Keyframe constraints (separate from boundary conditions)
            for(int key_ind = 0; key_ind <= numKeyCons-1; key_ind++){ // iterate over x,y,z,yaw
                // create the constraint
                int flat_ind = keys[key_ind].flat_ind;
                int ti = keys[key_ind].t_ind;
                for(int i = 0; i <= n; i++){  // iterate over degrees of polynomial segment
                    A_ind = (ti-1)*(n+1) + i + flat_ind*(n+1)*m + con*numVars;
                    A[A_ind] = pow(time[ti],i);
                }
                lbA[con] = keys[key_ind].value;
                ubA[con] = lbA[con];
                con += 1;
            }
                        
            // FOV ax,ay,az constraints (with slack variables)
            if (do_fov_constraint){
                // Interpolate theta from initial time to final time
                double x0 = bounds[0](0,0); double xm = bounds[0](1,0);
                double y0 = bounds[0](0,1); double ym = bounds[0](1,1);
                double z0 = bounds[0](0,2); double zm = bounds[0](1,2);
                double theta_x0 = atan2(l[0] - x0, z0 - l[2]); // initial theta for ax az constraint
                double theta_xm = atan2(l[0] - xm, zm - l[2]); // final theta for ax az constraint
                double theta_y0 = atan2(l[1] - y0, z0 - l[2]); // initial theta for ay az constraint
                double theta_ym = atan2(l[1] - ym, zm - l[2]); // final theta for ay az constraint
                
                vector<double> FOVtime = linspace(time[0], time.back(), numFOVtimes); // times to enforce FOV constraint
                vector<double> Theta_x = linspace(theta_x0, theta_xm, numFOVtimes);   // interpolation
                vector<double> Theta_y = linspace(theta_y0, theta_ym, numFOVtimes);
                
                // Define FOV in radians
                double alpha = min(alpha_x, alpha_y);
                // Check FOV constraint satisfied at initial and final times
                bool theta_x_BC_FOV = ( (alpha_x - abs(theta_x0) > 1e-3) && (alpha_x - abs(theta_xm) > 1e-3) );
                bool theta_y_BC_FOV = ( (alpha_y - abs(theta_y0) > 1e-3) && (alpha_y - abs(theta_ym) > 1e-3) );
                assert(theta_x_BC_FOV && theta_y_BC_FOV);
                
                double tanx, tany;
                double grav = 9.81;
                int p = 2;      // constraint is on acceleration
                int flat_ind_x = 0;
                int flat_ind_y = 1;
                int flat_ind_z = 2;
                int A_ind_x, A_ind_y, A_ind_z; 
                int A_ind_ax, A_ind_ay, A_ind_az;
                int A_ind_slack_x, A_ind_slack_y;
                double pos_coeff, acc_coeff;
                vector<int> seg_ind = polySegmentInds(FOVtime); // get vector of segment indices, given vecotr of times

                // Satisfying the FOV constraint requires a position and acceleration constraints
                for(int ti = 1; ti <= numFOVtimes - 2; ti++){
                    // Position constraints
                    for(int ang_ind = 0; ang_ind <= posConSides - 1; ang_ind++){ // iterate over approx position constraint cone
                        // create position cone constraint
                        double cos_term = cos((double)(2.0*M_PI*ang_ind)/(double)posConSides);
                        double sin_term = sin((double)(2.0*M_PI*ang_ind)/(double)posConSides);
                        for(int i = 0; i <= n; i++){  // iterate over degrees of polynomial segment
                            A_ind_x = seg_ind[ti]*(n+1) + i + flat_ind_x*(n+1)*m + con*numVars;
                            A_ind_y = seg_ind[ti]*(n+1) + i + flat_ind_y*(n+1)*m + con*numVars;
                            A_ind_z = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                            pos_coeff = pow(FOVtime[ti],i);
                            A[A_ind_x] = cos_term*atan(alpha)*pos_coeff;
                            A[A_ind_y] = sin_term*atan(alpha)*pos_coeff;
                            A[A_ind_z] = pos_coeff;
                        }
                        lbA[con] = l[2] + atan(alpha)*(l[0]*cos_term + l[1]*sin_term);
                        ubA[con] = INFINITY;
                        con += 1;
                    }
                    // create position constraint to agree with Theta_x values
                    for(int i = 0; i <= n; i++){  // iterate over degrees of polynomial segment
                        A_ind_x = seg_ind[ti]*(n+1) + i + flat_ind_x*(n+1)*m + con*numVars;
                        A_ind_z = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                        pos_coeff = pow(FOVtime[ti],i);
                        A[A_ind_x] = pos_coeff;
                        A[A_ind_z] = tan(Theta_x[ti])*pos_coeff;
                    }
                    lbA[con] = l[2]*tan(Theta_x[ti]) + l[0] - 0.25;
                    ubA[con] = l[2]*tan(Theta_x[ti]) + l[0] + 0.25;
                    con += 1;
                    // create position constraint to agree with Theta_y values
                    for(int i = 0; i <= n; i++){  // iterate over degrees of polynomial segment
                        A_ind_y = seg_ind[ti]*(n+1) + i + flat_ind_y*(n+1)*m + con*numVars;
                        A_ind_z = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                        pos_coeff = pow(FOVtime[ti],i);
                        A[A_ind_y] = pos_coeff;
                        A[A_ind_z] = tan(Theta_y[ti])*pos_coeff;
                    }
                    lbA[con] = l[2]*tan(Theta_y[ti]) + l[1] - 0.25;
                    ubA[con] = l[2]*tan(Theta_y[ti]) + l[1] + 0.25;
                    con += 1;
                    //*/
                    
                    // Acceleration constraints
                    tanx = tan(alpha_x - abs(Theta_x[ti]));
                    tany = tan(alpha_y - abs(Theta_y[ti]));
                    // create first ax az constraint
                    for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                        A_ind_ax = seg_ind[ti]*(n+1) + i + flat_ind_x*(n+1)*m + con*numVars;
                        A_ind_az = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                        acc_coeff = pow(FOVtime[ti],i-p)*fact(i)/fact(i-p);
                        A[A_ind_ax] = acc_coeff;
                        A[A_ind_az] = -acc_coeff*tanx;
                    }
                    A_ind_slack_x = (con + 1)*numVars - 2;
                    A[A_ind_slack_x] = -1;
                    lbA[con] = -1.0e4;
                    ubA[con] = grav * tanx;
                    con += 1;
                    
                    // create second ax az constraint
                    for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                        A_ind_ax = seg_ind[ti]*(n+1) + i + flat_ind_x*(n+1)*m + con*numVars;
                        A_ind_az = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                        acc_coeff = pow(FOVtime[ti],i-p)*fact(i)/fact(i-p);
                        A[A_ind_ax] = -acc_coeff;
                        A[A_ind_az] = -acc_coeff*tanx;
                    }
                    A_ind_slack_x = (con + 1)*numVars - 2;
                    A[A_ind_slack_x] = -1;
                    lbA[con] = -1.0e4;
                    ubA[con] = grav * tanx;
                    con += 1;
                    
                    // create first ay az constraint
                    for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                        A_ind_ay = seg_ind[ti]*(n+1) + i + flat_ind_y*(n+1)*m + con*numVars;
                        A_ind_az = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                        acc_coeff = pow(FOVtime[ti],i-p)*fact(i)/fact(i-p);
                        A[A_ind_ay] = acc_coeff;
                        A[A_ind_az] = -acc_coeff*tany;
                    }
                    A_ind_slack_y = (con + 1)*numVars - 1;
                    A[A_ind_slack_y] = -1;
                    lbA[con] = -1.0e4;
                    ubA[con] = grav * tany;
                    con += 1;
                    
                    // create second ay az constraint
                    for(int i = p; i <= n; i++){  // iterate over degrees of polynomial segment
                        A_ind_ay = seg_ind[ti]*(n+1) + i + flat_ind_y*(n+1)*m + con*numVars;
                        A_ind_az = seg_ind[ti]*(n+1) + i + flat_ind_z*(n+1)*m + con*numVars;
                        acc_coeff = pow(FOVtime[ti],i-p)*fact(i)/fact(i-p);
                        A[A_ind_ay] = -acc_coeff;
                        A[A_ind_az] = -acc_coeff*tany;
                    }
                    A_ind_slack_y = (con + 1)*numVars - 1;
                    A[A_ind_slack_y] = -1;
                    lbA[con] = -1.0e4;
                    ubA[con] = grav * tany;
                    con += 1;
                }
            }
        }
        
        double* buildGvector() {
            // -------- Graident Vector: Build g vector -------- //
            g = new double[numVars] (); // set g to zeros
            if(do_fov_constraint){
                g[numVars - 2] = 9000; // slack variable weight for ax az constraints
                g[numVars - 1] = 9000; // slack variable weight for ay az constraints
            }
            return g;
        }
        
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
        void printArrayAsMatrix(double *A, int rows, int cols) {
            int i, j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++)
                    printf("  % 9.2f", A[i*cols+j]);
                printf(",\n");
            }
            printf(",\n");
        }
        
        void printVectorOfVectors(const vector<vector<double>>& inputVector, int dec_places) {
            size_t innerVecSize = inputVector[0].size(); // Assumes all inner vectors have the same size
            // Print each inner vector as a row
            for (const auto& innerVec : inputVector) {
                for (size_t i = 0; i < innerVecSize; ++i) {
                    // Print the element with aligned decimals
                    cout << fixed << setprecision(dec_places) << setw(6 + dec_places) << innerVec[i];
                }
                cout << endl;
            }
            cout << endl;
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
        
        Eigen::Quaterniond acc2quaternion(const Eigen::Vector3d &vector_acc, double yaw) {
            Eigen::Quaterniond quat;
            Eigen::Vector3d zb_des, yb_des, xb_des, yc, grav;
            Eigen::Matrix3d rotmat;
            
            grav << 0.0, 0.0, 9.8;
            yc << -std::sin(yaw), std::cos(yaw), 0.0;

            zb_des = (vector_acc + grav) / (vector_acc + grav).norm();
            xb_des = yc.cross(zb_des) / (yc.cross(zb_des)).norm();
            yb_des = zb_des.cross(xb_des) / (zb_des.cross(xb_des)).norm();

            rotmat << xb_des(0), yb_des(0), zb_des(0), xb_des(1), yb_des(1), zb_des(1), xb_des(2), yb_des(2), zb_des(2);
            quat = Eigen::Quaterniond(rotmat);
            return quat; // access with quat.coeffs().w() etc.
        }
};

#endif // MIN_SNAP_TRAJ_H_INCLUDED
