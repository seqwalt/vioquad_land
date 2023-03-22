/*    rpg_quadrotor_mpc
 *    A model predictive control implementation for quadrotors.
 *    Copyright (C) 2017-2018 Philipp Foehn,
 *    Robotics and Perception Group, University of Zurich
 *
 *    Intended to be used with rpg_quadrotor_control and rpg_quadrotor_common.
 *    https://github.com/uzh-rpg/rpg_quadrotor_control
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <memory>
#include <acado_optimal_control.hpp>
#include <acado_code_generation.hpp>

// Standalone code generation for a parameter-free quadrotor model
// with thrust and rates input.

int main( ){
  // Use Acado
  USING_NAMESPACE_ACADO

  // System variables
  DifferentialState     p_x, p_y, p_z;
  DifferentialState     q_w, q_x, q_y, q_z;
  DifferentialState     v_x, v_y, v_z;
  Control               T, w_x, w_y, w_z;
  DifferentialEquation  f;
  Function              h, hN;

  // Parameters with exemplary values. These are set/overwritten at runtime.
  const double t_start = 0.0;     // Initial time [s]
  const double t_end = 5.0;       // Time horizon [s]
  const double dt = 0.1;          // Discretization time [s]
  const int N = round(t_end/dt);  // Number of nodes
  const double g_z = 9.8066;      // Gravity is everywhere [m/s^2]
  const double w_max_yaw = 1;     // Maximal yaw rate [rad/s]
  const double w_max_xy = 3;      // Maximal pitch and roll rate [rad/s]
  const double T_min = 2;         // Minimal thrust [N]
  const double T_max = 20;        // Maximal thrust [N]

  // System Dynamics
  f << dot(p_x) ==  v_x;
  f << dot(p_y) ==  v_y;
  f << dot(p_z) ==  v_z;
  f << dot(q_w) ==  0.5 * ( - w_x * q_x - w_y * q_y - w_z * q_z);
  f << dot(q_x) ==  0.5 * ( w_x * q_w + w_z * q_y - w_y * q_z);
  f << dot(q_y) ==  0.5 * ( w_y * q_w - w_z * q_x + w_x * q_z);
  f << dot(q_z) ==  0.5 * ( w_z * q_w + w_y * q_x - w_x * q_y);
  f << dot(v_x) ==  2 * ( q_w * q_y + q_x * q_z ) * T;
  f << dot(v_y) ==  2 * ( q_y * q_z - q_w * q_x ) * T;
  f << dot(v_z) ==  ( 1 - 2 * q_x * q_x - 2 * q_y * q_y ) * T - g_z;

  // Cost: Sum(i=0, ..., N-1){h_i' * Q * h_i} + h_N' * Q_N * h_N
  // Running cost vector consists of all states and inputs.
  h << p_x << p_y << p_z
    << q_w << q_x << q_y << q_z
    << v_x << v_y << v_z
    << T << w_x << w_y << w_z;

  // End cost vector consists of all states (no inputs at last state).
  hN << p_x << p_y << p_z
    << q_w << q_x << q_y << q_z
    << v_x << v_y << v_z;

  // Running cost weight matrix, set at run-time
  BMatrix Q(h.getDim(), h.getDim());
  Q.setIdentity();

  // End cost weight matrix, set at run-time
  BMatrix QN(hN.getDim(), hN.getDim());
  QN.setIdentity();

  // DEFINE AN OPTIMAL CONTROL PROBLEM:
  // ----------------------------------
  OCP ocp(t_start, t_end, N);
  ocp.minimizeLSQ(Q, h);
  ocp.minimizeLSQEndTerm(QN, hN);

  // Add system dynamics
  ocp.subjectTo( f );
  // Add constraints
  ocp.subjectTo(-w_max_xy <= w_x <= w_max_xy);
  ocp.subjectTo(-w_max_xy <= w_y <= w_max_xy);
  ocp.subjectTo(-w_max_yaw <= w_z <= w_max_yaw);
  ocp.subjectTo( T_min <= T <= T_max);

  // For code generation, we can set some properties.
  // The main reason for a setting is given as comment.
  OCPexport mpc(ocp);

  mpc.set(HESSIAN_APPROXIMATION,  GAUSS_NEWTON);        // is robust, stable
  mpc.set(DISCRETIZATION_TYPE,    MULTIPLE_SHOOTING);   // good convergence
  mpc.set(SPARSE_QP_SOLUTION,     FULL_CONDENSING_N2);  // due to qpOASES
  mpc.set(INTEGRATOR_TYPE,        INT_IRK_GL4);         // accurate
  //mpc.set(NUM_INTEGRATOR_STEPS,   N);			// original
  mpc.set(NUM_INTEGRATOR_STEPS,   20);
  mpc.set(QP_SOLVER,              QP_QPOASES);          // free, source code
  mpc.set(HOTSTART_QP,            YES);
  mpc.set(CG_USE_OPENMP,                    YES);       // paralellization
  mpc.set(CG_HARDCODE_CONSTRAINT_VALUES,    NO);        // set on runtime
  mpc.set(CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);       // time-varying costs
  mpc.set(USE_SINGLE_PRECISION,        YES);            // Single precision

  // Do not generate tests, makes or matlab-related interfaces.
  mpc.set(GENERATE_TEST_FILE,          NO);
  mpc.set(GENERATE_MAKE_FILE,          YES);
  mpc.set(GENERATE_MATLAB_INTERFACE,   NO);
  mpc.set(GENERATE_SIMULINK_INTERFACE, NO);

  // Finally, export everything.
  if(mpc.exportCode("autogen_code") != SUCCESSFUL_RETURN)
    exit( EXIT_FAILURE );
  mpc.printDimensionsQP( );


  return EXIT_SUCCESS;
}
