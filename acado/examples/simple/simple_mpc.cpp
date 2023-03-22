#include <acado_code_generation.hpp>

int main( ){
  // Use Acado
  USING_NAMESPACE_ACADO

  // System variables
  DifferentialState     p_x, p_y, p_z;
  DifferentialState     v_x, v_y, v_z;
  Control               T, a_x, a_y;
  DifferentialEquation  f;
  Function              h, hN;

  // Parameters with exemplary values. These are set/overwritten at runtime.
  const double t_start = 0.0;     // Initial time [s]
  const double t_end = 2.0;       // Time horizon [s]
  const double dt = 0.1;          // Discretization time [s]
  const int N = round(t_end/dt);  // Number of nodes
  const double g_z = 9.8066;      // Gravity is everywhere [m/s^2]
  const double T_min = 2;         // Minimal thrust [N]
  const double T_max = 20;        // Maximal thrust [N]

  // System Dynamics
  f << dot(p_x) ==  v_x;
  f << dot(p_y) ==  v_y;
  f << dot(p_z) ==  v_z;
  f << dot(v_x) ==  a_x;
  f << dot(v_y) ==  a_y;
  f << dot(v_z) ==  T - g_z;

  // Cost: Sum(i=0, ..., N-1){h_i' * Q * h_i} + h_N' * Q_N * h_N
  // Running cost vector consists of all states and inputs.
  h << p_x << p_y << p_z
    << v_x << v_y << v_z
    << T << a_x << a_y;

  // End cost vector consists of all states (no inputs at last state).
  hN << p_x << p_y << p_z
    << v_x << v_y << v_z;

  // Running cost weight matrix, defined at run-time
  BMatrix Q(h.getDim(), h.getDim());
  Q.setIdentity();

  // End cost weight matrix, defined at run-time
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
  ocp.subjectTo( T_min <= T <= T_max);
  ocp.subjectTo( -2.0 <= a_x <= 2.0);
  ocp.subjectTo( -2.0 <= a_y <= 2.0);

  // For code generation, we can set some properties.
  // The main reason for a setting is given as comment.
  OCPexport mpc(ocp);

  mpc.set(HESSIAN_APPROXIMATION,  GAUSS_NEWTON);        // is robust, stable
  mpc.set(DISCRETIZATION_TYPE,    MULTIPLE_SHOOTING);   // good convergence
  mpc.set(SPARSE_QP_SOLUTION,     FULL_CONDENSING_N2);  // due to qpOASES
  mpc.set(INTEGRATOR_TYPE,        INT_IRK_GL4);         // accurate
  mpc.set(NUM_INTEGRATOR_STEPS,   N);
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
