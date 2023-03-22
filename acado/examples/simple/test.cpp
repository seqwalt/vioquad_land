/*
 *    This file was auto-generated using the ACADO Toolkit.
 *
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */



/*

IMPORTANT: This file should serve as a starting point to develop the user
code for the OCP solver. The code below is for illustration purposes. Most
likely you will not get good results if you execute this code without any
modification(s).

Please read the examples in order to understand how to write user code how
to run the OCP solver. You can find more info on the website:
www.acadotoolkit.org

*/

#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#include <stdio.h>
#include <iostream>
#include <random>
#include <assert.h>

/* Some convenient definitions. */
#define NX          ACADO_NX  /* Number of differential state variables.  */
#define NXA         ACADO_NXA /* Number of algebraic variables. */
#define NU          ACADO_NU  /* Number of control inputs. */
#define NOD         ACADO_NOD  /* Number of online data values. */

#define NY          ACADO_NY  /* Number of measurements/references on nodes 0..N - 1. */
#define NYN         ACADO_NYN /* Number of measurements/references on node N. */

#define N           ACADO_N   /* Number of intervals in the horizon. */

#define NUM_STEPS   10        /* Number of real-time iterations. */
#define VERBOSE     1         /* Show iterations: 1, silent: 0.  */

/* Global variables used by the solver. */
ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

void init_acadoVariables();
float lerp(float a, float b, float t);
float lerp(float a, float b, float c, float t);
double rand_num(double min, double max);

/* A template for testing of the solver. */
int main( )
{
	/* Some temporary variables. */
	int    i, iter, status;
	acado_timer t;

	/* Initialize the solver. */
	init_acadoVariables();
	status = acado_initializeSolver();
	if (status != 0){
		std::cout << "Initialize solver problem! QP status: " << status << std::endl;
	}

	if( VERBOSE ) acado_printHeader();

	/* Prepare first step */
	status = acado_preparationStep();
	if (status != 0){
		std::cout << " Preparation step problem! QP status: " << status << std::endl;
	}

	/* Get the time before start of the loop. */
	acado_tic( &t );

	/* The "real-time iterations" loop. */
	for(iter = 0; iter < NUM_STEPS; ++iter)
	{
        /* Perform the feedback step. */
		status = acado_feedbackStep( );

		if (status != 0){
			std::cout << "Iteration:" << iter << ", QP problem! QP status: " << status << std::endl;
			break;
		}

		/* Apply the new control immediately to the process, first NU components. */

		std::cout << "Real-Time Iteration " << iter << ", KKT Tolerance = " << std::scientific << acado_getKKT()
		<< " curr objective = " << std::scientific << acado_getObjective() << std::endl;

		/* Optional: shift the initialization (look at acado_common.h). */
		acado_shiftStates(2, 0, 0);
		acado_shiftControls( 0 );

		/* Prepare for the next step. */
		acado_preparationStep();
	}
	/* Read the elapsed time. */
	real_t te = acado_toc( &t );

	if( VERBOSE ) printf("\n\nEnd of the RTI loop. \n\n\n");

	/* Eye-candy. */

	if( !VERBOSE ) printf("\n\n Average time of one real-time iteration:   %.3g microseconds\n\n", 1e6 * te / NUM_STEPS);

	if (status == 0){
		acado_printDifferentialVariables();
		acado_printControlVariables();
	}

    return 0;
}

void init_acadoVariables(){
    int i;
    float s;
    // Initialize the states
		float xi, yi, zi, xf, yf, zf;	// initial and final position
		xi = 1.0f; yi = 0.5f; zi = 1.0f;
		xf = 0.0f; yf = 0.0f; zf = 2.0f;
		float vxi, vyi, vzi, vxf, vyf, vzf; // initial and final velocity
		vxi = 1.0f; vyi = 0.3f; vzi = -0.2f;
		vxf = 0.0f; vyf = 0.0f; vzf = 0.0f;
    for (i = 0; i < N + 1; ++i){
        s = (float)i/(float)N; // interp parameter in range [0,1]
				acadoVariables.x[i * NX + 0] = lerp(xi, xf, s); // x; 2 if s==0, 0 if s==1
        acadoVariables.x[i * NX + 1] = lerp(yi, yf, s); // y
        acadoVariables.x[i * NX + 2] = lerp(zi, zf, s); // z
        acadoVariables.x[i * NX + 3] = lerp(vxi, xf-xi, vxf, s); // vx
        acadoVariables.x[i * NX + 4] = lerp(vyi, yf-yi, vyf, s); // vy
        acadoVariables.x[i * NX + 5] = lerp(vzi, zf-zi, vzf, s); // vz
    }
    // Initialize the controls and references
    for (i = 0; i < N; ++i){
        // Initize the controls
        acadoVariables.u[i * NU + 0] = 9.0; // Thrust (mass-normalized)
        acadoVariables.u[i * NU + 1] = 0.0; // a_x
        acadoVariables.u[i * NU + 2] = 0.0; // a_y

        // Initialize the references
				acadoVariables.y[i * NY + 0] = xf;  //  x
        acadoVariables.y[i * NY + 1] = yf;  //  y
        acadoVariables.y[i * NY + 2] = zf;  //  z
        acadoVariables.y[i * NY + 3] = vxf; // vx
        acadoVariables.y[i * NY + 4] = vyf; // vy
        acadoVariables.y[i * NY + 5] = vzf; // vz
    }
    // Initialize the final reference
    acadoVariables.yN[0] = xf;  //  x
    acadoVariables.yN[1] = yf;  //  y
    acadoVariables.yN[2] = zf;  //  z
    acadoVariables.yN[3] = vxf; // vx
    acadoVariables.yN[4] = vyf; // vy
    acadoVariables.yN[5] = vzf; // vz

    /* MPC: initialize the current state feedback. */
#if ACADO_INITIAL_STATE_FIXED
		acadoVariables.x0[0] = xi; // x
		acadoVariables.x0[1] = yi; // y
		acadoVariables.x0[2] = zi; // z
		acadoVariables.x0[3] = vxi; // vx
		acadoVariables.x0[4] = vyi; // vy
		acadoVariables.x0[5] = vzi; // vz
#endif

		int n_sc = 9; // number of states and controls combined
		int n_s = 6; // number of states
		int blk_size = n_sc*n_sc;
		//for(size_t i = 0; i < N*blk_size; ++i) acadoVariables.WN[i] = 0.0;
		for(size_t i = 0; i < N; ++i){
			acadoVariables.W[i*blk_size] = 100; // x gain
			acadoVariables.W[i*blk_size + n_sc + 1] = 100; // y gain
			acadoVariables.W[i*blk_size + n_sc*2 + 2] = 100; // z gain
			acadoVariables.W[i*blk_size + n_sc*3 + 3] = 10; // v_x gain
			acadoVariables.W[i*blk_size + n_sc*4 + 4] = 10; // v_y gain
			acadoVariables.W[i*blk_size + n_sc*5 + 5] = 10; // v_z gain
			acadoVariables.W[i*blk_size + n_sc*6 + 6] = 1; // T gain
			acadoVariables.W[i*blk_size + n_sc*7 + 7] = 1; // a_x gain
			acadoVariables.W[i*blk_size + n_sc*8 + 8] = 1; // a_y gain
		}

		for(size_t i = 0; i < n_s*n_s; ++i) acadoVariables.WN[i] = 0.0;
		acadoVariables.WN[0] = 10000; // x gain
		acadoVariables.WN[n_s+1] = 10000; // y gain
		acadoVariables.WN[n_s*2 + 2] = 10000; // z gain
		acadoVariables.WN[n_s*3 + 3] = 1000; // vx gain
		acadoVariables.WN[n_s*4 + 4] = 1000; // vy gain
		acadoVariables.WN[n_s*5 + 5] = 1000; // vz gain

		int n_c = 3; // number of box constraints
		for(size_t i = 0; i < N; ++i){
			acadoVariables.lbValues[i*n_c] = 2.0; // T
			acadoVariables.ubValues[i*n_c] = 20.0;
			acadoVariables.lbValues[i*n_c + 1] = -3.0; // a_x
			acadoVariables.ubValues[i*n_c + 1] = 3.0;
			acadoVariables.lbValues[i*n_c + 2] = -1.0; // a_y
			acadoVariables.ubValues[i*n_c + 2] = 1.0;
		}
}

float lerp(float a, float b, float t){
		// linear interpolation between a and b
		assert(t >= 0.0f);
		assert(t <= 1.0f);
    return a + t * (b - a);
}
float lerp(float a, float b, float c, float t){
		// linear interpolation between a, b and c
		assert(t >= 0.0f);
		assert(t <= 1.0f);
		if (t <= 0.5f) return a + t * (b - a);
		else return b + t * (c - b);
}
double rand_num(double min, double max) {
		// Making rng static ensures that it stays the same
		// between different invocations of the function
		static std::default_random_engine rng;
		std::uniform_real_distribution<double> dist(min, max);
		return dist(rng);
}
