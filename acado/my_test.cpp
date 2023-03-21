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

#define NUM_STEPS   20        /* Number of real-time iterations. */
#define VERBOSE     1         /* Show iterations: 1, silent: 0.  */

/* Global variables used by the solver. */
ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

void init_acadoVariables();
float lerp(float a, float b, float t);
double rand_num();
double rand_num(double min, double max);

/* A template for testing of the solver. */
int main( )
{
	/* Some temporary variables. */
	int    i, iter;
	acado_timer t;

	/* Initialize the solver. */
	init_acadoVariables();
	acado_initializeSolver();

	if( VERBOSE ) acado_printHeader();

	/* Prepare first step */
	acado_preparationStep();

	/* Get the time before start of the loop. */
	acado_tic( &t );

	/* The "real-time iterations" loop. */
	for(iter = 0; iter < NUM_STEPS; ++iter)
	{
        /* Perform the feedback step. */
		acado_feedbackStep( );

		/* Apply the new control immediately to the process, first NU components. */

		if( VERBOSE ) printf("\tReal-Time Iteration %d:  KKT Tolerance = %.3e\n", iter, acado_getKKT() );

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

	if( !VERBOSE )
	printf("\n\n Average time of one real-time iteration:   %.3g microseconds\n\n", 1e6 * te / NUM_STEPS);

	acado_printDifferentialVariables();
	acado_printControlVariables();

    return 0;
}

void init_acadoVariables(){
    int i;
    float s;
    // Initialize the states
    for (i = 0; i < N + 1; ++i){
        s = (float)i/(float)N; // linear interp parameter in range [0,1]
        acadoVariables.x[i * NX + 0] = lerp(2.0f, 0.0f, s); // x; 2 if s==0, 0 if s==1
        acadoVariables.x[i * NX + 1] = lerp(2.0f, 0.0f, s); // y
        acadoVariables.x[i * NX + 2] = lerp(1.5f, 1.0f, s); // z
        acadoVariables.x[i * NX + 3] = 1.0; // qw
        acadoVariables.x[i * NX + 4] = 0.0; // qx
        acadoVariables.x[i * NX + 5] = 0.0; // qy
        acadoVariables.x[i * NX + 6] = 0.0; // qz
        acadoVariables.x[i * NX + 7] = 0.0 + rand_num()*(double)(i>0); // vx
        acadoVariables.x[i * NX + 8] = 0.0 + rand_num()*(double)(i>0); // vy
        acadoVariables.x[i * NX + 9] = 0.0 + rand_num()*(double)(i>0); // vz
    }
    // Initialize the controls and references
    for (i = 0; i < N; ++i){
        // Initize the controls
        acadoVariables.u[i * NU + 0] = 9.8 + rand_num(-6.0, 10); // Thrust (mass-normalized)
        acadoVariables.u[i * NU + 1] = 0.0 + rand_num(-3.0, 3.0); // w_x
        acadoVariables.u[i * NU + 2] = 0.0 + rand_num(-3.0, 3.0); // w_y
        acadoVariables.u[i * NU + 3] = 0.0 + rand_num(); // w_z

        // Initialize the references
        acadoVariables.y[i * NY + 0] = 0.0 + rand_num(); // x
        acadoVariables.y[i * NY + 1] = 0.0 + rand_num(); // y
        acadoVariables.y[i * NY + 2] = 1.0 + rand_num(); // z
        acadoVariables.y[i * NY + 3] = 1.0 + rand_num(); // qw
        acadoVariables.y[i * NY + 4] = 0.0 + rand_num(); // qx
        acadoVariables.y[i * NY + 5] = 0.0 + rand_num(); // qy
        acadoVariables.y[i * NY + 6] = 0.0 + rand_num(); // qz
        acadoVariables.y[i * NY + 7] = 0.0 + rand_num(); // vx
        acadoVariables.y[i * NY + 8] = 0.0 + rand_num(); // vy
        acadoVariables.y[i * NY + 9] = 0.0 + rand_num(); // vz
    }
    // Initialize the final reference
    acadoVariables.yN[0] = 0.0; // x
    acadoVariables.yN[1] = 0.0; // y
    acadoVariables.yN[2] = 1.0; // z
    acadoVariables.yN[3] = 1.0; // qw
    acadoVariables.yN[4] = 0.0; // qx
    acadoVariables.yN[5] = 0.0; // qy
    acadoVariables.yN[6] = 0.0; // qz
    acadoVariables.yN[7] = 0.0; // vx
    acadoVariables.yN[8] = 0.0; // vy
    acadoVariables.yN[9] = 0.0; // vz

    /* MPC: initialize the current state feedback. */
#if ACADO_INITIAL_STATE_FIXED
		acadoVariables.x0[0] = 2.0f; // x
		acadoVariables.x0[1] = 2.0f; // y
		acadoVariables.x0[2] = 1.5f; // z
		acadoVariables.x0[3] = 1.0; // qw
		acadoVariables.x0[4] = 0.0; // qx
		acadoVariables.x0[5] = 0.0; // qy
		acadoVariables.x0[6] = 0.0; // qz
		acadoVariables.x0[7] = 0.0; // vx
		acadoVariables.x0[8] = 0.0; // vy
		acadoVariables.x0[9] = 0.0; // vz
#endif
}

float lerp(float a, float b, float t){
		assert(t >= 0.0f);
		assert(t <= 1.0f);
    return a + t * (b - a);
}
double rand_num() {
		// Making rng static ensures that it stays the same
		// between different invocations of the function
		static std::default_random_engine rng;
		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		return dist(rng);
}
double rand_num(double min, double max) {
		// Making rng static ensures that it stays the same
		// between different invocations of the function
		static std::default_random_engine rng;
		std::uniform_real_distribution<double> dist(min, max);
		return dist(rng);
}
