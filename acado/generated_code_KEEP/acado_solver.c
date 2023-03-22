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


#include "acado_common.h"




/******************************************************************************/
/*                                                                            */
/* ACADO code generation                                                      */
/*                                                                            */
/******************************************************************************/


/** Row vector of size: 63 */
real_t state[ 63 ];

int acado_modelSimulation(  )
{
int ret;

int lRun1;
ret = 0;
#pragma omp parallel for private(lRun1, state) shared(acadoWorkspace, acadoVariables)
for (lRun1 = 0; lRun1 < 20; ++lRun1)
{
state[0] = acadoVariables.x[lRun1 * 6];
state[1] = acadoVariables.x[lRun1 * 6 + 1];
state[2] = acadoVariables.x[lRun1 * 6 + 2];
state[3] = acadoVariables.x[lRun1 * 6 + 3];
state[4] = acadoVariables.x[lRun1 * 6 + 4];
state[5] = acadoVariables.x[lRun1 * 6 + 5];

state[60] = acadoVariables.u[lRun1 * 3];
state[61] = acadoVariables.u[lRun1 * 3 + 1];
state[62] = acadoVariables.u[lRun1 * 3 + 2];

ret = acado_integrate(state, 1);

acadoWorkspace.d[lRun1 * 6] = state[0] - acadoVariables.x[lRun1 * 6 + 6];
acadoWorkspace.d[lRun1 * 6 + 1] = state[1] - acadoVariables.x[lRun1 * 6 + 7];
acadoWorkspace.d[lRun1 * 6 + 2] = state[2] - acadoVariables.x[lRun1 * 6 + 8];
acadoWorkspace.d[lRun1 * 6 + 3] = state[3] - acadoVariables.x[lRun1 * 6 + 9];
acadoWorkspace.d[lRun1 * 6 + 4] = state[4] - acadoVariables.x[lRun1 * 6 + 10];
acadoWorkspace.d[lRun1 * 6 + 5] = state[5] - acadoVariables.x[lRun1 * 6 + 11];

acadoWorkspace.evGx[lRun1 * 36] = state[6];
acadoWorkspace.evGx[lRun1 * 36 + 1] = state[7];
acadoWorkspace.evGx[lRun1 * 36 + 2] = state[8];
acadoWorkspace.evGx[lRun1 * 36 + 3] = state[9];
acadoWorkspace.evGx[lRun1 * 36 + 4] = state[10];
acadoWorkspace.evGx[lRun1 * 36 + 5] = state[11];
acadoWorkspace.evGx[lRun1 * 36 + 6] = state[12];
acadoWorkspace.evGx[lRun1 * 36 + 7] = state[13];
acadoWorkspace.evGx[lRun1 * 36 + 8] = state[14];
acadoWorkspace.evGx[lRun1 * 36 + 9] = state[15];
acadoWorkspace.evGx[lRun1 * 36 + 10] = state[16];
acadoWorkspace.evGx[lRun1 * 36 + 11] = state[17];
acadoWorkspace.evGx[lRun1 * 36 + 12] = state[18];
acadoWorkspace.evGx[lRun1 * 36 + 13] = state[19];
acadoWorkspace.evGx[lRun1 * 36 + 14] = state[20];
acadoWorkspace.evGx[lRun1 * 36 + 15] = state[21];
acadoWorkspace.evGx[lRun1 * 36 + 16] = state[22];
acadoWorkspace.evGx[lRun1 * 36 + 17] = state[23];
acadoWorkspace.evGx[lRun1 * 36 + 18] = state[24];
acadoWorkspace.evGx[lRun1 * 36 + 19] = state[25];
acadoWorkspace.evGx[lRun1 * 36 + 20] = state[26];
acadoWorkspace.evGx[lRun1 * 36 + 21] = state[27];
acadoWorkspace.evGx[lRun1 * 36 + 22] = state[28];
acadoWorkspace.evGx[lRun1 * 36 + 23] = state[29];
acadoWorkspace.evGx[lRun1 * 36 + 24] = state[30];
acadoWorkspace.evGx[lRun1 * 36 + 25] = state[31];
acadoWorkspace.evGx[lRun1 * 36 + 26] = state[32];
acadoWorkspace.evGx[lRun1 * 36 + 27] = state[33];
acadoWorkspace.evGx[lRun1 * 36 + 28] = state[34];
acadoWorkspace.evGx[lRun1 * 36 + 29] = state[35];
acadoWorkspace.evGx[lRun1 * 36 + 30] = state[36];
acadoWorkspace.evGx[lRun1 * 36 + 31] = state[37];
acadoWorkspace.evGx[lRun1 * 36 + 32] = state[38];
acadoWorkspace.evGx[lRun1 * 36 + 33] = state[39];
acadoWorkspace.evGx[lRun1 * 36 + 34] = state[40];
acadoWorkspace.evGx[lRun1 * 36 + 35] = state[41];

acadoWorkspace.evGu[lRun1 * 18] = state[42];
acadoWorkspace.evGu[lRun1 * 18 + 1] = state[43];
acadoWorkspace.evGu[lRun1 * 18 + 2] = state[44];
acadoWorkspace.evGu[lRun1 * 18 + 3] = state[45];
acadoWorkspace.evGu[lRun1 * 18 + 4] = state[46];
acadoWorkspace.evGu[lRun1 * 18 + 5] = state[47];
acadoWorkspace.evGu[lRun1 * 18 + 6] = state[48];
acadoWorkspace.evGu[lRun1 * 18 + 7] = state[49];
acadoWorkspace.evGu[lRun1 * 18 + 8] = state[50];
acadoWorkspace.evGu[lRun1 * 18 + 9] = state[51];
acadoWorkspace.evGu[lRun1 * 18 + 10] = state[52];
acadoWorkspace.evGu[lRun1 * 18 + 11] = state[53];
acadoWorkspace.evGu[lRun1 * 18 + 12] = state[54];
acadoWorkspace.evGu[lRun1 * 18 + 13] = state[55];
acadoWorkspace.evGu[lRun1 * 18 + 14] = state[56];
acadoWorkspace.evGu[lRun1 * 18 + 15] = state[57];
acadoWorkspace.evGu[lRun1 * 18 + 16] = state[58];
acadoWorkspace.evGu[lRun1 * 18 + 17] = state[59];
}
return ret;
}

void acado_evaluateLSQ(const real_t* in, real_t* out)
{
const real_t* xd = in;
const real_t* u = in + 6;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = xd[4];
out[5] = xd[5];
out[6] = u[0];
out[7] = u[1];
out[8] = u[2];
}

void acado_evaluateLSQEndTerm(const real_t* in, real_t* out)
{
const real_t* xd = in;

/* Compute outputs: */
out[0] = xd[0];
out[1] = xd[1];
out[2] = xd[2];
out[3] = xd[3];
out[4] = xd[4];
out[5] = xd[5];
}

void acado_setObjQ1Q2( real_t* const tmpObjS, real_t* const tmpQ1, real_t* const tmpQ2 )
{
tmpQ2[0] = +tmpObjS[0];
tmpQ2[1] = +tmpObjS[1];
tmpQ2[2] = +tmpObjS[2];
tmpQ2[3] = +tmpObjS[3];
tmpQ2[4] = +tmpObjS[4];
tmpQ2[5] = +tmpObjS[5];
tmpQ2[6] = +tmpObjS[6];
tmpQ2[7] = +tmpObjS[7];
tmpQ2[8] = +tmpObjS[8];
tmpQ2[9] = +tmpObjS[9];
tmpQ2[10] = +tmpObjS[10];
tmpQ2[11] = +tmpObjS[11];
tmpQ2[12] = +tmpObjS[12];
tmpQ2[13] = +tmpObjS[13];
tmpQ2[14] = +tmpObjS[14];
tmpQ2[15] = +tmpObjS[15];
tmpQ2[16] = +tmpObjS[16];
tmpQ2[17] = +tmpObjS[17];
tmpQ2[18] = +tmpObjS[18];
tmpQ2[19] = +tmpObjS[19];
tmpQ2[20] = +tmpObjS[20];
tmpQ2[21] = +tmpObjS[21];
tmpQ2[22] = +tmpObjS[22];
tmpQ2[23] = +tmpObjS[23];
tmpQ2[24] = +tmpObjS[24];
tmpQ2[25] = +tmpObjS[25];
tmpQ2[26] = +tmpObjS[26];
tmpQ2[27] = +tmpObjS[27];
tmpQ2[28] = +tmpObjS[28];
tmpQ2[29] = +tmpObjS[29];
tmpQ2[30] = +tmpObjS[30];
tmpQ2[31] = +tmpObjS[31];
tmpQ2[32] = +tmpObjS[32];
tmpQ2[33] = +tmpObjS[33];
tmpQ2[34] = +tmpObjS[34];
tmpQ2[35] = +tmpObjS[35];
tmpQ2[36] = +tmpObjS[36];
tmpQ2[37] = +tmpObjS[37];
tmpQ2[38] = +tmpObjS[38];
tmpQ2[39] = +tmpObjS[39];
tmpQ2[40] = +tmpObjS[40];
tmpQ2[41] = +tmpObjS[41];
tmpQ2[42] = +tmpObjS[42];
tmpQ2[43] = +tmpObjS[43];
tmpQ2[44] = +tmpObjS[44];
tmpQ2[45] = +tmpObjS[45];
tmpQ2[46] = +tmpObjS[46];
tmpQ2[47] = +tmpObjS[47];
tmpQ2[48] = +tmpObjS[48];
tmpQ2[49] = +tmpObjS[49];
tmpQ2[50] = +tmpObjS[50];
tmpQ2[51] = +tmpObjS[51];
tmpQ2[52] = +tmpObjS[52];
tmpQ2[53] = +tmpObjS[53];
tmpQ1[0] = + tmpQ2[0];
tmpQ1[1] = + tmpQ2[1];
tmpQ1[2] = + tmpQ2[2];
tmpQ1[3] = + tmpQ2[3];
tmpQ1[4] = + tmpQ2[4];
tmpQ1[5] = + tmpQ2[5];
tmpQ1[6] = + tmpQ2[9];
tmpQ1[7] = + tmpQ2[10];
tmpQ1[8] = + tmpQ2[11];
tmpQ1[9] = + tmpQ2[12];
tmpQ1[10] = + tmpQ2[13];
tmpQ1[11] = + tmpQ2[14];
tmpQ1[12] = + tmpQ2[18];
tmpQ1[13] = + tmpQ2[19];
tmpQ1[14] = + tmpQ2[20];
tmpQ1[15] = + tmpQ2[21];
tmpQ1[16] = + tmpQ2[22];
tmpQ1[17] = + tmpQ2[23];
tmpQ1[18] = + tmpQ2[27];
tmpQ1[19] = + tmpQ2[28];
tmpQ1[20] = + tmpQ2[29];
tmpQ1[21] = + tmpQ2[30];
tmpQ1[22] = + tmpQ2[31];
tmpQ1[23] = + tmpQ2[32];
tmpQ1[24] = + tmpQ2[36];
tmpQ1[25] = + tmpQ2[37];
tmpQ1[26] = + tmpQ2[38];
tmpQ1[27] = + tmpQ2[39];
tmpQ1[28] = + tmpQ2[40];
tmpQ1[29] = + tmpQ2[41];
tmpQ1[30] = + tmpQ2[45];
tmpQ1[31] = + tmpQ2[46];
tmpQ1[32] = + tmpQ2[47];
tmpQ1[33] = + tmpQ2[48];
tmpQ1[34] = + tmpQ2[49];
tmpQ1[35] = + tmpQ2[50];
}

void acado_setObjR1R2( real_t* const tmpObjS, real_t* const tmpR1, real_t* const tmpR2 )
{
tmpR2[0] = +tmpObjS[54];
tmpR2[1] = +tmpObjS[55];
tmpR2[2] = +tmpObjS[56];
tmpR2[3] = +tmpObjS[57];
tmpR2[4] = +tmpObjS[58];
tmpR2[5] = +tmpObjS[59];
tmpR2[6] = +tmpObjS[60];
tmpR2[7] = +tmpObjS[61];
tmpR2[8] = +tmpObjS[62];
tmpR2[9] = +tmpObjS[63];
tmpR2[10] = +tmpObjS[64];
tmpR2[11] = +tmpObjS[65];
tmpR2[12] = +tmpObjS[66];
tmpR2[13] = +tmpObjS[67];
tmpR2[14] = +tmpObjS[68];
tmpR2[15] = +tmpObjS[69];
tmpR2[16] = +tmpObjS[70];
tmpR2[17] = +tmpObjS[71];
tmpR2[18] = +tmpObjS[72];
tmpR2[19] = +tmpObjS[73];
tmpR2[20] = +tmpObjS[74];
tmpR2[21] = +tmpObjS[75];
tmpR2[22] = +tmpObjS[76];
tmpR2[23] = +tmpObjS[77];
tmpR2[24] = +tmpObjS[78];
tmpR2[25] = +tmpObjS[79];
tmpR2[26] = +tmpObjS[80];
tmpR1[0] = + tmpR2[6];
tmpR1[1] = + tmpR2[7];
tmpR1[2] = + tmpR2[8];
tmpR1[3] = + tmpR2[15];
tmpR1[4] = + tmpR2[16];
tmpR1[5] = + tmpR2[17];
tmpR1[6] = + tmpR2[24];
tmpR1[7] = + tmpR2[25];
tmpR1[8] = + tmpR2[26];
}

void acado_setObjQN1QN2( real_t* const tmpObjSEndTerm, real_t* const tmpQN1, real_t* const tmpQN2 )
{
tmpQN2[0] = +tmpObjSEndTerm[0];
tmpQN2[1] = +tmpObjSEndTerm[1];
tmpQN2[2] = +tmpObjSEndTerm[2];
tmpQN2[3] = +tmpObjSEndTerm[3];
tmpQN2[4] = +tmpObjSEndTerm[4];
tmpQN2[5] = +tmpObjSEndTerm[5];
tmpQN2[6] = +tmpObjSEndTerm[6];
tmpQN2[7] = +tmpObjSEndTerm[7];
tmpQN2[8] = +tmpObjSEndTerm[8];
tmpQN2[9] = +tmpObjSEndTerm[9];
tmpQN2[10] = +tmpObjSEndTerm[10];
tmpQN2[11] = +tmpObjSEndTerm[11];
tmpQN2[12] = +tmpObjSEndTerm[12];
tmpQN2[13] = +tmpObjSEndTerm[13];
tmpQN2[14] = +tmpObjSEndTerm[14];
tmpQN2[15] = +tmpObjSEndTerm[15];
tmpQN2[16] = +tmpObjSEndTerm[16];
tmpQN2[17] = +tmpObjSEndTerm[17];
tmpQN2[18] = +tmpObjSEndTerm[18];
tmpQN2[19] = +tmpObjSEndTerm[19];
tmpQN2[20] = +tmpObjSEndTerm[20];
tmpQN2[21] = +tmpObjSEndTerm[21];
tmpQN2[22] = +tmpObjSEndTerm[22];
tmpQN2[23] = +tmpObjSEndTerm[23];
tmpQN2[24] = +tmpObjSEndTerm[24];
tmpQN2[25] = +tmpObjSEndTerm[25];
tmpQN2[26] = +tmpObjSEndTerm[26];
tmpQN2[27] = +tmpObjSEndTerm[27];
tmpQN2[28] = +tmpObjSEndTerm[28];
tmpQN2[29] = +tmpObjSEndTerm[29];
tmpQN2[30] = +tmpObjSEndTerm[30];
tmpQN2[31] = +tmpObjSEndTerm[31];
tmpQN2[32] = +tmpObjSEndTerm[32];
tmpQN2[33] = +tmpObjSEndTerm[33];
tmpQN2[34] = +tmpObjSEndTerm[34];
tmpQN2[35] = +tmpObjSEndTerm[35];
tmpQN1[0] = + tmpQN2[0];
tmpQN1[1] = + tmpQN2[1];
tmpQN1[2] = + tmpQN2[2];
tmpQN1[3] = + tmpQN2[3];
tmpQN1[4] = + tmpQN2[4];
tmpQN1[5] = + tmpQN2[5];
tmpQN1[6] = + tmpQN2[6];
tmpQN1[7] = + tmpQN2[7];
tmpQN1[8] = + tmpQN2[8];
tmpQN1[9] = + tmpQN2[9];
tmpQN1[10] = + tmpQN2[10];
tmpQN1[11] = + tmpQN2[11];
tmpQN1[12] = + tmpQN2[12];
tmpQN1[13] = + tmpQN2[13];
tmpQN1[14] = + tmpQN2[14];
tmpQN1[15] = + tmpQN2[15];
tmpQN1[16] = + tmpQN2[16];
tmpQN1[17] = + tmpQN2[17];
tmpQN1[18] = + tmpQN2[18];
tmpQN1[19] = + tmpQN2[19];
tmpQN1[20] = + tmpQN2[20];
tmpQN1[21] = + tmpQN2[21];
tmpQN1[22] = + tmpQN2[22];
tmpQN1[23] = + tmpQN2[23];
tmpQN1[24] = + tmpQN2[24];
tmpQN1[25] = + tmpQN2[25];
tmpQN1[26] = + tmpQN2[26];
tmpQN1[27] = + tmpQN2[27];
tmpQN1[28] = + tmpQN2[28];
tmpQN1[29] = + tmpQN2[29];
tmpQN1[30] = + tmpQN2[30];
tmpQN1[31] = + tmpQN2[31];
tmpQN1[32] = + tmpQN2[32];
tmpQN1[33] = + tmpQN2[33];
tmpQN1[34] = + tmpQN2[34];
tmpQN1[35] = + tmpQN2[35];
}

void acado_evaluateObjective(  )
{
int runObj;
for (runObj = 0; runObj < 20; ++runObj)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[runObj * 6];
acadoWorkspace.objValueIn[1] = acadoVariables.x[runObj * 6 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[runObj * 6 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[runObj * 6 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.x[runObj * 6 + 4];
acadoWorkspace.objValueIn[5] = acadoVariables.x[runObj * 6 + 5];
acadoWorkspace.objValueIn[6] = acadoVariables.u[runObj * 3];
acadoWorkspace.objValueIn[7] = acadoVariables.u[runObj * 3 + 1];
acadoWorkspace.objValueIn[8] = acadoVariables.u[runObj * 3 + 2];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[runObj * 9] = acadoWorkspace.objValueOut[0];
acadoWorkspace.Dy[runObj * 9 + 1] = acadoWorkspace.objValueOut[1];
acadoWorkspace.Dy[runObj * 9 + 2] = acadoWorkspace.objValueOut[2];
acadoWorkspace.Dy[runObj * 9 + 3] = acadoWorkspace.objValueOut[3];
acadoWorkspace.Dy[runObj * 9 + 4] = acadoWorkspace.objValueOut[4];
acadoWorkspace.Dy[runObj * 9 + 5] = acadoWorkspace.objValueOut[5];
acadoWorkspace.Dy[runObj * 9 + 6] = acadoWorkspace.objValueOut[6];
acadoWorkspace.Dy[runObj * 9 + 7] = acadoWorkspace.objValueOut[7];
acadoWorkspace.Dy[runObj * 9 + 8] = acadoWorkspace.objValueOut[8];

acado_setObjQ1Q2( &(acadoVariables.W[ runObj * 81 ]), &(acadoWorkspace.Q1[ runObj * 36 ]), &(acadoWorkspace.Q2[ runObj * 54 ]) );

acado_setObjR1R2( &(acadoVariables.W[ runObj * 81 ]), &(acadoWorkspace.R1[ runObj * 9 ]), &(acadoWorkspace.R2[ runObj * 27 ]) );

}
acadoWorkspace.objValueIn[0] = acadoVariables.x[120];
acadoWorkspace.objValueIn[1] = acadoVariables.x[121];
acadoWorkspace.objValueIn[2] = acadoVariables.x[122];
acadoWorkspace.objValueIn[3] = acadoVariables.x[123];
acadoWorkspace.objValueIn[4] = acadoVariables.x[124];
acadoWorkspace.objValueIn[5] = acadoVariables.x[125];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );

acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3];
acadoWorkspace.DyN[4] = acadoWorkspace.objValueOut[4];
acadoWorkspace.DyN[5] = acadoWorkspace.objValueOut[5];

acado_setObjQN1QN2( acadoVariables.WN, acadoWorkspace.QN1, acadoWorkspace.QN2 );

}

void acado_multGxd( real_t* const dOld, real_t* const Gx1, real_t* const dNew )
{
dNew[0] += + Gx1[0]*dOld[0] + Gx1[1]*dOld[1] + Gx1[2]*dOld[2] + Gx1[3]*dOld[3] + Gx1[4]*dOld[4] + Gx1[5]*dOld[5];
dNew[1] += + Gx1[6]*dOld[0] + Gx1[7]*dOld[1] + Gx1[8]*dOld[2] + Gx1[9]*dOld[3] + Gx1[10]*dOld[4] + Gx1[11]*dOld[5];
dNew[2] += + Gx1[12]*dOld[0] + Gx1[13]*dOld[1] + Gx1[14]*dOld[2] + Gx1[15]*dOld[3] + Gx1[16]*dOld[4] + Gx1[17]*dOld[5];
dNew[3] += + Gx1[18]*dOld[0] + Gx1[19]*dOld[1] + Gx1[20]*dOld[2] + Gx1[21]*dOld[3] + Gx1[22]*dOld[4] + Gx1[23]*dOld[5];
dNew[4] += + Gx1[24]*dOld[0] + Gx1[25]*dOld[1] + Gx1[26]*dOld[2] + Gx1[27]*dOld[3] + Gx1[28]*dOld[4] + Gx1[29]*dOld[5];
dNew[5] += + Gx1[30]*dOld[0] + Gx1[31]*dOld[1] + Gx1[32]*dOld[2] + Gx1[33]*dOld[3] + Gx1[34]*dOld[4] + Gx1[35]*dOld[5];
}

void acado_moveGxT( real_t* const Gx1, real_t* const Gx2 )
{
Gx2[0] = Gx1[0];
Gx2[1] = Gx1[1];
Gx2[2] = Gx1[2];
Gx2[3] = Gx1[3];
Gx2[4] = Gx1[4];
Gx2[5] = Gx1[5];
Gx2[6] = Gx1[6];
Gx2[7] = Gx1[7];
Gx2[8] = Gx1[8];
Gx2[9] = Gx1[9];
Gx2[10] = Gx1[10];
Gx2[11] = Gx1[11];
Gx2[12] = Gx1[12];
Gx2[13] = Gx1[13];
Gx2[14] = Gx1[14];
Gx2[15] = Gx1[15];
Gx2[16] = Gx1[16];
Gx2[17] = Gx1[17];
Gx2[18] = Gx1[18];
Gx2[19] = Gx1[19];
Gx2[20] = Gx1[20];
Gx2[21] = Gx1[21];
Gx2[22] = Gx1[22];
Gx2[23] = Gx1[23];
Gx2[24] = Gx1[24];
Gx2[25] = Gx1[25];
Gx2[26] = Gx1[26];
Gx2[27] = Gx1[27];
Gx2[28] = Gx1[28];
Gx2[29] = Gx1[29];
Gx2[30] = Gx1[30];
Gx2[31] = Gx1[31];
Gx2[32] = Gx1[32];
Gx2[33] = Gx1[33];
Gx2[34] = Gx1[34];
Gx2[35] = Gx1[35];
}

void acado_multGxGx( real_t* const Gx1, real_t* const Gx2, real_t* const Gx3 )
{
Gx3[0] = + Gx1[0]*Gx2[0] + Gx1[1]*Gx2[6] + Gx1[2]*Gx2[12] + Gx1[3]*Gx2[18] + Gx1[4]*Gx2[24] + Gx1[5]*Gx2[30];
Gx3[1] = + Gx1[0]*Gx2[1] + Gx1[1]*Gx2[7] + Gx1[2]*Gx2[13] + Gx1[3]*Gx2[19] + Gx1[4]*Gx2[25] + Gx1[5]*Gx2[31];
Gx3[2] = + Gx1[0]*Gx2[2] + Gx1[1]*Gx2[8] + Gx1[2]*Gx2[14] + Gx1[3]*Gx2[20] + Gx1[4]*Gx2[26] + Gx1[5]*Gx2[32];
Gx3[3] = + Gx1[0]*Gx2[3] + Gx1[1]*Gx2[9] + Gx1[2]*Gx2[15] + Gx1[3]*Gx2[21] + Gx1[4]*Gx2[27] + Gx1[5]*Gx2[33];
Gx3[4] = + Gx1[0]*Gx2[4] + Gx1[1]*Gx2[10] + Gx1[2]*Gx2[16] + Gx1[3]*Gx2[22] + Gx1[4]*Gx2[28] + Gx1[5]*Gx2[34];
Gx3[5] = + Gx1[0]*Gx2[5] + Gx1[1]*Gx2[11] + Gx1[2]*Gx2[17] + Gx1[3]*Gx2[23] + Gx1[4]*Gx2[29] + Gx1[5]*Gx2[35];
Gx3[6] = + Gx1[6]*Gx2[0] + Gx1[7]*Gx2[6] + Gx1[8]*Gx2[12] + Gx1[9]*Gx2[18] + Gx1[10]*Gx2[24] + Gx1[11]*Gx2[30];
Gx3[7] = + Gx1[6]*Gx2[1] + Gx1[7]*Gx2[7] + Gx1[8]*Gx2[13] + Gx1[9]*Gx2[19] + Gx1[10]*Gx2[25] + Gx1[11]*Gx2[31];
Gx3[8] = + Gx1[6]*Gx2[2] + Gx1[7]*Gx2[8] + Gx1[8]*Gx2[14] + Gx1[9]*Gx2[20] + Gx1[10]*Gx2[26] + Gx1[11]*Gx2[32];
Gx3[9] = + Gx1[6]*Gx2[3] + Gx1[7]*Gx2[9] + Gx1[8]*Gx2[15] + Gx1[9]*Gx2[21] + Gx1[10]*Gx2[27] + Gx1[11]*Gx2[33];
Gx3[10] = + Gx1[6]*Gx2[4] + Gx1[7]*Gx2[10] + Gx1[8]*Gx2[16] + Gx1[9]*Gx2[22] + Gx1[10]*Gx2[28] + Gx1[11]*Gx2[34];
Gx3[11] = + Gx1[6]*Gx2[5] + Gx1[7]*Gx2[11] + Gx1[8]*Gx2[17] + Gx1[9]*Gx2[23] + Gx1[10]*Gx2[29] + Gx1[11]*Gx2[35];
Gx3[12] = + Gx1[12]*Gx2[0] + Gx1[13]*Gx2[6] + Gx1[14]*Gx2[12] + Gx1[15]*Gx2[18] + Gx1[16]*Gx2[24] + Gx1[17]*Gx2[30];
Gx3[13] = + Gx1[12]*Gx2[1] + Gx1[13]*Gx2[7] + Gx1[14]*Gx2[13] + Gx1[15]*Gx2[19] + Gx1[16]*Gx2[25] + Gx1[17]*Gx2[31];
Gx3[14] = + Gx1[12]*Gx2[2] + Gx1[13]*Gx2[8] + Gx1[14]*Gx2[14] + Gx1[15]*Gx2[20] + Gx1[16]*Gx2[26] + Gx1[17]*Gx2[32];
Gx3[15] = + Gx1[12]*Gx2[3] + Gx1[13]*Gx2[9] + Gx1[14]*Gx2[15] + Gx1[15]*Gx2[21] + Gx1[16]*Gx2[27] + Gx1[17]*Gx2[33];
Gx3[16] = + Gx1[12]*Gx2[4] + Gx1[13]*Gx2[10] + Gx1[14]*Gx2[16] + Gx1[15]*Gx2[22] + Gx1[16]*Gx2[28] + Gx1[17]*Gx2[34];
Gx3[17] = + Gx1[12]*Gx2[5] + Gx1[13]*Gx2[11] + Gx1[14]*Gx2[17] + Gx1[15]*Gx2[23] + Gx1[16]*Gx2[29] + Gx1[17]*Gx2[35];
Gx3[18] = + Gx1[18]*Gx2[0] + Gx1[19]*Gx2[6] + Gx1[20]*Gx2[12] + Gx1[21]*Gx2[18] + Gx1[22]*Gx2[24] + Gx1[23]*Gx2[30];
Gx3[19] = + Gx1[18]*Gx2[1] + Gx1[19]*Gx2[7] + Gx1[20]*Gx2[13] + Gx1[21]*Gx2[19] + Gx1[22]*Gx2[25] + Gx1[23]*Gx2[31];
Gx3[20] = + Gx1[18]*Gx2[2] + Gx1[19]*Gx2[8] + Gx1[20]*Gx2[14] + Gx1[21]*Gx2[20] + Gx1[22]*Gx2[26] + Gx1[23]*Gx2[32];
Gx3[21] = + Gx1[18]*Gx2[3] + Gx1[19]*Gx2[9] + Gx1[20]*Gx2[15] + Gx1[21]*Gx2[21] + Gx1[22]*Gx2[27] + Gx1[23]*Gx2[33];
Gx3[22] = + Gx1[18]*Gx2[4] + Gx1[19]*Gx2[10] + Gx1[20]*Gx2[16] + Gx1[21]*Gx2[22] + Gx1[22]*Gx2[28] + Gx1[23]*Gx2[34];
Gx3[23] = + Gx1[18]*Gx2[5] + Gx1[19]*Gx2[11] + Gx1[20]*Gx2[17] + Gx1[21]*Gx2[23] + Gx1[22]*Gx2[29] + Gx1[23]*Gx2[35];
Gx3[24] = + Gx1[24]*Gx2[0] + Gx1[25]*Gx2[6] + Gx1[26]*Gx2[12] + Gx1[27]*Gx2[18] + Gx1[28]*Gx2[24] + Gx1[29]*Gx2[30];
Gx3[25] = + Gx1[24]*Gx2[1] + Gx1[25]*Gx2[7] + Gx1[26]*Gx2[13] + Gx1[27]*Gx2[19] + Gx1[28]*Gx2[25] + Gx1[29]*Gx2[31];
Gx3[26] = + Gx1[24]*Gx2[2] + Gx1[25]*Gx2[8] + Gx1[26]*Gx2[14] + Gx1[27]*Gx2[20] + Gx1[28]*Gx2[26] + Gx1[29]*Gx2[32];
Gx3[27] = + Gx1[24]*Gx2[3] + Gx1[25]*Gx2[9] + Gx1[26]*Gx2[15] + Gx1[27]*Gx2[21] + Gx1[28]*Gx2[27] + Gx1[29]*Gx2[33];
Gx3[28] = + Gx1[24]*Gx2[4] + Gx1[25]*Gx2[10] + Gx1[26]*Gx2[16] + Gx1[27]*Gx2[22] + Gx1[28]*Gx2[28] + Gx1[29]*Gx2[34];
Gx3[29] = + Gx1[24]*Gx2[5] + Gx1[25]*Gx2[11] + Gx1[26]*Gx2[17] + Gx1[27]*Gx2[23] + Gx1[28]*Gx2[29] + Gx1[29]*Gx2[35];
Gx3[30] = + Gx1[30]*Gx2[0] + Gx1[31]*Gx2[6] + Gx1[32]*Gx2[12] + Gx1[33]*Gx2[18] + Gx1[34]*Gx2[24] + Gx1[35]*Gx2[30];
Gx3[31] = + Gx1[30]*Gx2[1] + Gx1[31]*Gx2[7] + Gx1[32]*Gx2[13] + Gx1[33]*Gx2[19] + Gx1[34]*Gx2[25] + Gx1[35]*Gx2[31];
Gx3[32] = + Gx1[30]*Gx2[2] + Gx1[31]*Gx2[8] + Gx1[32]*Gx2[14] + Gx1[33]*Gx2[20] + Gx1[34]*Gx2[26] + Gx1[35]*Gx2[32];
Gx3[33] = + Gx1[30]*Gx2[3] + Gx1[31]*Gx2[9] + Gx1[32]*Gx2[15] + Gx1[33]*Gx2[21] + Gx1[34]*Gx2[27] + Gx1[35]*Gx2[33];
Gx3[34] = + Gx1[30]*Gx2[4] + Gx1[31]*Gx2[10] + Gx1[32]*Gx2[16] + Gx1[33]*Gx2[22] + Gx1[34]*Gx2[28] + Gx1[35]*Gx2[34];
Gx3[35] = + Gx1[30]*Gx2[5] + Gx1[31]*Gx2[11] + Gx1[32]*Gx2[17] + Gx1[33]*Gx2[23] + Gx1[34]*Gx2[29] + Gx1[35]*Gx2[35];
}

void acado_multGxGu( real_t* const Gx1, real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = + Gx1[0]*Gu1[0] + Gx1[1]*Gu1[3] + Gx1[2]*Gu1[6] + Gx1[3]*Gu1[9] + Gx1[4]*Gu1[12] + Gx1[5]*Gu1[15];
Gu2[1] = + Gx1[0]*Gu1[1] + Gx1[1]*Gu1[4] + Gx1[2]*Gu1[7] + Gx1[3]*Gu1[10] + Gx1[4]*Gu1[13] + Gx1[5]*Gu1[16];
Gu2[2] = + Gx1[0]*Gu1[2] + Gx1[1]*Gu1[5] + Gx1[2]*Gu1[8] + Gx1[3]*Gu1[11] + Gx1[4]*Gu1[14] + Gx1[5]*Gu1[17];
Gu2[3] = + Gx1[6]*Gu1[0] + Gx1[7]*Gu1[3] + Gx1[8]*Gu1[6] + Gx1[9]*Gu1[9] + Gx1[10]*Gu1[12] + Gx1[11]*Gu1[15];
Gu2[4] = + Gx1[6]*Gu1[1] + Gx1[7]*Gu1[4] + Gx1[8]*Gu1[7] + Gx1[9]*Gu1[10] + Gx1[10]*Gu1[13] + Gx1[11]*Gu1[16];
Gu2[5] = + Gx1[6]*Gu1[2] + Gx1[7]*Gu1[5] + Gx1[8]*Gu1[8] + Gx1[9]*Gu1[11] + Gx1[10]*Gu1[14] + Gx1[11]*Gu1[17];
Gu2[6] = + Gx1[12]*Gu1[0] + Gx1[13]*Gu1[3] + Gx1[14]*Gu1[6] + Gx1[15]*Gu1[9] + Gx1[16]*Gu1[12] + Gx1[17]*Gu1[15];
Gu2[7] = + Gx1[12]*Gu1[1] + Gx1[13]*Gu1[4] + Gx1[14]*Gu1[7] + Gx1[15]*Gu1[10] + Gx1[16]*Gu1[13] + Gx1[17]*Gu1[16];
Gu2[8] = + Gx1[12]*Gu1[2] + Gx1[13]*Gu1[5] + Gx1[14]*Gu1[8] + Gx1[15]*Gu1[11] + Gx1[16]*Gu1[14] + Gx1[17]*Gu1[17];
Gu2[9] = + Gx1[18]*Gu1[0] + Gx1[19]*Gu1[3] + Gx1[20]*Gu1[6] + Gx1[21]*Gu1[9] + Gx1[22]*Gu1[12] + Gx1[23]*Gu1[15];
Gu2[10] = + Gx1[18]*Gu1[1] + Gx1[19]*Gu1[4] + Gx1[20]*Gu1[7] + Gx1[21]*Gu1[10] + Gx1[22]*Gu1[13] + Gx1[23]*Gu1[16];
Gu2[11] = + Gx1[18]*Gu1[2] + Gx1[19]*Gu1[5] + Gx1[20]*Gu1[8] + Gx1[21]*Gu1[11] + Gx1[22]*Gu1[14] + Gx1[23]*Gu1[17];
Gu2[12] = + Gx1[24]*Gu1[0] + Gx1[25]*Gu1[3] + Gx1[26]*Gu1[6] + Gx1[27]*Gu1[9] + Gx1[28]*Gu1[12] + Gx1[29]*Gu1[15];
Gu2[13] = + Gx1[24]*Gu1[1] + Gx1[25]*Gu1[4] + Gx1[26]*Gu1[7] + Gx1[27]*Gu1[10] + Gx1[28]*Gu1[13] + Gx1[29]*Gu1[16];
Gu2[14] = + Gx1[24]*Gu1[2] + Gx1[25]*Gu1[5] + Gx1[26]*Gu1[8] + Gx1[27]*Gu1[11] + Gx1[28]*Gu1[14] + Gx1[29]*Gu1[17];
Gu2[15] = + Gx1[30]*Gu1[0] + Gx1[31]*Gu1[3] + Gx1[32]*Gu1[6] + Gx1[33]*Gu1[9] + Gx1[34]*Gu1[12] + Gx1[35]*Gu1[15];
Gu2[16] = + Gx1[30]*Gu1[1] + Gx1[31]*Gu1[4] + Gx1[32]*Gu1[7] + Gx1[33]*Gu1[10] + Gx1[34]*Gu1[13] + Gx1[35]*Gu1[16];
Gu2[17] = + Gx1[30]*Gu1[2] + Gx1[31]*Gu1[5] + Gx1[32]*Gu1[8] + Gx1[33]*Gu1[11] + Gx1[34]*Gu1[14] + Gx1[35]*Gu1[17];
}

void acado_moveGuE( real_t* const Gu1, real_t* const Gu2 )
{
Gu2[0] = Gu1[0];
Gu2[1] = Gu1[1];
Gu2[2] = Gu1[2];
Gu2[3] = Gu1[3];
Gu2[4] = Gu1[4];
Gu2[5] = Gu1[5];
Gu2[6] = Gu1[6];
Gu2[7] = Gu1[7];
Gu2[8] = Gu1[8];
Gu2[9] = Gu1[9];
Gu2[10] = Gu1[10];
Gu2[11] = Gu1[11];
Gu2[12] = Gu1[12];
Gu2[13] = Gu1[13];
Gu2[14] = Gu1[14];
Gu2[15] = Gu1[15];
Gu2[16] = Gu1[16];
Gu2[17] = Gu1[17];
}

void acado_setBlockH11( int iRow, int iCol, real_t* const Gu1, real_t* const Gu2 )
{
acadoWorkspace.H[(iRow * 180) + (iCol * 3)] += + Gu1[0]*Gu2[0] + Gu1[3]*Gu2[3] + Gu1[6]*Gu2[6] + Gu1[9]*Gu2[9] + Gu1[12]*Gu2[12] + Gu1[15]*Gu2[15];
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 1)] += + Gu1[0]*Gu2[1] + Gu1[3]*Gu2[4] + Gu1[6]*Gu2[7] + Gu1[9]*Gu2[10] + Gu1[12]*Gu2[13] + Gu1[15]*Gu2[16];
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 2)] += + Gu1[0]*Gu2[2] + Gu1[3]*Gu2[5] + Gu1[6]*Gu2[8] + Gu1[9]*Gu2[11] + Gu1[12]*Gu2[14] + Gu1[15]*Gu2[17];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3)] += + Gu1[1]*Gu2[0] + Gu1[4]*Gu2[3] + Gu1[7]*Gu2[6] + Gu1[10]*Gu2[9] + Gu1[13]*Gu2[12] + Gu1[16]*Gu2[15];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 1)] += + Gu1[1]*Gu2[1] + Gu1[4]*Gu2[4] + Gu1[7]*Gu2[7] + Gu1[10]*Gu2[10] + Gu1[13]*Gu2[13] + Gu1[16]*Gu2[16];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 2)] += + Gu1[1]*Gu2[2] + Gu1[4]*Gu2[5] + Gu1[7]*Gu2[8] + Gu1[10]*Gu2[11] + Gu1[13]*Gu2[14] + Gu1[16]*Gu2[17];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3)] += + Gu1[2]*Gu2[0] + Gu1[5]*Gu2[3] + Gu1[8]*Gu2[6] + Gu1[11]*Gu2[9] + Gu1[14]*Gu2[12] + Gu1[17]*Gu2[15];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 1)] += + Gu1[2]*Gu2[1] + Gu1[5]*Gu2[4] + Gu1[8]*Gu2[7] + Gu1[11]*Gu2[10] + Gu1[14]*Gu2[13] + Gu1[17]*Gu2[16];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 2)] += + Gu1[2]*Gu2[2] + Gu1[5]*Gu2[5] + Gu1[8]*Gu2[8] + Gu1[11]*Gu2[11] + Gu1[14]*Gu2[14] + Gu1[17]*Gu2[17];
}

void acado_setBlockH11_R1( int iRow, int iCol, real_t* const R11 )
{
acadoWorkspace.H[(iRow * 180) + (iCol * 3)] = R11[0];
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 1)] = R11[1];
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 2)] = R11[2];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3)] = R11[3];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 1)] = R11[4];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 2)] = R11[5];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3)] = R11[6];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 1)] = R11[7];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 2)] = R11[8];
}

void acado_zeroBlockH11( int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 180) + (iCol * 3)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 1)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 2)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 1)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 2)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 1)] = 0.0000000000000000e+00;
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 2)] = 0.0000000000000000e+00;
}

void acado_copyHTH( int iRow, int iCol )
{
acadoWorkspace.H[(iRow * 180) + (iCol * 3)] = acadoWorkspace.H[(iCol * 180) + (iRow * 3)];
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 1)] = acadoWorkspace.H[(iCol * 180 + 60) + (iRow * 3)];
acadoWorkspace.H[(iRow * 180) + (iCol * 3 + 2)] = acadoWorkspace.H[(iCol * 180 + 120) + (iRow * 3)];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3)] = acadoWorkspace.H[(iCol * 180) + (iRow * 3 + 1)];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 1)] = acadoWorkspace.H[(iCol * 180 + 60) + (iRow * 3 + 1)];
acadoWorkspace.H[(iRow * 180 + 60) + (iCol * 3 + 2)] = acadoWorkspace.H[(iCol * 180 + 120) + (iRow * 3 + 1)];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3)] = acadoWorkspace.H[(iCol * 180) + (iRow * 3 + 2)];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 1)] = acadoWorkspace.H[(iCol * 180 + 60) + (iRow * 3 + 2)];
acadoWorkspace.H[(iRow * 180 + 120) + (iCol * 3 + 2)] = acadoWorkspace.H[(iCol * 180 + 120) + (iRow * 3 + 2)];
}

void acado_multQ1d( real_t* const Gx1, real_t* const dOld, real_t* const dNew )
{
dNew[0] = + Gx1[0]*dOld[0] + Gx1[1]*dOld[1] + Gx1[2]*dOld[2] + Gx1[3]*dOld[3] + Gx1[4]*dOld[4] + Gx1[5]*dOld[5];
dNew[1] = + Gx1[6]*dOld[0] + Gx1[7]*dOld[1] + Gx1[8]*dOld[2] + Gx1[9]*dOld[3] + Gx1[10]*dOld[4] + Gx1[11]*dOld[5];
dNew[2] = + Gx1[12]*dOld[0] + Gx1[13]*dOld[1] + Gx1[14]*dOld[2] + Gx1[15]*dOld[3] + Gx1[16]*dOld[4] + Gx1[17]*dOld[5];
dNew[3] = + Gx1[18]*dOld[0] + Gx1[19]*dOld[1] + Gx1[20]*dOld[2] + Gx1[21]*dOld[3] + Gx1[22]*dOld[4] + Gx1[23]*dOld[5];
dNew[4] = + Gx1[24]*dOld[0] + Gx1[25]*dOld[1] + Gx1[26]*dOld[2] + Gx1[27]*dOld[3] + Gx1[28]*dOld[4] + Gx1[29]*dOld[5];
dNew[5] = + Gx1[30]*dOld[0] + Gx1[31]*dOld[1] + Gx1[32]*dOld[2] + Gx1[33]*dOld[3] + Gx1[34]*dOld[4] + Gx1[35]*dOld[5];
}

void acado_multQN1d( real_t* const QN1, real_t* const dOld, real_t* const dNew )
{
dNew[0] = + acadoWorkspace.QN1[0]*dOld[0] + acadoWorkspace.QN1[1]*dOld[1] + acadoWorkspace.QN1[2]*dOld[2] + acadoWorkspace.QN1[3]*dOld[3] + acadoWorkspace.QN1[4]*dOld[4] + acadoWorkspace.QN1[5]*dOld[5];
dNew[1] = + acadoWorkspace.QN1[6]*dOld[0] + acadoWorkspace.QN1[7]*dOld[1] + acadoWorkspace.QN1[8]*dOld[2] + acadoWorkspace.QN1[9]*dOld[3] + acadoWorkspace.QN1[10]*dOld[4] + acadoWorkspace.QN1[11]*dOld[5];
dNew[2] = + acadoWorkspace.QN1[12]*dOld[0] + acadoWorkspace.QN1[13]*dOld[1] + acadoWorkspace.QN1[14]*dOld[2] + acadoWorkspace.QN1[15]*dOld[3] + acadoWorkspace.QN1[16]*dOld[4] + acadoWorkspace.QN1[17]*dOld[5];
dNew[3] = + acadoWorkspace.QN1[18]*dOld[0] + acadoWorkspace.QN1[19]*dOld[1] + acadoWorkspace.QN1[20]*dOld[2] + acadoWorkspace.QN1[21]*dOld[3] + acadoWorkspace.QN1[22]*dOld[4] + acadoWorkspace.QN1[23]*dOld[5];
dNew[4] = + acadoWorkspace.QN1[24]*dOld[0] + acadoWorkspace.QN1[25]*dOld[1] + acadoWorkspace.QN1[26]*dOld[2] + acadoWorkspace.QN1[27]*dOld[3] + acadoWorkspace.QN1[28]*dOld[4] + acadoWorkspace.QN1[29]*dOld[5];
dNew[5] = + acadoWorkspace.QN1[30]*dOld[0] + acadoWorkspace.QN1[31]*dOld[1] + acadoWorkspace.QN1[32]*dOld[2] + acadoWorkspace.QN1[33]*dOld[3] + acadoWorkspace.QN1[34]*dOld[4] + acadoWorkspace.QN1[35]*dOld[5];
}

void acado_multRDy( real_t* const R2, real_t* const Dy1, real_t* const RDy1 )
{
RDy1[0] = + R2[0]*Dy1[0] + R2[1]*Dy1[1] + R2[2]*Dy1[2] + R2[3]*Dy1[3] + R2[4]*Dy1[4] + R2[5]*Dy1[5] + R2[6]*Dy1[6] + R2[7]*Dy1[7] + R2[8]*Dy1[8];
RDy1[1] = + R2[9]*Dy1[0] + R2[10]*Dy1[1] + R2[11]*Dy1[2] + R2[12]*Dy1[3] + R2[13]*Dy1[4] + R2[14]*Dy1[5] + R2[15]*Dy1[6] + R2[16]*Dy1[7] + R2[17]*Dy1[8];
RDy1[2] = + R2[18]*Dy1[0] + R2[19]*Dy1[1] + R2[20]*Dy1[2] + R2[21]*Dy1[3] + R2[22]*Dy1[4] + R2[23]*Dy1[5] + R2[24]*Dy1[6] + R2[25]*Dy1[7] + R2[26]*Dy1[8];
}

void acado_multQDy( real_t* const Q2, real_t* const Dy1, real_t* const QDy1 )
{
QDy1[0] = + Q2[0]*Dy1[0] + Q2[1]*Dy1[1] + Q2[2]*Dy1[2] + Q2[3]*Dy1[3] + Q2[4]*Dy1[4] + Q2[5]*Dy1[5] + Q2[6]*Dy1[6] + Q2[7]*Dy1[7] + Q2[8]*Dy1[8];
QDy1[1] = + Q2[9]*Dy1[0] + Q2[10]*Dy1[1] + Q2[11]*Dy1[2] + Q2[12]*Dy1[3] + Q2[13]*Dy1[4] + Q2[14]*Dy1[5] + Q2[15]*Dy1[6] + Q2[16]*Dy1[7] + Q2[17]*Dy1[8];
QDy1[2] = + Q2[18]*Dy1[0] + Q2[19]*Dy1[1] + Q2[20]*Dy1[2] + Q2[21]*Dy1[3] + Q2[22]*Dy1[4] + Q2[23]*Dy1[5] + Q2[24]*Dy1[6] + Q2[25]*Dy1[7] + Q2[26]*Dy1[8];
QDy1[3] = + Q2[27]*Dy1[0] + Q2[28]*Dy1[1] + Q2[29]*Dy1[2] + Q2[30]*Dy1[3] + Q2[31]*Dy1[4] + Q2[32]*Dy1[5] + Q2[33]*Dy1[6] + Q2[34]*Dy1[7] + Q2[35]*Dy1[8];
QDy1[4] = + Q2[36]*Dy1[0] + Q2[37]*Dy1[1] + Q2[38]*Dy1[2] + Q2[39]*Dy1[3] + Q2[40]*Dy1[4] + Q2[41]*Dy1[5] + Q2[42]*Dy1[6] + Q2[43]*Dy1[7] + Q2[44]*Dy1[8];
QDy1[5] = + Q2[45]*Dy1[0] + Q2[46]*Dy1[1] + Q2[47]*Dy1[2] + Q2[48]*Dy1[3] + Q2[49]*Dy1[4] + Q2[50]*Dy1[5] + Q2[51]*Dy1[6] + Q2[52]*Dy1[7] + Q2[53]*Dy1[8];
}

void acado_multEQDy( real_t* const E1, real_t* const QDy1, real_t* const U1 )
{
U1[0] += + E1[0]*QDy1[0] + E1[3]*QDy1[1] + E1[6]*QDy1[2] + E1[9]*QDy1[3] + E1[12]*QDy1[4] + E1[15]*QDy1[5];
U1[1] += + E1[1]*QDy1[0] + E1[4]*QDy1[1] + E1[7]*QDy1[2] + E1[10]*QDy1[3] + E1[13]*QDy1[4] + E1[16]*QDy1[5];
U1[2] += + E1[2]*QDy1[0] + E1[5]*QDy1[1] + E1[8]*QDy1[2] + E1[11]*QDy1[3] + E1[14]*QDy1[4] + E1[17]*QDy1[5];
}

void acado_multQETGx( real_t* const E1, real_t* const Gx1, real_t* const H101 )
{
H101[0] += + E1[0]*Gx1[0] + E1[3]*Gx1[6] + E1[6]*Gx1[12] + E1[9]*Gx1[18] + E1[12]*Gx1[24] + E1[15]*Gx1[30];
H101[1] += + E1[0]*Gx1[1] + E1[3]*Gx1[7] + E1[6]*Gx1[13] + E1[9]*Gx1[19] + E1[12]*Gx1[25] + E1[15]*Gx1[31];
H101[2] += + E1[0]*Gx1[2] + E1[3]*Gx1[8] + E1[6]*Gx1[14] + E1[9]*Gx1[20] + E1[12]*Gx1[26] + E1[15]*Gx1[32];
H101[3] += + E1[0]*Gx1[3] + E1[3]*Gx1[9] + E1[6]*Gx1[15] + E1[9]*Gx1[21] + E1[12]*Gx1[27] + E1[15]*Gx1[33];
H101[4] += + E1[0]*Gx1[4] + E1[3]*Gx1[10] + E1[6]*Gx1[16] + E1[9]*Gx1[22] + E1[12]*Gx1[28] + E1[15]*Gx1[34];
H101[5] += + E1[0]*Gx1[5] + E1[3]*Gx1[11] + E1[6]*Gx1[17] + E1[9]*Gx1[23] + E1[12]*Gx1[29] + E1[15]*Gx1[35];
H101[6] += + E1[1]*Gx1[0] + E1[4]*Gx1[6] + E1[7]*Gx1[12] + E1[10]*Gx1[18] + E1[13]*Gx1[24] + E1[16]*Gx1[30];
H101[7] += + E1[1]*Gx1[1] + E1[4]*Gx1[7] + E1[7]*Gx1[13] + E1[10]*Gx1[19] + E1[13]*Gx1[25] + E1[16]*Gx1[31];
H101[8] += + E1[1]*Gx1[2] + E1[4]*Gx1[8] + E1[7]*Gx1[14] + E1[10]*Gx1[20] + E1[13]*Gx1[26] + E1[16]*Gx1[32];
H101[9] += + E1[1]*Gx1[3] + E1[4]*Gx1[9] + E1[7]*Gx1[15] + E1[10]*Gx1[21] + E1[13]*Gx1[27] + E1[16]*Gx1[33];
H101[10] += + E1[1]*Gx1[4] + E1[4]*Gx1[10] + E1[7]*Gx1[16] + E1[10]*Gx1[22] + E1[13]*Gx1[28] + E1[16]*Gx1[34];
H101[11] += + E1[1]*Gx1[5] + E1[4]*Gx1[11] + E1[7]*Gx1[17] + E1[10]*Gx1[23] + E1[13]*Gx1[29] + E1[16]*Gx1[35];
H101[12] += + E1[2]*Gx1[0] + E1[5]*Gx1[6] + E1[8]*Gx1[12] + E1[11]*Gx1[18] + E1[14]*Gx1[24] + E1[17]*Gx1[30];
H101[13] += + E1[2]*Gx1[1] + E1[5]*Gx1[7] + E1[8]*Gx1[13] + E1[11]*Gx1[19] + E1[14]*Gx1[25] + E1[17]*Gx1[31];
H101[14] += + E1[2]*Gx1[2] + E1[5]*Gx1[8] + E1[8]*Gx1[14] + E1[11]*Gx1[20] + E1[14]*Gx1[26] + E1[17]*Gx1[32];
H101[15] += + E1[2]*Gx1[3] + E1[5]*Gx1[9] + E1[8]*Gx1[15] + E1[11]*Gx1[21] + E1[14]*Gx1[27] + E1[17]*Gx1[33];
H101[16] += + E1[2]*Gx1[4] + E1[5]*Gx1[10] + E1[8]*Gx1[16] + E1[11]*Gx1[22] + E1[14]*Gx1[28] + E1[17]*Gx1[34];
H101[17] += + E1[2]*Gx1[5] + E1[5]*Gx1[11] + E1[8]*Gx1[17] + E1[11]*Gx1[23] + E1[14]*Gx1[29] + E1[17]*Gx1[35];
}

void acado_zeroBlockH10( real_t* const H101 )
{
{ int lCopy; for (lCopy = 0; lCopy < 18; lCopy++) H101[ lCopy ] = 0; }
}

void acado_multEDu( real_t* const E1, real_t* const U1, real_t* const dNew )
{
dNew[0] += + E1[0]*U1[0] + E1[1]*U1[1] + E1[2]*U1[2];
dNew[1] += + E1[3]*U1[0] + E1[4]*U1[1] + E1[5]*U1[2];
dNew[2] += + E1[6]*U1[0] + E1[7]*U1[1] + E1[8]*U1[2];
dNew[3] += + E1[9]*U1[0] + E1[10]*U1[1] + E1[11]*U1[2];
dNew[4] += + E1[12]*U1[0] + E1[13]*U1[1] + E1[14]*U1[2];
dNew[5] += + E1[15]*U1[0] + E1[16]*U1[1] + E1[17]*U1[2];
}

void acado_macETSlu( real_t* const E0, real_t* const g1 )
{
g1[0] += 0.0;
;
g1[1] += 0.0;
;
g1[2] += 0.0;
;
}

void acado_condensePrep(  )
{
acado_moveGuE( acadoWorkspace.evGu, acadoWorkspace.E );
acado_moveGxT( &(acadoWorkspace.evGx[ 36 ]), acadoWorkspace.T );
acado_multGxd( acadoWorkspace.d, &(acadoWorkspace.evGx[ 36 ]), &(acadoWorkspace.d[ 6 ]) );
acado_multGxGx( acadoWorkspace.T, acadoWorkspace.evGx, &(acadoWorkspace.evGx[ 36 ]) );

acado_multGxGu( acadoWorkspace.T, acadoWorkspace.E, &(acadoWorkspace.E[ 18 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 18 ]), &(acadoWorkspace.E[ 36 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 72 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 6 ]), &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.d[ 12 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 36 ]), &(acadoWorkspace.evGx[ 72 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 18 ]), &(acadoWorkspace.E[ 54 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 36 ]), &(acadoWorkspace.E[ 72 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 36 ]), &(acadoWorkspace.E[ 90 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 108 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 12 ]), &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.d[ 18 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.evGx[ 108 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.E[ 108 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.E[ 126 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 90 ]), &(acadoWorkspace.E[ 144 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 54 ]), &(acadoWorkspace.E[ 162 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 18 ]), &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.d[ 24 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.evGx[ 144 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.E[ 180 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.E[ 198 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.E[ 216 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 162 ]), &(acadoWorkspace.E[ 234 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 72 ]), &(acadoWorkspace.E[ 252 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 24 ]), &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.d[ 30 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.evGx[ 180 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.E[ 270 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.E[ 288 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.E[ 306 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.E[ 324 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 252 ]), &(acadoWorkspace.E[ 342 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 90 ]), &(acadoWorkspace.E[ 360 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 30 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.d[ 36 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.evGx[ 216 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.E[ 378 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.E[ 396 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.E[ 414 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.E[ 432 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.E[ 450 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 360 ]), &(acadoWorkspace.E[ 468 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 108 ]), &(acadoWorkspace.E[ 486 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 36 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.d[ 42 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.evGx[ 252 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.E[ 504 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.E[ 522 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.E[ 540 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.E[ 558 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.E[ 576 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 468 ]), &(acadoWorkspace.E[ 594 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 486 ]), &(acadoWorkspace.E[ 612 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 126 ]), &(acadoWorkspace.E[ 630 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 42 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.d[ 48 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.evGx[ 288 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.E[ 648 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.E[ 666 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.E[ 684 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.E[ 702 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.E[ 720 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.E[ 738 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.E[ 756 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 630 ]), &(acadoWorkspace.E[ 774 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 144 ]), &(acadoWorkspace.E[ 792 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 48 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.d[ 54 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.evGx[ 324 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.E[ 810 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.E[ 828 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.E[ 846 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.E[ 864 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.E[ 882 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.E[ 900 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.E[ 918 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.E[ 936 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 792 ]), &(acadoWorkspace.E[ 954 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 162 ]), &(acadoWorkspace.E[ 972 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 360 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 54 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.d[ 60 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.evGx[ 360 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.E[ 990 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.E[ 1008 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.E[ 1026 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.E[ 1044 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.E[ 1062 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.E[ 1080 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.E[ 1098 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.E[ 1116 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 954 ]), &(acadoWorkspace.E[ 1134 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 972 ]), &(acadoWorkspace.E[ 1152 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 180 ]), &(acadoWorkspace.E[ 1170 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 396 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 60 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.d[ 66 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.evGx[ 396 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.E[ 1188 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.E[ 1206 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.E[ 1224 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.E[ 1242 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.E[ 1260 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.E[ 1278 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.E[ 1296 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.E[ 1314 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.E[ 1332 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1152 ]), &(acadoWorkspace.E[ 1350 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1170 ]), &(acadoWorkspace.E[ 1368 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 198 ]), &(acadoWorkspace.E[ 1386 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 432 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 66 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.d[ 72 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.evGx[ 432 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.E[ 1404 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.E[ 1422 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.E[ 1440 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.E[ 1458 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.E[ 1476 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.E[ 1494 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.E[ 1512 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.E[ 1530 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.E[ 1548 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.E[ 1566 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1368 ]), &(acadoWorkspace.E[ 1584 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1386 ]), &(acadoWorkspace.E[ 1602 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 216 ]), &(acadoWorkspace.E[ 1620 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 468 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 72 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.d[ 78 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.evGx[ 468 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.E[ 1638 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.E[ 1656 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.E[ 1674 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.E[ 1692 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.E[ 1710 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.E[ 1728 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.E[ 1746 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.E[ 1764 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.E[ 1782 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.E[ 1800 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.E[ 1818 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1602 ]), &(acadoWorkspace.E[ 1836 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1620 ]), &(acadoWorkspace.E[ 1854 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 234 ]), &(acadoWorkspace.E[ 1872 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 504 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 78 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.d[ 84 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.evGx[ 504 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.E[ 1890 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.E[ 1908 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.E[ 1926 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.E[ 1944 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.E[ 1962 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.E[ 1980 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.E[ 1998 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.E[ 2016 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.E[ 2034 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.E[ 2052 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.E[ 2070 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.E[ 2088 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1854 ]), &(acadoWorkspace.E[ 2106 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1872 ]), &(acadoWorkspace.E[ 2124 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 252 ]), &(acadoWorkspace.E[ 2142 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 540 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 84 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.d[ 90 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.evGx[ 540 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.E[ 2160 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.E[ 2178 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.E[ 2196 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.E[ 2214 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.E[ 2232 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.E[ 2250 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.E[ 2268 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.E[ 2286 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.E[ 2304 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.E[ 2322 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.E[ 2340 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.E[ 2358 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.E[ 2376 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2124 ]), &(acadoWorkspace.E[ 2394 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2142 ]), &(acadoWorkspace.E[ 2412 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 270 ]), &(acadoWorkspace.E[ 2430 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 576 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 90 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.d[ 96 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.evGx[ 576 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.E[ 2448 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.E[ 2466 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.E[ 2484 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.E[ 2502 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.E[ 2520 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.E[ 2538 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.E[ 2556 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.E[ 2574 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.E[ 2592 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.E[ 2610 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.E[ 2628 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.E[ 2646 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.E[ 2664 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.E[ 2682 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2412 ]), &(acadoWorkspace.E[ 2700 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2430 ]), &(acadoWorkspace.E[ 2718 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 288 ]), &(acadoWorkspace.E[ 2736 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 612 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 96 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.d[ 102 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.evGx[ 612 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.E[ 2754 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.E[ 2772 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.E[ 2790 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.E[ 2808 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.E[ 2826 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.E[ 2844 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.E[ 2862 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.E[ 2880 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.E[ 2898 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.E[ 2916 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.E[ 2934 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.E[ 2952 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.E[ 2970 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.E[ 2988 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.E[ 3006 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2718 ]), &(acadoWorkspace.E[ 3024 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2736 ]), &(acadoWorkspace.E[ 3042 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 306 ]), &(acadoWorkspace.E[ 3060 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 648 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 102 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.d[ 108 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.evGx[ 648 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.E[ 3078 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.E[ 3096 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.E[ 3114 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.E[ 3132 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.E[ 3150 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.E[ 3168 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.E[ 3186 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.E[ 3204 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.E[ 3222 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.E[ 3240 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.E[ 3258 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.E[ 3276 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.E[ 3294 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.E[ 3312 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.E[ 3330 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.E[ 3348 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3042 ]), &(acadoWorkspace.E[ 3366 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3060 ]), &(acadoWorkspace.E[ 3384 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 324 ]), &(acadoWorkspace.E[ 3402 ]) );

acado_moveGxT( &(acadoWorkspace.evGx[ 684 ]), acadoWorkspace.T );
acado_multGxd( &(acadoWorkspace.d[ 108 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.d[ 114 ]) );
acado_multGxGx( acadoWorkspace.T, &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.evGx[ 684 ]) );

acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.E[ 3420 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.E[ 3438 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.E[ 3456 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.E[ 3474 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.E[ 3492 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.E[ 3510 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.E[ 3528 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.E[ 3546 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.E[ 3564 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.E[ 3582 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.E[ 3600 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.E[ 3618 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.E[ 3636 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.E[ 3654 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.E[ 3672 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.E[ 3690 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.E[ 3708 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3384 ]), &(acadoWorkspace.E[ 3726 ]) );
acado_multGxGu( acadoWorkspace.T, &(acadoWorkspace.E[ 3402 ]), &(acadoWorkspace.E[ 3744 ]) );

acado_moveGuE( &(acadoWorkspace.evGu[ 342 ]), &(acadoWorkspace.E[ 3762 ]) );

acado_multGxGu( &(acadoWorkspace.Q1[ 36 ]), acadoWorkspace.E, acadoWorkspace.QE );
acado_multGxGu( &(acadoWorkspace.Q1[ 72 ]), &(acadoWorkspace.E[ 18 ]), &(acadoWorkspace.QE[ 18 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 72 ]), &(acadoWorkspace.E[ 36 ]), &(acadoWorkspace.QE[ 36 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.QE[ 54 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.QE[ 72 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.E[ 90 ]), &(acadoWorkspace.QE[ 90 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.QE[ 108 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.QE[ 126 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.QE[ 144 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.E[ 162 ]), &(acadoWorkspace.QE[ 162 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QE[ 180 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.QE[ 198 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.QE[ 216 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.QE[ 234 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.E[ 252 ]), &(acadoWorkspace.QE[ 252 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 270 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QE[ 288 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.QE[ 306 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.QE[ 324 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.QE[ 342 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.E[ 360 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 378 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 396 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QE[ 414 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.QE[ 432 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.QE[ 450 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 468 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.E[ 486 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 504 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 522 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 540 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QE[ 558 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.QE[ 576 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.E[ 630 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 648 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 666 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 684 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 702 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QE[ 720 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.E[ 792 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 810 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 828 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 846 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 864 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 882 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 954 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.E[ 972 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 990 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1008 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1026 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1044 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1062 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1152 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.E[ 1170 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1188 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1206 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1224 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1242 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1260 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1368 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.E[ 1386 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1404 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1422 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1440 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1458 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1476 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1602 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.E[ 1620 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1638 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1656 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1674 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1692 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1710 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1854 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.E[ 1872 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1890 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1908 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 1926 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 1944 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 1962 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2124 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.E[ 2142 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2160 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2178 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2196 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2214 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2232 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2412 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.E[ 2430 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2448 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2466 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2484 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2502 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2520 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2718 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.E[ 2736 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2754 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2772 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2790 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2808 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2826 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 3042 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.E[ 3060 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3078 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3096 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3114 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3132 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3150 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3384 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_multGxGu( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.E[ 3402 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3420 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3438 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3456 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3474 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3492 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3510 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3528 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3546 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3564 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3582 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3600 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3618 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3636 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3654 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3672 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QE[ 3690 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.QE[ 3708 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3726 ]), &(acadoWorkspace.QE[ 3726 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3744 ]), &(acadoWorkspace.QE[ 3744 ]) );
acado_multGxGu( acadoWorkspace.QN1, &(acadoWorkspace.E[ 3762 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_zeroBlockH10( acadoWorkspace.H10 );
acado_multQETGx( acadoWorkspace.QE, acadoWorkspace.evGx, acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 18 ]), &(acadoWorkspace.evGx[ 36 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 54 ]), &(acadoWorkspace.evGx[ 72 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 108 ]), &(acadoWorkspace.evGx[ 108 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 180 ]), &(acadoWorkspace.evGx[ 144 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 270 ]), &(acadoWorkspace.evGx[ 180 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 378 ]), &(acadoWorkspace.evGx[ 216 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 504 ]), &(acadoWorkspace.evGx[ 252 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 648 ]), &(acadoWorkspace.evGx[ 288 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 810 ]), &(acadoWorkspace.evGx[ 324 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 990 ]), &(acadoWorkspace.evGx[ 360 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 1188 ]), &(acadoWorkspace.evGx[ 396 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 1404 ]), &(acadoWorkspace.evGx[ 432 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 1638 ]), &(acadoWorkspace.evGx[ 468 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 1890 ]), &(acadoWorkspace.evGx[ 504 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 2160 ]), &(acadoWorkspace.evGx[ 540 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 2448 ]), &(acadoWorkspace.evGx[ 576 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 2754 ]), &(acadoWorkspace.evGx[ 612 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 3078 ]), &(acadoWorkspace.evGx[ 648 ]), acadoWorkspace.H10 );
acado_multQETGx( &(acadoWorkspace.QE[ 3420 ]), &(acadoWorkspace.evGx[ 684 ]), acadoWorkspace.H10 );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 36 ]), &(acadoWorkspace.evGx[ 36 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 72 ]), &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 126 ]), &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 198 ]), &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 288 ]), &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 396 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 522 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 666 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 828 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1008 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1206 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1422 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1656 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1908 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2178 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2466 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2772 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3096 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3438 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 18 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 90 ]), &(acadoWorkspace.evGx[ 72 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 144 ]), &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 216 ]), &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 306 ]), &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 414 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 540 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 684 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 846 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1026 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1224 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1440 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1674 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1926 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2196 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2484 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2790 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3114 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3456 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 36 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 162 ]), &(acadoWorkspace.evGx[ 108 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 234 ]), &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 324 ]), &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 432 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 558 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 702 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 864 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1044 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1242 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1458 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1692 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1944 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2214 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2502 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2808 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3132 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3474 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 54 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 252 ]), &(acadoWorkspace.evGx[ 144 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 342 ]), &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 450 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 576 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 720 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 882 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1062 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1260 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1476 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1710 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1962 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2232 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2520 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2826 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3150 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3492 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 72 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 360 ]), &(acadoWorkspace.evGx[ 180 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 468 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 594 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 738 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 900 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1080 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1278 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1494 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1728 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1980 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2250 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2538 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2844 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3168 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3510 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 90 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 486 ]), &(acadoWorkspace.evGx[ 216 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 612 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 756 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 918 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1098 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1296 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1512 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1746 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1998 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2268 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2556 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2862 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3186 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3528 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 108 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 630 ]), &(acadoWorkspace.evGx[ 252 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 774 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 936 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1116 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1314 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1530 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1764 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2016 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2286 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2574 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2880 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3204 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3546 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 126 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 792 ]), &(acadoWorkspace.evGx[ 288 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 954 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1134 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1332 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1548 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1782 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2034 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2304 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2592 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2898 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3222 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3564 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 144 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 972 ]), &(acadoWorkspace.evGx[ 324 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1152 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1350 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1566 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1800 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2052 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2322 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2610 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2916 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3240 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3582 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 162 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1170 ]), &(acadoWorkspace.evGx[ 360 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1368 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1584 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1818 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2070 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2340 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2628 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2934 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3258 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3600 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 180 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1386 ]), &(acadoWorkspace.evGx[ 396 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1602 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1836 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2088 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2358 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2646 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2952 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3276 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3618 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 198 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1620 ]), &(acadoWorkspace.evGx[ 432 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1854 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2106 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2376 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2664 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2970 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3294 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3636 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 216 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 1872 ]), &(acadoWorkspace.evGx[ 468 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2124 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2394 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2682 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2988 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3312 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3654 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 234 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 252 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2142 ]), &(acadoWorkspace.evGx[ 504 ]), &(acadoWorkspace.H10[ 252 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2412 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 252 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2700 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 252 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3006 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 252 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3330 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 252 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3672 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 252 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 270 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2430 ]), &(acadoWorkspace.evGx[ 540 ]), &(acadoWorkspace.H10[ 270 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2718 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 270 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3024 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 270 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3348 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 270 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3690 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 270 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 288 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 2736 ]), &(acadoWorkspace.evGx[ 576 ]), &(acadoWorkspace.H10[ 288 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3042 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 288 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3366 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 288 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3708 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 288 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 306 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3060 ]), &(acadoWorkspace.evGx[ 612 ]), &(acadoWorkspace.H10[ 306 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3384 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 306 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3726 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 306 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 324 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3402 ]), &(acadoWorkspace.evGx[ 648 ]), &(acadoWorkspace.H10[ 324 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3744 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 324 ]) );
acado_zeroBlockH10( &(acadoWorkspace.H10[ 342 ]) );
acado_multQETGx( &(acadoWorkspace.QE[ 3762 ]), &(acadoWorkspace.evGx[ 684 ]), &(acadoWorkspace.H10[ 342 ]) );

acado_setBlockH11_R1( 0, 0, acadoWorkspace.R1 );
acado_setBlockH11( 0, 0, acadoWorkspace.E, acadoWorkspace.QE );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 18 ]), &(acadoWorkspace.QE[ 18 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.QE[ 54 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.QE[ 108 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QE[ 180 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 270 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 378 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 504 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 648 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 810 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 990 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1188 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1404 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1638 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1890 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2160 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2448 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2754 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3078 ]) );
acado_setBlockH11( 0, 0, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3420 ]) );

acado_zeroBlockH11( 0, 1 );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 18 ]), &(acadoWorkspace.QE[ 36 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.QE[ 72 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.QE[ 126 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QE[ 198 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 288 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 396 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 522 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 666 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 828 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1008 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1206 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1422 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1656 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1908 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2178 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2466 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2772 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3096 ]) );
acado_setBlockH11( 0, 1, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3438 ]) );

acado_zeroBlockH11( 0, 2 );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.QE[ 90 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.QE[ 144 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QE[ 216 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 306 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 414 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 540 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 684 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 846 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1026 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1224 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1440 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1674 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1926 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2196 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2484 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2790 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3114 ]) );
acado_setBlockH11( 0, 2, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3456 ]) );

acado_zeroBlockH11( 0, 3 );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.QE[ 162 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QE[ 234 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 324 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 432 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 558 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 702 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 864 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1044 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1242 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1458 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1692 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1944 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2214 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2502 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2808 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3132 ]) );
acado_setBlockH11( 0, 3, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3474 ]) );

acado_zeroBlockH11( 0, 4 );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QE[ 252 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 342 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 450 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 576 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 720 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 882 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1062 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1260 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1476 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1710 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1962 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2232 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2520 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2826 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3150 ]) );
acado_setBlockH11( 0, 4, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3492 ]) );

acado_zeroBlockH11( 0, 5 );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_setBlockH11( 0, 5, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3510 ]) );

acado_zeroBlockH11( 0, 6 );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 0, 6, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 0, 7 );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 0, 7, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 0, 8 );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 0, 8, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 0, 9 );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 0, 9, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 0, 10 );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 0, 10, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 0, 11 );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 0, 11, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 0, 12 );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 0, 12, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 0, 13 );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 0, 13, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 0, 14 );
acado_setBlockH11( 0, 14, &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 0, 14, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 0, 14, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 0, 14, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 0, 14, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 0, 14, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 0, 15 );
acado_setBlockH11( 0, 15, &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 0, 15, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 0, 15, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 0, 15, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 0, 15, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 0, 16 );
acado_setBlockH11( 0, 16, &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 0, 16, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 0, 16, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 0, 16, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 0, 17 );
acado_setBlockH11( 0, 17, &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 0, 17, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 0, 17, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 0, 18 );
acado_setBlockH11( 0, 18, &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 0, 18, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 0, 19 );
acado_setBlockH11( 0, 19, &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 1, 1, &(acadoWorkspace.R1[ 9 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 36 ]), &(acadoWorkspace.QE[ 36 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.QE[ 72 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.QE[ 126 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.QE[ 198 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QE[ 288 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 396 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 522 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 666 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 828 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1008 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1206 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1422 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1656 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1908 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2178 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2466 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2772 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3096 ]) );
acado_setBlockH11( 1, 1, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3438 ]) );

acado_zeroBlockH11( 1, 2 );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.QE[ 90 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.QE[ 144 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.QE[ 216 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QE[ 306 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 414 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 540 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 684 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 846 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1026 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1224 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1440 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1674 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1926 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2196 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2484 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2790 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3114 ]) );
acado_setBlockH11( 1, 2, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3456 ]) );

acado_zeroBlockH11( 1, 3 );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.QE[ 162 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.QE[ 234 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QE[ 324 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 432 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 558 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 702 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 864 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1044 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1242 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1458 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1692 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1944 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2214 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2502 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2808 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3132 ]) );
acado_setBlockH11( 1, 3, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3474 ]) );

acado_zeroBlockH11( 1, 4 );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.QE[ 252 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QE[ 342 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 450 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 576 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 720 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 882 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1062 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1260 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1476 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1710 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1962 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2232 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2520 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2826 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3150 ]) );
acado_setBlockH11( 1, 4, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3492 ]) );

acado_zeroBlockH11( 1, 5 );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_setBlockH11( 1, 5, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3510 ]) );

acado_zeroBlockH11( 1, 6 );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 1, 6, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 1, 7 );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 1, 7, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 1, 8 );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 1, 8, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 1, 9 );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 1, 9, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 1, 10 );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 1, 10, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 1, 11 );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 1, 11, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 1, 12 );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 1, 12, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 1, 13 );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 1, 13, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 1, 14 );
acado_setBlockH11( 1, 14, &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 1, 14, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 1, 14, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 1, 14, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 1, 14, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 1, 14, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 1, 15 );
acado_setBlockH11( 1, 15, &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 1, 15, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 1, 15, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 1, 15, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 1, 15, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 1, 16 );
acado_setBlockH11( 1, 16, &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 1, 16, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 1, 16, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 1, 16, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 1, 17 );
acado_setBlockH11( 1, 17, &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 1, 17, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 1, 17, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 1, 18 );
acado_setBlockH11( 1, 18, &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 1, 18, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 1, 19 );
acado_setBlockH11( 1, 19, &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 2, 2, &(acadoWorkspace.R1[ 18 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 90 ]), &(acadoWorkspace.QE[ 90 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.QE[ 144 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.QE[ 216 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.QE[ 306 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QE[ 414 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 540 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 684 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 846 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1026 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1224 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1440 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1674 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 1926 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2196 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2484 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2790 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3114 ]) );
acado_setBlockH11( 2, 2, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3456 ]) );

acado_zeroBlockH11( 2, 3 );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.QE[ 162 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.QE[ 234 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.QE[ 324 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QE[ 432 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 558 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 702 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 864 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1044 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1242 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1458 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1692 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 1944 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2214 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2502 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2808 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3132 ]) );
acado_setBlockH11( 2, 3, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3474 ]) );

acado_zeroBlockH11( 2, 4 );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.QE[ 252 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.QE[ 342 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QE[ 450 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 576 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 720 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 882 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1062 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1260 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1476 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1710 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 1962 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2232 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2520 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2826 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3150 ]) );
acado_setBlockH11( 2, 4, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3492 ]) );

acado_zeroBlockH11( 2, 5 );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_setBlockH11( 2, 5, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3510 ]) );

acado_zeroBlockH11( 2, 6 );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 2, 6, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 2, 7 );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 2, 7, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 2, 8 );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 2, 8, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 2, 9 );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 2, 9, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 2, 10 );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 2, 10, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 2, 11 );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 2, 11, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 2, 12 );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 2, 12, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 2, 13 );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 2, 13, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 2, 14 );
acado_setBlockH11( 2, 14, &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 2, 14, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 2, 14, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 2, 14, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 2, 14, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 2, 14, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 2, 15 );
acado_setBlockH11( 2, 15, &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 2, 15, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 2, 15, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 2, 15, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 2, 15, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 2, 16 );
acado_setBlockH11( 2, 16, &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 2, 16, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 2, 16, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 2, 16, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 2, 17 );
acado_setBlockH11( 2, 17, &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 2, 17, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 2, 17, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 2, 18 );
acado_setBlockH11( 2, 18, &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 2, 18, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 2, 19 );
acado_setBlockH11( 2, 19, &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 3, 3, &(acadoWorkspace.R1[ 27 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 162 ]), &(acadoWorkspace.QE[ 162 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.QE[ 234 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.QE[ 324 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.QE[ 432 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QE[ 558 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 702 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 864 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1044 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1242 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1458 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1692 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 1944 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2214 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2502 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2808 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3132 ]) );
acado_setBlockH11( 3, 3, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3474 ]) );

acado_zeroBlockH11( 3, 4 );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.QE[ 252 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.QE[ 342 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.QE[ 450 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QE[ 576 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 720 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 882 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1062 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1260 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1476 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1710 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 1962 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2232 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2520 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2826 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3150 ]) );
acado_setBlockH11( 3, 4, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3492 ]) );

acado_zeroBlockH11( 3, 5 );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_setBlockH11( 3, 5, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3510 ]) );

acado_zeroBlockH11( 3, 6 );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 3, 6, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 3, 7 );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 3, 7, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 3, 8 );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 3, 8, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 3, 9 );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 3, 9, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 3, 10 );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 3, 10, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 3, 11 );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 3, 11, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 3, 12 );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 3, 12, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 3, 13 );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 3, 13, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 3, 14 );
acado_setBlockH11( 3, 14, &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 3, 14, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 3, 14, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 3, 14, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 3, 14, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 3, 14, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 3, 15 );
acado_setBlockH11( 3, 15, &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 3, 15, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 3, 15, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 3, 15, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 3, 15, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 3, 16 );
acado_setBlockH11( 3, 16, &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 3, 16, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 3, 16, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 3, 16, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 3, 17 );
acado_setBlockH11( 3, 17, &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 3, 17, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 3, 17, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 3, 18 );
acado_setBlockH11( 3, 18, &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 3, 18, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 3, 19 );
acado_setBlockH11( 3, 19, &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 4, 4, &(acadoWorkspace.R1[ 36 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 252 ]), &(acadoWorkspace.QE[ 252 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.QE[ 342 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.QE[ 450 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.QE[ 576 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QE[ 720 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 882 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1062 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1260 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1476 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1710 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 1962 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2232 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2520 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2826 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3150 ]) );
acado_setBlockH11( 4, 4, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3492 ]) );

acado_zeroBlockH11( 4, 5 );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_setBlockH11( 4, 5, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3510 ]) );

acado_zeroBlockH11( 4, 6 );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 4, 6, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 4, 7 );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 4, 7, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 4, 8 );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 4, 8, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 4, 9 );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 4, 9, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 4, 10 );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 4, 10, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 4, 11 );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 4, 11, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 4, 12 );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 4, 12, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 4, 13 );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 4, 13, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 4, 14 );
acado_setBlockH11( 4, 14, &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 4, 14, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 4, 14, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 4, 14, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 4, 14, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 4, 14, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 4, 15 );
acado_setBlockH11( 4, 15, &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 4, 15, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 4, 15, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 4, 15, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 4, 15, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 4, 16 );
acado_setBlockH11( 4, 16, &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 4, 16, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 4, 16, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 4, 16, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 4, 17 );
acado_setBlockH11( 4, 17, &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 4, 17, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 4, 17, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 4, 18 );
acado_setBlockH11( 4, 18, &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 4, 18, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 4, 19 );
acado_setBlockH11( 4, 19, &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 5, 5, &(acadoWorkspace.R1[ 45 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 360 ]), &(acadoWorkspace.QE[ 360 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 468 ]), &(acadoWorkspace.QE[ 468 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.QE[ 594 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.QE[ 738 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QE[ 900 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1080 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1278 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1494 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1728 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 1980 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2250 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2538 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2844 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3168 ]) );
acado_setBlockH11( 5, 5, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3510 ]) );

acado_zeroBlockH11( 5, 6 );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 468 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 5, 6, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 5, 7 );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 5, 7, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 5, 8 );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 5, 8, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 5, 9 );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 5, 9, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 5, 10 );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 5, 10, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 5, 11 );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 5, 11, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 5, 12 );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 5, 12, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 5, 13 );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 5, 13, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 5, 14 );
acado_setBlockH11( 5, 14, &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 5, 14, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 5, 14, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 5, 14, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 5, 14, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 5, 14, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 5, 15 );
acado_setBlockH11( 5, 15, &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 5, 15, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 5, 15, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 5, 15, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 5, 15, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 5, 16 );
acado_setBlockH11( 5, 16, &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 5, 16, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 5, 16, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 5, 16, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 5, 17 );
acado_setBlockH11( 5, 17, &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 5, 17, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 5, 17, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 5, 18 );
acado_setBlockH11( 5, 18, &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 5, 18, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 5, 19 );
acado_setBlockH11( 5, 19, &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 6, 6, &(acadoWorkspace.R1[ 54 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 486 ]), &(acadoWorkspace.QE[ 486 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.QE[ 612 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.QE[ 756 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.QE[ 918 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QE[ 1098 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1296 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1512 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1746 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 1998 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2268 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2556 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2862 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3186 ]) );
acado_setBlockH11( 6, 6, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3528 ]) );

acado_zeroBlockH11( 6, 7 );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 6, 7, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 6, 8 );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 6, 8, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 6, 9 );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 6, 9, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 6, 10 );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 6, 10, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 6, 11 );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 6, 11, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 6, 12 );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 6, 12, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 6, 13 );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 6, 13, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 6, 14 );
acado_setBlockH11( 6, 14, &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 6, 14, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 6, 14, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 6, 14, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 6, 14, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 6, 14, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 6, 15 );
acado_setBlockH11( 6, 15, &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 6, 15, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 6, 15, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 6, 15, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 6, 15, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 6, 16 );
acado_setBlockH11( 6, 16, &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 6, 16, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 6, 16, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 6, 16, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 6, 17 );
acado_setBlockH11( 6, 17, &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 6, 17, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 6, 17, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 6, 18 );
acado_setBlockH11( 6, 18, &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 6, 18, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 6, 19 );
acado_setBlockH11( 6, 19, &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 7, 7, &(acadoWorkspace.R1[ 63 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 630 ]), &(acadoWorkspace.QE[ 630 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.QE[ 774 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.QE[ 936 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.QE[ 1116 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QE[ 1314 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1530 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1764 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2016 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2286 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2574 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2880 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3204 ]) );
acado_setBlockH11( 7, 7, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3546 ]) );

acado_zeroBlockH11( 7, 8 );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 7, 8, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 7, 9 );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 7, 9, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 7, 10 );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 7, 10, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 7, 11 );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 7, 11, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 7, 12 );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 7, 12, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 7, 13 );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 7, 13, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 7, 14 );
acado_setBlockH11( 7, 14, &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 7, 14, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 7, 14, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 7, 14, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 7, 14, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 7, 14, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 7, 15 );
acado_setBlockH11( 7, 15, &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 7, 15, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 7, 15, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 7, 15, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 7, 15, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 7, 16 );
acado_setBlockH11( 7, 16, &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 7, 16, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 7, 16, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 7, 16, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 7, 17 );
acado_setBlockH11( 7, 17, &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 7, 17, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 7, 17, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 7, 18 );
acado_setBlockH11( 7, 18, &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 7, 18, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 7, 19 );
acado_setBlockH11( 7, 19, &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 8, 8, &(acadoWorkspace.R1[ 72 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 792 ]), &(acadoWorkspace.QE[ 792 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 954 ]), &(acadoWorkspace.QE[ 954 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.QE[ 1134 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.QE[ 1332 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QE[ 1548 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1782 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2034 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2304 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2592 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2898 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3222 ]) );
acado_setBlockH11( 8, 8, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3564 ]) );

acado_zeroBlockH11( 8, 9 );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 954 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 8, 9, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 8, 10 );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 8, 10, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 8, 11 );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 8, 11, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 8, 12 );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 8, 12, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 8, 13 );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 8, 13, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 8, 14 );
acado_setBlockH11( 8, 14, &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 8, 14, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 8, 14, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 8, 14, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 8, 14, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 8, 14, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 8, 15 );
acado_setBlockH11( 8, 15, &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 8, 15, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 8, 15, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 8, 15, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 8, 15, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 8, 16 );
acado_setBlockH11( 8, 16, &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 8, 16, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 8, 16, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 8, 16, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 8, 17 );
acado_setBlockH11( 8, 17, &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 8, 17, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 8, 17, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 8, 18 );
acado_setBlockH11( 8, 18, &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 8, 18, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 8, 19 );
acado_setBlockH11( 8, 19, &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 9, 9, &(acadoWorkspace.R1[ 81 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 972 ]), &(acadoWorkspace.QE[ 972 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 1152 ]), &(acadoWorkspace.QE[ 1152 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.QE[ 1350 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.QE[ 1566 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QE[ 1800 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2052 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2322 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2610 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 2916 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3240 ]) );
acado_setBlockH11( 9, 9, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3582 ]) );

acado_zeroBlockH11( 9, 10 );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 1152 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 9, 10, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 9, 11 );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 9, 11, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 9, 12 );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 9, 12, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 9, 13 );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 9, 13, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 9, 14 );
acado_setBlockH11( 9, 14, &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 9, 14, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 9, 14, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 9, 14, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 9, 14, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 9, 14, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 9, 15 );
acado_setBlockH11( 9, 15, &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 9, 15, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 9, 15, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 9, 15, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 9, 15, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 9, 16 );
acado_setBlockH11( 9, 16, &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 9, 16, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 9, 16, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 9, 16, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 9, 17 );
acado_setBlockH11( 9, 17, &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 9, 17, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 9, 17, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 9, 18 );
acado_setBlockH11( 9, 18, &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 9, 18, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 9, 19 );
acado_setBlockH11( 9, 19, &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 10, 10, &(acadoWorkspace.R1[ 90 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 1170 ]), &(acadoWorkspace.QE[ 1170 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 1368 ]), &(acadoWorkspace.QE[ 1368 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.QE[ 1584 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.QE[ 1818 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QE[ 2070 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2340 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2628 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 2934 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3258 ]) );
acado_setBlockH11( 10, 10, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3600 ]) );

acado_zeroBlockH11( 10, 11 );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 1368 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 10, 11, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 10, 12 );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 10, 12, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 10, 13 );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 10, 13, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 10, 14 );
acado_setBlockH11( 10, 14, &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 10, 14, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 10, 14, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 10, 14, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 10, 14, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 10, 14, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 10, 15 );
acado_setBlockH11( 10, 15, &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 10, 15, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 10, 15, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 10, 15, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 10, 15, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 10, 16 );
acado_setBlockH11( 10, 16, &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 10, 16, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 10, 16, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 10, 16, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 10, 17 );
acado_setBlockH11( 10, 17, &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 10, 17, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 10, 17, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 10, 18 );
acado_setBlockH11( 10, 18, &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 10, 18, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 10, 19 );
acado_setBlockH11( 10, 19, &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 11, 11, &(acadoWorkspace.R1[ 99 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 1386 ]), &(acadoWorkspace.QE[ 1386 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 1602 ]), &(acadoWorkspace.QE[ 1602 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.QE[ 1836 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.QE[ 2088 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QE[ 2358 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2646 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 2952 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3276 ]) );
acado_setBlockH11( 11, 11, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3618 ]) );

acado_zeroBlockH11( 11, 12 );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 1602 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 11, 12, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 11, 13 );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 11, 13, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 11, 14 );
acado_setBlockH11( 11, 14, &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 11, 14, &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 11, 14, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 11, 14, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 11, 14, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 11, 14, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 11, 15 );
acado_setBlockH11( 11, 15, &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 11, 15, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 11, 15, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 11, 15, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 11, 15, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 11, 16 );
acado_setBlockH11( 11, 16, &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 11, 16, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 11, 16, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 11, 16, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 11, 17 );
acado_setBlockH11( 11, 17, &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 11, 17, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 11, 17, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 11, 18 );
acado_setBlockH11( 11, 18, &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 11, 18, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 11, 19 );
acado_setBlockH11( 11, 19, &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 12, 12, &(acadoWorkspace.R1[ 108 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 1620 ]), &(acadoWorkspace.QE[ 1620 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 1854 ]), &(acadoWorkspace.QE[ 1854 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.QE[ 2106 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.QE[ 2376 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QE[ 2664 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 2970 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3294 ]) );
acado_setBlockH11( 12, 12, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3636 ]) );

acado_zeroBlockH11( 12, 13 );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 1854 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 12, 13, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 12, 14 );
acado_setBlockH11( 12, 14, &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 12, 14, &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 12, 14, &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 12, 14, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 12, 14, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 12, 14, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 12, 15 );
acado_setBlockH11( 12, 15, &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 12, 15, &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 12, 15, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 12, 15, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 12, 15, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 12, 16 );
acado_setBlockH11( 12, 16, &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 12, 16, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 12, 16, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 12, 16, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 12, 17 );
acado_setBlockH11( 12, 17, &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 12, 17, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 12, 17, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 12, 18 );
acado_setBlockH11( 12, 18, &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 12, 18, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 12, 19 );
acado_setBlockH11( 12, 19, &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 13, 13, &(acadoWorkspace.R1[ 117 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 1872 ]), &(acadoWorkspace.QE[ 1872 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 2124 ]), &(acadoWorkspace.QE[ 2124 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.QE[ 2394 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.QE[ 2682 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QE[ 2988 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3312 ]) );
acado_setBlockH11( 13, 13, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3654 ]) );

acado_zeroBlockH11( 13, 14 );
acado_setBlockH11( 13, 14, &(acadoWorkspace.E[ 2124 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 13, 14, &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 13, 14, &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 13, 14, &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 13, 14, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 13, 14, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 13, 15 );
acado_setBlockH11( 13, 15, &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 13, 15, &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 13, 15, &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 13, 15, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 13, 15, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 13, 16 );
acado_setBlockH11( 13, 16, &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 13, 16, &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 13, 16, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 13, 16, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 13, 17 );
acado_setBlockH11( 13, 17, &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 13, 17, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 13, 17, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 13, 18 );
acado_setBlockH11( 13, 18, &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 13, 18, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 13, 19 );
acado_setBlockH11( 13, 19, &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 14, 14, &(acadoWorkspace.R1[ 126 ]) );
acado_setBlockH11( 14, 14, &(acadoWorkspace.E[ 2142 ]), &(acadoWorkspace.QE[ 2142 ]) );
acado_setBlockH11( 14, 14, &(acadoWorkspace.E[ 2412 ]), &(acadoWorkspace.QE[ 2412 ]) );
acado_setBlockH11( 14, 14, &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.QE[ 2700 ]) );
acado_setBlockH11( 14, 14, &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.QE[ 3006 ]) );
acado_setBlockH11( 14, 14, &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QE[ 3330 ]) );
acado_setBlockH11( 14, 14, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3672 ]) );

acado_zeroBlockH11( 14, 15 );
acado_setBlockH11( 14, 15, &(acadoWorkspace.E[ 2412 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 14, 15, &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 14, 15, &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 14, 15, &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 14, 15, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 14, 16 );
acado_setBlockH11( 14, 16, &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 14, 16, &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 14, 16, &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 14, 16, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 14, 17 );
acado_setBlockH11( 14, 17, &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 14, 17, &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 14, 17, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 14, 18 );
acado_setBlockH11( 14, 18, &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 14, 18, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 14, 19 );
acado_setBlockH11( 14, 19, &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 15, 15, &(acadoWorkspace.R1[ 135 ]) );
acado_setBlockH11( 15, 15, &(acadoWorkspace.E[ 2430 ]), &(acadoWorkspace.QE[ 2430 ]) );
acado_setBlockH11( 15, 15, &(acadoWorkspace.E[ 2718 ]), &(acadoWorkspace.QE[ 2718 ]) );
acado_setBlockH11( 15, 15, &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.QE[ 3024 ]) );
acado_setBlockH11( 15, 15, &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.QE[ 3348 ]) );
acado_setBlockH11( 15, 15, &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QE[ 3690 ]) );

acado_zeroBlockH11( 15, 16 );
acado_setBlockH11( 15, 16, &(acadoWorkspace.E[ 2718 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 15, 16, &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 15, 16, &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 15, 16, &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 15, 17 );
acado_setBlockH11( 15, 17, &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 15, 17, &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 15, 17, &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 15, 18 );
acado_setBlockH11( 15, 18, &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 15, 18, &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 15, 19 );
acado_setBlockH11( 15, 19, &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 16, 16, &(acadoWorkspace.R1[ 144 ]) );
acado_setBlockH11( 16, 16, &(acadoWorkspace.E[ 2736 ]), &(acadoWorkspace.QE[ 2736 ]) );
acado_setBlockH11( 16, 16, &(acadoWorkspace.E[ 3042 ]), &(acadoWorkspace.QE[ 3042 ]) );
acado_setBlockH11( 16, 16, &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.QE[ 3366 ]) );
acado_setBlockH11( 16, 16, &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.QE[ 3708 ]) );

acado_zeroBlockH11( 16, 17 );
acado_setBlockH11( 16, 17, &(acadoWorkspace.E[ 3042 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 16, 17, &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 16, 17, &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 16, 18 );
acado_setBlockH11( 16, 18, &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 16, 18, &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 16, 19 );
acado_setBlockH11( 16, 19, &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 17, 17, &(acadoWorkspace.R1[ 153 ]) );
acado_setBlockH11( 17, 17, &(acadoWorkspace.E[ 3060 ]), &(acadoWorkspace.QE[ 3060 ]) );
acado_setBlockH11( 17, 17, &(acadoWorkspace.E[ 3384 ]), &(acadoWorkspace.QE[ 3384 ]) );
acado_setBlockH11( 17, 17, &(acadoWorkspace.E[ 3726 ]), &(acadoWorkspace.QE[ 3726 ]) );

acado_zeroBlockH11( 17, 18 );
acado_setBlockH11( 17, 18, &(acadoWorkspace.E[ 3384 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 17, 18, &(acadoWorkspace.E[ 3726 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 17, 19 );
acado_setBlockH11( 17, 19, &(acadoWorkspace.E[ 3726 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 18, 18, &(acadoWorkspace.R1[ 162 ]) );
acado_setBlockH11( 18, 18, &(acadoWorkspace.E[ 3402 ]), &(acadoWorkspace.QE[ 3402 ]) );
acado_setBlockH11( 18, 18, &(acadoWorkspace.E[ 3744 ]), &(acadoWorkspace.QE[ 3744 ]) );

acado_zeroBlockH11( 18, 19 );
acado_setBlockH11( 18, 19, &(acadoWorkspace.E[ 3744 ]), &(acadoWorkspace.QE[ 3762 ]) );

acado_setBlockH11_R1( 19, 19, &(acadoWorkspace.R1[ 171 ]) );
acado_setBlockH11( 19, 19, &(acadoWorkspace.E[ 3762 ]), &(acadoWorkspace.QE[ 3762 ]) );


acado_copyHTH( 1, 0 );
acado_copyHTH( 2, 0 );
acado_copyHTH( 2, 1 );
acado_copyHTH( 3, 0 );
acado_copyHTH( 3, 1 );
acado_copyHTH( 3, 2 );
acado_copyHTH( 4, 0 );
acado_copyHTH( 4, 1 );
acado_copyHTH( 4, 2 );
acado_copyHTH( 4, 3 );
acado_copyHTH( 5, 0 );
acado_copyHTH( 5, 1 );
acado_copyHTH( 5, 2 );
acado_copyHTH( 5, 3 );
acado_copyHTH( 5, 4 );
acado_copyHTH( 6, 0 );
acado_copyHTH( 6, 1 );
acado_copyHTH( 6, 2 );
acado_copyHTH( 6, 3 );
acado_copyHTH( 6, 4 );
acado_copyHTH( 6, 5 );
acado_copyHTH( 7, 0 );
acado_copyHTH( 7, 1 );
acado_copyHTH( 7, 2 );
acado_copyHTH( 7, 3 );
acado_copyHTH( 7, 4 );
acado_copyHTH( 7, 5 );
acado_copyHTH( 7, 6 );
acado_copyHTH( 8, 0 );
acado_copyHTH( 8, 1 );
acado_copyHTH( 8, 2 );
acado_copyHTH( 8, 3 );
acado_copyHTH( 8, 4 );
acado_copyHTH( 8, 5 );
acado_copyHTH( 8, 6 );
acado_copyHTH( 8, 7 );
acado_copyHTH( 9, 0 );
acado_copyHTH( 9, 1 );
acado_copyHTH( 9, 2 );
acado_copyHTH( 9, 3 );
acado_copyHTH( 9, 4 );
acado_copyHTH( 9, 5 );
acado_copyHTH( 9, 6 );
acado_copyHTH( 9, 7 );
acado_copyHTH( 9, 8 );
acado_copyHTH( 10, 0 );
acado_copyHTH( 10, 1 );
acado_copyHTH( 10, 2 );
acado_copyHTH( 10, 3 );
acado_copyHTH( 10, 4 );
acado_copyHTH( 10, 5 );
acado_copyHTH( 10, 6 );
acado_copyHTH( 10, 7 );
acado_copyHTH( 10, 8 );
acado_copyHTH( 10, 9 );
acado_copyHTH( 11, 0 );
acado_copyHTH( 11, 1 );
acado_copyHTH( 11, 2 );
acado_copyHTH( 11, 3 );
acado_copyHTH( 11, 4 );
acado_copyHTH( 11, 5 );
acado_copyHTH( 11, 6 );
acado_copyHTH( 11, 7 );
acado_copyHTH( 11, 8 );
acado_copyHTH( 11, 9 );
acado_copyHTH( 11, 10 );
acado_copyHTH( 12, 0 );
acado_copyHTH( 12, 1 );
acado_copyHTH( 12, 2 );
acado_copyHTH( 12, 3 );
acado_copyHTH( 12, 4 );
acado_copyHTH( 12, 5 );
acado_copyHTH( 12, 6 );
acado_copyHTH( 12, 7 );
acado_copyHTH( 12, 8 );
acado_copyHTH( 12, 9 );
acado_copyHTH( 12, 10 );
acado_copyHTH( 12, 11 );
acado_copyHTH( 13, 0 );
acado_copyHTH( 13, 1 );
acado_copyHTH( 13, 2 );
acado_copyHTH( 13, 3 );
acado_copyHTH( 13, 4 );
acado_copyHTH( 13, 5 );
acado_copyHTH( 13, 6 );
acado_copyHTH( 13, 7 );
acado_copyHTH( 13, 8 );
acado_copyHTH( 13, 9 );
acado_copyHTH( 13, 10 );
acado_copyHTH( 13, 11 );
acado_copyHTH( 13, 12 );
acado_copyHTH( 14, 0 );
acado_copyHTH( 14, 1 );
acado_copyHTH( 14, 2 );
acado_copyHTH( 14, 3 );
acado_copyHTH( 14, 4 );
acado_copyHTH( 14, 5 );
acado_copyHTH( 14, 6 );
acado_copyHTH( 14, 7 );
acado_copyHTH( 14, 8 );
acado_copyHTH( 14, 9 );
acado_copyHTH( 14, 10 );
acado_copyHTH( 14, 11 );
acado_copyHTH( 14, 12 );
acado_copyHTH( 14, 13 );
acado_copyHTH( 15, 0 );
acado_copyHTH( 15, 1 );
acado_copyHTH( 15, 2 );
acado_copyHTH( 15, 3 );
acado_copyHTH( 15, 4 );
acado_copyHTH( 15, 5 );
acado_copyHTH( 15, 6 );
acado_copyHTH( 15, 7 );
acado_copyHTH( 15, 8 );
acado_copyHTH( 15, 9 );
acado_copyHTH( 15, 10 );
acado_copyHTH( 15, 11 );
acado_copyHTH( 15, 12 );
acado_copyHTH( 15, 13 );
acado_copyHTH( 15, 14 );
acado_copyHTH( 16, 0 );
acado_copyHTH( 16, 1 );
acado_copyHTH( 16, 2 );
acado_copyHTH( 16, 3 );
acado_copyHTH( 16, 4 );
acado_copyHTH( 16, 5 );
acado_copyHTH( 16, 6 );
acado_copyHTH( 16, 7 );
acado_copyHTH( 16, 8 );
acado_copyHTH( 16, 9 );
acado_copyHTH( 16, 10 );
acado_copyHTH( 16, 11 );
acado_copyHTH( 16, 12 );
acado_copyHTH( 16, 13 );
acado_copyHTH( 16, 14 );
acado_copyHTH( 16, 15 );
acado_copyHTH( 17, 0 );
acado_copyHTH( 17, 1 );
acado_copyHTH( 17, 2 );
acado_copyHTH( 17, 3 );
acado_copyHTH( 17, 4 );
acado_copyHTH( 17, 5 );
acado_copyHTH( 17, 6 );
acado_copyHTH( 17, 7 );
acado_copyHTH( 17, 8 );
acado_copyHTH( 17, 9 );
acado_copyHTH( 17, 10 );
acado_copyHTH( 17, 11 );
acado_copyHTH( 17, 12 );
acado_copyHTH( 17, 13 );
acado_copyHTH( 17, 14 );
acado_copyHTH( 17, 15 );
acado_copyHTH( 17, 16 );
acado_copyHTH( 18, 0 );
acado_copyHTH( 18, 1 );
acado_copyHTH( 18, 2 );
acado_copyHTH( 18, 3 );
acado_copyHTH( 18, 4 );
acado_copyHTH( 18, 5 );
acado_copyHTH( 18, 6 );
acado_copyHTH( 18, 7 );
acado_copyHTH( 18, 8 );
acado_copyHTH( 18, 9 );
acado_copyHTH( 18, 10 );
acado_copyHTH( 18, 11 );
acado_copyHTH( 18, 12 );
acado_copyHTH( 18, 13 );
acado_copyHTH( 18, 14 );
acado_copyHTH( 18, 15 );
acado_copyHTH( 18, 16 );
acado_copyHTH( 18, 17 );
acado_copyHTH( 19, 0 );
acado_copyHTH( 19, 1 );
acado_copyHTH( 19, 2 );
acado_copyHTH( 19, 3 );
acado_copyHTH( 19, 4 );
acado_copyHTH( 19, 5 );
acado_copyHTH( 19, 6 );
acado_copyHTH( 19, 7 );
acado_copyHTH( 19, 8 );
acado_copyHTH( 19, 9 );
acado_copyHTH( 19, 10 );
acado_copyHTH( 19, 11 );
acado_copyHTH( 19, 12 );
acado_copyHTH( 19, 13 );
acado_copyHTH( 19, 14 );
acado_copyHTH( 19, 15 );
acado_copyHTH( 19, 16 );
acado_copyHTH( 19, 17 );
acado_copyHTH( 19, 18 );

acado_multQ1d( &(acadoWorkspace.Q1[ 36 ]), acadoWorkspace.d, acadoWorkspace.Qd );
acado_multQ1d( &(acadoWorkspace.Q1[ 72 ]), &(acadoWorkspace.d[ 6 ]), &(acadoWorkspace.Qd[ 6 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 108 ]), &(acadoWorkspace.d[ 12 ]), &(acadoWorkspace.Qd[ 12 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 144 ]), &(acadoWorkspace.d[ 18 ]), &(acadoWorkspace.Qd[ 18 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 180 ]), &(acadoWorkspace.d[ 24 ]), &(acadoWorkspace.Qd[ 24 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 216 ]), &(acadoWorkspace.d[ 30 ]), &(acadoWorkspace.Qd[ 30 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 252 ]), &(acadoWorkspace.d[ 36 ]), &(acadoWorkspace.Qd[ 36 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 288 ]), &(acadoWorkspace.d[ 42 ]), &(acadoWorkspace.Qd[ 42 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 324 ]), &(acadoWorkspace.d[ 48 ]), &(acadoWorkspace.Qd[ 48 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 360 ]), &(acadoWorkspace.d[ 54 ]), &(acadoWorkspace.Qd[ 54 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 396 ]), &(acadoWorkspace.d[ 60 ]), &(acadoWorkspace.Qd[ 60 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 432 ]), &(acadoWorkspace.d[ 66 ]), &(acadoWorkspace.Qd[ 66 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 468 ]), &(acadoWorkspace.d[ 72 ]), &(acadoWorkspace.Qd[ 72 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 504 ]), &(acadoWorkspace.d[ 78 ]), &(acadoWorkspace.Qd[ 78 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 540 ]), &(acadoWorkspace.d[ 84 ]), &(acadoWorkspace.Qd[ 84 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 576 ]), &(acadoWorkspace.d[ 90 ]), &(acadoWorkspace.Qd[ 90 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 612 ]), &(acadoWorkspace.d[ 96 ]), &(acadoWorkspace.Qd[ 96 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 648 ]), &(acadoWorkspace.d[ 102 ]), &(acadoWorkspace.Qd[ 102 ]) );
acado_multQ1d( &(acadoWorkspace.Q1[ 684 ]), &(acadoWorkspace.d[ 108 ]), &(acadoWorkspace.Qd[ 108 ]) );
acado_multQN1d( acadoWorkspace.QN1, &(acadoWorkspace.d[ 114 ]), &(acadoWorkspace.Qd[ 114 ]) );

acado_macETSlu( acadoWorkspace.QE, acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 18 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 54 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 108 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 180 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 270 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 378 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 504 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 648 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 810 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 990 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 1188 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 1404 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 1638 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 1890 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 2160 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 2448 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 2754 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 3078 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 3420 ]), acadoWorkspace.g );
acado_macETSlu( &(acadoWorkspace.QE[ 36 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 72 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 126 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 198 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 288 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 396 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 522 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 666 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 828 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1008 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1206 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1422 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1656 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1908 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2178 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2466 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2772 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3096 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3438 ]), &(acadoWorkspace.g[ 3 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 90 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 144 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 216 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 306 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 414 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 540 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 684 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 846 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1026 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1224 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1440 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1674 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1926 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2196 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2484 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2790 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3114 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3456 ]), &(acadoWorkspace.g[ 6 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 162 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 234 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 324 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 432 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 558 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 702 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 864 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1044 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1242 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1458 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1692 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1944 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2214 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2502 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2808 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3132 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3474 ]), &(acadoWorkspace.g[ 9 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 252 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 342 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 450 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 576 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 720 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 882 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1062 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1260 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1476 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1710 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1962 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2232 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2520 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2826 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3150 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3492 ]), &(acadoWorkspace.g[ 12 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 360 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 468 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 594 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 738 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 900 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1080 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1278 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1494 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1728 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1980 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2250 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2538 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2844 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3168 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3510 ]), &(acadoWorkspace.g[ 15 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 486 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 612 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 756 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 918 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1098 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1296 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1512 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1746 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1998 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2268 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2556 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2862 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3186 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3528 ]), &(acadoWorkspace.g[ 18 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 630 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 774 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 936 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1116 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1314 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1530 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1764 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2016 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2286 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2574 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2880 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3204 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3546 ]), &(acadoWorkspace.g[ 21 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 792 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 954 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1134 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1332 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1548 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1782 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2034 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2304 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2592 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2898 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3222 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3564 ]), &(acadoWorkspace.g[ 24 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 972 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1152 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1350 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1566 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1800 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2052 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2322 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2610 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2916 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3240 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3582 ]), &(acadoWorkspace.g[ 27 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1170 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1368 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1584 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1818 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2070 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2340 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2628 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2934 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3258 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3600 ]), &(acadoWorkspace.g[ 30 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1386 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1602 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1836 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2088 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2358 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2646 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2952 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3276 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3618 ]), &(acadoWorkspace.g[ 33 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1620 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1854 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2106 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2376 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2664 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2970 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3294 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3636 ]), &(acadoWorkspace.g[ 36 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 1872 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2124 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2394 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2682 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2988 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3312 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3654 ]), &(acadoWorkspace.g[ 39 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2142 ]), &(acadoWorkspace.g[ 42 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2412 ]), &(acadoWorkspace.g[ 42 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2700 ]), &(acadoWorkspace.g[ 42 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3006 ]), &(acadoWorkspace.g[ 42 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3330 ]), &(acadoWorkspace.g[ 42 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3672 ]), &(acadoWorkspace.g[ 42 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2430 ]), &(acadoWorkspace.g[ 45 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2718 ]), &(acadoWorkspace.g[ 45 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3024 ]), &(acadoWorkspace.g[ 45 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3348 ]), &(acadoWorkspace.g[ 45 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3690 ]), &(acadoWorkspace.g[ 45 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 2736 ]), &(acadoWorkspace.g[ 48 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3042 ]), &(acadoWorkspace.g[ 48 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3366 ]), &(acadoWorkspace.g[ 48 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3708 ]), &(acadoWorkspace.g[ 48 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3060 ]), &(acadoWorkspace.g[ 51 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3384 ]), &(acadoWorkspace.g[ 51 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3726 ]), &(acadoWorkspace.g[ 51 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3402 ]), &(acadoWorkspace.g[ 54 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3744 ]), &(acadoWorkspace.g[ 54 ]) );
acado_macETSlu( &(acadoWorkspace.QE[ 3762 ]), &(acadoWorkspace.g[ 57 ]) );
}

void acado_condenseFdb(  )
{
int lRun1;
acadoWorkspace.Dx0[0] = acadoVariables.x0[0] - acadoVariables.x[0];
acadoWorkspace.Dx0[1] = acadoVariables.x0[1] - acadoVariables.x[1];
acadoWorkspace.Dx0[2] = acadoVariables.x0[2] - acadoVariables.x[2];
acadoWorkspace.Dx0[3] = acadoVariables.x0[3] - acadoVariables.x[3];
acadoWorkspace.Dx0[4] = acadoVariables.x0[4] - acadoVariables.x[4];
acadoWorkspace.Dx0[5] = acadoVariables.x0[5] - acadoVariables.x[5];

for (lRun1 = 0; lRun1 < 180; ++lRun1)
acadoWorkspace.Dy[lRun1] -= acadoVariables.y[lRun1];

acadoWorkspace.DyN[0] -= acadoVariables.yN[0];
acadoWorkspace.DyN[1] -= acadoVariables.yN[1];
acadoWorkspace.DyN[2] -= acadoVariables.yN[2];
acadoWorkspace.DyN[3] -= acadoVariables.yN[3];
acadoWorkspace.DyN[4] -= acadoVariables.yN[4];
acadoWorkspace.DyN[5] -= acadoVariables.yN[5];

acado_multRDy( acadoWorkspace.R2, acadoWorkspace.Dy, acadoWorkspace.g );
acado_multRDy( &(acadoWorkspace.R2[ 27 ]), &(acadoWorkspace.Dy[ 9 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 54 ]), &(acadoWorkspace.Dy[ 18 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 81 ]), &(acadoWorkspace.Dy[ 27 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 108 ]), &(acadoWorkspace.Dy[ 36 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 135 ]), &(acadoWorkspace.Dy[ 45 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 162 ]), &(acadoWorkspace.Dy[ 54 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 189 ]), &(acadoWorkspace.Dy[ 63 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 216 ]), &(acadoWorkspace.Dy[ 72 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 243 ]), &(acadoWorkspace.Dy[ 81 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 270 ]), &(acadoWorkspace.Dy[ 90 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 297 ]), &(acadoWorkspace.Dy[ 99 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 324 ]), &(acadoWorkspace.Dy[ 108 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 351 ]), &(acadoWorkspace.Dy[ 117 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 378 ]), &(acadoWorkspace.Dy[ 126 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 405 ]), &(acadoWorkspace.Dy[ 135 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 432 ]), &(acadoWorkspace.Dy[ 144 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 459 ]), &(acadoWorkspace.Dy[ 153 ]), &(acadoWorkspace.g[ 51 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 486 ]), &(acadoWorkspace.Dy[ 162 ]), &(acadoWorkspace.g[ 54 ]) );
acado_multRDy( &(acadoWorkspace.R2[ 513 ]), &(acadoWorkspace.Dy[ 171 ]), &(acadoWorkspace.g[ 57 ]) );

acado_multQDy( acadoWorkspace.Q2, acadoWorkspace.Dy, acadoWorkspace.QDy );
acado_multQDy( &(acadoWorkspace.Q2[ 54 ]), &(acadoWorkspace.Dy[ 9 ]), &(acadoWorkspace.QDy[ 6 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 108 ]), &(acadoWorkspace.Dy[ 18 ]), &(acadoWorkspace.QDy[ 12 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 162 ]), &(acadoWorkspace.Dy[ 27 ]), &(acadoWorkspace.QDy[ 18 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 216 ]), &(acadoWorkspace.Dy[ 36 ]), &(acadoWorkspace.QDy[ 24 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 270 ]), &(acadoWorkspace.Dy[ 45 ]), &(acadoWorkspace.QDy[ 30 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 324 ]), &(acadoWorkspace.Dy[ 54 ]), &(acadoWorkspace.QDy[ 36 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 378 ]), &(acadoWorkspace.Dy[ 63 ]), &(acadoWorkspace.QDy[ 42 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 432 ]), &(acadoWorkspace.Dy[ 72 ]), &(acadoWorkspace.QDy[ 48 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 486 ]), &(acadoWorkspace.Dy[ 81 ]), &(acadoWorkspace.QDy[ 54 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 540 ]), &(acadoWorkspace.Dy[ 90 ]), &(acadoWorkspace.QDy[ 60 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 594 ]), &(acadoWorkspace.Dy[ 99 ]), &(acadoWorkspace.QDy[ 66 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 648 ]), &(acadoWorkspace.Dy[ 108 ]), &(acadoWorkspace.QDy[ 72 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 702 ]), &(acadoWorkspace.Dy[ 117 ]), &(acadoWorkspace.QDy[ 78 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 756 ]), &(acadoWorkspace.Dy[ 126 ]), &(acadoWorkspace.QDy[ 84 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 810 ]), &(acadoWorkspace.Dy[ 135 ]), &(acadoWorkspace.QDy[ 90 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 864 ]), &(acadoWorkspace.Dy[ 144 ]), &(acadoWorkspace.QDy[ 96 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 918 ]), &(acadoWorkspace.Dy[ 153 ]), &(acadoWorkspace.QDy[ 102 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 972 ]), &(acadoWorkspace.Dy[ 162 ]), &(acadoWorkspace.QDy[ 108 ]) );
acado_multQDy( &(acadoWorkspace.Q2[ 1026 ]), &(acadoWorkspace.Dy[ 171 ]), &(acadoWorkspace.QDy[ 114 ]) );

acadoWorkspace.QDy[120] = + acadoWorkspace.QN2[0]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[1]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[2]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[3]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[4]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[5]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[121] = + acadoWorkspace.QN2[6]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[7]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[8]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[9]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[10]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[11]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[122] = + acadoWorkspace.QN2[12]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[13]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[14]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[15]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[16]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[17]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[123] = + acadoWorkspace.QN2[18]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[19]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[20]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[21]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[22]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[23]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[124] = + acadoWorkspace.QN2[24]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[25]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[26]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[27]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[28]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[29]*acadoWorkspace.DyN[5];
acadoWorkspace.QDy[125] = + acadoWorkspace.QN2[30]*acadoWorkspace.DyN[0] + acadoWorkspace.QN2[31]*acadoWorkspace.DyN[1] + acadoWorkspace.QN2[32]*acadoWorkspace.DyN[2] + acadoWorkspace.QN2[33]*acadoWorkspace.DyN[3] + acadoWorkspace.QN2[34]*acadoWorkspace.DyN[4] + acadoWorkspace.QN2[35]*acadoWorkspace.DyN[5];

acadoWorkspace.QDy[6] += acadoWorkspace.Qd[0];
acadoWorkspace.QDy[7] += acadoWorkspace.Qd[1];
acadoWorkspace.QDy[8] += acadoWorkspace.Qd[2];
acadoWorkspace.QDy[9] += acadoWorkspace.Qd[3];
acadoWorkspace.QDy[10] += acadoWorkspace.Qd[4];
acadoWorkspace.QDy[11] += acadoWorkspace.Qd[5];
acadoWorkspace.QDy[12] += acadoWorkspace.Qd[6];
acadoWorkspace.QDy[13] += acadoWorkspace.Qd[7];
acadoWorkspace.QDy[14] += acadoWorkspace.Qd[8];
acadoWorkspace.QDy[15] += acadoWorkspace.Qd[9];
acadoWorkspace.QDy[16] += acadoWorkspace.Qd[10];
acadoWorkspace.QDy[17] += acadoWorkspace.Qd[11];
acadoWorkspace.QDy[18] += acadoWorkspace.Qd[12];
acadoWorkspace.QDy[19] += acadoWorkspace.Qd[13];
acadoWorkspace.QDy[20] += acadoWorkspace.Qd[14];
acadoWorkspace.QDy[21] += acadoWorkspace.Qd[15];
acadoWorkspace.QDy[22] += acadoWorkspace.Qd[16];
acadoWorkspace.QDy[23] += acadoWorkspace.Qd[17];
acadoWorkspace.QDy[24] += acadoWorkspace.Qd[18];
acadoWorkspace.QDy[25] += acadoWorkspace.Qd[19];
acadoWorkspace.QDy[26] += acadoWorkspace.Qd[20];
acadoWorkspace.QDy[27] += acadoWorkspace.Qd[21];
acadoWorkspace.QDy[28] += acadoWorkspace.Qd[22];
acadoWorkspace.QDy[29] += acadoWorkspace.Qd[23];
acadoWorkspace.QDy[30] += acadoWorkspace.Qd[24];
acadoWorkspace.QDy[31] += acadoWorkspace.Qd[25];
acadoWorkspace.QDy[32] += acadoWorkspace.Qd[26];
acadoWorkspace.QDy[33] += acadoWorkspace.Qd[27];
acadoWorkspace.QDy[34] += acadoWorkspace.Qd[28];
acadoWorkspace.QDy[35] += acadoWorkspace.Qd[29];
acadoWorkspace.QDy[36] += acadoWorkspace.Qd[30];
acadoWorkspace.QDy[37] += acadoWorkspace.Qd[31];
acadoWorkspace.QDy[38] += acadoWorkspace.Qd[32];
acadoWorkspace.QDy[39] += acadoWorkspace.Qd[33];
acadoWorkspace.QDy[40] += acadoWorkspace.Qd[34];
acadoWorkspace.QDy[41] += acadoWorkspace.Qd[35];
acadoWorkspace.QDy[42] += acadoWorkspace.Qd[36];
acadoWorkspace.QDy[43] += acadoWorkspace.Qd[37];
acadoWorkspace.QDy[44] += acadoWorkspace.Qd[38];
acadoWorkspace.QDy[45] += acadoWorkspace.Qd[39];
acadoWorkspace.QDy[46] += acadoWorkspace.Qd[40];
acadoWorkspace.QDy[47] += acadoWorkspace.Qd[41];
acadoWorkspace.QDy[48] += acadoWorkspace.Qd[42];
acadoWorkspace.QDy[49] += acadoWorkspace.Qd[43];
acadoWorkspace.QDy[50] += acadoWorkspace.Qd[44];
acadoWorkspace.QDy[51] += acadoWorkspace.Qd[45];
acadoWorkspace.QDy[52] += acadoWorkspace.Qd[46];
acadoWorkspace.QDy[53] += acadoWorkspace.Qd[47];
acadoWorkspace.QDy[54] += acadoWorkspace.Qd[48];
acadoWorkspace.QDy[55] += acadoWorkspace.Qd[49];
acadoWorkspace.QDy[56] += acadoWorkspace.Qd[50];
acadoWorkspace.QDy[57] += acadoWorkspace.Qd[51];
acadoWorkspace.QDy[58] += acadoWorkspace.Qd[52];
acadoWorkspace.QDy[59] += acadoWorkspace.Qd[53];
acadoWorkspace.QDy[60] += acadoWorkspace.Qd[54];
acadoWorkspace.QDy[61] += acadoWorkspace.Qd[55];
acadoWorkspace.QDy[62] += acadoWorkspace.Qd[56];
acadoWorkspace.QDy[63] += acadoWorkspace.Qd[57];
acadoWorkspace.QDy[64] += acadoWorkspace.Qd[58];
acadoWorkspace.QDy[65] += acadoWorkspace.Qd[59];
acadoWorkspace.QDy[66] += acadoWorkspace.Qd[60];
acadoWorkspace.QDy[67] += acadoWorkspace.Qd[61];
acadoWorkspace.QDy[68] += acadoWorkspace.Qd[62];
acadoWorkspace.QDy[69] += acadoWorkspace.Qd[63];
acadoWorkspace.QDy[70] += acadoWorkspace.Qd[64];
acadoWorkspace.QDy[71] += acadoWorkspace.Qd[65];
acadoWorkspace.QDy[72] += acadoWorkspace.Qd[66];
acadoWorkspace.QDy[73] += acadoWorkspace.Qd[67];
acadoWorkspace.QDy[74] += acadoWorkspace.Qd[68];
acadoWorkspace.QDy[75] += acadoWorkspace.Qd[69];
acadoWorkspace.QDy[76] += acadoWorkspace.Qd[70];
acadoWorkspace.QDy[77] += acadoWorkspace.Qd[71];
acadoWorkspace.QDy[78] += acadoWorkspace.Qd[72];
acadoWorkspace.QDy[79] += acadoWorkspace.Qd[73];
acadoWorkspace.QDy[80] += acadoWorkspace.Qd[74];
acadoWorkspace.QDy[81] += acadoWorkspace.Qd[75];
acadoWorkspace.QDy[82] += acadoWorkspace.Qd[76];
acadoWorkspace.QDy[83] += acadoWorkspace.Qd[77];
acadoWorkspace.QDy[84] += acadoWorkspace.Qd[78];
acadoWorkspace.QDy[85] += acadoWorkspace.Qd[79];
acadoWorkspace.QDy[86] += acadoWorkspace.Qd[80];
acadoWorkspace.QDy[87] += acadoWorkspace.Qd[81];
acadoWorkspace.QDy[88] += acadoWorkspace.Qd[82];
acadoWorkspace.QDy[89] += acadoWorkspace.Qd[83];
acadoWorkspace.QDy[90] += acadoWorkspace.Qd[84];
acadoWorkspace.QDy[91] += acadoWorkspace.Qd[85];
acadoWorkspace.QDy[92] += acadoWorkspace.Qd[86];
acadoWorkspace.QDy[93] += acadoWorkspace.Qd[87];
acadoWorkspace.QDy[94] += acadoWorkspace.Qd[88];
acadoWorkspace.QDy[95] += acadoWorkspace.Qd[89];
acadoWorkspace.QDy[96] += acadoWorkspace.Qd[90];
acadoWorkspace.QDy[97] += acadoWorkspace.Qd[91];
acadoWorkspace.QDy[98] += acadoWorkspace.Qd[92];
acadoWorkspace.QDy[99] += acadoWorkspace.Qd[93];
acadoWorkspace.QDy[100] += acadoWorkspace.Qd[94];
acadoWorkspace.QDy[101] += acadoWorkspace.Qd[95];
acadoWorkspace.QDy[102] += acadoWorkspace.Qd[96];
acadoWorkspace.QDy[103] += acadoWorkspace.Qd[97];
acadoWorkspace.QDy[104] += acadoWorkspace.Qd[98];
acadoWorkspace.QDy[105] += acadoWorkspace.Qd[99];
acadoWorkspace.QDy[106] += acadoWorkspace.Qd[100];
acadoWorkspace.QDy[107] += acadoWorkspace.Qd[101];
acadoWorkspace.QDy[108] += acadoWorkspace.Qd[102];
acadoWorkspace.QDy[109] += acadoWorkspace.Qd[103];
acadoWorkspace.QDy[110] += acadoWorkspace.Qd[104];
acadoWorkspace.QDy[111] += acadoWorkspace.Qd[105];
acadoWorkspace.QDy[112] += acadoWorkspace.Qd[106];
acadoWorkspace.QDy[113] += acadoWorkspace.Qd[107];
acadoWorkspace.QDy[114] += acadoWorkspace.Qd[108];
acadoWorkspace.QDy[115] += acadoWorkspace.Qd[109];
acadoWorkspace.QDy[116] += acadoWorkspace.Qd[110];
acadoWorkspace.QDy[117] += acadoWorkspace.Qd[111];
acadoWorkspace.QDy[118] += acadoWorkspace.Qd[112];
acadoWorkspace.QDy[119] += acadoWorkspace.Qd[113];
acadoWorkspace.QDy[120] += acadoWorkspace.Qd[114];
acadoWorkspace.QDy[121] += acadoWorkspace.Qd[115];
acadoWorkspace.QDy[122] += acadoWorkspace.Qd[116];
acadoWorkspace.QDy[123] += acadoWorkspace.Qd[117];
acadoWorkspace.QDy[124] += acadoWorkspace.Qd[118];
acadoWorkspace.QDy[125] += acadoWorkspace.Qd[119];

acado_multEQDy( acadoWorkspace.E, &(acadoWorkspace.QDy[ 6 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 18 ]), &(acadoWorkspace.QDy[ 12 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 54 ]), &(acadoWorkspace.QDy[ 18 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 108 ]), &(acadoWorkspace.QDy[ 24 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 180 ]), &(acadoWorkspace.QDy[ 30 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 270 ]), &(acadoWorkspace.QDy[ 36 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 378 ]), &(acadoWorkspace.QDy[ 42 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 504 ]), &(acadoWorkspace.QDy[ 48 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 648 ]), &(acadoWorkspace.QDy[ 54 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 810 ]), &(acadoWorkspace.QDy[ 60 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 990 ]), &(acadoWorkspace.QDy[ 66 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 1188 ]), &(acadoWorkspace.QDy[ 72 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 1404 ]), &(acadoWorkspace.QDy[ 78 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 1638 ]), &(acadoWorkspace.QDy[ 84 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 1890 ]), &(acadoWorkspace.QDy[ 90 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 2160 ]), &(acadoWorkspace.QDy[ 96 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 2448 ]), &(acadoWorkspace.QDy[ 102 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 2754 ]), &(acadoWorkspace.QDy[ 108 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 3078 ]), &(acadoWorkspace.QDy[ 114 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 3420 ]), &(acadoWorkspace.QDy[ 120 ]), acadoWorkspace.g );
acado_multEQDy( &(acadoWorkspace.E[ 36 ]), &(acadoWorkspace.QDy[ 12 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.QDy[ 18 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.QDy[ 24 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.QDy[ 30 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.QDy[ 36 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.QDy[ 42 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 3 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 90 ]), &(acadoWorkspace.QDy[ 18 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.QDy[ 24 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.QDy[ 30 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.QDy[ 36 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.QDy[ 42 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 6 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 162 ]), &(acadoWorkspace.QDy[ 24 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.QDy[ 30 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.QDy[ 36 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.QDy[ 42 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 9 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 252 ]), &(acadoWorkspace.QDy[ 30 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.QDy[ 36 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.QDy[ 42 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 12 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 360 ]), &(acadoWorkspace.QDy[ 36 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 468 ]), &(acadoWorkspace.QDy[ 42 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 15 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 486 ]), &(acadoWorkspace.QDy[ 42 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 18 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 630 ]), &(acadoWorkspace.QDy[ 48 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 21 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 792 ]), &(acadoWorkspace.QDy[ 54 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 954 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 24 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 972 ]), &(acadoWorkspace.QDy[ 60 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1152 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 27 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1170 ]), &(acadoWorkspace.QDy[ 66 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1368 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 30 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1386 ]), &(acadoWorkspace.QDy[ 72 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1602 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 33 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1620 ]), &(acadoWorkspace.QDy[ 78 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1854 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 36 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 1872 ]), &(acadoWorkspace.QDy[ 84 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2124 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 39 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2142 ]), &(acadoWorkspace.QDy[ 90 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2412 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 42 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2430 ]), &(acadoWorkspace.QDy[ 96 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2718 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 45 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 2736 ]), &(acadoWorkspace.QDy[ 102 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3042 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 48 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3060 ]), &(acadoWorkspace.QDy[ 108 ]), &(acadoWorkspace.g[ 51 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3384 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 51 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3726 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 51 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3402 ]), &(acadoWorkspace.QDy[ 114 ]), &(acadoWorkspace.g[ 54 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3744 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 54 ]) );
acado_multEQDy( &(acadoWorkspace.E[ 3762 ]), &(acadoWorkspace.QDy[ 120 ]), &(acadoWorkspace.g[ 57 ]) );

acadoWorkspace.g[0] += + acadoWorkspace.H10[0]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[1]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[2]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[3]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[4]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[5]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[1] += + acadoWorkspace.H10[6]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[7]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[8]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[9]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[10]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[11]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[2] += + acadoWorkspace.H10[12]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[13]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[14]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[15]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[16]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[17]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[3] += + acadoWorkspace.H10[18]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[19]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[20]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[21]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[22]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[23]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[4] += + acadoWorkspace.H10[24]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[25]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[26]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[27]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[28]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[29]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[5] += + acadoWorkspace.H10[30]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[31]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[32]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[33]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[34]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[35]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[6] += + acadoWorkspace.H10[36]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[37]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[38]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[39]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[40]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[41]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[7] += + acadoWorkspace.H10[42]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[43]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[44]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[45]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[46]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[47]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[8] += + acadoWorkspace.H10[48]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[49]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[50]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[51]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[52]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[53]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[9] += + acadoWorkspace.H10[54]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[55]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[56]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[57]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[58]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[59]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[10] += + acadoWorkspace.H10[60]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[61]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[62]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[63]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[64]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[65]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[11] += + acadoWorkspace.H10[66]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[67]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[68]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[69]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[70]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[71]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[12] += + acadoWorkspace.H10[72]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[73]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[74]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[75]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[76]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[77]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[13] += + acadoWorkspace.H10[78]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[79]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[80]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[81]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[82]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[83]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[14] += + acadoWorkspace.H10[84]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[85]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[86]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[87]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[88]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[89]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[15] += + acadoWorkspace.H10[90]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[91]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[92]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[93]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[94]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[95]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[16] += + acadoWorkspace.H10[96]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[97]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[98]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[99]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[100]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[101]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[17] += + acadoWorkspace.H10[102]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[103]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[104]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[105]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[106]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[107]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[18] += + acadoWorkspace.H10[108]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[109]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[110]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[111]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[112]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[113]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[19] += + acadoWorkspace.H10[114]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[115]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[116]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[117]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[118]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[119]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[20] += + acadoWorkspace.H10[120]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[121]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[122]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[123]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[124]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[125]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[21] += + acadoWorkspace.H10[126]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[127]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[128]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[129]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[130]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[131]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[22] += + acadoWorkspace.H10[132]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[133]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[134]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[135]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[136]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[137]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[23] += + acadoWorkspace.H10[138]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[139]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[140]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[141]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[142]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[143]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[24] += + acadoWorkspace.H10[144]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[145]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[146]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[147]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[148]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[149]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[25] += + acadoWorkspace.H10[150]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[151]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[152]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[153]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[154]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[155]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[26] += + acadoWorkspace.H10[156]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[157]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[158]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[159]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[160]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[161]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[27] += + acadoWorkspace.H10[162]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[163]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[164]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[165]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[166]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[167]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[28] += + acadoWorkspace.H10[168]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[169]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[170]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[171]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[172]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[173]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[29] += + acadoWorkspace.H10[174]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[175]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[176]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[177]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[178]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[179]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[30] += + acadoWorkspace.H10[180]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[181]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[182]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[183]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[184]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[185]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[31] += + acadoWorkspace.H10[186]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[187]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[188]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[189]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[190]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[191]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[32] += + acadoWorkspace.H10[192]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[193]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[194]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[195]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[196]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[197]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[33] += + acadoWorkspace.H10[198]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[199]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[200]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[201]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[202]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[203]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[34] += + acadoWorkspace.H10[204]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[205]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[206]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[207]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[208]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[209]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[35] += + acadoWorkspace.H10[210]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[211]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[212]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[213]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[214]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[215]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[36] += + acadoWorkspace.H10[216]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[217]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[218]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[219]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[220]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[221]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[37] += + acadoWorkspace.H10[222]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[223]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[224]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[225]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[226]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[227]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[38] += + acadoWorkspace.H10[228]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[229]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[230]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[231]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[232]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[233]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[39] += + acadoWorkspace.H10[234]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[235]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[236]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[237]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[238]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[239]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[40] += + acadoWorkspace.H10[240]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[241]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[242]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[243]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[244]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[245]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[41] += + acadoWorkspace.H10[246]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[247]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[248]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[249]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[250]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[251]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[42] += + acadoWorkspace.H10[252]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[253]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[254]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[255]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[256]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[257]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[43] += + acadoWorkspace.H10[258]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[259]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[260]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[261]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[262]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[263]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[44] += + acadoWorkspace.H10[264]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[265]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[266]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[267]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[268]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[269]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[45] += + acadoWorkspace.H10[270]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[271]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[272]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[273]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[274]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[275]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[46] += + acadoWorkspace.H10[276]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[277]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[278]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[279]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[280]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[281]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[47] += + acadoWorkspace.H10[282]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[283]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[284]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[285]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[286]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[287]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[48] += + acadoWorkspace.H10[288]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[289]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[290]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[291]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[292]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[293]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[49] += + acadoWorkspace.H10[294]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[295]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[296]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[297]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[298]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[299]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[50] += + acadoWorkspace.H10[300]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[301]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[302]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[303]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[304]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[305]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[51] += + acadoWorkspace.H10[306]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[307]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[308]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[309]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[310]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[311]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[52] += + acadoWorkspace.H10[312]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[313]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[314]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[315]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[316]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[317]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[53] += + acadoWorkspace.H10[318]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[319]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[320]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[321]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[322]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[323]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[54] += + acadoWorkspace.H10[324]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[325]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[326]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[327]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[328]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[329]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[55] += + acadoWorkspace.H10[330]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[331]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[332]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[333]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[334]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[335]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[56] += + acadoWorkspace.H10[336]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[337]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[338]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[339]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[340]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[341]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[57] += + acadoWorkspace.H10[342]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[343]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[344]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[345]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[346]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[347]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[58] += + acadoWorkspace.H10[348]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[349]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[350]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[351]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[352]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[353]*acadoWorkspace.Dx0[5];
acadoWorkspace.g[59] += + acadoWorkspace.H10[354]*acadoWorkspace.Dx0[0] + acadoWorkspace.H10[355]*acadoWorkspace.Dx0[1] + acadoWorkspace.H10[356]*acadoWorkspace.Dx0[2] + acadoWorkspace.H10[357]*acadoWorkspace.Dx0[3] + acadoWorkspace.H10[358]*acadoWorkspace.Dx0[4] + acadoWorkspace.H10[359]*acadoWorkspace.Dx0[5];

acadoWorkspace.lb[0] = acadoVariables.lbValues[0] - acadoVariables.u[0];
acadoWorkspace.lb[1] = acadoVariables.lbValues[1] - acadoVariables.u[1];
acadoWorkspace.lb[2] = acadoVariables.lbValues[2] - acadoVariables.u[2];
acadoWorkspace.lb[3] = acadoVariables.lbValues[3] - acadoVariables.u[3];
acadoWorkspace.lb[4] = acadoVariables.lbValues[4] - acadoVariables.u[4];
acadoWorkspace.lb[5] = acadoVariables.lbValues[5] - acadoVariables.u[5];
acadoWorkspace.lb[6] = acadoVariables.lbValues[6] - acadoVariables.u[6];
acadoWorkspace.lb[7] = acadoVariables.lbValues[7] - acadoVariables.u[7];
acadoWorkspace.lb[8] = acadoVariables.lbValues[8] - acadoVariables.u[8];
acadoWorkspace.lb[9] = acadoVariables.lbValues[9] - acadoVariables.u[9];
acadoWorkspace.lb[10] = acadoVariables.lbValues[10] - acadoVariables.u[10];
acadoWorkspace.lb[11] = acadoVariables.lbValues[11] - acadoVariables.u[11];
acadoWorkspace.lb[12] = acadoVariables.lbValues[12] - acadoVariables.u[12];
acadoWorkspace.lb[13] = acadoVariables.lbValues[13] - acadoVariables.u[13];
acadoWorkspace.lb[14] = acadoVariables.lbValues[14] - acadoVariables.u[14];
acadoWorkspace.lb[15] = acadoVariables.lbValues[15] - acadoVariables.u[15];
acadoWorkspace.lb[16] = acadoVariables.lbValues[16] - acadoVariables.u[16];
acadoWorkspace.lb[17] = acadoVariables.lbValues[17] - acadoVariables.u[17];
acadoWorkspace.lb[18] = acadoVariables.lbValues[18] - acadoVariables.u[18];
acadoWorkspace.lb[19] = acadoVariables.lbValues[19] - acadoVariables.u[19];
acadoWorkspace.lb[20] = acadoVariables.lbValues[20] - acadoVariables.u[20];
acadoWorkspace.lb[21] = acadoVariables.lbValues[21] - acadoVariables.u[21];
acadoWorkspace.lb[22] = acadoVariables.lbValues[22] - acadoVariables.u[22];
acadoWorkspace.lb[23] = acadoVariables.lbValues[23] - acadoVariables.u[23];
acadoWorkspace.lb[24] = acadoVariables.lbValues[24] - acadoVariables.u[24];
acadoWorkspace.lb[25] = acadoVariables.lbValues[25] - acadoVariables.u[25];
acadoWorkspace.lb[26] = acadoVariables.lbValues[26] - acadoVariables.u[26];
acadoWorkspace.lb[27] = acadoVariables.lbValues[27] - acadoVariables.u[27];
acadoWorkspace.lb[28] = acadoVariables.lbValues[28] - acadoVariables.u[28];
acadoWorkspace.lb[29] = acadoVariables.lbValues[29] - acadoVariables.u[29];
acadoWorkspace.lb[30] = acadoVariables.lbValues[30] - acadoVariables.u[30];
acadoWorkspace.lb[31] = acadoVariables.lbValues[31] - acadoVariables.u[31];
acadoWorkspace.lb[32] = acadoVariables.lbValues[32] - acadoVariables.u[32];
acadoWorkspace.lb[33] = acadoVariables.lbValues[33] - acadoVariables.u[33];
acadoWorkspace.lb[34] = acadoVariables.lbValues[34] - acadoVariables.u[34];
acadoWorkspace.lb[35] = acadoVariables.lbValues[35] - acadoVariables.u[35];
acadoWorkspace.lb[36] = acadoVariables.lbValues[36] - acadoVariables.u[36];
acadoWorkspace.lb[37] = acadoVariables.lbValues[37] - acadoVariables.u[37];
acadoWorkspace.lb[38] = acadoVariables.lbValues[38] - acadoVariables.u[38];
acadoWorkspace.lb[39] = acadoVariables.lbValues[39] - acadoVariables.u[39];
acadoWorkspace.lb[40] = acadoVariables.lbValues[40] - acadoVariables.u[40];
acadoWorkspace.lb[41] = acadoVariables.lbValues[41] - acadoVariables.u[41];
acadoWorkspace.lb[42] = acadoVariables.lbValues[42] - acadoVariables.u[42];
acadoWorkspace.lb[43] = acadoVariables.lbValues[43] - acadoVariables.u[43];
acadoWorkspace.lb[44] = acadoVariables.lbValues[44] - acadoVariables.u[44];
acadoWorkspace.lb[45] = acadoVariables.lbValues[45] - acadoVariables.u[45];
acadoWorkspace.lb[46] = acadoVariables.lbValues[46] - acadoVariables.u[46];
acadoWorkspace.lb[47] = acadoVariables.lbValues[47] - acadoVariables.u[47];
acadoWorkspace.lb[48] = acadoVariables.lbValues[48] - acadoVariables.u[48];
acadoWorkspace.lb[49] = acadoVariables.lbValues[49] - acadoVariables.u[49];
acadoWorkspace.lb[50] = acadoVariables.lbValues[50] - acadoVariables.u[50];
acadoWorkspace.lb[51] = acadoVariables.lbValues[51] - acadoVariables.u[51];
acadoWorkspace.lb[52] = acadoVariables.lbValues[52] - acadoVariables.u[52];
acadoWorkspace.lb[53] = acadoVariables.lbValues[53] - acadoVariables.u[53];
acadoWorkspace.lb[54] = acadoVariables.lbValues[54] - acadoVariables.u[54];
acadoWorkspace.lb[55] = acadoVariables.lbValues[55] - acadoVariables.u[55];
acadoWorkspace.lb[56] = acadoVariables.lbValues[56] - acadoVariables.u[56];
acadoWorkspace.lb[57] = acadoVariables.lbValues[57] - acadoVariables.u[57];
acadoWorkspace.lb[58] = acadoVariables.lbValues[58] - acadoVariables.u[58];
acadoWorkspace.lb[59] = acadoVariables.lbValues[59] - acadoVariables.u[59];
acadoWorkspace.ub[0] = acadoVariables.ubValues[0] - acadoVariables.u[0];
acadoWorkspace.ub[1] = acadoVariables.ubValues[1] - acadoVariables.u[1];
acadoWorkspace.ub[2] = acadoVariables.ubValues[2] - acadoVariables.u[2];
acadoWorkspace.ub[3] = acadoVariables.ubValues[3] - acadoVariables.u[3];
acadoWorkspace.ub[4] = acadoVariables.ubValues[4] - acadoVariables.u[4];
acadoWorkspace.ub[5] = acadoVariables.ubValues[5] - acadoVariables.u[5];
acadoWorkspace.ub[6] = acadoVariables.ubValues[6] - acadoVariables.u[6];
acadoWorkspace.ub[7] = acadoVariables.ubValues[7] - acadoVariables.u[7];
acadoWorkspace.ub[8] = acadoVariables.ubValues[8] - acadoVariables.u[8];
acadoWorkspace.ub[9] = acadoVariables.ubValues[9] - acadoVariables.u[9];
acadoWorkspace.ub[10] = acadoVariables.ubValues[10] - acadoVariables.u[10];
acadoWorkspace.ub[11] = acadoVariables.ubValues[11] - acadoVariables.u[11];
acadoWorkspace.ub[12] = acadoVariables.ubValues[12] - acadoVariables.u[12];
acadoWorkspace.ub[13] = acadoVariables.ubValues[13] - acadoVariables.u[13];
acadoWorkspace.ub[14] = acadoVariables.ubValues[14] - acadoVariables.u[14];
acadoWorkspace.ub[15] = acadoVariables.ubValues[15] - acadoVariables.u[15];
acadoWorkspace.ub[16] = acadoVariables.ubValues[16] - acadoVariables.u[16];
acadoWorkspace.ub[17] = acadoVariables.ubValues[17] - acadoVariables.u[17];
acadoWorkspace.ub[18] = acadoVariables.ubValues[18] - acadoVariables.u[18];
acadoWorkspace.ub[19] = acadoVariables.ubValues[19] - acadoVariables.u[19];
acadoWorkspace.ub[20] = acadoVariables.ubValues[20] - acadoVariables.u[20];
acadoWorkspace.ub[21] = acadoVariables.ubValues[21] - acadoVariables.u[21];
acadoWorkspace.ub[22] = acadoVariables.ubValues[22] - acadoVariables.u[22];
acadoWorkspace.ub[23] = acadoVariables.ubValues[23] - acadoVariables.u[23];
acadoWorkspace.ub[24] = acadoVariables.ubValues[24] - acadoVariables.u[24];
acadoWorkspace.ub[25] = acadoVariables.ubValues[25] - acadoVariables.u[25];
acadoWorkspace.ub[26] = acadoVariables.ubValues[26] - acadoVariables.u[26];
acadoWorkspace.ub[27] = acadoVariables.ubValues[27] - acadoVariables.u[27];
acadoWorkspace.ub[28] = acadoVariables.ubValues[28] - acadoVariables.u[28];
acadoWorkspace.ub[29] = acadoVariables.ubValues[29] - acadoVariables.u[29];
acadoWorkspace.ub[30] = acadoVariables.ubValues[30] - acadoVariables.u[30];
acadoWorkspace.ub[31] = acadoVariables.ubValues[31] - acadoVariables.u[31];
acadoWorkspace.ub[32] = acadoVariables.ubValues[32] - acadoVariables.u[32];
acadoWorkspace.ub[33] = acadoVariables.ubValues[33] - acadoVariables.u[33];
acadoWorkspace.ub[34] = acadoVariables.ubValues[34] - acadoVariables.u[34];
acadoWorkspace.ub[35] = acadoVariables.ubValues[35] - acadoVariables.u[35];
acadoWorkspace.ub[36] = acadoVariables.ubValues[36] - acadoVariables.u[36];
acadoWorkspace.ub[37] = acadoVariables.ubValues[37] - acadoVariables.u[37];
acadoWorkspace.ub[38] = acadoVariables.ubValues[38] - acadoVariables.u[38];
acadoWorkspace.ub[39] = acadoVariables.ubValues[39] - acadoVariables.u[39];
acadoWorkspace.ub[40] = acadoVariables.ubValues[40] - acadoVariables.u[40];
acadoWorkspace.ub[41] = acadoVariables.ubValues[41] - acadoVariables.u[41];
acadoWorkspace.ub[42] = acadoVariables.ubValues[42] - acadoVariables.u[42];
acadoWorkspace.ub[43] = acadoVariables.ubValues[43] - acadoVariables.u[43];
acadoWorkspace.ub[44] = acadoVariables.ubValues[44] - acadoVariables.u[44];
acadoWorkspace.ub[45] = acadoVariables.ubValues[45] - acadoVariables.u[45];
acadoWorkspace.ub[46] = acadoVariables.ubValues[46] - acadoVariables.u[46];
acadoWorkspace.ub[47] = acadoVariables.ubValues[47] - acadoVariables.u[47];
acadoWorkspace.ub[48] = acadoVariables.ubValues[48] - acadoVariables.u[48];
acadoWorkspace.ub[49] = acadoVariables.ubValues[49] - acadoVariables.u[49];
acadoWorkspace.ub[50] = acadoVariables.ubValues[50] - acadoVariables.u[50];
acadoWorkspace.ub[51] = acadoVariables.ubValues[51] - acadoVariables.u[51];
acadoWorkspace.ub[52] = acadoVariables.ubValues[52] - acadoVariables.u[52];
acadoWorkspace.ub[53] = acadoVariables.ubValues[53] - acadoVariables.u[53];
acadoWorkspace.ub[54] = acadoVariables.ubValues[54] - acadoVariables.u[54];
acadoWorkspace.ub[55] = acadoVariables.ubValues[55] - acadoVariables.u[55];
acadoWorkspace.ub[56] = acadoVariables.ubValues[56] - acadoVariables.u[56];
acadoWorkspace.ub[57] = acadoVariables.ubValues[57] - acadoVariables.u[57];
acadoWorkspace.ub[58] = acadoVariables.ubValues[58] - acadoVariables.u[58];
acadoWorkspace.ub[59] = acadoVariables.ubValues[59] - acadoVariables.u[59];

}

void acado_expand(  )
{
acadoVariables.u[0] += acadoWorkspace.x[0];
acadoVariables.u[1] += acadoWorkspace.x[1];
acadoVariables.u[2] += acadoWorkspace.x[2];
acadoVariables.u[3] += acadoWorkspace.x[3];
acadoVariables.u[4] += acadoWorkspace.x[4];
acadoVariables.u[5] += acadoWorkspace.x[5];
acadoVariables.u[6] += acadoWorkspace.x[6];
acadoVariables.u[7] += acadoWorkspace.x[7];
acadoVariables.u[8] += acadoWorkspace.x[8];
acadoVariables.u[9] += acadoWorkspace.x[9];
acadoVariables.u[10] += acadoWorkspace.x[10];
acadoVariables.u[11] += acadoWorkspace.x[11];
acadoVariables.u[12] += acadoWorkspace.x[12];
acadoVariables.u[13] += acadoWorkspace.x[13];
acadoVariables.u[14] += acadoWorkspace.x[14];
acadoVariables.u[15] += acadoWorkspace.x[15];
acadoVariables.u[16] += acadoWorkspace.x[16];
acadoVariables.u[17] += acadoWorkspace.x[17];
acadoVariables.u[18] += acadoWorkspace.x[18];
acadoVariables.u[19] += acadoWorkspace.x[19];
acadoVariables.u[20] += acadoWorkspace.x[20];
acadoVariables.u[21] += acadoWorkspace.x[21];
acadoVariables.u[22] += acadoWorkspace.x[22];
acadoVariables.u[23] += acadoWorkspace.x[23];
acadoVariables.u[24] += acadoWorkspace.x[24];
acadoVariables.u[25] += acadoWorkspace.x[25];
acadoVariables.u[26] += acadoWorkspace.x[26];
acadoVariables.u[27] += acadoWorkspace.x[27];
acadoVariables.u[28] += acadoWorkspace.x[28];
acadoVariables.u[29] += acadoWorkspace.x[29];
acadoVariables.u[30] += acadoWorkspace.x[30];
acadoVariables.u[31] += acadoWorkspace.x[31];
acadoVariables.u[32] += acadoWorkspace.x[32];
acadoVariables.u[33] += acadoWorkspace.x[33];
acadoVariables.u[34] += acadoWorkspace.x[34];
acadoVariables.u[35] += acadoWorkspace.x[35];
acadoVariables.u[36] += acadoWorkspace.x[36];
acadoVariables.u[37] += acadoWorkspace.x[37];
acadoVariables.u[38] += acadoWorkspace.x[38];
acadoVariables.u[39] += acadoWorkspace.x[39];
acadoVariables.u[40] += acadoWorkspace.x[40];
acadoVariables.u[41] += acadoWorkspace.x[41];
acadoVariables.u[42] += acadoWorkspace.x[42];
acadoVariables.u[43] += acadoWorkspace.x[43];
acadoVariables.u[44] += acadoWorkspace.x[44];
acadoVariables.u[45] += acadoWorkspace.x[45];
acadoVariables.u[46] += acadoWorkspace.x[46];
acadoVariables.u[47] += acadoWorkspace.x[47];
acadoVariables.u[48] += acadoWorkspace.x[48];
acadoVariables.u[49] += acadoWorkspace.x[49];
acadoVariables.u[50] += acadoWorkspace.x[50];
acadoVariables.u[51] += acadoWorkspace.x[51];
acadoVariables.u[52] += acadoWorkspace.x[52];
acadoVariables.u[53] += acadoWorkspace.x[53];
acadoVariables.u[54] += acadoWorkspace.x[54];
acadoVariables.u[55] += acadoWorkspace.x[55];
acadoVariables.u[56] += acadoWorkspace.x[56];
acadoVariables.u[57] += acadoWorkspace.x[57];
acadoVariables.u[58] += acadoWorkspace.x[58];
acadoVariables.u[59] += acadoWorkspace.x[59];

acadoVariables.x[0] += acadoWorkspace.Dx0[0];
acadoVariables.x[1] += acadoWorkspace.Dx0[1];
acadoVariables.x[2] += acadoWorkspace.Dx0[2];
acadoVariables.x[3] += acadoWorkspace.Dx0[3];
acadoVariables.x[4] += acadoWorkspace.Dx0[4];
acadoVariables.x[5] += acadoWorkspace.Dx0[5];

acadoVariables.x[6] += + acadoWorkspace.evGx[0]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[1]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[2]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[3]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[4]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[5]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[0];
acadoVariables.x[7] += + acadoWorkspace.evGx[6]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[7]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[8]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[9]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[10]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[11]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[1];
acadoVariables.x[8] += + acadoWorkspace.evGx[12]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[13]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[14]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[15]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[16]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[17]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[2];
acadoVariables.x[9] += + acadoWorkspace.evGx[18]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[19]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[20]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[21]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[22]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[23]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[3];
acadoVariables.x[10] += + acadoWorkspace.evGx[24]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[25]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[26]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[27]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[28]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[29]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[4];
acadoVariables.x[11] += + acadoWorkspace.evGx[30]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[31]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[32]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[33]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[34]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[35]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[5];
acadoVariables.x[12] += + acadoWorkspace.evGx[36]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[37]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[38]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[39]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[40]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[41]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[6];
acadoVariables.x[13] += + acadoWorkspace.evGx[42]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[43]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[44]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[45]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[46]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[47]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[7];
acadoVariables.x[14] += + acadoWorkspace.evGx[48]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[49]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[50]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[51]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[52]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[53]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[8];
acadoVariables.x[15] += + acadoWorkspace.evGx[54]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[55]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[56]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[57]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[58]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[59]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[9];
acadoVariables.x[16] += + acadoWorkspace.evGx[60]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[61]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[62]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[63]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[64]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[65]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[10];
acadoVariables.x[17] += + acadoWorkspace.evGx[66]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[67]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[68]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[69]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[70]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[71]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[11];
acadoVariables.x[18] += + acadoWorkspace.evGx[72]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[73]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[74]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[75]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[76]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[77]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[12];
acadoVariables.x[19] += + acadoWorkspace.evGx[78]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[79]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[80]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[81]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[82]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[83]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[13];
acadoVariables.x[20] += + acadoWorkspace.evGx[84]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[85]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[86]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[87]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[88]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[89]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[14];
acadoVariables.x[21] += + acadoWorkspace.evGx[90]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[91]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[92]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[93]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[94]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[95]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[15];
acadoVariables.x[22] += + acadoWorkspace.evGx[96]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[97]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[98]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[99]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[100]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[101]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[16];
acadoVariables.x[23] += + acadoWorkspace.evGx[102]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[103]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[104]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[105]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[106]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[107]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[17];
acadoVariables.x[24] += + acadoWorkspace.evGx[108]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[109]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[110]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[111]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[112]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[113]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[18];
acadoVariables.x[25] += + acadoWorkspace.evGx[114]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[115]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[116]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[117]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[118]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[119]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[19];
acadoVariables.x[26] += + acadoWorkspace.evGx[120]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[121]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[122]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[123]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[124]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[125]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[20];
acadoVariables.x[27] += + acadoWorkspace.evGx[126]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[127]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[128]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[129]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[130]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[131]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[21];
acadoVariables.x[28] += + acadoWorkspace.evGx[132]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[133]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[134]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[135]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[136]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[137]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[22];
acadoVariables.x[29] += + acadoWorkspace.evGx[138]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[139]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[140]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[141]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[142]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[143]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[23];
acadoVariables.x[30] += + acadoWorkspace.evGx[144]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[145]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[146]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[147]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[148]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[149]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[24];
acadoVariables.x[31] += + acadoWorkspace.evGx[150]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[151]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[152]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[153]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[154]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[155]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[25];
acadoVariables.x[32] += + acadoWorkspace.evGx[156]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[157]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[158]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[159]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[160]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[161]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[26];
acadoVariables.x[33] += + acadoWorkspace.evGx[162]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[163]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[164]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[165]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[166]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[167]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[27];
acadoVariables.x[34] += + acadoWorkspace.evGx[168]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[169]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[170]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[171]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[172]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[173]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[28];
acadoVariables.x[35] += + acadoWorkspace.evGx[174]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[175]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[176]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[177]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[178]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[179]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[29];
acadoVariables.x[36] += + acadoWorkspace.evGx[180]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[181]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[182]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[183]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[184]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[185]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[30];
acadoVariables.x[37] += + acadoWorkspace.evGx[186]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[187]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[188]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[189]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[190]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[191]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[31];
acadoVariables.x[38] += + acadoWorkspace.evGx[192]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[193]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[194]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[195]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[196]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[197]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[32];
acadoVariables.x[39] += + acadoWorkspace.evGx[198]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[199]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[200]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[201]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[202]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[203]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[33];
acadoVariables.x[40] += + acadoWorkspace.evGx[204]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[205]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[206]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[207]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[208]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[209]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[34];
acadoVariables.x[41] += + acadoWorkspace.evGx[210]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[211]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[212]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[213]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[214]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[215]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[35];
acadoVariables.x[42] += + acadoWorkspace.evGx[216]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[217]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[218]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[219]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[220]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[221]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[36];
acadoVariables.x[43] += + acadoWorkspace.evGx[222]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[223]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[224]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[225]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[226]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[227]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[37];
acadoVariables.x[44] += + acadoWorkspace.evGx[228]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[229]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[230]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[231]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[232]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[233]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[38];
acadoVariables.x[45] += + acadoWorkspace.evGx[234]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[235]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[236]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[237]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[238]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[239]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[39];
acadoVariables.x[46] += + acadoWorkspace.evGx[240]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[241]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[242]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[243]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[244]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[245]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[40];
acadoVariables.x[47] += + acadoWorkspace.evGx[246]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[247]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[248]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[249]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[250]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[251]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[41];
acadoVariables.x[48] += + acadoWorkspace.evGx[252]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[253]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[254]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[255]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[256]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[257]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[42];
acadoVariables.x[49] += + acadoWorkspace.evGx[258]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[259]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[260]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[261]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[262]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[263]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[43];
acadoVariables.x[50] += + acadoWorkspace.evGx[264]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[265]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[266]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[267]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[268]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[269]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[44];
acadoVariables.x[51] += + acadoWorkspace.evGx[270]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[271]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[272]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[273]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[274]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[275]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[45];
acadoVariables.x[52] += + acadoWorkspace.evGx[276]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[277]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[278]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[279]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[280]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[281]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[46];
acadoVariables.x[53] += + acadoWorkspace.evGx[282]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[283]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[284]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[285]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[286]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[287]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[47];
acadoVariables.x[54] += + acadoWorkspace.evGx[288]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[289]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[290]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[291]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[292]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[293]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[48];
acadoVariables.x[55] += + acadoWorkspace.evGx[294]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[295]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[296]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[297]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[298]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[299]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[49];
acadoVariables.x[56] += + acadoWorkspace.evGx[300]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[301]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[302]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[303]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[304]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[305]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[50];
acadoVariables.x[57] += + acadoWorkspace.evGx[306]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[307]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[308]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[309]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[310]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[311]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[51];
acadoVariables.x[58] += + acadoWorkspace.evGx[312]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[313]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[314]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[315]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[316]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[317]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[52];
acadoVariables.x[59] += + acadoWorkspace.evGx[318]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[319]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[320]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[321]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[322]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[323]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[53];
acadoVariables.x[60] += + acadoWorkspace.evGx[324]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[325]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[326]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[327]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[328]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[329]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[54];
acadoVariables.x[61] += + acadoWorkspace.evGx[330]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[331]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[332]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[333]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[334]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[335]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[55];
acadoVariables.x[62] += + acadoWorkspace.evGx[336]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[337]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[338]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[339]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[340]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[341]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[56];
acadoVariables.x[63] += + acadoWorkspace.evGx[342]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[343]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[344]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[345]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[346]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[347]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[57];
acadoVariables.x[64] += + acadoWorkspace.evGx[348]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[349]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[350]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[351]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[352]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[353]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[58];
acadoVariables.x[65] += + acadoWorkspace.evGx[354]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[355]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[356]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[357]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[358]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[359]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[59];
acadoVariables.x[66] += + acadoWorkspace.evGx[360]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[361]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[362]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[363]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[364]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[365]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[60];
acadoVariables.x[67] += + acadoWorkspace.evGx[366]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[367]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[368]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[369]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[370]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[371]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[61];
acadoVariables.x[68] += + acadoWorkspace.evGx[372]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[373]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[374]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[375]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[376]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[377]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[62];
acadoVariables.x[69] += + acadoWorkspace.evGx[378]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[379]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[380]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[381]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[382]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[383]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[63];
acadoVariables.x[70] += + acadoWorkspace.evGx[384]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[385]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[386]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[387]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[388]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[389]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[64];
acadoVariables.x[71] += + acadoWorkspace.evGx[390]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[391]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[392]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[393]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[394]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[395]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[65];
acadoVariables.x[72] += + acadoWorkspace.evGx[396]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[397]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[398]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[399]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[400]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[401]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[66];
acadoVariables.x[73] += + acadoWorkspace.evGx[402]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[403]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[404]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[405]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[406]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[407]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[67];
acadoVariables.x[74] += + acadoWorkspace.evGx[408]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[409]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[410]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[411]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[412]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[413]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[68];
acadoVariables.x[75] += + acadoWorkspace.evGx[414]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[415]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[416]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[417]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[418]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[419]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[69];
acadoVariables.x[76] += + acadoWorkspace.evGx[420]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[421]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[422]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[423]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[424]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[425]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[70];
acadoVariables.x[77] += + acadoWorkspace.evGx[426]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[427]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[428]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[429]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[430]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[431]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[71];
acadoVariables.x[78] += + acadoWorkspace.evGx[432]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[433]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[434]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[435]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[436]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[437]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[72];
acadoVariables.x[79] += + acadoWorkspace.evGx[438]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[439]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[440]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[441]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[442]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[443]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[73];
acadoVariables.x[80] += + acadoWorkspace.evGx[444]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[445]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[446]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[447]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[448]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[449]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[74];
acadoVariables.x[81] += + acadoWorkspace.evGx[450]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[451]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[452]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[453]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[454]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[455]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[75];
acadoVariables.x[82] += + acadoWorkspace.evGx[456]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[457]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[458]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[459]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[460]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[461]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[76];
acadoVariables.x[83] += + acadoWorkspace.evGx[462]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[463]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[464]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[465]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[466]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[467]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[77];
acadoVariables.x[84] += + acadoWorkspace.evGx[468]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[469]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[470]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[471]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[472]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[473]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[78];
acadoVariables.x[85] += + acadoWorkspace.evGx[474]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[475]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[476]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[477]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[478]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[479]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[79];
acadoVariables.x[86] += + acadoWorkspace.evGx[480]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[481]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[482]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[483]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[484]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[485]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[80];
acadoVariables.x[87] += + acadoWorkspace.evGx[486]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[487]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[488]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[489]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[490]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[491]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[81];
acadoVariables.x[88] += + acadoWorkspace.evGx[492]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[493]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[494]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[495]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[496]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[497]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[82];
acadoVariables.x[89] += + acadoWorkspace.evGx[498]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[499]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[500]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[501]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[502]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[503]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[83];
acadoVariables.x[90] += + acadoWorkspace.evGx[504]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[505]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[506]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[507]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[508]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[509]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[84];
acadoVariables.x[91] += + acadoWorkspace.evGx[510]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[511]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[512]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[513]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[514]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[515]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[85];
acadoVariables.x[92] += + acadoWorkspace.evGx[516]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[517]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[518]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[519]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[520]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[521]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[86];
acadoVariables.x[93] += + acadoWorkspace.evGx[522]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[523]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[524]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[525]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[526]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[527]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[87];
acadoVariables.x[94] += + acadoWorkspace.evGx[528]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[529]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[530]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[531]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[532]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[533]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[88];
acadoVariables.x[95] += + acadoWorkspace.evGx[534]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[535]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[536]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[537]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[538]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[539]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[89];
acadoVariables.x[96] += + acadoWorkspace.evGx[540]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[541]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[542]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[543]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[544]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[545]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[90];
acadoVariables.x[97] += + acadoWorkspace.evGx[546]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[547]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[548]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[549]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[550]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[551]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[91];
acadoVariables.x[98] += + acadoWorkspace.evGx[552]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[553]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[554]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[555]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[556]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[557]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[92];
acadoVariables.x[99] += + acadoWorkspace.evGx[558]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[559]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[560]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[561]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[562]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[563]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[93];
acadoVariables.x[100] += + acadoWorkspace.evGx[564]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[565]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[566]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[567]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[568]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[569]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[94];
acadoVariables.x[101] += + acadoWorkspace.evGx[570]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[571]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[572]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[573]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[574]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[575]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[95];
acadoVariables.x[102] += + acadoWorkspace.evGx[576]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[577]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[578]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[579]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[580]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[581]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[96];
acadoVariables.x[103] += + acadoWorkspace.evGx[582]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[583]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[584]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[585]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[586]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[587]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[97];
acadoVariables.x[104] += + acadoWorkspace.evGx[588]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[589]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[590]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[591]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[592]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[593]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[98];
acadoVariables.x[105] += + acadoWorkspace.evGx[594]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[595]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[596]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[597]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[598]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[599]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[99];
acadoVariables.x[106] += + acadoWorkspace.evGx[600]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[601]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[602]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[603]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[604]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[605]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[100];
acadoVariables.x[107] += + acadoWorkspace.evGx[606]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[607]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[608]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[609]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[610]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[611]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[101];
acadoVariables.x[108] += + acadoWorkspace.evGx[612]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[613]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[614]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[615]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[616]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[617]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[102];
acadoVariables.x[109] += + acadoWorkspace.evGx[618]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[619]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[620]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[621]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[622]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[623]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[103];
acadoVariables.x[110] += + acadoWorkspace.evGx[624]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[625]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[626]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[627]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[628]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[629]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[104];
acadoVariables.x[111] += + acadoWorkspace.evGx[630]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[631]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[632]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[633]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[634]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[635]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[105];
acadoVariables.x[112] += + acadoWorkspace.evGx[636]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[637]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[638]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[639]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[640]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[641]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[106];
acadoVariables.x[113] += + acadoWorkspace.evGx[642]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[643]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[644]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[645]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[646]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[647]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[107];
acadoVariables.x[114] += + acadoWorkspace.evGx[648]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[649]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[650]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[651]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[652]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[653]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[108];
acadoVariables.x[115] += + acadoWorkspace.evGx[654]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[655]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[656]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[657]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[658]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[659]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[109];
acadoVariables.x[116] += + acadoWorkspace.evGx[660]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[661]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[662]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[663]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[664]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[665]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[110];
acadoVariables.x[117] += + acadoWorkspace.evGx[666]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[667]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[668]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[669]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[670]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[671]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[111];
acadoVariables.x[118] += + acadoWorkspace.evGx[672]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[673]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[674]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[675]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[676]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[677]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[112];
acadoVariables.x[119] += + acadoWorkspace.evGx[678]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[679]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[680]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[681]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[682]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[683]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[113];
acadoVariables.x[120] += + acadoWorkspace.evGx[684]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[685]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[686]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[687]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[688]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[689]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[114];
acadoVariables.x[121] += + acadoWorkspace.evGx[690]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[691]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[692]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[693]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[694]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[695]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[115];
acadoVariables.x[122] += + acadoWorkspace.evGx[696]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[697]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[698]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[699]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[700]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[701]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[116];
acadoVariables.x[123] += + acadoWorkspace.evGx[702]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[703]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[704]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[705]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[706]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[707]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[117];
acadoVariables.x[124] += + acadoWorkspace.evGx[708]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[709]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[710]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[711]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[712]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[713]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[118];
acadoVariables.x[125] += + acadoWorkspace.evGx[714]*acadoWorkspace.Dx0[0] + acadoWorkspace.evGx[715]*acadoWorkspace.Dx0[1] + acadoWorkspace.evGx[716]*acadoWorkspace.Dx0[2] + acadoWorkspace.evGx[717]*acadoWorkspace.Dx0[3] + acadoWorkspace.evGx[718]*acadoWorkspace.Dx0[4] + acadoWorkspace.evGx[719]*acadoWorkspace.Dx0[5] + acadoWorkspace.d[119];

acado_multEDu( acadoWorkspace.E, acadoWorkspace.x, &(acadoVariables.x[ 6 ]) );
acado_multEDu( &(acadoWorkspace.E[ 18 ]), acadoWorkspace.x, &(acadoVariables.x[ 12 ]) );
acado_multEDu( &(acadoWorkspace.E[ 36 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 12 ]) );
acado_multEDu( &(acadoWorkspace.E[ 54 ]), acadoWorkspace.x, &(acadoVariables.x[ 18 ]) );
acado_multEDu( &(acadoWorkspace.E[ 72 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 18 ]) );
acado_multEDu( &(acadoWorkspace.E[ 90 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 18 ]) );
acado_multEDu( &(acadoWorkspace.E[ 108 ]), acadoWorkspace.x, &(acadoVariables.x[ 24 ]) );
acado_multEDu( &(acadoWorkspace.E[ 126 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 24 ]) );
acado_multEDu( &(acadoWorkspace.E[ 144 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 24 ]) );
acado_multEDu( &(acadoWorkspace.E[ 162 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 24 ]) );
acado_multEDu( &(acadoWorkspace.E[ 180 ]), acadoWorkspace.x, &(acadoVariables.x[ 30 ]) );
acado_multEDu( &(acadoWorkspace.E[ 198 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 30 ]) );
acado_multEDu( &(acadoWorkspace.E[ 216 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 30 ]) );
acado_multEDu( &(acadoWorkspace.E[ 234 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 30 ]) );
acado_multEDu( &(acadoWorkspace.E[ 252 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 30 ]) );
acado_multEDu( &(acadoWorkspace.E[ 270 ]), acadoWorkspace.x, &(acadoVariables.x[ 36 ]) );
acado_multEDu( &(acadoWorkspace.E[ 288 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 36 ]) );
acado_multEDu( &(acadoWorkspace.E[ 306 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 36 ]) );
acado_multEDu( &(acadoWorkspace.E[ 324 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 36 ]) );
acado_multEDu( &(acadoWorkspace.E[ 342 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 36 ]) );
acado_multEDu( &(acadoWorkspace.E[ 360 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 36 ]) );
acado_multEDu( &(acadoWorkspace.E[ 378 ]), acadoWorkspace.x, &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 396 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 414 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 432 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 450 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 468 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 486 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 42 ]) );
acado_multEDu( &(acadoWorkspace.E[ 504 ]), acadoWorkspace.x, &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 522 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 540 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 558 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 576 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 594 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 612 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 630 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 48 ]) );
acado_multEDu( &(acadoWorkspace.E[ 648 ]), acadoWorkspace.x, &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 666 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 684 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 702 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 720 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 738 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 756 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 774 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 792 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 54 ]) );
acado_multEDu( &(acadoWorkspace.E[ 810 ]), acadoWorkspace.x, &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 828 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 846 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 864 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 882 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 900 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 918 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 936 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 954 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 972 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 60 ]) );
acado_multEDu( &(acadoWorkspace.E[ 990 ]), acadoWorkspace.x, &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1008 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1026 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1044 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1062 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1080 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1098 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1116 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1134 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1152 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1170 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 66 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1188 ]), acadoWorkspace.x, &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1206 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1224 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1242 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1260 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1278 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1296 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1314 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1332 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1350 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1368 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1386 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 72 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1404 ]), acadoWorkspace.x, &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1422 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1440 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1458 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1476 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1494 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1512 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1530 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1548 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1566 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1584 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1602 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1620 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 78 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1638 ]), acadoWorkspace.x, &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1656 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1674 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1692 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1710 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1728 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1746 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1764 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1782 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1800 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1818 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1836 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1854 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1872 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 84 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1890 ]), acadoWorkspace.x, &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1908 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1926 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1944 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1962 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1980 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 1998 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2016 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2034 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2052 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2070 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2088 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2106 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2124 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2142 ]), &(acadoWorkspace.x[ 42 ]), &(acadoVariables.x[ 90 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2160 ]), acadoWorkspace.x, &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2178 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2196 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2214 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2232 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2250 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2268 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2286 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2304 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2322 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2340 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2358 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2376 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2394 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2412 ]), &(acadoWorkspace.x[ 42 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2430 ]), &(acadoWorkspace.x[ 45 ]), &(acadoVariables.x[ 96 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2448 ]), acadoWorkspace.x, &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2466 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2484 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2502 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2520 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2538 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2556 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2574 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2592 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2610 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2628 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2646 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2664 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2682 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2700 ]), &(acadoWorkspace.x[ 42 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2718 ]), &(acadoWorkspace.x[ 45 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2736 ]), &(acadoWorkspace.x[ 48 ]), &(acadoVariables.x[ 102 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2754 ]), acadoWorkspace.x, &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2772 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2790 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2808 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2826 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2844 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2862 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2880 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2898 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2916 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2934 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2952 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2970 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 2988 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3006 ]), &(acadoWorkspace.x[ 42 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3024 ]), &(acadoWorkspace.x[ 45 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3042 ]), &(acadoWorkspace.x[ 48 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3060 ]), &(acadoWorkspace.x[ 51 ]), &(acadoVariables.x[ 108 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3078 ]), acadoWorkspace.x, &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3096 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3114 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3132 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3150 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3168 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3186 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3204 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3222 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3240 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3258 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3276 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3294 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3312 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3330 ]), &(acadoWorkspace.x[ 42 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3348 ]), &(acadoWorkspace.x[ 45 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3366 ]), &(acadoWorkspace.x[ 48 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3384 ]), &(acadoWorkspace.x[ 51 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3402 ]), &(acadoWorkspace.x[ 54 ]), &(acadoVariables.x[ 114 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3420 ]), acadoWorkspace.x, &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3438 ]), &(acadoWorkspace.x[ 3 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3456 ]), &(acadoWorkspace.x[ 6 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3474 ]), &(acadoWorkspace.x[ 9 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3492 ]), &(acadoWorkspace.x[ 12 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3510 ]), &(acadoWorkspace.x[ 15 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3528 ]), &(acadoWorkspace.x[ 18 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3546 ]), &(acadoWorkspace.x[ 21 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3564 ]), &(acadoWorkspace.x[ 24 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3582 ]), &(acadoWorkspace.x[ 27 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3600 ]), &(acadoWorkspace.x[ 30 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3618 ]), &(acadoWorkspace.x[ 33 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3636 ]), &(acadoWorkspace.x[ 36 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3654 ]), &(acadoWorkspace.x[ 39 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3672 ]), &(acadoWorkspace.x[ 42 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3690 ]), &(acadoWorkspace.x[ 45 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3708 ]), &(acadoWorkspace.x[ 48 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3726 ]), &(acadoWorkspace.x[ 51 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3744 ]), &(acadoWorkspace.x[ 54 ]), &(acadoVariables.x[ 120 ]) );
acado_multEDu( &(acadoWorkspace.E[ 3762 ]), &(acadoWorkspace.x[ 57 ]), &(acadoVariables.x[ 120 ]) );
}

int acado_preparationStep(  )
{
int ret;

ret = acado_modelSimulation();
acado_evaluateObjective(  );
acado_condensePrep(  );
return ret;
}

int acado_feedbackStep(  )
{
int tmp;

acado_condenseFdb(  );

tmp = acado_solve( );

acado_expand(  );
return tmp;
}

int acado_initializeSolver(  )
{
int ret;

/* This is a function which must be called once before any other function call! */


ret = 0;

memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
acadoVariables.lbValues[0] = 2.0000000000000000e+00;
acadoVariables.lbValues[1] = -2.0000000000000000e+00;
acadoVariables.lbValues[2] = -2.0000000000000000e+00;
acadoVariables.lbValues[3] = 2.0000000000000000e+00;
acadoVariables.lbValues[4] = -2.0000000000000000e+00;
acadoVariables.lbValues[5] = -2.0000000000000000e+00;
acadoVariables.lbValues[6] = 2.0000000000000000e+00;
acadoVariables.lbValues[7] = -2.0000000000000000e+00;
acadoVariables.lbValues[8] = -2.0000000000000000e+00;
acadoVariables.lbValues[9] = 2.0000000000000000e+00;
acadoVariables.lbValues[10] = -2.0000000000000000e+00;
acadoVariables.lbValues[11] = -2.0000000000000000e+00;
acadoVariables.lbValues[12] = 2.0000000000000000e+00;
acadoVariables.lbValues[13] = -2.0000000000000000e+00;
acadoVariables.lbValues[14] = -2.0000000000000000e+00;
acadoVariables.lbValues[15] = 2.0000000000000000e+00;
acadoVariables.lbValues[16] = -2.0000000000000000e+00;
acadoVariables.lbValues[17] = -2.0000000000000000e+00;
acadoVariables.lbValues[18] = 2.0000000000000000e+00;
acadoVariables.lbValues[19] = -2.0000000000000000e+00;
acadoVariables.lbValues[20] = -2.0000000000000000e+00;
acadoVariables.lbValues[21] = 2.0000000000000000e+00;
acadoVariables.lbValues[22] = -2.0000000000000000e+00;
acadoVariables.lbValues[23] = -2.0000000000000000e+00;
acadoVariables.lbValues[24] = 2.0000000000000000e+00;
acadoVariables.lbValues[25] = -2.0000000000000000e+00;
acadoVariables.lbValues[26] = -2.0000000000000000e+00;
acadoVariables.lbValues[27] = 2.0000000000000000e+00;
acadoVariables.lbValues[28] = -2.0000000000000000e+00;
acadoVariables.lbValues[29] = -2.0000000000000000e+00;
acadoVariables.lbValues[30] = 2.0000000000000000e+00;
acadoVariables.lbValues[31] = -2.0000000000000000e+00;
acadoVariables.lbValues[32] = -2.0000000000000000e+00;
acadoVariables.lbValues[33] = 2.0000000000000000e+00;
acadoVariables.lbValues[34] = -2.0000000000000000e+00;
acadoVariables.lbValues[35] = -2.0000000000000000e+00;
acadoVariables.lbValues[36] = 2.0000000000000000e+00;
acadoVariables.lbValues[37] = -2.0000000000000000e+00;
acadoVariables.lbValues[38] = -2.0000000000000000e+00;
acadoVariables.lbValues[39] = 2.0000000000000000e+00;
acadoVariables.lbValues[40] = -2.0000000000000000e+00;
acadoVariables.lbValues[41] = -2.0000000000000000e+00;
acadoVariables.lbValues[42] = 2.0000000000000000e+00;
acadoVariables.lbValues[43] = -2.0000000000000000e+00;
acadoVariables.lbValues[44] = -2.0000000000000000e+00;
acadoVariables.lbValues[45] = 2.0000000000000000e+00;
acadoVariables.lbValues[46] = -2.0000000000000000e+00;
acadoVariables.lbValues[47] = -2.0000000000000000e+00;
acadoVariables.lbValues[48] = 2.0000000000000000e+00;
acadoVariables.lbValues[49] = -2.0000000000000000e+00;
acadoVariables.lbValues[50] = -2.0000000000000000e+00;
acadoVariables.lbValues[51] = 2.0000000000000000e+00;
acadoVariables.lbValues[52] = -2.0000000000000000e+00;
acadoVariables.lbValues[53] = -2.0000000000000000e+00;
acadoVariables.lbValues[54] = 2.0000000000000000e+00;
acadoVariables.lbValues[55] = -2.0000000000000000e+00;
acadoVariables.lbValues[56] = -2.0000000000000000e+00;
acadoVariables.lbValues[57] = 2.0000000000000000e+00;
acadoVariables.lbValues[58] = -2.0000000000000000e+00;
acadoVariables.lbValues[59] = -2.0000000000000000e+00;
acadoVariables.ubValues[0] = 2.0000000000000000e+01;
acadoVariables.ubValues[1] = 2.0000000000000000e+00;
acadoVariables.ubValues[2] = 2.0000000000000000e+00;
acadoVariables.ubValues[3] = 2.0000000000000000e+01;
acadoVariables.ubValues[4] = 2.0000000000000000e+00;
acadoVariables.ubValues[5] = 2.0000000000000000e+00;
acadoVariables.ubValues[6] = 2.0000000000000000e+01;
acadoVariables.ubValues[7] = 2.0000000000000000e+00;
acadoVariables.ubValues[8] = 2.0000000000000000e+00;
acadoVariables.ubValues[9] = 2.0000000000000000e+01;
acadoVariables.ubValues[10] = 2.0000000000000000e+00;
acadoVariables.ubValues[11] = 2.0000000000000000e+00;
acadoVariables.ubValues[12] = 2.0000000000000000e+01;
acadoVariables.ubValues[13] = 2.0000000000000000e+00;
acadoVariables.ubValues[14] = 2.0000000000000000e+00;
acadoVariables.ubValues[15] = 2.0000000000000000e+01;
acadoVariables.ubValues[16] = 2.0000000000000000e+00;
acadoVariables.ubValues[17] = 2.0000000000000000e+00;
acadoVariables.ubValues[18] = 2.0000000000000000e+01;
acadoVariables.ubValues[19] = 2.0000000000000000e+00;
acadoVariables.ubValues[20] = 2.0000000000000000e+00;
acadoVariables.ubValues[21] = 2.0000000000000000e+01;
acadoVariables.ubValues[22] = 2.0000000000000000e+00;
acadoVariables.ubValues[23] = 2.0000000000000000e+00;
acadoVariables.ubValues[24] = 2.0000000000000000e+01;
acadoVariables.ubValues[25] = 2.0000000000000000e+00;
acadoVariables.ubValues[26] = 2.0000000000000000e+00;
acadoVariables.ubValues[27] = 2.0000000000000000e+01;
acadoVariables.ubValues[28] = 2.0000000000000000e+00;
acadoVariables.ubValues[29] = 2.0000000000000000e+00;
acadoVariables.ubValues[30] = 2.0000000000000000e+01;
acadoVariables.ubValues[31] = 2.0000000000000000e+00;
acadoVariables.ubValues[32] = 2.0000000000000000e+00;
acadoVariables.ubValues[33] = 2.0000000000000000e+01;
acadoVariables.ubValues[34] = 2.0000000000000000e+00;
acadoVariables.ubValues[35] = 2.0000000000000000e+00;
acadoVariables.ubValues[36] = 2.0000000000000000e+01;
acadoVariables.ubValues[37] = 2.0000000000000000e+00;
acadoVariables.ubValues[38] = 2.0000000000000000e+00;
acadoVariables.ubValues[39] = 2.0000000000000000e+01;
acadoVariables.ubValues[40] = 2.0000000000000000e+00;
acadoVariables.ubValues[41] = 2.0000000000000000e+00;
acadoVariables.ubValues[42] = 2.0000000000000000e+01;
acadoVariables.ubValues[43] = 2.0000000000000000e+00;
acadoVariables.ubValues[44] = 2.0000000000000000e+00;
acadoVariables.ubValues[45] = 2.0000000000000000e+01;
acadoVariables.ubValues[46] = 2.0000000000000000e+00;
acadoVariables.ubValues[47] = 2.0000000000000000e+00;
acadoVariables.ubValues[48] = 2.0000000000000000e+01;
acadoVariables.ubValues[49] = 2.0000000000000000e+00;
acadoVariables.ubValues[50] = 2.0000000000000000e+00;
acadoVariables.ubValues[51] = 2.0000000000000000e+01;
acadoVariables.ubValues[52] = 2.0000000000000000e+00;
acadoVariables.ubValues[53] = 2.0000000000000000e+00;
acadoVariables.ubValues[54] = 2.0000000000000000e+01;
acadoVariables.ubValues[55] = 2.0000000000000000e+00;
acadoVariables.ubValues[56] = 2.0000000000000000e+00;
acadoVariables.ubValues[57] = 2.0000000000000000e+01;
acadoVariables.ubValues[58] = 2.0000000000000000e+00;
acadoVariables.ubValues[59] = 2.0000000000000000e+00;
return ret;
}

void acado_initializeNodesByForwardSimulation(  )
{
int index;
for (index = 0; index < 20; ++index)
{
state[0] = acadoVariables.x[index * 6];
state[1] = acadoVariables.x[index * 6 + 1];
state[2] = acadoVariables.x[index * 6 + 2];
state[3] = acadoVariables.x[index * 6 + 3];
state[4] = acadoVariables.x[index * 6 + 4];
state[5] = acadoVariables.x[index * 6 + 5];
state[60] = acadoVariables.u[index * 3];
state[61] = acadoVariables.u[index * 3 + 1];
state[62] = acadoVariables.u[index * 3 + 2];

acado_integrate(state, index == 0);

acadoVariables.x[index * 6 + 6] = state[0];
acadoVariables.x[index * 6 + 7] = state[1];
acadoVariables.x[index * 6 + 8] = state[2];
acadoVariables.x[index * 6 + 9] = state[3];
acadoVariables.x[index * 6 + 10] = state[4];
acadoVariables.x[index * 6 + 11] = state[5];
}
}

void acado_shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd )
{
int index;
for (index = 0; index < 20; ++index)
{
acadoVariables.x[index * 6] = acadoVariables.x[index * 6 + 6];
acadoVariables.x[index * 6 + 1] = acadoVariables.x[index * 6 + 7];
acadoVariables.x[index * 6 + 2] = acadoVariables.x[index * 6 + 8];
acadoVariables.x[index * 6 + 3] = acadoVariables.x[index * 6 + 9];
acadoVariables.x[index * 6 + 4] = acadoVariables.x[index * 6 + 10];
acadoVariables.x[index * 6 + 5] = acadoVariables.x[index * 6 + 11];
}

if (strategy == 1 && xEnd != 0)
{
acadoVariables.x[120] = xEnd[0];
acadoVariables.x[121] = xEnd[1];
acadoVariables.x[122] = xEnd[2];
acadoVariables.x[123] = xEnd[3];
acadoVariables.x[124] = xEnd[4];
acadoVariables.x[125] = xEnd[5];
}
else if (strategy == 2) 
{
state[0] = acadoVariables.x[120];
state[1] = acadoVariables.x[121];
state[2] = acadoVariables.x[122];
state[3] = acadoVariables.x[123];
state[4] = acadoVariables.x[124];
state[5] = acadoVariables.x[125];
if (uEnd != 0)
{
state[60] = uEnd[0];
state[61] = uEnd[1];
state[62] = uEnd[2];
}
else
{
state[60] = acadoVariables.u[57];
state[61] = acadoVariables.u[58];
state[62] = acadoVariables.u[59];
}

acado_integrate(state, 1);

acadoVariables.x[120] = state[0];
acadoVariables.x[121] = state[1];
acadoVariables.x[122] = state[2];
acadoVariables.x[123] = state[3];
acadoVariables.x[124] = state[4];
acadoVariables.x[125] = state[5];
}
}

void acado_shiftControls( real_t* const uEnd )
{
int index;
for (index = 0; index < 19; ++index)
{
acadoVariables.u[index * 3] = acadoVariables.u[index * 3 + 3];
acadoVariables.u[index * 3 + 1] = acadoVariables.u[index * 3 + 4];
acadoVariables.u[index * 3 + 2] = acadoVariables.u[index * 3 + 5];
}

if (uEnd != 0)
{
acadoVariables.u[57] = uEnd[0];
acadoVariables.u[58] = uEnd[1];
acadoVariables.u[59] = uEnd[2];
}
}

real_t acado_getKKT(  )
{
real_t kkt;

int index;
real_t prd;

kkt = + acadoWorkspace.g[0]*acadoWorkspace.x[0] + acadoWorkspace.g[1]*acadoWorkspace.x[1] + acadoWorkspace.g[2]*acadoWorkspace.x[2] + acadoWorkspace.g[3]*acadoWorkspace.x[3] + acadoWorkspace.g[4]*acadoWorkspace.x[4] + acadoWorkspace.g[5]*acadoWorkspace.x[5] + acadoWorkspace.g[6]*acadoWorkspace.x[6] + acadoWorkspace.g[7]*acadoWorkspace.x[7] + acadoWorkspace.g[8]*acadoWorkspace.x[8] + acadoWorkspace.g[9]*acadoWorkspace.x[9] + acadoWorkspace.g[10]*acadoWorkspace.x[10] + acadoWorkspace.g[11]*acadoWorkspace.x[11] + acadoWorkspace.g[12]*acadoWorkspace.x[12] + acadoWorkspace.g[13]*acadoWorkspace.x[13] + acadoWorkspace.g[14]*acadoWorkspace.x[14] + acadoWorkspace.g[15]*acadoWorkspace.x[15] + acadoWorkspace.g[16]*acadoWorkspace.x[16] + acadoWorkspace.g[17]*acadoWorkspace.x[17] + acadoWorkspace.g[18]*acadoWorkspace.x[18] + acadoWorkspace.g[19]*acadoWorkspace.x[19] + acadoWorkspace.g[20]*acadoWorkspace.x[20] + acadoWorkspace.g[21]*acadoWorkspace.x[21] + acadoWorkspace.g[22]*acadoWorkspace.x[22] + acadoWorkspace.g[23]*acadoWorkspace.x[23] + acadoWorkspace.g[24]*acadoWorkspace.x[24] + acadoWorkspace.g[25]*acadoWorkspace.x[25] + acadoWorkspace.g[26]*acadoWorkspace.x[26] + acadoWorkspace.g[27]*acadoWorkspace.x[27] + acadoWorkspace.g[28]*acadoWorkspace.x[28] + acadoWorkspace.g[29]*acadoWorkspace.x[29] + acadoWorkspace.g[30]*acadoWorkspace.x[30] + acadoWorkspace.g[31]*acadoWorkspace.x[31] + acadoWorkspace.g[32]*acadoWorkspace.x[32] + acadoWorkspace.g[33]*acadoWorkspace.x[33] + acadoWorkspace.g[34]*acadoWorkspace.x[34] + acadoWorkspace.g[35]*acadoWorkspace.x[35] + acadoWorkspace.g[36]*acadoWorkspace.x[36] + acadoWorkspace.g[37]*acadoWorkspace.x[37] + acadoWorkspace.g[38]*acadoWorkspace.x[38] + acadoWorkspace.g[39]*acadoWorkspace.x[39] + acadoWorkspace.g[40]*acadoWorkspace.x[40] + acadoWorkspace.g[41]*acadoWorkspace.x[41] + acadoWorkspace.g[42]*acadoWorkspace.x[42] + acadoWorkspace.g[43]*acadoWorkspace.x[43] + acadoWorkspace.g[44]*acadoWorkspace.x[44] + acadoWorkspace.g[45]*acadoWorkspace.x[45] + acadoWorkspace.g[46]*acadoWorkspace.x[46] + acadoWorkspace.g[47]*acadoWorkspace.x[47] + acadoWorkspace.g[48]*acadoWorkspace.x[48] + acadoWorkspace.g[49]*acadoWorkspace.x[49] + acadoWorkspace.g[50]*acadoWorkspace.x[50] + acadoWorkspace.g[51]*acadoWorkspace.x[51] + acadoWorkspace.g[52]*acadoWorkspace.x[52] + acadoWorkspace.g[53]*acadoWorkspace.x[53] + acadoWorkspace.g[54]*acadoWorkspace.x[54] + acadoWorkspace.g[55]*acadoWorkspace.x[55] + acadoWorkspace.g[56]*acadoWorkspace.x[56] + acadoWorkspace.g[57]*acadoWorkspace.x[57] + acadoWorkspace.g[58]*acadoWorkspace.x[58] + acadoWorkspace.g[59]*acadoWorkspace.x[59];
kkt = fabs( kkt );
for (index = 0; index < 60; ++index)
{
prd = acadoWorkspace.y[index];
if (prd > 1e-12)
kkt += fabs(acadoWorkspace.lb[index] * prd);
else if (prd < -1e-12)
kkt += fabs(acadoWorkspace.ub[index] * prd);
}
return kkt;
}

real_t acado_getObjective(  )
{
real_t objVal;

int lRun1;
/** Row vector of size: 9 */
real_t tmpDy[ 9 ];

/** Row vector of size: 6 */
real_t tmpDyN[ 6 ];

for (lRun1 = 0; lRun1 < 20; ++lRun1)
{
acadoWorkspace.objValueIn[0] = acadoVariables.x[lRun1 * 6];
acadoWorkspace.objValueIn[1] = acadoVariables.x[lRun1 * 6 + 1];
acadoWorkspace.objValueIn[2] = acadoVariables.x[lRun1 * 6 + 2];
acadoWorkspace.objValueIn[3] = acadoVariables.x[lRun1 * 6 + 3];
acadoWorkspace.objValueIn[4] = acadoVariables.x[lRun1 * 6 + 4];
acadoWorkspace.objValueIn[5] = acadoVariables.x[lRun1 * 6 + 5];
acadoWorkspace.objValueIn[6] = acadoVariables.u[lRun1 * 3];
acadoWorkspace.objValueIn[7] = acadoVariables.u[lRun1 * 3 + 1];
acadoWorkspace.objValueIn[8] = acadoVariables.u[lRun1 * 3 + 2];

acado_evaluateLSQ( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.Dy[lRun1 * 9] = acadoWorkspace.objValueOut[0] - acadoVariables.y[lRun1 * 9];
acadoWorkspace.Dy[lRun1 * 9 + 1] = acadoWorkspace.objValueOut[1] - acadoVariables.y[lRun1 * 9 + 1];
acadoWorkspace.Dy[lRun1 * 9 + 2] = acadoWorkspace.objValueOut[2] - acadoVariables.y[lRun1 * 9 + 2];
acadoWorkspace.Dy[lRun1 * 9 + 3] = acadoWorkspace.objValueOut[3] - acadoVariables.y[lRun1 * 9 + 3];
acadoWorkspace.Dy[lRun1 * 9 + 4] = acadoWorkspace.objValueOut[4] - acadoVariables.y[lRun1 * 9 + 4];
acadoWorkspace.Dy[lRun1 * 9 + 5] = acadoWorkspace.objValueOut[5] - acadoVariables.y[lRun1 * 9 + 5];
acadoWorkspace.Dy[lRun1 * 9 + 6] = acadoWorkspace.objValueOut[6] - acadoVariables.y[lRun1 * 9 + 6];
acadoWorkspace.Dy[lRun1 * 9 + 7] = acadoWorkspace.objValueOut[7] - acadoVariables.y[lRun1 * 9 + 7];
acadoWorkspace.Dy[lRun1 * 9 + 8] = acadoWorkspace.objValueOut[8] - acadoVariables.y[lRun1 * 9 + 8];
}
acadoWorkspace.objValueIn[0] = acadoVariables.x[120];
acadoWorkspace.objValueIn[1] = acadoVariables.x[121];
acadoWorkspace.objValueIn[2] = acadoVariables.x[122];
acadoWorkspace.objValueIn[3] = acadoVariables.x[123];
acadoWorkspace.objValueIn[4] = acadoVariables.x[124];
acadoWorkspace.objValueIn[5] = acadoVariables.x[125];
acado_evaluateLSQEndTerm( acadoWorkspace.objValueIn, acadoWorkspace.objValueOut );
acadoWorkspace.DyN[0] = acadoWorkspace.objValueOut[0] - acadoVariables.yN[0];
acadoWorkspace.DyN[1] = acadoWorkspace.objValueOut[1] - acadoVariables.yN[1];
acadoWorkspace.DyN[2] = acadoWorkspace.objValueOut[2] - acadoVariables.yN[2];
acadoWorkspace.DyN[3] = acadoWorkspace.objValueOut[3] - acadoVariables.yN[3];
acadoWorkspace.DyN[4] = acadoWorkspace.objValueOut[4] - acadoVariables.yN[4];
acadoWorkspace.DyN[5] = acadoWorkspace.objValueOut[5] - acadoVariables.yN[5];
objVal = 0.0000000000000000e+00;
for (lRun1 = 0; lRun1 < 20; ++lRun1)
{
tmpDy[0] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 9] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 18] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 27] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 36] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 45] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 54] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 63] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 72];
tmpDy[1] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 1] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 10] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 19] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 28] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 37] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 46] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 55] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 64] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 73];
tmpDy[2] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 2] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 11] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 20] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 29] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 38] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 47] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 56] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 65] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 74];
tmpDy[3] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 3] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 12] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 21] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 30] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 39] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 48] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 57] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 66] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 75];
tmpDy[4] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 4] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 13] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 22] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 31] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 40] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 49] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 58] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 67] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 76];
tmpDy[5] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 5] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 14] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 23] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 32] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 41] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 50] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 59] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 68] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 77];
tmpDy[6] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 6] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 15] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 24] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 33] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 42] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 51] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 60] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 69] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 78];
tmpDy[7] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 7] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 16] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 25] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 34] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 43] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 52] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 61] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 70] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 79];
tmpDy[8] = + acadoWorkspace.Dy[lRun1 * 9]*acadoVariables.W[lRun1 * 81 + 8] + acadoWorkspace.Dy[lRun1 * 9 + 1]*acadoVariables.W[lRun1 * 81 + 17] + acadoWorkspace.Dy[lRun1 * 9 + 2]*acadoVariables.W[lRun1 * 81 + 26] + acadoWorkspace.Dy[lRun1 * 9 + 3]*acadoVariables.W[lRun1 * 81 + 35] + acadoWorkspace.Dy[lRun1 * 9 + 4]*acadoVariables.W[lRun1 * 81 + 44] + acadoWorkspace.Dy[lRun1 * 9 + 5]*acadoVariables.W[lRun1 * 81 + 53] + acadoWorkspace.Dy[lRun1 * 9 + 6]*acadoVariables.W[lRun1 * 81 + 62] + acadoWorkspace.Dy[lRun1 * 9 + 7]*acadoVariables.W[lRun1 * 81 + 71] + acadoWorkspace.Dy[lRun1 * 9 + 8]*acadoVariables.W[lRun1 * 81 + 80];
objVal += + acadoWorkspace.Dy[lRun1 * 9]*tmpDy[0] + acadoWorkspace.Dy[lRun1 * 9 + 1]*tmpDy[1] + acadoWorkspace.Dy[lRun1 * 9 + 2]*tmpDy[2] + acadoWorkspace.Dy[lRun1 * 9 + 3]*tmpDy[3] + acadoWorkspace.Dy[lRun1 * 9 + 4]*tmpDy[4] + acadoWorkspace.Dy[lRun1 * 9 + 5]*tmpDy[5] + acadoWorkspace.Dy[lRun1 * 9 + 6]*tmpDy[6] + acadoWorkspace.Dy[lRun1 * 9 + 7]*tmpDy[7] + acadoWorkspace.Dy[lRun1 * 9 + 8]*tmpDy[8];
}

tmpDyN[0] = + acadoWorkspace.DyN[0]*acadoVariables.WN[0];
tmpDyN[1] = + acadoWorkspace.DyN[1]*acadoVariables.WN[7];
tmpDyN[2] = + acadoWorkspace.DyN[2]*acadoVariables.WN[14];
tmpDyN[3] = + acadoWorkspace.DyN[3]*acadoVariables.WN[21];
tmpDyN[4] = + acadoWorkspace.DyN[4]*acadoVariables.WN[28];
tmpDyN[5] = + acadoWorkspace.DyN[5]*acadoVariables.WN[35];
objVal += + acadoWorkspace.DyN[0]*tmpDyN[0] + acadoWorkspace.DyN[1]*tmpDyN[1] + acadoWorkspace.DyN[2]*tmpDyN[2] + acadoWorkspace.DyN[3]*tmpDyN[3] + acadoWorkspace.DyN[4]*tmpDyN[4] + acadoWorkspace.DyN[5]*tmpDyN[5];

objVal *= 0.5;
return objVal;
}

