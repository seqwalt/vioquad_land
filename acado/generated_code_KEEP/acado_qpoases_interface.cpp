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


extern "C"
{
#include "acado_common.h"
}

#include "INCLUDE/QProblemB.hpp"
//#include <iostream>
//#include <assert.h>

#if ACADO_COMPUTE_COVARIANCE_MATRIX == 1
#include "INCLUDE/EXTRAS/SolutionAnalysis.hpp"
#endif /* ACADO_COMPUTE_COVARIANCE_MATRIX */

static int acado_nWSR;



#if ACADO_COMPUTE_COVARIANCE_MATRIX == 1
static SolutionAnalysis acado_sa;
#endif /* ACADO_COMPUTE_COVARIANCE_MATRIX */

int acado_solve( void )
{
	acado_nWSR = QPOASES_NWSRMAX;

	QProblemB qp( 60 );

	returnValue retVal = qp.init(acadoWorkspace.H, acadoWorkspace.g, acadoWorkspace.lb, acadoWorkspace.ub, acado_nWSR, acadoWorkspace.y);
	// std::cout << "retVal = " << (int)retVal << std::endl<< std::endl;
	// std::cout << "H = "<< std::endl;
	// for(std::size_t row = 0; row != 60; ++row){
	// 	for(std::size_t col = 0; col != 60; ++col){
	// 		int el = row*60 + col;
	// 		assert(el <= 3600);
	// 		std::cout << acadoWorkspace.H[el] << " ";
	// 	}
	// 	std::cout << std::endl;
	// }
	// std::cout << std::endl<< "g = "<< std::endl;
	// for(std::size_t i = 0; i != 60; ++i) std::cout << acadoWorkspace.g[i] << " ";
	// std::cout << std::endl<< std::endl<< "lb = "<< std::endl;
	// for(std::size_t i = 0; i != 60; ++i) std::cout << acadoWorkspace.lb[i] << " ";
	// std::cout << std::endl<< std::endl<< "ub = "<< std::endl;
	// for(std::size_t i = 0; i != 60; ++i) std::cout << acadoWorkspace.ub[i] << " ";
	// std::cout << std::endl<< std::endl<< "y = "<< std::endl;
	// for(std::size_t i = 0; i != 60; ++i) std::cout << acadoWorkspace.y[i] << " ";
	// std::cout << std::endl<< std::endl<< "nWSR = " << acado_nWSR << std::endl<< std::endl;

  qp.getPrimalSolution( acadoWorkspace.x );
  qp.getDualSolution( acadoWorkspace.y );

#if ACADO_COMPUTE_COVARIANCE_MATRIX == 1

	if (retVal != SUCCESSFUL_RETURN)
		return (int)retVal;

	retVal = acado_sa.getHessianInverse( &qp,var );

#endif /* ACADO_COMPUTE_COVARIANCE_MATRIX */

	return (int)retVal;
}

int acado_getNWSR( void )
{
	return acado_nWSR;
}

const char* acado_getErrorString( int error )
{
	return MessageHandling::getErrorString( error );
}
