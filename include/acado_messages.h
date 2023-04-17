#pragma once

#include <iostream>

class ACADO_MSG {
    public:
        static void showStatus(int status){
            
            // Defines pairs of acado return values and messages.
            AcadoReturnValueList returnValueList[] =
            {
            /* miscellaneous */
            { SUCCESSFUL_RETURN, "Successful return" },
            { RET_DIV_BY_ZERO, "Division by zero" },
            { RET_INDEX_OUT_OF_BOUNDS, "Index out of bounds" },
            { RET_INVALID_ARGUMENTS, "At least one of the arguments is invalid" },
            { RET_ERROR_UNDEFINED, "Error number undefined" },
            { RET_WARNING_UNDEFINED, "Warning number undefined" },
            { RET_INFO_UNDEFINED, "Info number undefined" },
            { RET_EWI_UNDEFINED, "Error/warning/info number undefined" },
            { RET_AVAILABLE_WITH_LINUX_ONLY, "This function is available under Linux only" },
            { RET_UNKNOWN_BUG, "The error occured is not yet known" },
            { RET_PRINTLEVEL_CHANGED, "Print level changed" },
            { RET_NOT_YET_IMPLEMENTED, "Requested function is not yet implemented." },
            /* Indexlist */
            { RET_INDEXLIST_MUST_BE_REORDERD, "Index list has to be reordered" },
            { RET_INDEXLIST_EXCEEDS_MAX_LENGTH, "Index list exceeds its maximal physical length" },
            { RET_INDEXLIST_CORRUPTED, "Index list corrupted" },
            { RET_INDEXLIST_OUTOFBOUNDS, "Physical index is out of bounds" },
            { RET_INDEXLIST_ADD_FAILED, "Adding indices from another index set failed" },
            { RET_INDEXLIST_INTERSECT_FAILED, "Intersection with another index set failed" },
            /* SubjectTo / Bounds / Constraints */
            { RET_INDEX_ALREADY_OF_DESIRED_STATUS, "Index is already of desired status" },
            { RET_SWAPINDEX_FAILED, "Cannot swap between different indexsets" },
            { RET_ADDINDEX_FAILED, "Adding index to index set failed" },
            { RET_NOTHING_TO_DO, "Nothing to do" },
            { RET_SETUP_BOUND_FAILED, "Setting up bound index failed" },
            { RET_SETUP_CONSTRAINT_FAILED, "Setting up constraint index failed" },
            { RET_MOVING_BOUND_FAILED, "Moving bound between index sets failed" },
            { RET_MOVING_CONSTRAINT_FAILED, "Moving constraint between index sets failed" },
            /* QProblem */
            { RET_QP_ALREADY_INITIALISED, "QProblem has already been initialised" },
            { RET_NO_INIT_WITH_STANDARD_SOLVER, "Initialisation via extern QP solver is not yet implemented" },
            { RET_RESET_FAILED, "Reset failed" },
            { RET_INIT_FAILED, "Initialisation failed" },
            { RET_INIT_FAILED_TQ, "Initialisation failed due to TQ factorisation" },
            { RET_INIT_FAILED_CHOLESKY, "Initialisation failed due to Cholesky decomposition" },
            { RET_INIT_FAILED_HOTSTART, "Initialisation failed! QP could not be solved!" },
            { RET_INIT_FAILED_INFEASIBILITY, "Initial QP could not be solved due to infeasibility!" },
            { RET_INIT_FAILED_UNBOUNDEDNESS, "Initial QP could not be solved due to unboundedness!" },
            { RET_INIT_SUCCESSFUL, "Initialisation done" },
            { RET_OBTAINING_WORKINGSET_FAILED, "Failed to obtain working set for auxiliary QP" },
            { RET_SETUP_WORKINGSET_FAILED, "Failed to setup working set for auxiliary QP" },
            { RET_SETUP_AUXILIARYQP_FAILED, "Failed to setup auxiliary QP for initialised homotopy" },
            { RET_NO_EXTERN_SOLVER, "No extern QP solver available" },
            { RET_QP_UNBOUNDED, "QP is unbounded" },
            { RET_QP_INFEASIBLE, "QP is infeasible" },
            { RET_QP_NOT_SOLVED, "Problems occured while solving QP with standard solver" },
            { RET_QP_SOLVED, "QP successfully solved" },
            { RET_UNABLE_TO_SOLVE_QP, "Problems occured while solving QP" },
            { RET_INITIALISATION_STARTED, "Starting problem initialisation..." },
            { RET_HOTSTART_FAILED, "Unable to perform homotopy due to internal error" },
            { RET_HOTSTART_FAILED_TO_INIT, "Unable to initialise problem" },
            { RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED, "Unable to perform homotopy as previous QP is not solved" },
            { RET_ITERATION_STARTED, "Iteration" },
            { RET_SHIFT_DETERMINATION_FAILED, "Determination of shift of the QP data failed" },
            { RET_STEPDIRECTION_DETERMINATION_FAILED, "Determination of step direction failed" },
            { RET_STEPLENGTH_DETERMINATION_FAILED, "Determination of step direction failed" },
            { RET_OPTIMAL_SOLUTION_FOUND, "Optimal solution of neighbouring QP found" },
            { RET_HOMOTOPY_STEP_FAILED, "Unable to perform homotopy step" },
            { RET_HOTSTART_STOPPED_INFEASIBILITY, "Premature homotopy termination because QP is infeasible" },
            { RET_HOTSTART_STOPPED_UNBOUNDEDNESS, "Premature homotopy termination because QP is unbounded" },
            { RET_WORKINGSET_UPDATE_FAILED, "Unable to update working sets according to initial guesses" },
            { RET_MAX_NWSR_REACHED, "Maximum number of working set recalculations performed" },
            { RET_CONSTRAINTS_NOT_SPECIFIED, "Problem does comprise constraints! You have to specify new constraints' bounds" },
            { RET_INVALID_FACTORISATION_FLAG, "Invalid factorisation flag" },
            { RET_UNABLE_TO_SAVE_QPDATA, "Unable to save QP data" },
            { RET_STEPDIRECTION_FAILED_TQ, "Abnormal termination due to TQ factorisation" },
            { RET_STEPDIRECTION_FAILED_CHOLESKY, "Abnormal termination due to Cholesky factorisation" },
            { RET_CYCLING_DETECTED, "Cycling detected" },
            { RET_CYCLING_NOT_RESOLVED, "Cycling cannot be resolved, QP is probably infeasible" },
            { RET_CYCLING_RESOLVED, "Cycling probably resolved" },
            { RET_STEPSIZE, "" },
            { RET_STEPSIZE_NONPOSITIVE, "" },
            { RET_SETUPSUBJECTTOTYPE_FAILED, "Setup of SubjectToTypes failed" },
            { RET_ADDCONSTRAINT_FAILED, "Addition of constraint to working set failed" },
            { RET_ADDCONSTRAINT_FAILED_INFEASIBILITY, "Addition of constraint to working set failed" },
            { RET_ADDBOUND_FAILED, "Addition of bound to working set failed" },
            { RET_ADDBOUND_FAILED_INFEASIBILITY, "Addition of bound to working set failed" },
            { RET_REMOVECONSTRAINT_FAILED, "Removal of constraint from working set failed" },
            { RET_REMOVEBOUND_FAILED, "Removal of bound from working set failed" },
            { RET_REMOVE_FROM_ACTIVESET, "Removing from active set:" },
            { RET_ADD_TO_ACTIVESET, "Adding to active set:" },
            { RET_REMOVE_FROM_ACTIVESET_FAILED, "Removing from active set failed" },
            { RET_ADD_TO_ACTIVESET_FAILED, "Adding to active set failed" },
            { RET_CONSTRAINT_ALREADY_ACTIVE, "Constraint is already active" },
            { RET_ALL_CONSTRAINTS_ACTIVE, "All constraints are active, no further constraint can be added" },
            { RET_LINEARLY_DEPENDENT, "New bound/constraint is linearly dependent" },
            { RET_LINEARLY_INDEPENDENT, "New bound/constraint is linearly independent" },
            { RET_LI_RESOLVED, "Linear independence of active contraint matrix successfully resolved" },
            { RET_ENSURELI_FAILED, "Failed to ensure linear indepence of active contraint matrix" },
            { RET_ENSURELI_FAILED_TQ, "Abnormal termination due to TQ factorisation" },
            { RET_ENSURELI_FAILED_NOINDEX, "No index found, QP is probably infeasible" },
            { RET_ENSURELI_FAILED_CYCLING, "Cycling detected, QP is probably infeasible" },
            { RET_BOUND_ALREADY_ACTIVE, "Bound is already active" },
            { RET_ALL_BOUNDS_ACTIVE, "All bounds are active, no further bound can be added" },
            { RET_CONSTRAINT_NOT_ACTIVE, "Constraint is not active" },
            { RET_BOUND_NOT_ACTIVE, "Bound is not active" },
            { RET_HESSIAN_NOT_SPD, "Projected Hessian matrix not positive definite" },
            { RET_MATRIX_SHIFT_FAILED, "Unable to update matrices or to transform vectors" },
            { RET_MATRIX_FACTORISATION_FAILED, "Unable to calculate new matrix factorisations" },
            { RET_PRINT_ITERATION_FAILED, "Unable to print information on current iteration" },
            { RET_NO_GLOBAL_MESSAGE_OUTPUTFILE, "No global message output file initialised" },
            /* Utils */
            { RET_UNABLE_TO_OPEN_FILE, "Unable to open file" },
            { RET_UNABLE_TO_WRITE_FILE, "Unable to write into file" },
            { RET_UNABLE_TO_READ_FILE, "Unable to read from file" },
            { RET_FILEDATA_INCONSISTENT, "File contains inconsistent data" },
            /* SolutionAnalysis */
            { RET_NO_SOLUTION, "QP solution does not satisfy KKT optimality conditions" },
            { RET_INACCURATE_SOLUTION, "KKT optimality conditions not satisfied to sufficient accuracy" },
            };
            
            std::cout << std::endl << returnValueList[ status ].data << std::endl << std::endl;
        }

    private:
        // Defines symbols for acado return values.
        // Important: All return values are assumed to be nonnegative!
        enum returnValue
        {
        /* miscellaneous */
        SUCCESSFUL_RETURN = 0,							/**< Successful return. */
        RET_DIV_BY_ZERO,		   						/**< Division by zero. */
        RET_INDEX_OUT_OF_BOUNDS,						/**< Index out of bounds. */
        RET_INVALID_ARGUMENTS,							/**< At least one of the arguments is invalid. */
        RET_ERROR_UNDEFINED,							/**< Error number undefined. */
        RET_WARNING_UNDEFINED,							/**< Warning number undefined. */
        RET_INFO_UNDEFINED,								/**< Info number undefined. */
        RET_EWI_UNDEFINED,								/**< Error/warning/info number undefined. */
        RET_AVAILABLE_WITH_LINUX_ONLY,					/**< This function is available under Linux only. */
        RET_UNKNOWN_BUG,								/**< The error occured is not yet known. */
        RET_PRINTLEVEL_CHANGED,							/**< 10 Print level changed. */
        RET_NOT_YET_IMPLEMENTED,						/**< Requested function is not yet implemented in this version of qpOASES. */
        /* Indexlist */
        RET_INDEXLIST_MUST_BE_REORDERD,					/**< Index list has to be reordered. */
        RET_INDEXLIST_EXCEEDS_MAX_LENGTH,				/**< Index list exceeds its maximal physical length. */
        RET_INDEXLIST_CORRUPTED,						/**< Index list corrupted. */
        RET_INDEXLIST_OUTOFBOUNDS,						/**< Physical index is out of bounds. */
        RET_INDEXLIST_ADD_FAILED,						/**< Adding indices from another index set failed. */
        RET_INDEXLIST_INTERSECT_FAILED,					/**< Intersection with another index set failed. */
        /* SubjectTo / Bounds / Constraints */
        RET_INDEX_ALREADY_OF_DESIRED_STATUS,			/**< Index is already of desired status. */
        RET_ADDINDEX_FAILED,							/**< Cannot swap between different indexsets. */
        RET_SWAPINDEX_FAILED,							/**< 20 Adding index to index set failed. */
        RET_NOTHING_TO_DO,								/**< Nothing to do. */
        RET_SETUP_BOUND_FAILED,							/**< Setting up bound index failed. */
        RET_SETUP_CONSTRAINT_FAILED,					/**< Setting up constraint index failed. */
        RET_MOVING_BOUND_FAILED,						/**< Moving bound between index sets failed. */
        RET_MOVING_CONSTRAINT_FAILED,					/**< Moving constraint between index sets failed. */
        /* QProblem */
        RET_QP_ALREADY_INITIALISED,						/**< QProblem has already been initialised. */
        RET_NO_INIT_WITH_STANDARD_SOLVER,				/**< Initialisation via extern QP solver is not yet implemented. */
        RET_RESET_FAILED,								/**< Reset failed. */
        RET_INIT_FAILED,								/**< Initialisation failed. */
        RET_INIT_FAILED_TQ,								/**< 30 Initialisation failed due to TQ factorisation. */
        RET_INIT_FAILED_CHOLESKY,						/**< Initialisation failed due to Cholesky decomposition. */
        RET_INIT_FAILED_HOTSTART,						/**< Initialisation failed! QP could not be solved! */
        RET_INIT_FAILED_INFEASIBILITY,					/**< Initial QP could not be solved due to infeasibility! */
        RET_INIT_FAILED_UNBOUNDEDNESS,					/**< Initial QP could not be solved due to unboundedness! */
        RET_INIT_SUCCESSFUL,							/**< Initialisation done. */
        RET_OBTAINING_WORKINGSET_FAILED,				/**< Failed to obtain working set for auxiliary QP. */
        RET_SETUP_WORKINGSET_FAILED,					/**< Failed to setup working set for auxiliary QP. */
        RET_SETUP_AUXILIARYQP_FAILED,					/**< Failed to setup auxiliary QP for initialised homotopy. */
        RET_NO_EXTERN_SOLVER,							/**< No extern QP solver available. */
        RET_QP_UNBOUNDED,								/**< 40 QP is unbounded. */
        RET_QP_INFEASIBLE,								/**< QP is infeasible. */
        RET_QP_NOT_SOLVED,								/**< Problems occured while solving QP with standard solver. */
        RET_QP_SOLVED,									/**< QP successfully solved. */
        RET_UNABLE_TO_SOLVE_QP,							/**< Problems occured while solving QP. */
        RET_INITIALISATION_STARTED,						/**< Starting problem initialisation. */
        RET_HOTSTART_FAILED,							/**< Unable to perform homotopy due to internal error. */
        RET_HOTSTART_FAILED_TO_INIT,					/**< Unable to initialise problem. */
        RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED,		/**< Unable to perform homotopy as previous QP is not solved. */
        RET_ITERATION_STARTED,							/**< Iteration... */
        RET_SHIFT_DETERMINATION_FAILED,					/**< 50 Determination of shift of the QP data failed. */
        RET_STEPDIRECTION_DETERMINATION_FAILED,			/**< Determination of step direction failed. */
        RET_STEPLENGTH_DETERMINATION_FAILED,			/**< Determination of step direction failed. */
        RET_OPTIMAL_SOLUTION_FOUND,						/**< Optimal solution of neighbouring QP found. */
        RET_HOMOTOPY_STEP_FAILED,						/**< Unable to perform homotopy step. */
        RET_HOTSTART_STOPPED_INFEASIBILITY,				/**< Premature homotopy termination because QP is infeasible. */
        RET_HOTSTART_STOPPED_UNBOUNDEDNESS,				/**< Premature homotopy termination because QP is unbounded. */
        RET_WORKINGSET_UPDATE_FAILED,					/**< Unable to update working sets according to initial guesses. */
        RET_MAX_NWSR_REACHED,							/**< Maximum number of working set recalculations performed. */
        RET_CONSTRAINTS_NOT_SPECIFIED,					/**< Problem does comprise constraints! You also have to specify new constraints' bounds. */
        RET_INVALID_FACTORISATION_FLAG,					/**< 60 Invalid factorisation flag. */
        RET_UNABLE_TO_SAVE_QPDATA,						/**< Unable to save QP data. */
        RET_STEPDIRECTION_FAILED_TQ,					/**< Abnormal termination due to TQ factorisation. */
        RET_STEPDIRECTION_FAILED_CHOLESKY,				/**< Abnormal termination due to Cholesky factorisation. */
        RET_CYCLING_DETECTED,							/**< Cycling detected. */
        RET_CYCLING_NOT_RESOLVED,						/**< Cycling cannot be resolved, QP probably infeasible. */
        RET_CYCLING_RESOLVED,							/**< Cycling probably resolved. */
        RET_STEPSIZE,									/**< For displaying performed stepsize. */
        RET_STEPSIZE_NONPOSITIVE,						/**< For displaying non-positive stepsize. */
        RET_SETUPSUBJECTTOTYPE_FAILED,					/**< Setup of SubjectToTypes failed. */
        RET_ADDCONSTRAINT_FAILED,						/**< 70 Addition of constraint to working set failed. */
        RET_ADDCONSTRAINT_FAILED_INFEASIBILITY,			/**< Addition of constraint to working set failed (due to QP infeasibility). */
        RET_ADDBOUND_FAILED,							/**< Addition of bound to working set failed. */
        RET_ADDBOUND_FAILED_INFEASIBILITY,				/**< Addition of bound to working set failed (due to QP infeasibility). */
        RET_REMOVECONSTRAINT_FAILED,					/**< Removal of constraint from working set failed. */
        RET_REMOVEBOUND_FAILED,							/**< Removal of bound from working set failed. */
        RET_REMOVE_FROM_ACTIVESET,						/**< Removing from active set... */
        RET_ADD_TO_ACTIVESET,							/**< Adding to active set... */
        RET_REMOVE_FROM_ACTIVESET_FAILED,				/**< Removing from active set failed. */
        RET_ADD_TO_ACTIVESET_FAILED,					/**< Adding to active set failed. */
        RET_CONSTRAINT_ALREADY_ACTIVE,					/**< 80 Constraint is already active. */
        RET_ALL_CONSTRAINTS_ACTIVE,						/**< All constraints are active, no further constraint can be added. */
        RET_LINEARLY_DEPENDENT,							/**< New bound/constraint is linearly dependent. */
        RET_LINEARLY_INDEPENDENT,						/**< New bound/constraint is linearly independent. */
        RET_LI_RESOLVED,								/**< Linear independence of active contraint matrix successfully resolved. */
        RET_ENSURELI_FAILED,							/**< Failed to ensure linear indepence of active contraint matrix. */
        RET_ENSURELI_FAILED_TQ,							/**< Abnormal termination due to TQ factorisation. */
        RET_ENSURELI_FAILED_NOINDEX,					/**< No index found, QP probably infeasible. */
        RET_ENSURELI_FAILED_CYCLING,					/**< Cycling detected, QP probably infeasible. */
        RET_BOUND_ALREADY_ACTIVE,						/**< Bound is already active. */
        RET_ALL_BOUNDS_ACTIVE,							/**< 90 All bounds are active, no further bound can be added. */
        RET_CONSTRAINT_NOT_ACTIVE,						/**< Constraint is not active. */
        RET_BOUND_NOT_ACTIVE,							/**< Bound is not active. */
        RET_HESSIAN_NOT_SPD,							/**< Projected Hessian matrix not positive definite. */
        RET_MATRIX_SHIFT_FAILED,						/**< Unable to update matrices or to transform vectors. */
        RET_MATRIX_FACTORISATION_FAILED,				/**< Unable to calculate new matrix factorisations. */
        RET_PRINT_ITERATION_FAILED,						/**< Unable to print information on current iteration. */
        RET_NO_GLOBAL_MESSAGE_OUTPUTFILE,				/**< No global message output file initialised. */
        /* Utils */
        RET_UNABLE_TO_OPEN_FILE,						/**< Unable to open file. */
        RET_UNABLE_TO_WRITE_FILE,						/**< Unable to write into file. */
        RET_UNABLE_TO_READ_FILE,						/**< 100 Unable to read from file. */
        RET_FILEDATA_INCONSISTENT,						/**< File contains inconsistent data. */
        /* SolutionAnalysis */
        RET_NO_SOLUTION, 								/**< QP solution does not satisfy KKT optimality conditions. */
        RET_INACCURATE_SOLUTION							/**< KKT optimality conditions not satisfied to sufficient accuracy. */
        };

        // Struct for message info
        typedef struct {
            returnValue key;							/**< Global return value. */
            const char* data;							/**< Corresponding message. */
        } AcadoReturnValueList;

};

    // ADDED by SW    

// ADDED by SW
// Print array as matrix
// #include <cstdio>
// #include <iostream>
// #include <string>
// #include <vector>
// void printArrayAsMatrix(float *A, int rows, int cols) {
// 		int i, j;
// 		for (i = 0; i < rows; i++) {
// 				for (j = 0; j < cols; j++)
// 						printf("  % 9.2f", A[i*cols+j]);
// 				printf(",\n");
// 		}
// 		printf(",\n");
// }

//     std::cout << "H: " << std::endl;
//     //(sizeof(acadoWorkspace.H)/sizeof(*acadoWorkspace.H)) << std::endl;
// 	printArrayAsMatrix(acadoWorkspace.H, 80, 80);
//     std::vector<std::string> str {"g: ", "lb: ", "ub: ", "y: "};
//     std::vector<float *> ptrs {acadoWorkspace.g, acadoWorkspace.lb, acadoWorkspace.ub, acadoWorkspace.y};
//     for (int i = 0; i <= 3; i++){
//         std::cout << str[i] << std::endl;
//         //(sizeof(ptrs[i])/sizeof(*ptrs[i])) << std::endl;
//         printArrayAsMatrix(ptrs[i], 1, 80);
//     }
//     std::cout << "retVal: " << retVal << std::endl;
