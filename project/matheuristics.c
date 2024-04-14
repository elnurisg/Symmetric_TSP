#include "matheuristics.h"

double second();  

int hard_fixing(instance *inst, double pfix){

    if (inst->timelimit == CPX_INFBOUND)
        print_error("Time limit not specified or invalid. \nFor Hard-fixing method you should specify the time limit as the goal is to minimize the primal integral in given time limit.\n");

    	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);

    inst->ncols = CPXgetnumcols(env, lp);

	if ( greedy_heuristic(inst, 2, 0) ) print_error(" error within greedy_heuristic()");
	if ( two_opt_refining_heuristic(inst, inst->best_sol, 0) ) print_error(" error within two_opt_refining_heuristic()");
	calculate_best_val(inst);
    
    double initial_best_value = inst->best_val;
    while (second() - inst->tstart < inst->timelimit)
    {
        fixing(env, lp, inst, pfix);
        if (branch_and_cut(inst, env, lp)) print_error("Error in branch_and_cut()");
        unfixing(env, lp, inst);
        if(VERBOSE >= 100 && initial_best_value > inst->best_val)    printf("Found better solution with value: %f\n", inst->best_val);
    }
    

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return 0; // or an appropriate nonzero error code

}

void fixing(CPXENVptr env, CPXLPptr lp, instance *inst, double pfix){

    double bd = 1; int index; char bound = 'L';
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    store_succ(inst, succ, inst->best_sol);
    // double lb_value, ub_value;

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
        {
            if((j == succ[i] || i == succ[j]) && random01() < pfix)
            {
                index = xpos(i,j, inst);
                if (CPXchgbds(env, lp, 1, &index, &bound, &bd)) print_error("Error in CPXchgbds() fixing");
                // if (CPXgetlb(env, lp, &lb_value, index, index)) print_error("Error in CPXgetlb()");
                // if (CPXgetub(env, lp, &ub_value, index, index)) print_error("Error in CPXgetub()");

                // printf("Variable x(%d,%d) lower bound: %f, upper bound: %f and pos is %d\n", i+1, j+1, lb_value, ub_value, index);
            }
        }
    }
    
    free(succ);
}

void unfixing(CPXENVptr env, CPXLPptr lp, instance *inst){

    double bd = 0; int index; char bound = 'L';
    // double lb_value, ub_value;

    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            index = xpos(i,j, inst);
            if (CPXchgbds(env, lp, 1, &index, &bound, &bd)) print_error("Error in CPXchgbds() unfixing");

            // if (CPXgetlb(env, lp, &lb_value, index, index)) print_error("Error in CPXgetlb()");
            // if (CPXgetub(env, lp, &ub_value, index, index)) print_error("Error in CPXgetub()");

            // printf("UNFIX: Variable x(%d,%d) lower bound: %f, upper bound: %f and pos is %d\n", i+1, j+1, lb_value, ub_value, index);

        }
    }

}