#include "matheuristics.h"


/////////////////////////////////// Hard-fixing /////////////////////////////////////////

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
    
    double absolute_best_value = inst->best_val;
    while (!time_limit_expired(inst))
    {
        fixing(env, lp, inst, pfix);
        if (branch_and_cut(inst, env, lp)) print_error("Error in branch_and_cut()");
        unfixing(env, lp, inst);
        if(VERBOSE >= 100 && absolute_best_value > inst->best_val)    printf("Found better solution with value: %f\n", inst->best_val);
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

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
        {
            if((j == succ[i] || i == succ[j]) && random01() < pfix)
            {
                index = xpos(i,j, inst);
                if (CPXchgbds(env, lp, 1, &index, &bound, &bd)) print_error("Error in CPXchgbds() fixing");
            }
        }
    }
    
    free(succ);
}

void unfixing(CPXENVptr env, CPXLPptr lp, instance *inst){

    double bd = 0; int index; char bound = 'L';

    for ( int i = 0; i < inst->nnodes; i++ )
    {
        for ( int j = i+1; j < inst->nnodes; j++ )
        {
            index = xpos(i,j, inst);
            if (CPXchgbds(env, lp, 1, &index, &bound, &bd)) print_error("Error in CPXchgbds() unfixing");
        }
    }

}


/////////////////////////////////// Local branching /////////////////////////////////////////

int local_branching(instance *inst, int k){
    
    if (inst->timelimit == CPX_INFBOUND)
        print_error("Time limit not specified or invalid. \nFor Local Branching method you should specify the time limit as the goal is to minimize the primal integral in given time limit.\n");

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
    
    double absolute_best_value = inst->best_val;
    int k_initial = k;

    while (!time_limit_expired(inst))
    {
        add_local_branching_constraints(env, lp, inst, k);
        if (branch_and_cut(inst, env, lp)) print_error("Error in branch_and_cut()");
        remove_last_constraints(env, lp);

        if(absolute_best_value <= inst->best_val){
            k += k_initial; // if stuck in a local minimum, increase the size of the neighborhood
            if(VERBOSE >= 100)    printf("Stuck in a local minimum with value: %f\n\t[!]So k is increased +%d (k=%d) for the next iteration\n", inst->best_val, k_initial, k);
        }
        else {
            k = k_initial; // if found a better solution, reset the size of the neighborhood
            if(VERBOSE >= 100)    printf("Found better solution with value: %f\n", inst->best_val);
            absolute_best_value = inst->best_val;
        }
    }
    

    // free and close cplex model   
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env); 

    return 0; // or an appropriate nonzero error code

}


void add_local_branching_constraints(CPXENVptr env, CPXLPptr lp, instance *inst, int k){
    
    int *succ = (int *) calloc(inst->nnodes, sizeof(int));
    store_succ(inst, succ, inst->best_sol);

	int *index = (int *) calloc(inst->ncols, sizeof(int));
	double *value = (double *) calloc(inst->ncols, sizeof(double));
	double rhs = inst->nnodes - k;
	char sense = 'G';                            // 'G' for greater or equal
	int nnz = 0;

	for ( int i = 0; i < inst->nnodes; i++ )
	{
        for (int j = i+1; j < inst->nnodes; j++)
        {
            if (j == succ[i] || i == succ[j]) // if Xij = 1
            {
                index[nnz] = xpos(i,j, inst);
                value[nnz] = 1.0;
                nnz++;
            }
        }
	}
	int izero = 0;

	if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, NULL) ) print_error("CPXaddrows(): error");
	
    free(succ);
	free(value);
	free(index);

}


void remove_last_constraints(CPXENVptr env, CPXLPptr lp){
        
    int begin = CPXgetnumrows(env, lp) - 1;
    int end = CPXgetnumrows(env, lp) - 1;

    if (CPXdelrows(env, lp, begin, end)) print_error("Error in CPXdelrows()");

}