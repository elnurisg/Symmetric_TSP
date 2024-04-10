#include "tsp.h"
// #include <time.h>
// #include "convex_hull.h"

// void print_error(const char *err);
double second();       
int time_limit_expired(instance *inst);

int is_fractional(double x) 						// it works for x in [0,1] only
{
	return ( (x > XSMALL) && (x < 1-XSMALL) );
}    

int is_all_integer(int n, const double *x) 			// it works for x_j in [0,1] only
{
	for ( int j = 0; j < n; j++ ) 
	{
		if ( is_fractional(x[j]) ) return 0; 
	}
	return 1;
}                                                                                                                               
                         

int time_limit_expired(instance *inst)	 
{
	double tspan = second() - inst->tstart;
	if (  tspan > inst->timelimit ) 
	{
		if ( VERBOSE >= 100 ) printf("\n\n$$$ time limit of %10.1lf sec.s expired after %10.1lf sec.s $$$\n\n", inst->timelimit, tspan);
		//exit(0); 
		return 1;
	}  
	return 0;
}


int benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp)
{

	double *xstar = (double *) calloc(inst->ncols, sizeof(double));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int *best_succ = (int *) calloc(inst->nnodes, sizeof(int));
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int *) calloc(1, sizeof(int));

	double LB = -CPX_INFBOUND; // lower bound
	if ( greedy_heuristic(inst, 0, 0) ) print_error(" error within greedy_heuristic()");
	if ( two_opt_refining_heuristic(inst, inst->best_sol, 0) ) print_error(" error within two_opt_refining_heuristic()");
	calculate_best_val(inst);
	double UB = inst->best_val; // upper bound
	double objval; double incumbent_value; int iteration = 0;

	do
	{
		CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);
		CPXsetintparam(env, CPX_PARAM_NODELIM, 100);
		if (CPXmipopt(env,lp)) print_error("CPXmipopt() error"); 
		if ( CPXgetbestobjval(env, lp, &objval) ) print_error("CPXgetbestobjval() error");	
		if ( CPXgetx(env, lp, xstar, 0, inst->ncols-1) ) continue;//print_error("CPXgetx() error");
		
		// new LB is assumed to increase as having more constraints
		LB = (LB > objval) ? LB : objval;

		build_sol(xstar, inst, succ, comp, ncomp);

		if (*ncomp <= 1){ // optimal solution is found by CPLEX so we can break the loop
			incumbent_value = calc_succ_value(succ, inst);
			LB = incumbent_value; // due to the type of cost(int and happens even using double),
			//when distance between nodes is very large, there can be numerical issues and
			// they can be not equal and we ensure that they are equal
		}	// as solution is already optimal
		else
		{
			for (int component_num = 1; component_num <= *ncomp; component_num++)
			{
				add_subtour_constraint(NULL, env, lp, inst, comp, component_num, inst->ncols);
				// add a subtour elimination constraint
			}	
			patching_heuristic(NULL, env, lp, inst->ncols, inst, succ, comp, ncomp);
			incumbent_value = calc_succ_value(succ, inst);

		}

		CPXsetdblparam(env, CPX_PARAM_TILIM, (inst->tstart + inst->timelimit - second())); 

		if (incumbent_value < UB)
		{
			UB = incumbent_value;
			best_succ = copy_array(succ, inst->nnodes);
			// store the best succ solution found so far
		}		
		iteration++;
		printf("\nITERATION: %d\t", iteration);
		printf("LB: %f\t",LB); printf("UB: %f\t",UB); (LB == UB) ? printf("||| Solution is Optimal\n") : printf("||| %f%% gap\n", (UB-LB)/LB*100);
	} while ((LB < (1-EPSILON) * UB) && (second() - inst->tstart < inst->timelimit));


	if (VERBOSE >= 100){
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = i+1; j < inst->nnodes; j++ )
			{
				if ( xstar[xpos(i,j,inst)] > 0.5 ) printf("  ... x(%3d,%3d) = 1\n", i+1,j+1);
			}
		}
	}

	double initialized_UB = inst->best_val; // Greedy heuristic + 2-OPT is used as intial upper bound
	store_solution(inst, best_succ, inst->best_sol);
	calculate_best_val(inst);
	if (initialized_UB < inst->best_val) printf("ATTENTION:	Initial Upper Bound (using Greedy heuristic + 2-OPT) is better than the best solution found by CPLEX, please increase the time_limit!\n");

	// CPXwriteprob(env, lp, "model.lp", NULL);  

	free(best_succ);
	free(xstar);
	free(succ);
	free(comp);
	free(ncomp);

	return 0;
}


int branch_and_cut(instance *inst, CPXENVptr env, CPXLPptr lp)
{

	double *xstar = (double *) calloc(inst->ncols, sizeof(double));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int *) calloc(1, sizeof(int));

	double LB = -CPX_INFBOUND; // lower bound
	if ( greedy_heuristic(inst, 0, 0) ) print_error(" error within greedy_heuristic()");
	if ( two_opt_refining_heuristic(inst, inst->best_sol, 0) ) print_error(" error within two_opt_refining_heuristic()");
	calculate_best_val(inst);
	double UB = inst->best_val; // upper bound
	
	CPXsetdblparam(env, CPX_PARAM_CUTUP, UB);
	CPXsetdblparam(env, CPX_PARAM_TILIM, (inst->tstart + inst->timelimit - second())); 
	// CPXsetintparam(env, CPX_PARAM_THREADS, 1); // trying with 1 thread to debug
	// CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
	// CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9) //abort Cplex when gap is below the given %

	// installing a lazyconstraint callback to cut infeasible integer solution
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // means lazy constraints
	if(CPXcallbacksetfunc(env, lp, contextid, my_cut_callback, inst)) print_error("CPXcallbacksetfunc() error");
	
// adding mipstart solution, this is different than PARAM_CUT, it will make sure that Cplex start from this solution
	heur_sol_to_mipstart(env, lp, inst);

	if (CPXmipopt(env, lp)) print_error("CPXmipopt() error"); 

	CPXgetbestobjval(env, lp, &LB);
	CPXgetobjval(env, lp, &UB);
	if(VERBOSE >= 60) {
		printf("\nLB: %f\t",LB); 
		printf("UB: %f\t",UB); 
		(LB == UB) ? printf("||| Solution is Optimal\n") : printf("||| %f%% gap\n", (UB-LB)/LB*100);
	}

	if ( CPXgetx(env, lp, xstar, 0, inst->ncols-1) ) printf("\nCplex could not find any better solution than the given Upper Bound (initialized by using Greedy heuristic + 2-OPT) in the given time limit.\n");
	else{
		build_sol(xstar, inst, succ, comp, ncomp);
		store_solution(inst, succ, inst->best_sol);
		calculate_best_val(inst);
	}

	free(xstar);
	free(succ);
	free(comp);
	free(ncomp);	

	return 0;
}

static int CPXPUBLIC my_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle )
{
	instance* inst = (instance *) userhandle;
	double* xstar = (double *) calloc(inst->ncols, sizeof(double));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int *) calloc(1, sizeof(int));

	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols-1, &objval)) print_error("CPXcallbackgetcandidatepoint() error");

	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
	if(VERBOSE >= 100) printf(".. Thread %2d, Node %5d, Incumbent: %10.2lf, Candidate value: %10.2lf\n", mythread, mynode, incumbent, objval);

	build_sol(xstar, inst, succ, comp, ncomp);
	
	if(*ncomp > 1)
	{
		for (int component_num = 1; component_num <= *ncomp; component_num++)
			add_subtour_constraint(context, NULL, NULL, inst, comp, component_num, inst->ncols);

		patching_heuristic(context, NULL, NULL, inst->ncols, inst, succ, comp, ncomp);
	}

	double succ_value = calc_succ_value(succ, inst);
	if (incumbent > succ_value) //if the succ solution from patching heuristic is better than the incumbent
		post_heuristic(context, inst, succ, succ_value);

	free(xstar);
	free(succ);
	free(comp);
	free(ncomp);	

	return 0;
}

void post_heuristic(CPXCALLBACKCONTEXTptr context, instance *inst, int *succ, double succ_value)
{
	//convert succ to xstar
	double *xheu = (double *) calloc(inst->ncols, sizeof(double)); //all zeros initially
	for (int i=0; i<inst->nnodes; i++)
		xheu[xpos(i,succ[i],inst)] = 1.0;

	int *indices = (int *) malloc(inst->ncols * sizeof(int));
	for (int i = 0; i < inst->ncols; i++) indices[i] = i;

	if(CPXcallbackpostheursoln(context, inst->ncols, indices, xheu, succ_value, CPXCALLBACKSOLUTION_NOCHECK)) print_error("CPXcallbackpostheursoln() error");
	
	free(xheu);
	free(indices);

}

int heur_sol_to_mipstart(CPXENVptr env, CPXLPptr lp, instance *inst)
{
	//convert best_sol to xstar
	double *xheu = (double *) calloc(inst->ncols, sizeof(double)); //all zeros initially
	for (int i=0; i<inst->nnodes; i++)
		xheu[xpos(inst->best_sol[i],inst->best_sol[i+1],inst)] = 1.0;

	int *indices = (int *) malloc(inst->ncols * sizeof(int));
	for (int i = 0; i < inst->ncols; i++) indices[i] = i;
	int effortlevel = CPX_MIPSTART_NOCHECK;
	int beg = 0;
	if(CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, indices, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");

	free(xheu);
	free(indices);
	
	return 0;
}

int TSPopt(instance *inst, int model_type)
{  
	
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	// CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	// if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	// CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	
	// CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0); 

	inst->ncols = CPXgetnumcols(env, lp);

	if (model_type == 0)
	{
		if (benders_loop(inst, env, lp)) print_error("Error in benders_loop()");
	}
	else if (model_type == 1)
	{
		if (branch_and_cut(inst, env, lp)) print_error("Error in branch_and_cut()");
	}
	else
	{
		print_error("Error! Model type is not defined correctly. Please choose 0 for [Benders' loop] and 1 for [Branch and Cut]");
	}	

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return 0; // or an appropriate nonzero error code

}

double delta_cost_patching(int a, int b, instance *inst, int *succ){
	
	double delta_cost = 0;
	double cost_of_old_edges = dist(a, succ[a], inst) + dist(b, succ[b], inst);
	double cost_of_new_edges = dist(a, succ[b], inst) + dist(b, succ[a], inst);

	delta_cost = cost_of_new_edges - cost_of_old_edges;

	return delta_cost;
}

void patching_heuristic(void *context_pointer, void *environment, void* linear_program, int ncols, instance *inst, int *succ, int *comp, int *ncomp){

	double min_delta_cost; double delta_cost;
	int min_a; int min_b;
	int *sol = (int *) calloc(inst->nnodes+1, sizeof(int));

	while (*ncomp != 1)
	{
		min_delta_cost = INFINITY;
		for (int a = 0; a < inst->nnodes; a++)
		{
			for (int b = 0; b < inst->nnodes; b++)
			{
				if (comp[a] > comp[b]) // if a and b is in different component,
				{                     //  '>' is trick to not to check the same pairs twice
					delta_cost = delta_cost_patching(a, b, inst, succ);
					if (delta_cost < min_delta_cost)
					{
						min_delta_cost = delta_cost;
						min_a = a; min_b = b;
					}
					
				}
				
			}
			
		}
		update_succ_and_comp(inst, min_a, min_b, succ, comp);
		--*ncomp;
		if(*ncomp != 1) add_subtour_constraint(context_pointer, environment, linear_program, inst, comp, comp[min_b], ncols); // as new component has component number of min_b
	}

	store_solution(inst, succ, sol);
	two_opt_refining_heuristic(inst, sol, 1);
	store_succ(inst, succ, sol);

	free(sol);
}

void store_succ(instance *inst, int *succ, int *sol){
	
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[sol[i]] = sol[i+1];
	}
	
}


void update_succ_and_comp(instance *inst, int min_a, int min_b, int *succ, int *comp){

	int temp = succ[min_a];
	succ[min_a] = succ[min_b];
	succ[min_b] = temp;

	int comp_a = comp[min_a];
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (comp[i] == comp_a) comp[i] = comp[min_b];
	}

}

double calc_succ_value(int *succ, instance *inst){

	int *sol = (int *) calloc(inst->nnodes+1, sizeof(int));
	double val = 0;

	store_solution(inst, succ, sol);
	for (int i = 0; i < inst->nnodes; i++)
	{
		val += dist(sol[i], sol[i+1], inst);
	}

	free(sol);
	return val;
}

/* Transforms the solution stored in the successor array format into a custom solution format for the Traveling Salesman Problem (TSP).
   Parameters:
   - inst: Pointer to the TSP instance structure.
   - succ: Array storing the successor of each node in the solution.
   - sol: Array to store the custom solution format.
*/
void store_solution(instance *inst, int *succ, int *sol){
	
	sol[0] = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		sol[i+1] = succ[sol[i]];
	}
	
}

void add_subtour_constraint(void *context_pointer, void *environment, void* linear_program, instance *inst, int *comp, int component_num, int ncols)
{

	int *index = (int *) calloc(ncols, sizeof(int));
	double *value = (double *) calloc(ncols, sizeof(double));
	double rhs = -1;
	char sense = 'L';                            // 'L' for less or equal
	// char *cname = (char *) calloc(100, sizeof(char));
	int nnz = 0;

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( comp[i] == component_num )
		{
			rhs++;
			for (int j = i+1; j < inst->nnodes; j++)
			{
				if (comp[j] == component_num)
				{
					index[nnz] = xpos(i,j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
		}
	}
	int izero = 0;
	// sprintf(cname, "SEC(%d)", *ncomp);

	if (context_pointer == NULL)
	{
		CPXENVptr env = (CPXENVptr) environment;
		CPXLPptr lp = (CPXLPptr) linear_program;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, NULL) ) print_error("CPXaddrows(): error");
	}
	else
	{
		CPXCALLBACKCONTEXTptr context = (CPXCALLBACKCONTEXTptr) context_pointer;
		if ( CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)) print_error("CPXcallbackrejectcandidate() error");
	}
	
	free(value);
	free(index);
	// free(cname);		
}

int xpos(int i, int j, instance *inst)      // to be verified                                           
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}
	

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
{    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = dist(i,j,inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

// add the degree constraints 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

}



void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
{   

#ifdef DEBUG
	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0)) > EPS ) print_error(" wrong xstar in build_sol()");
			if ( xstar[k] > 0.5 ) 
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);
#endif

	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
}