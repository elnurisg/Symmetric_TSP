#ifndef TSP_H_  

#define TSP_H_

#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  
#include <time.h>

#include "utilities.h"
#include "heuristics.h"
#include "metaheuristics.h"

#include <cplex.h>  
#include <pthread.h>  


//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 


int TSPopt(instance *inst, int model_type);
int xpos(int i, int j, instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
void add_subtour_constraint(void *context_pointer, void *environment, void* linear_program, instance *inst, int *comp, int component_num, int ncols);
void store_solution(instance *inst, int *succ, int *sol);
double calc_succ_value(int *succ, instance *inst);
void patching_heuristic(void *context_pointer, void *environment, void* linear_program, int ncols, instance *inst, int *succ, int *comp, int *ncomp);
double delta_cost_patching(int a, int b, instance *inst, int *succ);
void update_succ_and_comp(instance *inst, int min_a, int min_b, int *succ, int *comp);
void store_succ(instance *inst, int *succ, int *sol);
int benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp);
int branch_and_cut(instance *inst, CPXENVptr env, CPXLPptr lp);
static int CPXPUBLIC my_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);
int heur_sol_to_mipstart(CPXENVptr env, CPXLPptr lp, instance *inst);
void post_heuristic(CPXCALLBACKCONTEXTptr context, instance *inst, int *succ, double succ_value);

#endif   /* TSP_H_ */ 
 