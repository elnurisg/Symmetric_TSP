#include "../include/tsp.h"

/////////////////////////////////// Hard-fixing /////////////////////////////////////////

/**
 * @brief Applies hard fixing heuristic to the TSP instance to fix a portion of the tour with a given probability.
 *
 * @param[in, out] inst Pointer to the instance of the TSP problem.
 * @param[in] pfix Probability of fixing. For example, if pfix is 0.1, it fixes ~10% of the tour.
 * 
 * @return Returns 0 if the hard fixing heuristic is successfully applied.
**/
int hard_fixing(instance *inst, double pfix);

/**
 * @brief Fixes a portion of the tour in the CPLEX model with a given probability.
 *
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] pfix Probability of fixing.
**/
void fixing(CPXENVptr env, CPXLPptr lp, instance *inst, double pfix);

/**
 * @brief Unfixes all the previously fixed variables in the CPLEX model.
 *
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
 * @param[in] inst Input instance of the TSP problem.
**/
void unfixing(CPXENVptr env, CPXLPptr lp, instance *inst);


/////////////////////////////////// Local branching /////////////////////////////////////////

/**
 * @brief Applies the local branching metaheuristic to the TSP instance with a specified local branching constraint (k).
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] k Local branching constraint: the maximum number of cities allowed to be excluded from the tour.
 * @return Returns 0 if the local branching heuristic is successfully applied.
**/
int local_branching(instance *inst, int k);

/**
 * @brief Adds constraints to the CPLEX model to enforce a local branching heuristic with the specified constraint (k).
 *
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] k Local branching constraint: the maximum number of cities allowed to be excluded from the tour.
**/
void add_local_branching_constraints(CPXENVptr env, CPXLPptr lp, instance *inst, int k);

/**
 * @brief Removes the last set of constraints added to the CPLEX model.
 *
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
**/
void remove_last_constraints(CPXENVptr env, CPXLPptr lp);