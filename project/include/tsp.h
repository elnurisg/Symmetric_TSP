#ifndef TSP_H_  

#define TSP_H_

#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  
#include <time.h>

#include "../include/utilities.h"
#include "../include/heuristics.h"
#include "../include/metaheuristics.h"

#include <cplex.h>  
#include <pthread.h>  

/**
 * @brief Solves the Traveling Salesman Problem (TSP) instance using optimization techniques with the help of Cplex callbacks.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] model_type Integer indicating the type of optimization model to use:
 *                      - 0: Benders Decomposition.
 *                      - 1: Branch and Cut.
 * 
 * @return Returns 0 if it is successfully applied without encountering any errors.
**/
int TSPopt(instance *inst, int model_type);

/**
 * @brief Calculates and returns the position of the variable x(i,j) within the CPLEX model.
 *
 * @param[in] i Starting node index.
 * @param[in] j Ending node index.
 * @param[in] inst Input instance of the TSP problem.
 * 
 * @return Position of the variable x(i,j) in the CPLEX model.
**/
int xpos(int i, int j, instance *inst);

/**
 * @brief Constructs the optimization model for the Traveling Salesman Problem (TSP) instance using CPLEX with degree constraints.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
**/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);

/**
 * @brief Converts the CPLEX solution format to successor style and determines the connected components of nodes.
 *
 * @param[in] xstar Pointer to the CPLEX solution.
 * @param[in] inst Input instance of the TSP problem.
 * @param[out] succ Array to store the tour in successor style.
 * @param[out] comp Array to store the connected component of each node.
 * @param[out] ncomp Pointer to store the number of connected components.
**/
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

/**
 * @brief Adds a subtour elimination constraint to the linear program for the specified connected component.
 * The constraint aims to eliminate subtours in the TSP solution by restricting the solution space for the given component.
 *
 * @param[in] context_pointer Pointer (CPXCALLBACKCONTEXTptr) CPLEX callback context. Set to NULL if not used.
 * @param[in] environment Pointer (CPXENVptr) to the CPLEX environment. Set to NULL if not used.
 * @param[in] linear_program Pointer (CPXLPptr) to the CPLEX linear program. Set to NULL if not used.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] comp Array representing the connected components of nodes.
 * @param[in] component_num Index of the connected component for which the subtour elimination constraint is added.
 * @param[in] ncols Number of columns (variables) in the linear program.
 * 
 * @note If context_pointer is provided, CPXcallbackrejectcandidate is used; otherwise, CPXaddrows is used.
**/
void add_subtour_constraint(void *context_pointer, void *environment, void* linear_program, instance *inst, int *comp, int component_num, int ncols);

/**
 * @brief Converts a successor-style tour representation to a solution-style tour representation.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] succ Array representing the successor-style tour.
 * @param[out] sol Array to store the solution-style tour.
 * 
 * @note  successor-style tour representation: each node points to its successor in the tour. solution-style tour representation: connected nodes are stored consecutively in an array.
**/
void store_solution(instance *inst, int *succ, int *sol);

/**
 * @brief Calculates the total cost of a tour represented in successor-style.
 *
 * @param[in] succ Array representing the successor-style tour.
 * @param[in] inst Input instance of the TSP problem.
 * 
 * @return Total cost of the tour.
 * 
 * @note  successor-style tour representation: each node points to its successor in the tour.
**/
double calc_succ_value(int *succ, instance *inst);

/**
 * @brief Applies the patching heuristic to connect subtours until there is only one tour, and adds subtour elimination constraints for all subtours (original and newly formed subtours).
 *
 * @param[in] context_pointer Pointer (CPXCALLBACKCONTEXTptr) CPLEX callback context. Set to NULL if not used.
 * @param[in] environment Pointer (CPXENVptr) to the CPLEX environment. Set to NULL if not used.
 * @param[in] linear_program Pointer (CPXLPptr) to the CPLEX linear program. Set to NULL if not used.
 * @param[in] ncols Number of columns (variables) in the CPLEX linear program.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] succ Array representing the successor-style tour.
 * @param[in, out] comp Array to store the connected component of each node.
 * @param[in, out] ncomp Pointer to the number of connected components.
**/
void patching_heuristic(void *context_pointer, void *environment, void* linear_program, int ncols, instance *inst, int *succ, int *comp, int *ncomp);

/**
 * @brief Calculates the change in total cost (delta cost) by connecting two nodes (a and b) in the tour represented by the successor-style tour representation (succ).
 *
 * @param[in] a Index of the first node to connect.
 * @param[in] b Index of the second node to connect.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] succ Array representing the successor-style tour.
 * 
 * @return Change in total cost (delta cost) by connecting the two nodes.
**/
double delta_cost_patching(int a, int b, instance *inst, int *succ);

/**
 * @brief Updates the successor array and connected component array after connecting two subtours (from given nodes) in the tour represented by the successor-style tour representation.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] a Index of the node from the first subtour to connect.
 * @param[in] b Index of the node from the second subtour to connect.
 * @param[in, out] succ Array representing the successor-style tour.
 * @param[in, out] comp Array representing the connected components of nodes.
**/
void update_succ_and_comp(instance *inst, int a, int b, int *succ, int *comp);

/**
 * @brief Stores the solution-style tour as a successor-style tour representation.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] succ Array representing the successor-style tour.
 * @param[in] sol Array representing the solution-style tour.
 * 
 * @note  successor-style tour representation: each node points to its successor in the tour. solution-style tour representation: connected nodes are stored consecutively in an array.
**/
void store_succ(instance *inst, int *succ, int *sol);

/**
 * @brief Implements the Benders decomposition method to solve the TSP instance using Cplex.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
 * 
 * @return Returns 0 if the Benders decomposition method is successfully applied without encountering any errors.
**/
int benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp);

/**
 * @brief Implements the branch-and-cut method to solve the TSP instance using Cplex.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
 * 
 * @return Returns 0 if the branch-and-cut method is successfully applied without encountering any errors.
**/
int branch_and_cut(instance *inst, CPXENVptr env, CPXLPptr lp);

/**
 * @brief Callback function used to handle the addition of cutting plane constraints during the branch-and-cut process in CPLEX.
 *
 * This function is invoked by CPLEX during the branch-and-cut process to handle candidate solutions. It applies a patching heuristic to address subtours and posts improved solutions if found. 
 *
 * @param[in] context Pointer (CPXCALLBACKCONTEXTptr) CPLEX callback context.
 * @param[in] contextid Identifier for the context type.
 * @param[in] userhandle Pointer to the user-defined data structure containing problem-specific information.
 * 
 * @return Returns 0 to indicate successful execution.
**/
static int CPXPUBLIC my_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);

/**
 * @brief Adds the current best solution to the MIP start in the CPLEX environment.
 *
 * @param[in] env Pointer (CPXENVptr) to the CPLEX environment.
 * @param[in] lp Pointer (CPXLPptr) to the CPLEX linear program.
 * @param[in] inst Input instance of the TSP problem.
 * 
 * @return Returns 0 if the current best solution is successfully added to the MIP start.
**/
int best_sol_to_mipstart(CPXENVptr env, CPXLPptr lp, instance *inst);

/**
 * @brief Posts a successor-style solution to the CPLEX callback context.
 *
 * @param[in] context Pointer (CPXCALLBACKCONTEXTptr) CPLEX callback context.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] succ Array representing the successor-style tour.
 * @param[in] succ_value Total cost of the successor-style tour.
 * 
 * @note  successor-style tour representation: each node points to its successor in the tour.
**/
void post_heuristic(CPXCALLBACKCONTEXTptr context, instance *inst, int *succ, double succ_value);

#endif   /* TSP_H_ */ 
 