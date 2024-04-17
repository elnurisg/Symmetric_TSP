#include "utilities.h"
#include "convex_hull.h"

///////////////////////////////////Greedy Heuristic/////////////////////////////////////////

/**
 * Greeady Heuristic
 * @param[in, out] inst Input instance of the tsp problem.
 * @param[in] starting_mode The starting mode of the greedy heuristic, 0 for starting from position zero, 
 * 1 for starting from the random position, 2 is to try all the starting positions and choose the best one.
 * @param[in] grasp Switch for the greedy heuristic to use grasp or not. 0 is for standard greedy heuristic, 1 is for grasp.
 * When grasp is activated, it chooses randomly the minimum 40% of the time, the second minimum 30% of the time and the third minimum 30% of the time.
 * @return Returns 0 if the greedy heuristic is successfully applied without encountering any errors.
**/
int greedy_heuristic(instance *inst, int starting_mode, int grasp);
void calculate_greedy_steps(instance *inst, int starting_node_pos, int grasp);
int greedy_step(instance *inst, int current_node, int *uncovered_nodes, int current_length, int grasp);

///////////////////////////////////Extra Mileage Heuristic/////////////////////////////////////////

/**
 * Extra Mileage Heuristic
 * @param[in, out] inst Input instance of the tsp problem.
 * @param[in] starting_mode The starting mode for the extra mileage heuristic, 0 for largest distance nodes, 1 for two different random nodes,
 * 2 is to try all couple of nodes in the convexHull which are different and choose the best one.
 * @return Returns 0 if the extra mileage heuristic is successfully applied without encountering any errors.
**/
int extra_mileage_heuristic(instance *inst, int starting_mode);
void calculate_extra_mileage_heuristics(instance *inst, int *nodes_hierarchy);
void extra_mileage_step(instance *inst, int *uncovered_nodes, int current_length, int *nodes_hierarchy, int *best_values);
double delta_cost_extra_mileage(instance *inst, int i, int j, int h);
void calculate_best_val(instance *inst);
int maximum_cost_pos(instance *inst);
int * find_nodes(instance *inst, int pos);

///////////////////////////////////Two Opt Refining Heuristic/////////////////////////////////////////

int two_opt_refining_heuristic(instance *inst, int *tsp_sol, int is_instance);
double delta_cost_two_opt(int a, int b, instance *inst, int *tsp_sol);
int update_tour(int i, int j, instance *inst, int *tsp_sol, int best_val_update_switch);
