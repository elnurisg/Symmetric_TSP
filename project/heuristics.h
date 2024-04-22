#include "utilities.h"
#include "convex_hull.h"


///////////////////////////////////Greedy Heuristic/////////////////////////////////////////

/**
 * @brief Greedy Heuristic: Constructs a solution for the Traveling Salesman Problem (TSP) using a greedy approach.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] starting_mode The starting mode of the greedy heuristic.
 *                          - 0: Start from position zero.
 *                          - 1: Start from a random position.
 *                          - 2: Try all starting positions and choose the best one.
 * @param[in] grasp Switch for the greedy heuristic to use GRASP (Greedy Randomized Adaptive Search Procedure) or not:
 *                 - 0: Standard greedy heuristic.
 *                 - 1: GRASP activated. Chooses randomly the minimum 40% of the time, the second minimum 30% of the time,
 *                      and the third minimum 30% of the time.
 * @return Returns 0 if the greedy heuristic is successfully applied without encountering any errors.
**/
int greedy_heuristic(instance *inst, int starting_mode, int grasp);

/**
 * @brief Calculates the steps of the greedy heuristic based on the starting node and GRASP setting.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] starting_node_pos The position of the starting node.
 * @param[in] grasp Switch for the greedy heuristic to use GRASP (Greedy Randomized Adaptive Search Procedure) or not.
**/
void calculate_greedy_steps(instance *inst, int starting_node_pos, int grasp);

/**
 * @brief Performs a single step of the greedy heuristic.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] current_node The current node being considered.
 * @param[in] uncovered_nodes Array of indices representing the uncovered nodes.
 * @param[in] current_length The current length of uncovered nodes.
 * @param[in] grasp Switch for the greedy heuristic to use GRASP (Greedy Randomized Adaptive Search Procedure) or not.
 * @return Returns the index of the next node to be added to the tour.
**/
int greedy_step(instance *inst, int current_node, int *uncovered_nodes, int current_length, int grasp);


///////////////////////////////////Extra Mileage Heuristic/////////////////////////////////////////

/**
 * @brief Extra Mileage Heuristic: A method to construct a solution for the Traveling Salesman Problem (TSP) by considering extra mileage.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] starting_mode The starting mode for the extra mileage heuristic.
 *                          - 0: Start with nodes with the largest distance between them.
 *                          - 1: Start with two different random nodes.
 *                          - 2: Try all pairs of nodes in the convex hull that are different and choose the best one.
 * @return Returns 0 if the extra mileage heuristic is successfully applied without encountering any errors.
**/
int extra_mileage_heuristic(instance *inst, int starting_mode);

/**
 * @brief Calculates the steps of the extra mileage heuristic based on the starting mode and store the tour in nodes_hierarchy.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in, out] nodes_hierarchy An incomplete tour comprising only the starting node couple.
**/
void calculate_extra_mileage_heuristics(instance *inst, int *nodes_hierarchy);

/**
 * @brief Performs a single step of the extra mileage heuristic and add the new node in nodes_hierarchy.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] uncovered_nodes Array of indices representing the uncovered nodes.
 * @param[in] uncovered_length The number of uncovered nodes.
 * @param[in, out] nodes_hierarchy An incomplete tour at the given step.
**/
void extra_mileage_step(instance *inst, int *uncovered_nodes, int uncovered_length, int *nodes_hierarchy);

/**
 * @brief Calculates the change in cost by swapping two nodes in the extra mileage heuristic.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] i Index of the first node.
 * @param[in] j Index of the second node.
 * @param[in] h Index of the node to be replaced.
 * @return Returns the change in cost.
**/
double delta_cost_extra_mileage(instance *inst, int i, int j, int h);

/**
 * @brief Finds the position with the maximum cost in the extra mileage heuristic.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @return Returns the position with the maximum cost.
**/
int maximum_cost_pos(instance *inst);

/**
 * @brief Finds the nodes at a given position in the extra mileage heuristic.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] pos Position to find the nodes.
 * @return Returns an incomplete tour array containing only the found node couple.
**/
int * find_nodes(instance *inst, int pos);


///////////////////////////////////Two Opt Refining Heuristic/////////////////////////////////////////

/**
 * @brief Two-Opt Refining Heuristic: A method to refine a solution for the Traveling Salesman Problem (TSP) using the Two-Opt technique.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] tsp_sol Initial solution to be refined by the Two-Opt heuristic.
 * @param[in] is_instance Flag indicating whether the input solution is part of the instance structure or an external solution.
 *                        - 0: The solution is part of the instance structure.
 *                        - 1: The solution is external to the instance structure.
 * @return Returns 0 if the Two-Opt refining heuristic is successfully applied without encountering any errors.
**/
int two_opt_refining_heuristic(instance *inst, int *tsp_sol, int is_instance);

/**
 * @brief Calculates the change in cost by performing a Two-Opt move.
 *
 * @param[in] a Index of the first node.
 * @param[in] b Index of the second node.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] tsp_sol Current solution of the TSP.
 * @return Returns the change in cost.
**/
double delta_cost_two_opt(int a, int b, instance *inst, int *tsp_sol);

/**
 * @brief Updates the tour by performing a Two-Opt move.
 *
 * @param[in] i Index of the first node.
 * @param[in] j Index of the second node.
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] tsp_sol Current solution of the TSP.
 * @return Returns 1 if the tour is updated, 0 otherwise.
**/
int update_tour(int i, int j, int *tsp_sol);
