#ifndef TSP_H_  

#define TSP_H_

#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  
#include <time.h>
#include "utilities.h"

#include <cplex.h>  
#include <pthread.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-5		// 1e-5		// very small numerical tolerance 
// #define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
                                 
//data structures  

typedef struct {   
	
	//input data
	int nnodes; 	
	double *xcoord;
	double *ycoord;
	double *cost;    // c_{ij}
	// int count_of_random_points;

	// parameters 
	// double timelimit;						// overall time limit, in sec.s
	int random_seed;
	char input_file[1000];		  			// input file
	double timelimit;

	//global data
	double best_val;
	double	tstart;								
	// double zbest;							// best sol. available  
	// double tbest;							// time for the best sol. available  
	int *best_sol;						// best sol. available    
	int *tabu_list;						// if instance has a tabu list   
	int ncols;
	// double	best_lb;						// best lower bound available  
	// double *load_min;						// minimum load when leaving a node
	// double *load_max;						// maximum load when leaving a node

} instance;        

typedef struct{
	int *genes; // in our case it is tsp solution
	double fitness; // in our case it is cost of tsp and therefore smallest is the fittest
	double probability; // survival probability
						// its fitness divide by fitness of the champion in its generation
} Individual;

//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 
static inline double cost(int i, int j, instance *inst);

void plot_tsp_tour(instance *inst, int writing_to_file);
int dist(int i, int j, instance *inst);
int random_node_with_time_seed(int length);
int random_0_to_length(instance *inst, int length);
int verify_tour(instance *inst);

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

/**
 * Extra Mileage Heuristic
 * @param[in, out] inst Input instance of the tsp problem.
 * @param[in] starting_mode The starting mode for the extra mileage heuristic, 0 for largest distance nodes, 1 for two different random nodes,
 * 2 is to try all couple of nodes in the convexHull which are different and choose the best one.
 * @return Returns 0 if the extra mileage heuristic is successfully applied without encountering any errors.
**/
int extra_mileage_heuristic(instance *inst, int starting_mode);
void calculate_extra_mileage_heuristics(instance *inst, int *nodes_hierarchy);
int * extra_mileage_step(instance *inst, int *uncovered_nodes, int current_length, int *nodes_hierarchy);
double delta_cost_extra_mileage(instance *inst, int i, int j, int h);
void calculate_best_val(instance *inst);

int two_opt_refining_heuristic(instance *inst, int *tsp_sol, int is_instance);
double delta_cost_two_opt(int a, int b, instance *inst, int *tsp_sol);
int update_tour(int i, int j, instance *inst, int *tsp_sol, int best_val_update_switch);

/**
 * Tabu Search
 * @param[in, out] inst Input instance of the tsp problem. 
 * @param[in] tenure_mode The tenure mode of the tabu search, 0 for reactive step tenure, 1 for reactive line tenure and 2 is for random tenure.
 * @return Returns 0 if the tabu search is successfully applied without encountering any errors.
**/
int tabu_search(instance *inst, int tenure_mode);
int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode);
int update_tour_and_tabu_list(int a, int b, instance *inst);

/**
 * Variable Neighborhood Search
 * @param[in, out] inst Input instance of the tsp problem. 
 * @param[in] kick_neighborhood The neighborhood type of the kick, 3 for 3-opt, 5 for 5-opt and etc.
 * @return Returns 0 if the variable neighborhood search is successfully applied without encountering any errors.
 * @note The 1-OPT kick doesn't update the tour and also 2-OPT kick as 2-OPT refining is used.
 * So, try to use 3-OPT kick or more.
**/
int variable_neighborhood_search(instance *inst, int kick_neighborhood);
int n_opt_kick(instance *inst, int n);
int new_tour_from_break_positions(instance *inst, int *break_positions, int arr_size);
int copy_segment(instance *inst, int *old_solution, int starting_pos, int ending_pos, int into_pos);
int copy_segment_in_reverse_order(instance *inst, int *old_solution, int starting_pos, int ending_pos, int into_pos);

/**
 * Simulated Annealing
 * @param[in, out] inst Input instance of the tsp problem.
 * @return Returns 0 if the simulated annealing is successfully applied without encountering any errors.
**/
int simulated_annealing(instance *inst);
double metropolis_formula(double delta_cost, double Temprature, int scaler);
int annealing_process(instance *inst, int scaler);
double average_delta_cost_between_two_edges(instance *inst);

/**
 * Genetic Algorithm
 * @param[in, out] inst Input instance of the tsp problem.
 * @param[in] repair Switch for the genetic algorithm to use repair or not. 0 is to punish their fitness with penalty without repairing,
 * 1 is to repair the bad genes, 2 is to repair the bad genes but also use two opt refinining while repairing the genes.
 * @param[in] cutting_type The cutting type of the crossover, 0 for cutting from the middle, 1 for cutting from the random position
 * and 2 is for cutting from the one forth of the tour.
 * @return Returns 0 if the genetic algorithm is successfully applied without encountering any errors.
**/
int genetic_algorithm(instance *inst, int repair, int cutting_type);
int initialize_population_randomly(instance *inst, int population_size, Individual *population);
void initialize_individual(instance *inst, Individual *individual);
void print_population(instance *inst, int population_size, Individual *population);
void free_population(Individual *population, int population_size);
void free_Individual(Individual *individual);
Individual * find_champion_individual(instance *inst, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size);
void calculate_individual_fitness(instance *inst, Individual *individual);
void calculate_population_fitness(instance *inst, Individual *population, int population_size);
void avoid_bad_genes(instance *inst, Individual *children, int children_size);
void print_parents_and_child(instance *inst, Individual *parent1, Individual *parent2, Individual *child);
void crossover(instance *inst, Individual *population, int population_size, Individual *children, int children_size, int cutting_type );
void mutate_population(instance *inst, Individual *population, int population_size, Individual *mutations, int mutants_size);
void alteration_of_genes(instance *inst, int *genes);
Individual * survival_probabilities_of_generation(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion);
Individual * elitism(instance *inst, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion);
int compare_individuals(const void *a, const void *b);
void repair_bad_genes(instance *inst, Individual *children, int children_size, int apply_two_opt);
void eliminate_multiple_visits(instance *inst, Individual *indiviual);
void repair_extra_mileage(instance *inst, Individual *indiviual);

int TSPopt(instance *inst, int model_type);
int xpos(int i, int j, instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
void add_subtour_constraint(void *context_pointer, void *environment, void* linear_program, instance *inst, int *comp, int component_num, int ncols);
void store_solution(instance *inst, int *succ, int *sol);
double calc_incumbent_value(int *succ, instance *inst);
void patching_heuristic(CPXENVptr env, CPXLPptr lp, int ncols, instance *inst, int *succ, int *comp, int *ncomp);
double delta_cost_patching(int a, int b, instance *inst, int *succ);
void update_succ_and_comp(instance *inst, int min_a, int min_b, int *succ, int *comp);
void store_succ(instance *inst, int *succ, int *sol);
int benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp);
int branch_and_cut(instance *inst, CPXENVptr env, CPXLPptr lp);
static int CPXPUBLIC my_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);

#endif   /* TSP_H_ */ 
 