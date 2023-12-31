#ifndef TSP_H_  

#define TSP_H_

#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  
#include <time.h>
#include "utilities.h"

// #include <cplex.h>  
#include <pthread.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
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
	// double	best_lb;						// best lower bound available  
	// double *load_min;						// minimum load when leaving a node
	// double *load_max;						// maximum load when leaving a node

} instance;        

//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 
static inline double cost(int i, int j, instance *inst);

int dist(int i, int j, instance *inst);
int random_node_with_time_seed(int length);
int random_0_to_length(instance *inst, int length);
int verify_tour(instance *inst);

int greedy_heuristic(instance *inst, int starting_mode, int grasp);
void calculate_greedy_steps(instance *inst, int starting_node_pos, int grasp);
int greedy_step(instance *inst, int current_node, int *uncovered_nodes, int current_length, int grasp);

int extra_mileage_heuristic(instance *inst, int starting_mode);
void calculate_extra_mileage_heuristics(instance *inst, int *nodes_hierarchy);
int * extra_mileage_step(instance *inst, int *uncovered_nodes, int current_length, int *nodes_hierarchy);
double delta_cost_extra_mileage(instance *inst, int i, int j, int h);
void calculate_best_val(instance *inst);

int two_opt_refining_heuristic(instance *inst);
double delta_cost_two_opt(int i, int j, instance *inst);
int update_tour(int i, int j, instance *inst);

int tabu_search(instance *inst, int tenure_mode);
int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode);
int update_tour_and_tabu_list(int a, int b, instance *inst);

int variable_neighborhood_search(instance *inst, int kick_neighborhood);
int n_opt_kick(instance *inst, int n);
int new_tour_from_break_positions(instance *inst, int *break_positions, int arr_size);
int copy_segment(instance *inst, int *old_solution, int starting_pos, int ending_pos, int into_pos);
int copy_segment_in_reverse_order(instance *inst, int *old_solution, int starting_pos, int ending_pos, int into_pos);

int simulated_annealing(instance *inst);
double metropolis_formula(double delta_cost, double Temprature, int scaler);
int annealing_process(instance *inst, int scaler);
double average_delta_cost_between_two_edges(instance *inst);
#endif   /* TSP_H_ */ 
