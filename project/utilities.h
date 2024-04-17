#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
// #include "convex_hull.h"


#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)
//////////////////////////// add VERBOSE also for cplex output
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

void plot_tsp_tour(instance *inst, int writing_to_file);
int dist(int i, int j, instance *inst);
int random_node_with_time_seed(int length);
int random_0_to_length(instance *inst, int length);
int random_0_to_length_but_different_than_previous(instance *inst, int length, int previous_random_value);
int verify_tour(instance *inst, int *tour);
double random01();
double cost(int i, int j, instance *inst);
void print_error(const char *err);  

int * copy_to_new_array(int *arr, int size);
void copy_array(int *arr, int size, int *arr_to_copy);
void add_to_array(int starting_pos, int new_element, int *arr, int size);
void print_array(int *arr, int size);
int compare(const void* num1, const void* num2);
void swap(int *a, int *b);
void shuffle_array_for_kick(int arr[], int n);
void shuffle_tsp_sol(int *arr, int size);
void combine_two_tours_from_pos(int *arr1, int *arr2, int size, int position, int *result_arr);
void remove_from_array(int pos, int *arr, int size, int *new_arr);
void write_cost_to_file(double cost, const char *filename, int append);


#endif // UTILITIES_H