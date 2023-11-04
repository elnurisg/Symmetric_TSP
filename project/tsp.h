#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  

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
	double *uncovered_nodes;
	double best_val;
	// double	tstart;								
	// double zbest;							// best sol. available  
	// double tbest;							// time for the best sol. available  
	double *best_sol;						// best sol. available    
	// double	best_lb;						// best lower bound available  
	// double *load_min;						// minimum load when leaving a node
	// double *load_max;						// maximum load when leaving a node

} instance;        

//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 
int dist(int i, int j, instance *inst);
int greedy_heuristic(instance *inst);
#endif   /* TSP_H_ */ 
