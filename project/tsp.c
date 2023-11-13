#include "tsp.h"
#include <time.h>

void print_error(const char *err);
double second();       
void debug(const char *err);       
double random01() { return ((double) rand() / RAND_MAX); } // return a random value in range 0.0-1.0
int time_limit_expired(instance *inst);

inline double cost(int i, int j, instance *inst)
{
	return inst->cost[i*inst->nnodes + j];
}

int dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	int dis = sqrt(dx*dx+dy*dy) + 0.5; // nearest integer 
	return dis;
}        

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

// deterministic nearest_neighbor algorithm
// which starts from node 0 or random or try all

int random_node(int length){
	// double random_position = ((double)rand() / (double)RAND_MAX * (double)length);
	srand((unsigned int)time(NULL));
	int random_position = rand() % length;

	if (random_position<0 && random_position>=length)
	{
		printf("\n Error! Position exceeds the size of array");
		return -1;
	}
	
	return random_position;
}


// when the starting_node is known, it applies the nearest neighbor algorithm
void calculate_greedy_steps(instance *inst, int starting_node_pos, int grasp){
	int uncovered_nodes[inst->nnodes];
	int current_length = inst->nnodes;

	inst->best_sol = (double *) calloc(inst->nnodes, sizeof(double));

	inst->best_val = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		uncovered_nodes[i] = i; 
	}

	int best_node_pos = starting_node_pos;
	int best_node = uncovered_nodes[starting_node_pos]; // the first best_node, 
	//node of the best_sol is the starting position
	inst->best_sol[0] = best_node;
	int current_node = best_node;

	// starting_node is not uncovered anymore,
	// therefore, we substitute it with the last element of uncovered_nodes
	current_length--;
	uncovered_nodes[starting_node_pos] = uncovered_nodes[current_length];

	while (current_length != 0) // each step we consider two nodes 
	{ //so that we will have n-1 steps
		// update node with the next one which is the min dist among the uncovered_nodes
		// and update best_val and best_sol
		best_node_pos = greedy_step(inst, current_node, uncovered_nodes, current_length, grasp);
		best_node = uncovered_nodes[best_node_pos];
		inst->best_val += inst->cost[current_node*(inst->nnodes-1) + best_node];
		inst->best_sol[(inst->nnodes) - current_length]= best_node;
		current_node = best_node;
		current_length--;
		uncovered_nodes[best_node_pos] = uncovered_nodes[current_length];

	}
	
}

// get the current node as input, find the closest node and then 
// add the distance to best_val
// add the best_val into best_sol
// int greedy_step(instance *inst, int current_node, int *uncovered_nodes, int current_length){

// 	    if (current_length == 0) {
//         printf("Error: Empty array\n");
//         return -1;
//     	}

// 		double min = inst->cost[current_node*(inst->nnodes-1) + uncovered_nodes[0]];
// 		double d2 = min;
// 		int best_node_pos = -1;
// 		for (int i = 0; i < current_length; i++)
// 		{
// 			d2 = inst->cost[current_node*(inst->nnodes-1) + uncovered_nodes[i]];
// 			if (min > d2)
// 			{
// 				min = d2;
// 				best_node_pos = i;
// 			}
// 		}
// 		return best_node_pos;
// }

int greedy_step(instance *inst, int current_node, int *uncovered_nodes, int current_length, int grasp) {
    if (current_length == 0) {
        printf("Error: Empty array\n");
        return -1;
    }

    double min = inst->cost[current_node * (inst->nnodes - 1) + uncovered_nodes[0]];
    double d2 = min;
    int best_node_pos = -1;

    if (grasp) {
        // GRASP mode: choose 4 times the min and 3 times the second minimum
		// 3 times the third minimum randomly
        for (int i = 0; i < current_length; i++) {
            d2 = inst->cost[current_node * (inst->nnodes - 1) + uncovered_nodes[i]];
            if (min > d2) {
                min = d2;
                best_node_pos = i;
            }
        }

        // Randomly decide whether to choose the second minimum
        double random = random01();
        if (random < 0.4) {
            return best_node_pos;
        } else {
            int second_best_node_pos = -1;
            double second_min = min;

            for (int i = 0; i < current_length; i++) {
                if (i != best_node_pos) {
                    d2 = inst->cost[current_node * (inst->nnodes - 1) + uncovered_nodes[i]];
                    if (second_min > d2) {
                        second_min = d2;
                        second_best_node_pos = i;
                    }
                }
            }
			if (random < 0.7) {
				printf("second_best_node\n");
				return second_best_node_pos;
			} else {
				int third_best_node_pos = -1;
				double third_min = min;

				for (int i = 0; i < current_length; i++) {
					if (i != best_node_pos && i!= second_best_node_pos) {
						d2 = inst->cost[current_node * (inst->nnodes - 1) + uncovered_nodes[i]];
						if (third_min > d2) {
							third_min = d2;
							third_best_node_pos = i;
						}
					}
				}
				printf("third_best_node\n");
				return third_best_node_pos;
			}

        }
    } else {
        // Standard mode: choose the minimum
        for (int i = 0; i < current_length; i++) {
            d2 = inst->cost[current_node * (inst->nnodes - 1) + uncovered_nodes[i]];
            if (min > d2) {
                min = d2;
                best_node_pos = i;
            }
        }
        return best_node_pos;
    }
}


// we initialize the starting_node parameter 
// and then call calculate_greedy_steps algorithm
int greedy_heuristic(instance *inst, int starting_node_pos, int grasp)
{
	int starting_pos;
	int random_node_pos;
	int min_cost;
	int best_starting_node;

	switch (starting_node_pos)
	{
	case 0:  // node 0
		starting_pos = 0;
		calculate_greedy_steps(inst, starting_pos, grasp);
		break;
	case 1:  // random
		random_node_pos = random_node(inst->nnodes);
		printf("\n random_nod_pos is %d \n", random_node_pos);
		calculate_greedy_steps(inst, random_node_pos, grasp);
		break;
	case 2: // try all
		calculate_greedy_steps(inst, 0, grasp);
		min_cost = inst->best_val;
		best_starting_node = 0;

		for (int i = 1; i < inst->nnodes; i++)
		{
			calculate_greedy_steps(inst, i, grasp);
			if (min_cost > inst->best_val)
			{
				min_cost = inst->best_val;
				best_starting_node = i;
			}
		}		 
		printf("\n The best starting node is %d \n", best_starting_node);
		calculate_greedy_steps(inst, best_starting_node, grasp);
		break;
	default:
            printf("Error! \n starting_node is not in correct format\n");
	}

	return 0;


}