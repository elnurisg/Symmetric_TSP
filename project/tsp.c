#include "tsp.h"
#include <time.h>
#include "convex_hull.h"

void print_error(const char *err);
double second();       
void debug(const char *err);       
double random01() { return ((double) rand() / RAND_MAX); } // return a random value in range 0.0-1.0
int time_limit_expired(instance *inst);

static inline double cost(int i, int j, instance *inst)
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

int random_node_with_time_seed(int length){
	srand((unsigned int)time(NULL));
	int random_position = rand() % length;

	if (random_position<0 && random_position>=length)
	{
		printf("\n Error! Position exceeds the size of array");
		return -1;
	}
	
	return random_position;
}

int random_0_to_length(instance *inst, int length){

	if (inst->random_seed !=0){
		srand(2635623+abs(inst->random_seed));
		for (size_t i = 0; i < 1000; i++) random();
	}

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
	int *uncovered_nodes = (int *) calloc(inst->nnodes, sizeof(int));//[inst->nnodes];
	int current_length = inst->nnodes;

	inst->best_sol = (int *) calloc(inst->nnodes+1, sizeof(int)); 
	// last one to close the tour
	inst->best_val = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		uncovered_nodes[i] = i; 
	}

	int best_node_pos = starting_node_pos;
	int best_node = uncovered_nodes[starting_node_pos]; // the first best_node, 
	//first node of the best_sol is the starting position
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
		

		inst->best_val += inst->cost[current_node*inst->nnodes + best_node];
		inst->best_sol[(inst->nnodes) - current_length]= best_node;
		current_node = best_node;

		current_length--;
		uncovered_nodes[best_node_pos] = uncovered_nodes[current_length];

	}
		// close the tour
	inst->best_val += inst->cost[current_node*inst->nnodes + starting_node_pos];
	inst->best_sol[inst->nnodes]= starting_node_pos;

}

int greedy_step(instance *inst, int current_node, int *uncovered_nodes, int current_length, int grasp) {
    if (current_length == 0) {
        printf("Error: Empty array\n");
        return -1;
    }

    double min = __DBL_MAX__;
    double d2 = min;
    int best_node_pos = -1;

    int second_best_node_pos = -1;
    double second_min = __DBL_MAX__;

	int third_best_node_pos = -1;
	double third_min = __DBL_MAX__;

    if (grasp) {
        // GRASP mode: choose 4 times the min and 3 times the second minimum
		// 3 times the third minimum randomly
        for (int i = 0; i < current_length; i++) {
            d2 = inst->cost[current_node * inst->nnodes + uncovered_nodes[i]];
            if (min > d2) {
                min = d2;
                best_node_pos = i;
            }
        }

        // Randomly decide whether to choose the second minimum
        double random = random01();
        if (random < 0.4 || current_length < 3) { 
            return best_node_pos;
        } else {

            for (int i = 0; i < current_length; i++) {
                if (i != best_node_pos) {
                    d2 = inst->cost[current_node * inst->nnodes + uncovered_nodes[i]];
                    if (second_min > d2) {
                        second_min = d2;
                        second_best_node_pos = i;
                    }
                }
            }
			// return second_best_node_pos;
			if (random < 0.7) {
				return second_best_node_pos;
			} else {

				for (int i = 0; i < current_length; i++) {
					if (i != best_node_pos && i!= second_best_node_pos) {
						d2 = inst->cost[current_node * inst->nnodes + uncovered_nodes[i]];
						if (third_min > d2) {
							third_min = d2;
							third_best_node_pos = i;
						}
					}
				}
				return third_best_node_pos;
			}

        }
    } else {
        // Standard mode: choose the minimum
        for (int i = 0; i < current_length; i++) {
            d2 = inst->cost[current_node * inst->nnodes + uncovered_nodes[i]];
            if (min > d2) {
                min = d2;
                best_node_pos = i;
            }
        }
        return best_node_pos;
    }
}

// deterministic nearest_neighbor algorithm
// which starts from node 0 or random or try all

// we initialize the starting_mode parameter 
// and then call calculate_greedy_steps algorithm
int greedy_heuristic(instance *inst, int starting_mode, int grasp)
{
	printf("\n_________________________________________________________\nGreedy Heuristic:\n");
	int starting_pos;
	int random_node_with_time_seed_pos;
	int min_cost;
	int best_starting_node;

	switch (starting_mode)
	{
	case 0:  // node 0
		starting_pos = 0;
		calculate_greedy_steps(inst, starting_pos, grasp);
		break;
	case 1:  // random
		random_node_with_time_seed_pos = random_node_with_time_seed(inst->nnodes);
		printf("\n random_nod_pos is %d \n", random_node_with_time_seed_pos);
		calculate_greedy_steps(inst, random_node_with_time_seed_pos, grasp);
		break;
	case 2: // try all
		calculate_greedy_steps(inst, 0, grasp);
		min_cost = inst->best_val;
		best_starting_node = 0;

		for (int i = 0; i < inst->nnodes; i++)
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

// finding the position of edge with maximum cost between two nodes
int maximum_cost_pos(instance *inst){
	int max = inst->cost[1];
	int max_pos = 1;

	for (int i = 0; i < inst->nnodes*inst->nnodes; i++) 
	{ 
		if (max < inst->cost[i])
		{
			max = inst->cost[i];
			max_pos = i;
		}
		
	}

	return max_pos;
	
}

int * find_nodes(instance *inst, int pos){
	int *nodes_hierarchy = (int *) calloc(inst->nnodes, sizeof(int));;

	int first_node = pos / inst->nnodes;
	int second_node = pos % inst->nnodes;
	nodes_hierarchy[0] = first_node;
	nodes_hierarchy[1] = second_node;

	return nodes_hierarchy;
}

double delta_cost_extra_mileage(instance *inst, int i, int j, int h){
	double delta_cost = cost(i,h,inst) + cost(h,j,inst) - cost(i,j,inst);
	return delta_cost;
}

int * extra_mileage_step(instance *inst, int *uncovered_nodes, int current_length, int *nodes_hierarchy){

	int best_node;
	int best_edge;
	int *best_values = (int *) calloc(2, sizeof(int));
	double min_cost = __DBL_MAX__;
	double delta_cost = 0;

	for (int i = 0; i < inst->nnodes - current_length - 1; i++) // for each edge
	{	
		
		for (int j = 0; j < current_length; j++) // for each uncovered_node
		{
			delta_cost = delta_cost_extra_mileage(inst, nodes_hierarchy[i], nodes_hierarchy[i+1], uncovered_nodes[j]);
			if (min_cost > delta_cost) // choose the delta cost which is the minimumum
			{
				min_cost = delta_cost;
				best_node = j;
				best_edge = i;
			}
			
		}
		
	}
	best_values[0] = best_edge;
	best_values[1] = best_node;

	return best_values;

}

void calculate_best_val(instance *inst){
	double total_cost = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		total_cost += cost(inst->best_sol[i], inst->best_sol[i+1], inst);
	}

	inst->best_val = total_cost;
	
}
void calculate_extra_mileage_heuristics(instance *inst, int *nodes_hierarchy){
	int *uncovered_nodes = (int *) calloc(inst->nnodes, sizeof(int));
	int current_length = inst->nnodes;
	inst->best_sol = (int *) calloc(inst->nnodes+1, sizeof(int));
		// in order to close the tour +1 node in best_sol
	inst->best_val = (inst->cost[nodes_hierarchy[0]*inst->nnodes + nodes_hierarchy[1]]);

	for (int i = 0; i < inst->nnodes; i++)
	{
		uncovered_nodes[i] = i; 
	}
	int best_node_pos; int best_edge_pos; int *best_values_index; int node_to_new_successor;

	current_length --;
	uncovered_nodes[nodes_hierarchy[0]] = uncovered_nodes[current_length];
	current_length --;
	uncovered_nodes[nodes_hierarchy[1]] = uncovered_nodes[current_length];

	while (current_length != 0) 
	{ 
		best_values_index = extra_mileage_step(inst, uncovered_nodes, current_length, nodes_hierarchy);
		best_edge_pos = best_values_index[0]; 
		best_node_pos = best_values_index[1];

		// update nodes_hierarchy
		nodes_hierarchy = add_to_array(best_edge_pos, uncovered_nodes[best_node_pos], nodes_hierarchy, inst->nnodes);
		
		current_length--; // one node is already covered
		// update the uncovered_nodes
		uncovered_nodes[best_node_pos] = uncovered_nodes[current_length];
	}

	inst->best_sol = nodes_hierarchy;
	inst->best_val += inst->cost[inst->best_sol[inst->nnodes]*inst->nnodes + inst->best_sol[0]];
	inst->best_sol[inst->nnodes]= inst->best_sol[0]; // close the tour

	calculate_best_val(inst);

}

int extra_mileage_heuristic(instance *inst, int starting_mode){
	printf("\n_________________________________________________________\nInsertion Heuristic:\n");

	int *nodes_hierarchy = (int *) calloc(inst->nnodes, sizeof(int));
	double min_cost; int random_node_with_time_seed_pos;
	int *nodes_hierarchy_with_best_starting_couple = (int *) calloc(inst->nnodes, sizeof(int));
	int length_of_cost = inst->nnodes*inst->nnodes;
	int hullSize; Point* convexHull;

	switch (starting_mode)
	{
	case 0:  // distance A, B is max
		nodes_hierarchy = find_nodes(inst, maximum_cost_pos(inst));
		calculate_extra_mileage_heuristics(inst, nodes_hierarchy);
		break;
	case 1:  // A, B is random (A!=B)
		random_node_with_time_seed_pos = random_node_with_time_seed(length_of_cost);
		nodes_hierarchy = find_nodes(inst, random_node_with_time_seed_pos);
		// printf("\n random_nod_pos is %d \n", random_node_with_time_seed_pos);

		while (random_node_with_time_seed_pos % inst->nnodes == random_node_with_time_seed_pos / inst->nnodes) // A must be different from B
		{
			random_node_with_time_seed_pos = random_node_with_time_seed(length_of_cost);
			nodes_hierarchy = find_nodes(inst, random_node_with_time_seed_pos);
		}
		
		calculate_extra_mileage_heuristics(inst, nodes_hierarchy);
		break;

	case 2: // try all A, B but A!=B in the convexHull
		// because trying literally all takes more than 10 times slower
		// however, with convexHull we get best_value close to the other one
		convexHull = grahamScan(inst, &hullSize);
		
		nodes_hierarchy[0] = convexHull[0].id;
		nodes_hierarchy[1] = convexHull[1].id;
		calculate_extra_mileage_heuristics(inst, nodes_hierarchy);
		min_cost = inst->best_val;

		for (int i = 0; i < hullSize; i++)
		{
			for (int j = i+1; j < hullSize; j++)
			{ // cost index which A!=B and and check each couple just once
				nodes_hierarchy = find_nodes(inst, (convexHull[i].id*inst->nnodes)+convexHull[j].id); 
				calculate_extra_mileage_heuristics(inst, nodes_hierarchy);
				if (min_cost > inst->best_val)
				{
					min_cost = inst->best_val;
					nodes_hierarchy_with_best_starting_couple = copy_array(nodes_hierarchy, inst->nnodes);
				}
			}
			
		}
		
		printf("\n The best starting nodes are %d and %d \n"
		, nodes_hierarchy_with_best_starting_couple[0]
		, nodes_hierarchy_with_best_starting_couple[1]);
		
		free(convexHull);

		calculate_extra_mileage_heuristics(inst, nodes_hierarchy_with_best_starting_couple);
		break;
	default:
            printf("Error! \n starting_node is not in correct format\n");
	}

	return 0;

}

double delta_cost_two_opt(int a, int b, instance *inst){
	// if (a == b || a-b == 1 || b-a == 1 ) return 0;

	int a_node = inst->best_sol[a]; int b_node = inst->best_sol[b];
	int a_succ = inst->best_sol[a+1]; int b_succ = inst->best_sol[b+1];

	double old_cost = (cost(a_node,a_succ,inst) + cost(b_node,b_succ,inst));
	double new_cost = (cost(a_node,b_node,inst) + cost(a_succ,b_succ,inst));
	double delta_cost = new_cost - old_cost;

	return delta_cost;
}

int update_tour(int a, int b, instance *inst){

	int succ_a = a + 1;int succ_b = b+1;int counter = 1;int tmp;
	if (a < b)
	{	
		while( (a+counter) < (succ_b-counter) )
		{
			tmp = inst->best_sol[a+counter];
			inst->best_sol[a+counter] = inst->best_sol[succ_b-counter];
			inst->best_sol[succ_b-counter] = tmp;
			counter++;
		}
	}
	else if (a > b)
	{
		while( (b+counter) < (succ_a-counter) )
		{
			tmp = inst->best_sol[b+counter];
			inst->best_sol[b+counter] = inst->best_sol[succ_a-counter];
			inst->best_sol[succ_a-counter] = tmp;
			counter++;
		}
	}
	
	// printf("Improved the tour from %f", inst->best_val);
	calculate_best_val(inst);
	// printf(" to %f\n\n", inst->best_val);

	return 0;	
}

// we consider node a and node b
// and then the edge will be a, its successor(a+1), b and its successor(b+1)
int two_opt_refining_heuristic(instance *inst){
	printf("\n_________________________________________________________\n2-OPT Refining Heuristic:\n");

	double delta_cost;
	int a_with_min_delta_cost; int b_with_min_delta_cost;
	int update_switch = 0;
	double min_delta_cost;
	do
	{
		min_delta_cost = 0;
		update_switch = 0;
		for (int i = 0; i < inst->nnodes; i++) // nodes in 0 and 280 would have been be the same
		{
			for (int j = i+2; j < inst->nnodes; j++)
			{
				delta_cost = delta_cost_two_opt(i, j, inst);
				if (delta_cost < min_delta_cost)
				{
					min_delta_cost = delta_cost;
					a_with_min_delta_cost = i;
					b_with_min_delta_cost = j;
					update_switch = 1;
				}
				
			}
			
		}
		if (update_switch == 1)
		{
			// printf("Possible update detected, delta_cost: %f\n",min_delta_cost);
			update_tour(a_with_min_delta_cost, b_with_min_delta_cost, inst);
		}
		
	} while (update_switch == 1);
	printf("\n \tupdate in best_val after 2-OPT refining is %f\n", inst->best_val);

	return 0;
}

int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode){
	int lower_bound_tenure = upper_bound_tenure/5;
	int iteration_period = upper_bound_tenure - lower_bound_tenure;
	//it is hyperparameter, tenure shows the increase or decrease
			// becomes lower_bound or upper bound for subsequent period. 
	//f.e (upper_bound_tenure - lower_bound_tenure) iterations has been chosen 
	//so that both mode 0 and 1 has the same upper bounds
	int tenure;
	
	switch (tenure_mode)
		{
		case 0: // Tabu search with reactive step tenure
			// like _-_
			tenure = (iteration/iteration_period)%2 ? lower_bound_tenure : upper_bound_tenure;
			break;
		case 1:// Tabu search with reactive line tenure
			// tenure is steadily increasing and decreasing. like /\/
			tenure = (iteration/iteration_period)%2 ? --current_tenure : ++current_tenure;
			break;
		case 2: //random tenure length between upper and lower bound
			tenure = random_0_to_length(inst, upper_bound_tenure+1-lower_bound_tenure) + lower_bound_tenure;
			break;
		default:
			break;
		}
	
	return tenure;
}

int tabu_search(instance *inst, int tenure_mode){
	printf("\n_________________________________________________________\nTabu Search:\n");

	double t1 = second();

	int tenure; int upper_bound_tenure=100; tenure = upper_bound_tenure/5;//initialization

	double delta_cost; double min_delta_cost;
	int a_with_min_delta_cost; int b_with_min_delta_cost;
  
	int iteration = 0;
	int update_switch = 0;

	int* optimal_solution = copy_array(inst->best_sol, inst->nnodes+1);
	double optimal_value = inst->best_val;
	inst->tabu_list = (int *) calloc(inst->nnodes, sizeof(int));

	do
	{
		min_delta_cost = INFINITY;
		update_switch = 0;
		tenure = tenure_length_update(inst, tenure, iteration, upper_bound_tenure, tenure_mode);
		
		for (int i = 0; i < inst->nnodes; i++) // nodes in 0 and 280 would have been be the same
		{
			for (int j = i+2; j < inst->nnodes; j++)
			{
				delta_cost = delta_cost_two_opt(i, j, inst);

				if (iteration - inst->tabu_list[inst->best_sol[i]] <= tenure || 
								iteration - inst->tabu_list[inst->best_sol[j]] <= tenure)
				{  // do not consider tabu nodes
					if (delta_cost >= 0) // if it is negative then 
					{	// having an aspiration criteria that it would improve solution 
						//even it is in tabu list
						continue;
					}
					
				}
				
				if (delta_cost < min_delta_cost)
				{
					min_delta_cost = delta_cost;
					a_with_min_delta_cost = i;
					b_with_min_delta_cost = j;
					update_switch = 1;
				}
				
			}
			
		}
		if (update_switch == 1)
		{
			inst->tabu_list[inst->best_sol[a_with_min_delta_cost]] = iteration;
			inst->tabu_list[inst->best_sol[b_with_min_delta_cost]] = iteration;

			update_tour(a_with_min_delta_cost, b_with_min_delta_cost, inst);

			if (optimal_value > inst->best_val)
			{
				optimal_solution = copy_array(inst->best_sol, inst->nnodes);
				optimal_value = inst->best_val;
			}
			
		}
		iteration++;
	} while (second() - t1 < inst->timelimit);
	
	inst->best_sol = copy_array(optimal_solution, inst->nnodes);
	inst->best_val = optimal_value;

	free(optimal_solution);
	free(inst->tabu_list);

	return 0;
}
int copy_segment(instance *inst, int *old_solution, int starting_pos, int ending_pos, int into_pos){

	for (int i = starting_pos; i <= ending_pos; i++)
	{
		inst->best_sol[into_pos] = old_solution[i];
		into_pos++;
	}

	return into_pos;
}
int copy_segment_in_reverse_order(instance *inst, int *old_solution, int starting_pos, int ending_pos, int into_pos){

	for (int i = ending_pos; i >= starting_pos; i--)
	{
		inst->best_sol[into_pos] = old_solution[i];
		into_pos++;
	}

	return into_pos;

}

int new_tour_from_break_positions(instance *inst, int *break_positions, int arr_size){
	
	int* optimal_solution = copy_array(inst->best_sol, inst->nnodes+1);
	int into_pos = 0;

	for (int i = 0; i < arr_size; i+=2)
	{
		if(break_positions[i] > break_positions[i+1]) 
		{
			into_pos = copy_segment_in_reverse_order(inst, optimal_solution, break_positions[i+1], break_positions[i], into_pos);
			// into_pos += (break_positions[i] - break_positions[i+1]);
		}
		else 
		{
			into_pos = copy_segment(inst, optimal_solution, break_positions[i], break_positions[i+1], into_pos);
			// into_pos += (break_positions[i+1] - break_positions[i]);
		}
		// printf("\n%d hereeeeee %d\n", inst->best_sol[into_pos], into_pos);
		// print_array(break_positions,arr_size);
	}
	
	inst->best_sol[inst->nnodes] = inst->best_sol[0]; //close the tour

	calculate_best_val(inst);
	return 0;
}

int n_opt_kick(instance *inst, int n){
	printf("\n_________________________________________________________\n%d-OPT kick:\n", n);
	// we divide the tour into 5 pieces and then reconnect them randomly and
	int arr_size = n*2;
	int *break_positions = (int *) calloc(arr_size, sizeof(int));
	break_positions[0] = 0;
	int n_candidates = 1;
	while (n_candidates < arr_size -1) // bcz the last one is always known
	{
		break_positions[n_candidates] = random_0_to_length(inst, inst->nnodes-2) + 1;
		for (int j = 0; j < n_candidates; j++) // checking that they are different
		{									// and not neighbour
			if((break_positions[n_candidates] == break_positions[j]) ||
			 (break_positions[n_candidates] == (break_positions[j]+1)) || 
			 (break_positions[n_candidates]+1 == (break_positions[j]))) continue;
		}
		                                                                          
		n_candidates++;

		// successor in the previous solution
		break_positions[n_candidates] = break_positions[n_candidates-1] + 1; 
		n_candidates++;
	}
	break_positions[arr_size-1] = 279;
// ordering the elements of array bcz the random numbers we found were the positions
// in the tour, not the nodes. We consider ordered version as of original sol and then
// change the order in random way. We do like this bcz we divide the tour into n parts.
// if we use kick_positions elements as nodes then we would have to find initial nodes
// of those parts 
//but like this, it will be like from zero to a, from a to b, b to c, c to d, d to e, e to 0

	qsort(break_positions, arr_size, sizeof(int), compare);  
	shuffleArray(break_positions, arr_size);

	new_tour_from_break_positions(inst, break_positions, arr_size);
	// if(verify_tour(inst)==0) printf("\tIt is a tour!\n");

	printf("\n \tupdate in best_val after kick is %f\n", inst->best_val);
	return 0;
	//// verify if it's tour
}

int verify_tour(instance *inst){
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			if(i!=j && inst->best_sol[i]==inst->best_sol[j]) {
				printf("indexes %d and %d are the same nodes: %d\n", i,j, inst->best_sol[i]);
				return 1;
			}

		}

	}
	if (inst->best_sol[0] != inst->best_sol[inst->nnodes]){
		printf("first and the last nodes are not the same\n");
		return 1;
	}
	
	return 0;
}

int variable_neighborhood_search(instance *inst, int kick_neighborhood){
	printf("\n_________________________________________________________\nVariable Neighborhood Search:\n");

	double t1 = second();

	two_opt_refining_heuristic(inst);
	int* optimal_solution = copy_array(inst->best_sol, inst->nnodes+1);
	double optimal_value = inst->best_val;
	
	do
	{
		n_opt_kick(inst, kick_neighborhood);
		two_opt_refining_heuristic(inst);

		if (optimal_value > inst->best_val)
		{
			optimal_solution = copy_array(inst->best_sol, inst->nnodes+1);
			optimal_value = inst->best_val;
			printf("\n \tupdate in best_val %f\n", inst->best_val);
		}
		
	} while (second() - t1 < inst->timelimit);
	
	inst->best_sol = copy_array(optimal_solution, inst->nnodes+1);
	inst->best_val = optimal_value;

	free(optimal_solution);

	if(verify_tour(inst)==0) printf("\tIt is a tour!\n");
	return 0;
}