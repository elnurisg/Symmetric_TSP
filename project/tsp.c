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

// plot the solution in commands.txt file 
//and writes the solution to the file if writing_to_file is true
void plot_tsp_tour(instance *inst, int writing_to_file){

	if (writing_to_file == 1)
	{
		FILE *f = fopen("plot/data.dat", "w");
		if (f == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		for (int i = 0; i < inst->nnodes+1; i++)
		{
			fprintf(f, "%f %f\n", inst->xcoord[inst->best_sol[i]], inst->ycoord[inst->best_sol[i]]);
		}
		fclose(f);	
	}
	

    system("gnuplot ./plot/commands.txt");

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

int random_0_to_length_but_different_than_previous(instance *inst, int length, int previous_random_value){

	int random_value = random_0_to_length(inst, length);

	while (random_value == previous_random_value)
	{
		random_value = random_0_to_length(inst, length);
	}
	
	return random_value;
}

// when the starting_node is known, it applies the nearest neighbor algorithm
void calculate_greedy_steps(instance *inst, int starting_node_pos, int grasp){
	int *uncovered_nodes = (int *) calloc(inst->nnodes, sizeof(int));//[inst->nnodes];
	int current_length = inst->nnodes;

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
		

		inst->best_sol[(inst->nnodes) - current_length] = best_node;
		current_node = best_node;

		current_length--;
		uncovered_nodes[best_node_pos] = uncovered_nodes[current_length];

	}
		// close the tour
	inst->best_sol[inst->nnodes]= starting_node_pos;

	calculate_best_val(inst);
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
	// inst->best_sol = (int *) calloc(inst->nnodes+1, sizeof(int));
		// in order to close the tour +1 node in best_sol
	inst->best_val = (inst->cost[nodes_hierarchy[0]*inst->nnodes + nodes_hierarchy[1]]);

	for (int i = 0; i < inst->nnodes; i++)
	{
		uncovered_nodes[i] = i; 
	}
	int best_node_pos; int best_edge_pos; int *best_values_index;

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

double delta_cost_two_opt(int a, int b, instance *inst, int *tsp_sol){
	// if (a == b || a-b == 1 || b-a == 1 ) return 0;

	int a_node = tsp_sol[a]; int b_node = tsp_sol[b];
	int a_succ = tsp_sol[a+1]; int b_succ = tsp_sol[b+1];

	double old_cost = (cost(a_node,a_succ,inst) + cost(b_node,b_succ,inst));
	double new_cost = (cost(a_node,b_node,inst) + cost(a_succ,b_succ,inst));
	double delta_cost = new_cost - old_cost;

	return delta_cost;
}

int update_tour(int a, int b, instance *inst, int *tsp_sol, int is_instance){

	int succ_a = a + 1;int succ_b = b+1;int counter = 1;int tmp;
	if (a < b)
	{	
		while( (a+counter) < (succ_b-counter) )
		{
			tmp = tsp_sol[a+counter];
			tsp_sol[a+counter] = tsp_sol[succ_b-counter];
			tsp_sol[succ_b-counter] = tmp;
			counter++;
		}
	}
	else if (a > b)
	{
		while( (b+counter) < (succ_a-counter) )
		{
			tmp = tsp_sol[b+counter];
			tsp_sol[b+counter] = tsp_sol[succ_a-counter];
			tsp_sol[succ_a-counter] = tmp;
			counter++;
		}
	}
	
	if (is_instance == 0) // it is on
	{
		// printf("Improved the tour from %f", inst->best_val);
		calculate_best_val(inst);
		// printf(" to %f\n\n", inst->best_val);
	}
	

	return 0;	
}

// we consider node a and node b
// and then the edge will be a, its successor(a+1), b and its successor(b+1)
int two_opt_refining_heuristic(instance *inst, int *tsp_sol, int is_instance){
	if (is_instance == 0)	
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
				delta_cost = delta_cost_two_opt(i, j, inst, tsp_sol);
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
			update_tour(a_with_min_delta_cost, b_with_min_delta_cost, inst, tsp_sol, is_instance);
		}
		
	} while (update_switch == 1);
	
	tsp_sol[inst->nnodes] = tsp_sol[0];

	if (is_instance == 0)
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
				delta_cost = delta_cost_two_opt(i, j, inst, inst->best_sol);

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

			update_tour(a_with_min_delta_cost, b_with_min_delta_cost, inst, inst->best_sol, 0);

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
	inst->best_sol[inst->nnodes] = inst->best_sol[0];

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
	// we divide the tour into n pieces and then reconnect them randomly and
	int arr_size = n*2;
	int *break_positions = (int *) calloc(arr_size, sizeof(int));
	break_positions[0] = 0;
	int n_candidates = 1;
	int choose_candidate_again = 0;
	while (n_candidates < arr_size -1) // bcz the last one is always known
	{
		choose_candidate_again = 0;
		break_positions[n_candidates] = random_0_to_length(inst, inst->nnodes-2) + 1;
		for (int j = 0; j < n_candidates; j++) // checking that they are different
		{									// and not neighbour
			if((break_positions[n_candidates] == break_positions[j]) ||
			 (break_positions[n_candidates] == (break_positions[j]+1)) || 
			 (break_positions[n_candidates]+1 == (break_positions[j]))) choose_candidate_again = 1;
		}

		if (choose_candidate_again == 1) continue;

		n_candidates++;

		// successor in the previous solution
		break_positions[n_candidates] = break_positions[n_candidates-1] + 1; 
		n_candidates++;
	}
	break_positions[arr_size-1] = inst->nnodes-1;
// ordering the elements of array bcz the random numbers we found were the positions
// in the tour, not the nodes. We consider ordered version as of original sol and then
// change the order in random way. We do like this bcz we divide the tour into n parts.
// if we use kick_positions elements as nodes then we would have to find initial nodes
// of those parts 
//but like this, it will be like from zero to a, from a to b, b to c, c to d, d to e, e to 0

	qsort(break_positions, arr_size, sizeof(int), compare);  
	shuffle_array_for_kick(break_positions, arr_size);

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
	
	if (kick_neighborhood == 0 || kick_neighborhood > (inst->nnodes/3)) print_error("Error! Kick neighborhood can not be 0 or more than the one third of the number of nodes\n");
	// 0-OPT kick is not meaningful
	
	// char filename[50];
    // snprintf(filename, sizeof(filename), "cost_plot/costs_VNS_%d-OPT.txt", kick_neighborhood);

	double t1 = second();

	two_opt_refining_heuristic(inst, inst->best_sol, 0);
	// write_cost_to_file(inst->best_val, filename, 1);

	int* optimal_solution = copy_array(inst->best_sol, inst->nnodes+1);
	double optimal_value = inst->best_val;
	
	do
	{	
		n_opt_kick(inst, kick_neighborhood);
		// write_cost_to_file(inst->best_val, filename, 0);
	
		two_opt_refining_heuristic(inst, inst->best_sol, 0);
		// write_cost_to_file(inst->best_val, filename, 0);

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

double metropolis_formula(double delta_cost, double Temprature, int scaler){
	return exp(-(delta_cost / scaler) / Temprature);
}

int annealing_process(instance *inst, int scaler){
	
	// choosing random nodes for 2-OPT
	int a; int b; double delta_cost;
	int T_max = inst->best_val; //we use the proxy of initial sol

	for (int T = T_max; T > 1; T--)
	{
		do  //be sure b is after a and not successor of a
		{
			a = random_0_to_length(inst, inst->nnodes);
			b = random_0_to_length(inst, inst->nnodes);
		} while ( a+1 > b );
		delta_cost = delta_cost_two_opt(a, b, inst, inst->best_sol);

		if (random01() <= metropolis_formula(delta_cost, T, scaler))
			update_tour(a, b, inst, inst->best_sol, 0);

	}

	// two_opt_refining_heuristic(inst, inst->best_sol, 0); //it is 0 degree now and so,
	// all negative costs should be applied
	return 0;
}

double average_delta_cost_between_two_edges(instance *inst){

	double sum_delta_cost = 0; double numbers = 0;
	for (int i = 0; i < inst->nnodes; i++) // nodes in 0 and 280 would have been be the same
	{
		for (int j = i+2; j < inst->nnodes; j++)
		{
			sum_delta_cost += delta_cost_two_opt(i, j, inst, inst->best_sol);
			numbers++;
		}
	}
	double average_delta_cost = sum_delta_cost/numbers;

	return average_delta_cost;
	
}
int simulated_annealing(instance *inst){
	printf("\n_________________________________________________________\nSimulated Annealing:\n");

	double t1 = second();


	int scaler = average_delta_cost_between_two_edges(inst);
	double optimal_value = inst->best_val;
	int* optimal_solution = copy_array(inst->best_sol, inst->nnodes+1);
	
	do
	{
		//apply annealing twice (heat up right after cooling)
		annealing_process(inst, scaler);annealing_process(inst, scaler);
		two_opt_refining_heuristic(inst, inst->best_sol, 0); //it is 0 degree now and so,
								// all negative costs should be applied

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

	// if(verify_tour(inst)==0) printf("\tIt is a tour!\n");
	return 0;

}
void initialize_individual(instance *inst, Individual *individual){

	individual->genes = (int *) calloc((inst->nnodes+1), sizeof(int));

	for (int j = 0; j < inst->nnodes; j++) {
		individual->genes[j] = j;
	}

	// tour is closed inside of this function
	shuffle_tsp_sol(individual->genes, inst->nnodes); 

}
int initialize_population_randomly(instance *inst, int population_size, Individual *population){
	
	for (int i = 0; i < population_size; i++) 
	{
		initialize_individual(inst, &population[i]);

    }

	calculate_population_fitness(inst, population, population_size);

	return 0;
}

void print_population(instance *inst, int population_size, Individual *population){
	for (int i = 0; i < population_size; i++)
	{
		print_array(population[i].genes, inst->nnodes+1);
		printf("\n");
	}
	
}
void free_population(Individual *population, int population_size) {
    for (int i = 0; i < population_size; i++) {
        free_Individual(&population[i]);
    }
}

void free_Individual(Individual *individual) {
	if (individual->genes == NULL) // to avoid double free
	{
    	free(individual->genes);
		individual->genes = NULL;
	}
}


// calculate fitness values of all population and return to the ondex with best fitness value
Individual * find_champion_individual(instance *inst, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size){
	
	Individual *champion_individual = (Individual *)malloc(sizeof(Individual));
	champion_individual->fitness = INFINITY;

	for (int i = 0; i < population_size; i++)
	{
		if (population[i].fitness < champion_individual->fitness)
		{
			champion_individual = &population[i];
		}	
	}
	for (int i = 0; i < children_size; i++)
	{
		if (children[i].fitness < champion_individual->fitness)
		{
			champion_individual = &children[i];
		}	
	}
	for (int i = 0; i < mutants_size; i++)
	{
		if (mutations[i].fitness < champion_individual->fitness)
		{
			champion_individual = &mutations[i];
		}	
	}

	return champion_individual;
}

void calculate_individual_fitness(instance *inst, Individual *individual){

	double fitness = 0;
	for (int j = 0; j < (inst->nnodes+1); j++)
	{
		fitness += cost(individual->genes[j], individual->genes[j+1], inst);			
	}

	individual->fitness = fitness;

}

void calculate_population_fitness(instance *inst, Individual *population, int population_size){


	for (int i = 0; i < population_size; i++)
	{
		calculate_individual_fitness(inst, &population[i]);
	}

}
void crossover(instance *inst, Individual *population, int population_size, Individual *children, int children_size, int cutting_type ){
	
	int parent1; int parent2; int cutting_position;

	switch (cutting_type)
	{
	case 0: // close to middle
		cutting_position = inst->nnodes / 2;
		break;
	case 1: // random
		cutting_position = random_0_to_length(inst, inst->nnodes);
		break;
	case 2: // inst->nnodes/4
		cutting_position = inst->nnodes / 4;
		break;
	default:
		break;
	}

	for (int i = 0; i < children_size; i++)
	{
		parent1 = random_0_to_length(inst, inst->nnodes);
		parent2 = random_0_to_length_but_different_than_previous(inst, inst->nnodes, parent1);
		children[i].genes = combine_two_tours_from_pos(population[parent1].genes,
					 population[parent2].genes, inst->nnodes, cutting_position);
	}
	calculate_population_fitness(inst, children, children_size);
	
}

void print_parents_and_child(instance *inst, Individual *parent1, Individual *parent2, Individual *child){
	print_array(parent1->genes, inst->nnodes+1);
	print_array(parent2->genes, inst->nnodes+1);
	print_array(child->genes, inst->nnodes+1);
}

// detect deffects in genes, punish wrong genes with penalty and increase cost(fitness)
// in case of deffects, increase cost with the power of deffects
// and if there are more than 1 occurence of node means at least one missing node,
// so value for deffects is >=2 if not 0
void avoid_bad_genes(instance *inst, Individual *children, int children_size){
	
	int deffects; int count_of_node;

	for (int z = 0; z < children_size; z++)
	{
		deffects = 0;

		for (int i = 0; i < inst->nnodes; i++)
		{
			count_of_node = 0;
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (i == children[z].genes[j]) // there are subtours or missing elements
				{
					count_of_node++;
				}
			}

			if (count_of_node != 1)
			{
				deffects++;
			}
		}

		if (deffects != 0)
				children[z].fitness = pow(children[z].fitness, deffects);

	}
	
}

void alteration_of_genes(instance *inst, int *genes){
// 10% of mutation_rate has been decided
	int alternation_upper_bound = inst->nnodes / 10; // max amount of nodes can be alternated
	int alternation_size = random_0_to_length(inst, alternation_upper_bound);
	int swap_pos1; int swap_pos2; int temp;

	for (int i = 0; i < alternation_size; i++)
	{
		swap_pos1 = random_0_to_length(inst, inst->nnodes);
		swap_pos2 = random_0_to_length_but_different_than_previous(inst, inst->nnodes, swap_pos1);
		
		temp = genes[swap_pos1];
        genes[swap_pos1] = genes[swap_pos2];
        genes[swap_pos2] = temp;

	}
    // close the tour
    genes[inst->nnodes] = genes[0];
}

void mutate_population(instance *inst, Individual *population, int population_size, Individual *mutations, int mutants_size){

	int mutant;

	for (int i = 0; i < mutants_size; i++)
	{
		mutant = random_0_to_length(inst, population_size);
		mutations[i].genes = copy_array(population[mutant].genes, inst->nnodes+1);
		alteration_of_genes(inst, mutations[i].genes);
	}
	calculate_population_fitness(inst, mutations, mutants_size);
	
}

double probability_of_individual(Individual *individual, Individual *champion){
	return (champion->fitness / individual->fitness);
}

// Function to compare individuals based on their probabilities (for sorting)
int compare_individuals(const void *a, const void *b) {
    const Individual *individualA = (const Individual *)a;
    const Individual *individualB = (const Individual *)b;

    if (individualA->probability < individualB->probability) {
        return 1;  // sorting in descending order
    } else if (individualA->probability > individualB->probability) {
        return -1;
    } else {
        return 0;
    }
}

// calculate survival_probabilities with respect to its fitness and add into its structure
// gather all generation and sort them in the descending order of their survival rate
Individual * survival_probabilities_of_generation(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion){
	
	int generation_size = population_size + children_size + mutants_size;
    Individual *Generation = (Individual *)malloc(generation_size * sizeof(Individual));
	int index = 0;

	for (int j = 0; j < children_size; j++)
	{
		children[j].probability = probability_of_individual(&children[j], champion);
		Generation[index++] = children[j];
	}

	for (int z = 0; z < population_size; z++)
	{
		population[z].probability = probability_of_individual(&population[z], champion);
		Generation[index++] = population[z];
	}

	for (int h = 0; h < mutants_size; h++)
	{
		mutations[h].probability = probability_of_individual(&mutations[h], champion);
		Generation[index++] = mutations[h];
	}

    // sort the combined array based on probabilities
    qsort(Generation, generation_size, sizeof(Individual), compare_individuals);

    return Generation;
}

// kill population with bad genes and choose the fittest population for the next generation
// return new_population
Individual * elitism(instance *inst, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion){
	
	Individual *new_population = (Individual *) calloc(population_size, sizeof(Individual));
	Individual *Generation;

	int new_population_count = 0;
	Generation = survival_probabilities_of_generation(population, population_size, children, children_size, mutations, mutants_size, champion);

	while (new_population_count < population_size)
	{
		for (int i = 0; i < children_size + population_size + mutants_size; i++)
		{
			if (random01() <= Generation[i].probability){
				new_population[new_population_count] = Generation[i];
				new_population_count++;
			}

		}
				
	}

	return new_population;
}

void eliminate_multiple_visits(instance *inst, Individual *individual){

	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i+1; j < inst->nnodes; j++)
		{
				if ((individual->genes[i] == individual->genes[j]) && (individual->genes[i] != -1)) 
				{
					individual->genes = remove_from_array(j, individual->genes, inst->nnodes+1);
					j--; // bcz new pos j is different node now
					
				}
		}
	}
}

void repair_bad_genes(instance *inst, Individual *children, int children_size, int apply_two_opt){
		
	for (int z = 0; z < children_size; z++)
	{
		eliminate_multiple_visits(inst, &children[z]);
		repair_extra_mileage(inst, &children[z]);
		if(apply_two_opt == 2)
			two_opt_refining_heuristic(inst, children[z].genes, 1);
		calculate_individual_fitness(inst, &children[z]);
	}
}

void repair_extra_mileage(instance *inst, Individual *individual){

	int *uncovered_nodes = (int *) calloc(inst->nnodes, sizeof(int));
	int uncovered_length = 0; int missing_node = 0;
	
	// finding uncovered nodes
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			if (i == individual->genes[j] ) {
				missing_node = 1;
				break;
			}  
		}		
		if (missing_node == 0)
		{
			uncovered_nodes[uncovered_length] = i;
			uncovered_length++;		
		}
		missing_node = 0;
	}

	// and applying extra-mileage
	int best_node_pos; int best_edge_pos; int *best_values_index;
	while (uncovered_length != 0) 
	{ 
		best_values_index = extra_mileage_step(inst, uncovered_nodes, uncovered_length, individual->genes);
		best_edge_pos = best_values_index[0]; 
		best_node_pos = best_values_index[1];

		// update individual->genes
		individual->genes = add_to_array(best_edge_pos, uncovered_nodes[best_node_pos], individual->genes, inst->nnodes);
		
		uncovered_length--; // one node is already covered
		// update the uncovered_nodes
		uncovered_nodes[best_node_pos] = uncovered_nodes[uncovered_length];
	}

	// close the tour
	individual->genes[inst->nnodes] = individual->genes[0];

	// free the uncovered_nodes as it is not going to be used later
	free(uncovered_nodes);
}


int genetic_algorithm(instance *inst, int repair, int cutting_type){
	printf("\n_________________________________________________________\nGenetic Alghorithm:\n");
	double t1 = second();

	int population_size = 1000;
	int children_size = population_size/2;
	int mutants_size = population_size/10;
	Individual *champion; int count_generations = 0;

	Individual *population = (Individual *) calloc(population_size, sizeof(Individual));
	Individual *children = (Individual *) calloc(children_size, sizeof(Individual));
	Individual *mutations = (Individual *) calloc(mutants_size, sizeof(Individual));
	
	initialize_population_randomly(inst, population_size, population);
	count_generations++;

	do
	{
		crossover(inst, population, population_size, children, children_size, cutting_type);
		
		if (repair == 0) // OFF, punish their fitness with penalty
			avoid_bad_genes(inst, children, children_size);
		else  // ON 
			repair_bad_genes(inst, children, children_size, repair);
		
		mutate_population(inst, population, population_size, mutations, mutants_size);
		
		champion = find_champion_individual(inst, population, population_size, children, children_size, mutations, mutants_size);
		printf("Champion fitness of generation%d is:%f\n", count_generations, champion->fitness);

		population = elitism(inst, population, population_size, children, children_size, mutations, mutants_size, champion);
		
		count_generations++;
	} while (second() - t1 < inst->timelimit);
	
	// find the champion of the last survived generation of elitism
	champion = find_champion_individual(inst, population, population_size, children, children_size, mutations, mutants_size);

	// update the best_sol in case it is better
	if (champion->fitness < inst->best_val)
	{
		inst->best_sol = copy_array(champion->genes, inst->nnodes+1);
		inst->best_val = champion->fitness;
	}
 
	// free memory
	free_population(population, population_size);
	free_population(children, children_size);
	free_population(mutations, mutants_size);
	free_Individual(champion);

	return 0;
}

void fill_best_sol(int *edges, instance *inst){
	print_array(edges, inst->nnodes*2);

	int temp1; int temp2;
	// find another node which is the same as first_node and move it to the end with its pair
	// in order to have the same node in first and last position
	// so that it can be a tour
	move_node_in_edges(edges, inst->nnodes*2, 0, inst->nnodes*2-1);
	// make every same node to appear next to each other (except first and the last)
	// move_node_in_edges(edges, inst->nnodes*2, 1, 1+1);
	// move_node_in_edges(edges, inst->nnodes*2, 3, 3+1);
	for (int i = 1; i < inst->nnodes*2-5; i+=2)
	{
		move_node_in_edges(edges, inst->nnodes*2, edges[i], i+1);
	}
	print_array(edges, inst->nnodes*2);
	
}

int TSPopt(instance *inst) // enable the timelimit and check if it is optimal solution, mip function? and if not feasible, how to deal with it? mb printerror
/**************************************************************************************************************************/
{  
	// nodes in odd and in even positions are the nodes which are connected. we store them here
	// int *edges = (int *) calloc(inst->nnodes*2, sizeof(int));
	// int index_edge=0;
	
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if ( error ) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if ( error ) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	// CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	// if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	// CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);	
	// CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0); 

	error = CPXmipopt(env,lp);
	if ( error ) 
	{
		printf("CPX error code %d\n", error);
		print_error("CPXmipopt() error"); 
	}


	int ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(ncols, sizeof(double));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int *best_succ = (int *) calloc(inst->nnodes, sizeof(int));
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int *) calloc(1, sizeof(int));

	double LB = -CPX_INFBOUND; // lower bound
	if ( greedy_heuristic(inst, 0, 0) ) print_error(" error within greedy_heuristic()");
	if ( two_opt_refining_heuristic(inst, inst->best_sol, 0) ) print_error(" error within two_opt_refining_heuristic()");
	calculate_best_val(inst);
	double UB = inst->best_val; // upper bound
	double objval; double incumbent_value;

	do
	{
		if (CPXmipopt(env,lp)) print_error("CPXmipopt() error"); 
		if ( CPXgetbestobjval(env, lp, &objval) ) print_error("CPXgetbestobjval() error");	
		if ( CPXgetx(env, lp, xstar, 0, ncols-1) ) print_error("CPXgetx() error");
		
		// new LB is assumed to increase as having more constraints
		LB = (LB > objval) ? LB : objval;

		build_sol(xstar, inst, succ, comp, ncomp);
		printf("ncomp: %d\n", *ncomp);

		if (*ncomp <= 1){ // optimal solution is found by CPLEX so we can break the loop
			incumbent_value = calc_incumbent_value(succ, inst);
			LB = incumbent_value; // due to the type of cost(int and happens even using double),
			//when distance between nodes is very large,
			// they can be not equal and we ensure that they are equal
		}	// as solution is optimal
		else
		{
			for (int component_num = 1; component_num < *ncomp; component_num++)
			{
				add_subtour_constraint(env, lp, inst, comp, component_num, ncols);
				// add a subtour elimination constraint
			}	
			patching_heuristic(inst, succ, comp, ncomp);
			incumbent_value = calc_incumbent_value(succ, inst);
		}

		CPXsetdblparam(env, CPX_PARAM_TILIM, (inst->timelimit - second() - inst->tstart)); 
		
		if (incumbent_value < UB)
		{
			UB = incumbent_value;
			best_succ = copy_array(succ, inst->nnodes);
			// store the best succ solution found so far
		}		
		printf("LB: %f\t",LB); printf("UB: %f\t",UB); (LB == UB) ? printf("||| Solution is Optimal\n") : printf("||| %f%% gap\n", (UB-LB)/LB*100);
	} while ((LB < (1-EPSILON) * UB) && (second() - inst->tstart < inst->timelimit));




	if (VERBOSE >= 100){
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = i+1; j < inst->nnodes; j++ )
			{
				if ( xstar[xpos(i,j,inst)] > 0.5 ){
					printf("  ... x(%3d,%3d) = 1\n", i+1,j+1);
					// if (index_edge >= inst->nnodes*2)
					// {
					// 	print_error("index_edge is out of range: TSPopt");
					// 	return -1;
					// }
					
					// edges[index_edge] = i;
					// edges[index_edge+1] = j;
					// index_edge+=2;
				}
			}
		}
	}

	double initialized_UB = inst->best_val; // Greedy heuristic + 2-OPT is used as intial upper bound
	store_solution(inst, best_succ, inst->best_sol);
	calculate_best_val(inst);
	if (initialized_UB < inst->best_val) printf("ATTENTION:	Initial Upper Bound (using Greedy heuristic + 2-OPT) is better than the best solution found by CPLEX, please increase the time_limit!\n");

	// CPXwriteprob(env, lp, "model.lp", NULL);  

	free(best_succ);
	free(xstar);
	free(succ);
	free(comp);
	free(ncomp);
	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 


	// fill_best_sol(edges, inst);

	return 0; // or an appropriate nonzero error code

}

double delta_cost_patching(int a, int b, instance *inst, int *succ){
	
	double delta_cost = 0;
	double cost_of_old_edges = dist(a, succ[a], inst) + dist(b, succ[b], inst);
	double cost_of_new_edges = dist(a, succ[b], inst) + dist(b, succ[a], inst);

	delta_cost = cost_of_new_edges - cost_of_old_edges;

	return delta_cost;
}

void patching_heuristic(instance *inst, int *succ, int *comp, int *ncomp){

	double min_delta_cost; double delta_cost;
	int min_a; int min_b;
	int *sol = (int *) calloc(inst->nnodes+1, sizeof(int));

	while (*ncomp != 1)
	{
		min_delta_cost = INFINITY;
		for (int a = 0; a < inst->nnodes; a++)
		{
			for (int b = 0; b < inst->nnodes; b++)
			{
				if (comp[a] > comp[b]) // if a and b is in different component,
				{                     //  '>' is trick to not to check the same pairs twice
					delta_cost = delta_cost_patching(a, b, inst, succ);
					if (delta_cost < min_delta_cost)
					{
						min_delta_cost = delta_cost;
						min_a = a; min_b = b;
					}
					
				}
				
			}
			
		}
		update_succ_and_comp(inst, min_a, min_b, succ, comp);
		--*ncomp;
	}

	store_solution(inst, succ, sol);
	two_opt_refining_heuristic(inst, sol, 1);
	store_succ(inst, succ, sol);

	free(sol);
}

void store_succ(instance *inst, int *succ, int *sol){
	
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[sol[i]] = sol[i+1];
	}
	
}


void update_succ_and_comp(instance *inst, int min_a, int min_b, int *succ, int *comp){

	int temp = succ[min_a];
	succ[min_a] = succ[min_b];
	succ[min_b] = temp;

	int comp_a = comp[min_a];
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (comp[i] == comp_a) comp[i] = comp[min_b];
	}

}

double calc_incumbent_value(int *succ, instance *inst){

	int *sol = (int *) calloc(inst->nnodes+1, sizeof(int));
	double val = 0;

	store_solution(inst, succ, sol);
	for (int i = 0; i < inst->nnodes; i++)
	{
		val += dist(sol[i], sol[i+1], inst);
	}

	free(sol);
	return val;
}

/* Transforms the solution stored in the successor array format into a custom solution format for the Traveling Salesman Problem (TSP).
   Parameters:
   - inst: Pointer to the TSP instance structure.
   - succ: Array storing the successor of each node in the solution.
   - sol: Array to store the custom solution format.
*/
void store_solution(instance *inst, int *succ, int *sol){
	
	sol[0] = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		sol[i+1] = succ[sol[i]];
	}
	
}

void add_subtour_constraint(CPXENVptr env, CPXLPptr lp, instance *inst, int *comp, int component_num, int ncols)
{
	int *index = (int *) calloc(ncols, sizeof(int));
	double *value = (double *) calloc(ncols, sizeof(double));
	double rhs = -1;
	char sense = 'L';                            // 'L' for less or equal
	// char *cname = (char *) calloc(100, sizeof(char));
	int nnz = 0;

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( comp[i] == component_num )
		{
			rhs++;
			for (int j = i+1; j < inst->nnodes; j++)
			{
				if (comp[j] == component_num)
				{
					index[nnz] = xpos(i,j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
		}
	}
	int izero = 0;
	// sprintf(cname, "SEC(%d)", *ncomp);
	if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, NULL) ) print_error("CPXaddrows(): error 2");

	free(value);
	free(index);
	// free(cname);		
}
/***************************************************************************************************************************/
int xpos(int i, int j, instance *inst)      // to be verified                                           
/***************************************************************************************************************************/
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}
	

/***************************************************************************************************************************/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = dist(i,j,inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

// add the degree constraints 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

}



/*********************************************************************************************************************************/
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{   

#ifdef DEBUG
	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPS && fabs(xstar[k]-1.0)) > EPS ) print_error(" wrong xstar in build_sol()");
			if ( xstar[k] > 0.5 ) 
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);
#endif

	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
}