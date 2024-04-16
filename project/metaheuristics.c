#include "metaheuristics.h"

double second();       


///////////////////////////////////Tabu Search/////////////////////////////////////////

int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode){
	int lower_bound_tenure = upper_bound_tenure/5;
	int iteration_period = upper_bound_tenure - lower_bound_tenure;
	//it is hyperparameter, tenure shows the increase or decrease
			// becomes lower_bound or upper bound for subsequent period. 
	//f.e (upper_bound_tenure - lower_bound_tenure) iterations has been chosen 
	//so that both mode 0 and 1 has the same upper bounds
	int tenure = 0;
	
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
			print_error("Wrong tenure_mode\n");
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


///////////////////////////////////Variable Neighborhood Search/////////////////////////////////////////

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


///////////////////////////////////Simulated Annealing/////////////////////////////////////////

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


///////////////////////////////////Genetic Algorithm/////////////////////////////////////////

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
Individual * find_champion_individual(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size){
	
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
void crossover(instance *inst, Individual *population, Individual *children, int children_size, int cutting_type ){
	
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
Individual * elitism(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion){
	
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
		crossover(inst, population, children, children_size, cutting_type);
		
		if (repair == 0) // OFF, punish their fitness with penalty
			avoid_bad_genes(inst, children, children_size);
		else  // ON 
			repair_bad_genes(inst, children, children_size, repair);
		
		mutate_population(inst, population, population_size, mutations, mutants_size);
		
		champion = find_champion_individual(population, population_size, children, children_size, mutations, mutants_size);
		printf("Champion fitness of generation%d is:%f\n", count_generations, champion->fitness);

		population = elitism(population, population_size, children, children_size, mutations, mutants_size, champion);
		
		count_generations++;
	} while (second() - t1 < inst->timelimit);
	
	// find the champion of the last survived generation of elitism
	champion = find_champion_individual(population, population_size, children, children_size, mutations, mutants_size);

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
