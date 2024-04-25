#include "../include/metaheuristics.h"


///////////////////////////////////Tabu Search/////////////////////////////////////////

int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode){
	
	int lower_bound_tenure = upper_bound_tenure/5;
	int iteration_period = upper_bound_tenure - lower_bound_tenure; int tenure = 0;
	//it is hyperparameter, tenure shows the increase or decrease
			// becomes lower_bound or upper bound for subsequent period. 
	//f.e (upper_bound_tenure - lower_bound_tenure) iterations has been chosen 
	//so that both mode 0 and 1 has the same upper bounds
	
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


int tabu_search(instance *inst, int tenure_mode, int aspiration_criteria){

	printf("\n_________________________________________________________\nTabu Search:\n\n");
	
	if (inst->heur_flag == 0) print_error("Greedy or Insertion Heuristic should be applied before Metaheuristic (tabu search)\n");

	int upper_bound_tenure = (inst->nnodes/10 > 100) ? 100 : inst->nnodes/10;
	int tenure = upper_bound_tenure/5;//initialization
	double delta_cost; double min_delta_cost;int a_with_min_delta_cost; int b_with_min_delta_cost;
	int iteration = 0;
	int* curr_solution = copy_to_new_array(inst->best_sol, inst->nnodes+1);double curr_value = inst->best_val;
	
	inst->tabu_list = (int *) calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
		inst->tabu_list[i] = - LARGE_INT_NUMBER;

	do
	{
		min_delta_cost = INFINITY;
		tenure = tenure_length_update(inst, tenure, iteration, upper_bound_tenure, tenure_mode);
		
		for (int i = 0; i < inst->nnodes; i++) // nodes in 0 and 280 would have been be the same
		{
			for (int j = i+1; j < inst->nnodes; j++)
			{
				if(curr_solution[j] == curr_solution[i+1] || curr_solution[i] == curr_solution[j+1]) continue;

				delta_cost = delta_cost_two_opt(i, j, inst, curr_solution);

				if (iteration - inst->tabu_list[i] <= tenure || 
								iteration - inst->tabu_list[j] <= tenure)
				{  // do not consider tabu nodes
					if(aspiration_criteria == 1) continue; // it is off
					else if (aspiration_criteria == 0 && delta_cost >= 0) continue;
					//if it is negative then having an aspiration criteria that it would improve solution
						//even it is in tabu list
				}

				if (delta_cost < min_delta_cost)
				{
					min_delta_cost = delta_cost;
					a_with_min_delta_cost = i;
					b_with_min_delta_cost = j;
				}
				
			}
			
		}

		inst->tabu_list[a_with_min_delta_cost] = iteration;
		inst->tabu_list[b_with_min_delta_cost] = iteration;

		update_tour(a_with_min_delta_cost, b_with_min_delta_cost, curr_solution);
		curr_value = curr_value + min_delta_cost;

		if (curr_value < inst->best_val)
		{
			copy_array(curr_solution, inst->nnodes, inst->best_sol);
			inst->best_val = curr_value;
			printf("[Tabu Search] Update in best_val curr value is %f, iteration is %d tenure %d\n", inst->best_val, iteration, tenure);

		}
		printf("[Tabu Search] Current value is %f, iteration is %d tenure %d\n", curr_value, iteration, tenure);
		iteration++;
	} while (!time_limit_expired(inst));
	
	free(curr_solution);
	free(inst->tabu_list);

	return 0;

}


///////////////////////////////////Variable Neighborhood Search/////////////////////////////////////////

int copy_segment(instance *inst, int *current_solution, int starting_pos, int ending_pos, int from_pos){

	for (int i = starting_pos; i <= ending_pos; i++)
	{
		current_solution[from_pos] = inst->best_sol[i];
		from_pos++;
	}

	return from_pos;
}


int copy_segment_in_reverse_order(instance *inst, int *current_solution, int starting_pos, int ending_pos, int from_pos){

	for (int i = ending_pos; i >= starting_pos; i--)
	{
		current_solution[from_pos] = inst->best_sol[i];
		from_pos++;
	}

	return from_pos;
}


void new_tour_from_break_positions(instance *inst, int *break_positions, int arr_size, int *current_solution){
	
	int from_pos = 0;

	for (int i = 0; i < arr_size; i+=2) // for each segment
	{
		if(break_positions[i] > break_positions[i+1]) 
			from_pos = copy_segment_in_reverse_order(inst, current_solution, break_positions[i+1], break_positions[i], from_pos);
		else 	
			from_pos = copy_segment(inst, current_solution, break_positions[i], break_positions[i+1], from_pos);		
	}
	
	current_solution[inst->nnodes] = current_solution[0]; //close the tour

}


void n_opt_kick(instance *inst, int n, int *current_solution){

	// we divide the tour into n pieces and then reconnect them randomly and
	int arr_size = n*2;
	int *break_positions = (int *) calloc(arr_size, sizeof(int)); break_positions[0] = 0;
	int n_candidates = 1;int choose_candidate_again = 0;

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

	new_tour_from_break_positions(inst, break_positions, arr_size, current_solution);

	free(break_positions);

}


int variable_neighborhood_search(instance *inst, int kick_neighborhood){
	printf("\n_________________________________________________________\nVariable Neighborhood Search:\n\n");
	
	if (inst->heur_flag == 0) print_error("Greedy or Insertion Heuristic should be applied before Metaheuristic (VNS)\n");
	if (kick_neighborhood == 0 || kick_neighborhood > (inst->nnodes/3)) print_error("Kick neighborhood can not be 0 or more than the one third of the number of nodes\n");
	// 0-OPT kick is not meaningful
	
	// char filename[50];
    // snprintf(filename, sizeof(filename), "cost_plot/costs_VNS_%d-OPT.txt", kick_neighborhood);

	two_opt_refining_heuristic(inst, inst->best_sol, 0);
	// write_cost_to_file(inst->best_val, filename, 1);

	int* current_solution = copy_to_new_array(inst->best_sol, inst->nnodes+1);
	double current_value = inst->best_val;
	
	do
	{	
		n_opt_kick(inst, kick_neighborhood, current_solution);
		current_value = calculate_total_cost(inst, current_solution);
		printf("\n[%d-OPT kick] Total cost after kick is %f\n", kick_neighborhood, current_value);
		// write_cost_to_file(current_value, filename, 0);
	
		two_opt_refining_heuristic(inst, current_solution, 1);
		current_value = calculate_total_cost(inst, current_solution);
		printf("\n[2-OPT Refining] Total cost after refining is %f\n", current_value);
		// write_cost_to_file(current_value, filename, 0);

		if (current_value < inst->best_val)
		{
			copy_array(current_solution, inst->nnodes+1, inst->best_sol);
			inst->best_val = current_value;
			printf("\n \t[VNS] Update in best_val %f\n", inst->best_val);
		}
		
	} while (!time_limit_expired(inst));
	
	free(current_solution);

	return 0;

}


///////////////////////////////////Simulated Annealing/////////////////////////////////////////

double metropolis_formula(double delta_cost, double Temprature, int scaler){
	return exp(-(delta_cost / scaler) / Temprature);
}


void annealing_process(instance *inst, int scaler, int *tsp_sol){
	
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
		
		delta_cost = delta_cost_two_opt(a, b, inst, tsp_sol);

		if (random01() <= metropolis_formula(delta_cost, T, scaler))
			update_tour(a, b, tsp_sol);
	}

}


double average_delta_cost_between_two_edges(instance *inst, int *tsp_sol){

	double sum_delta_cost = 0; double numbers = 0; double average_delta_cost = 0;

	for (int i = 0; i < inst->nnodes; i++) // nodes in 0 and 280 would have been be the same
	{
		for (int j = i+1; j < inst->nnodes; j++)
		{
			sum_delta_cost += delta_cost_two_opt(i, j, inst, tsp_sol);
			numbers++;
		}
	}

	average_delta_cost = sum_delta_cost/numbers;

	return average_delta_cost;
	
}


int simulated_annealing(instance *inst, int annealing_iterations){
	
	printf("\n_________________________________________________________\nSimulated Annealing:\n\n");
	
	if (inst->heur_flag == 0) print_error("Greedy or Insertion Heuristic should be applied before Metaheuristic (Simulated Annealing)\n");

	int scaler = average_delta_cost_between_two_edges(inst, inst->best_sol);
	double optimal_value = inst->best_val;
	int* optimal_solution = copy_to_new_array(inst->best_sol, inst->nnodes+1);
	
	do
	{
		//apply annealing annealing_iterations times (heat up right after cooling)
		for (int i = 0; i < annealing_iterations; i++)
			annealing_process(inst, scaler, optimal_solution);		

		two_opt_refining_heuristic(inst, optimal_solution, 0); //it is 0 degree now and so,
													// all negative costs should be applied

		optimal_value = calculate_total_cost(inst, optimal_solution);

		printf("[Simulated Annealing] Current value is %f\n", optimal_value);

		if (optimal_value < inst->best_val)
		{
			copy_array(optimal_solution, inst->nnodes+1, inst->best_sol);
			inst->best_val = optimal_value;
			printf("[Simulated Annealing] Update in best_val %f\n", inst->best_val);
		}
		
	} while (!time_limit_expired(inst));

	free(optimal_solution);

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


void initialize_population_randomly(instance *inst, int population_size, Individual *population){
	
	for (int i = 0; i < population_size; i++) 
	{
		initialize_individual(inst, &population[i]);

    }

	calculate_population_fitness(inst, population, population_size);

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


void calculate_individual_fitness(instance *inst, Individual *individual){

	double fitness = 0;

	for (int j = 0; j < inst->nnodes; j++)
	{
		fitness += cost(individual->genes[j], individual->genes[j+1], inst);			
	}

	individual->fitness = fitness;

}


void calculate_population_fitness(instance *inst, Individual *population, int population_size){
	
	for (int i = 0; i < population_size; i++)
		calculate_individual_fitness(inst, &population[i]);

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
		parent1 = random_0_to_length(inst, population_size);
		parent2 = random_0_to_length_but_different_than_previous(inst, population_size, parent1);
		combine_two_tours_from_pos(population[parent1].genes,
					 population[parent2].genes, inst->nnodes, cutting_position, children[i].genes);
	}
	calculate_population_fitness(inst, children, children_size);
	
}


void print_parents_and_child(instance *inst, Individual *parent1, Individual *parent2, Individual *child){
	
	print_array(parent1->genes, inst->nnodes+1);
	print_array(parent2->genes, inst->nnodes+1);
	print_array(child->genes, inst->nnodes+1);

}


// detect deffects in genes, punish wrong genes with penalty and increase cost(fitness)
// in case of deffects, increase cost with multiplying by number of the deffects
// and if there are more than 1 occurence of node means at least one missing node,
// so value for deffects is >=2 if not 0

void avoid_bad_genes(instance *inst, Individual *population, int population_size){
	
	int deffects; int count_of_node;

	for (int z = 0; z < population_size; z++)
	{
		deffects = 0;

		for (int i = 0; i < inst->nnodes; i++)
		{
			count_of_node = 0;

			for (int j = 0; j < inst->nnodes; j++)
			{
				if (i == population[z].genes[j])	count_of_node++;
			}

			if (count_of_node != 1)		deffects++; // there are subtours or missing elements
		}

		// if it is not tour
		if (population[z].genes[inst->nnodes] != population[z].genes[0])	deffects++;

		if (deffects != 0)	population[z].fitness += population[z].fitness * deffects;

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
		copy_array(population[mutant].genes, inst->nnodes+1, mutations[i].genes);
		alteration_of_genes(inst, mutations[i].genes);
	}

	calculate_population_fitness(inst, mutations, mutants_size);
	
}


double probability_of_individual(Individual *individual, Individual *champion){
	return (champion->fitness / individual->fitness);
}


int compare_individuals(const void *a, const void *b) {	// Function to compare individuals based on their fitness (for sorting)
  
    const Individual *individualA = (const Individual *)a;
    const Individual *individualB = (const Individual *)b;

    if (individualA->fitness > individualB->fitness) {
        return 1;  // sorting in ascending order
    } else if (individualA->fitness < individualB->fitness) {
        return -1;
    } else {
        return 0;
    }

}


// calculate survival_probabilities with respect to its fitness and add into its structure
// gather all generation and sort them in the descending order of their survival rate

void evaluate_generation(instance *inst, Individual *Generation, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion){
	
	int generation_size = population_size + children_size + mutants_size;
	int index = 0;

	for (int j = 0; j < children_size; j++)
	{
		copy_array(children[j].genes, inst->nnodes+1, Generation[index].genes);
		calculate_individual_fitness(inst, &Generation[index]);
		index++;
	}

	for (int z = 0; z < population_size; z++)
	{
		copy_array(population[z].genes, inst->nnodes+1, Generation[index].genes);
		calculate_individual_fitness(inst, &Generation[index]);
		index++;
	}

	for (int h = 0; h < mutants_size; h++)
	{
		copy_array(mutations[h].genes, inst->nnodes+1, Generation[index].genes);
		calculate_individual_fitness(inst, &Generation[index]);
		index++;
	}

    // sort the combined array based on their fitness
    qsort(Generation, generation_size, sizeof(Individual), compare_individuals);
	
	copy_array(Generation[0].genes, inst->nnodes+1, champion->genes);
	champion->fitness = Generation[0].fitness;

	probability_of_population(Generation, generation_size, champion);

}


void probability_of_population(Individual *Generation, int generation_size, Individual *champion){
	
	for (int i = 0; i < generation_size; i++)
		Generation[i].probability = probability_of_individual(&Generation[i], champion);

}


// kill population with bad genes and choose the fittest population for the next generation

void elitism(instance *inst, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion){
		
	int generation_size = population_size + children_size + mutants_size;
    Individual *Generation = allocate_population(inst, generation_size);

	int new_population_count = 0;
	evaluate_generation(inst, Generation, population, population_size, children, children_size, mutations, mutants_size, champion);

	while (new_population_count < population_size)
	{
		for (int i = 0; i < generation_size; i++)
		{
			if (random01() <= Generation[i].probability && new_population_count < population_size){
				copy_array(Generation[i].genes, inst->nnodes+1, population[new_population_count].genes);
				new_population_count++;
			}
		}		
	}

	calculate_population_fitness(inst, population, population_size);
	
	free_population(Generation, generation_size);

}

void eliminate_multiple_visits(instance *inst, Individual *individual){

	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i+1; j < inst->nnodes; j++)
		{
			if ((individual->genes[i] == individual->genes[j]) && (individual->genes[i] != -1)) 
			{
				remove_from_array(j, individual->genes, inst->nnodes+1, individual->genes);
				j--; // bcz new pos j is different node now	
			}
		}
	}

}


void repair_bad_genes(instance *inst, Individual *population, int population_size, int apply_two_opt){
		
	for (int z = 0; z < population_size; z++)
	{
		eliminate_multiple_visits(inst, &population[z]);
		repair_extra_mileage(inst, &population[z]);
		if(apply_two_opt == 2)
			two_opt_refining_heuristic(inst, population[z].genes, 1);
		calculate_individual_fitness(inst, &population[z]);
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
	while (uncovered_length != 0) 
	{ 
		extra_mileage_step(inst, uncovered_nodes, uncovered_length, individual->genes);
		uncovered_length--; // one node is already covered
	}

	// close the tour
	individual->genes[inst->nnodes] = individual->genes[0];

	// free the uncovered_nodes as it is not going to be used later
	free(uncovered_nodes);

}


Individual * allocate_population(instance *inst, int population_size){

	Individual *population = (Individual *) calloc(population_size, sizeof(Individual));
	
	for (int i = 0; i < population_size; i++) 
		population[i].genes = (int *) calloc((inst->nnodes+1), sizeof(int));
	
	return population;

}


int genetic_algorithm(instance *inst, int repair, int cutting_type){
	
	printf("\n_________________________________________________________\nGenetic Alghorithm:\n\n");

	int population_size = 1000; int children_size = population_size/2; int mutants_size = population_size/10;
	Individual *champion = allocate_population(inst, 1); int count_generations = 0;
	Individual *population = (Individual *) calloc(population_size, sizeof(Individual));
	Individual *children = allocate_population(inst, children_size);
	Individual *mutations = allocate_population(inst, mutants_size);
	
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
		
		elitism(inst, population, population_size, children, children_size, mutations, mutants_size, champion);
		printf("[Genetic Algorithm] Champion fitness of generation%d is:%f\n", count_generations, champion->fitness);
		
		count_generations++;
		
	} while (!time_limit_expired(inst));

	// if repair is off, then it might not evolve to be a tour, so it has been decided to repair the last survivor champion
	// and as it is applied for repair=0, in order to be fair it is also applied for repair=1
	// even if repair (or it is on but without 2-OPT) is off, apply repairing with 2-OPT to champion
	if(repair == 0 || repair == 1) repair_bad_genes(inst, champion, 1, 2);

	printf("[Genetic Algorithm] Fitness of the Champion of Natural Selection over %d generations is:%f\n", --count_generations, champion->fitness);
	
	// update the best_sol in case champion is better
	if (champion->fitness < inst->best_val && verify_tour(inst, champion->genes)==0)
	{
		copy_array(champion->genes, inst->nnodes+1, inst->best_sol);
		inst->best_val = champion->fitness;
	}
 
	// free memory
	free_population(population, population_size);
	free_population(children, children_size);
	free_population(mutations, mutants_size);
	free_Individual(champion);

	return 0;

}