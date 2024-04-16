#include "utilities.h"
#include "heuristics.h"

typedef struct{
	int *genes; // in our case it is tsp solution
	double fitness; // in our case it is cost of tsp and therefore smallest is the fittest
	double probability; // survival probability
						// its fitness divide by fitness of the champion in its generation
} Individual;


///////////////////////////////////Tabu Search/////////////////////////////////////////

/**
 * Tabu Search
 * @param[in, out] inst Input instance of the tsp problem. 
 * @param[in] tenure_mode The tenure mode of the tabu search, 0 for reactive step tenure, 1 for reactive line tenure and 2 is for random tenure.
 * @return Returns 0 if the tabu search is successfully applied without encountering any errors.
**/
int tabu_search(instance *inst, int tenure_mode);
int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode);
int update_tour_and_tabu_list(int a, int b, instance *inst);

///////////////////////////////////Variable Neighborhood Search/////////////////////////////////////////	
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

///////////////////////////////////Simulated Annealing/////////////////////////////////////////
/**
 * Simulated Annealing
 * @param[in, out] inst Input instance of the tsp problem.
 * @return Returns 0 if the simulated annealing is successfully applied without encountering any errors.
**/
int simulated_annealing(instance *inst);
double metropolis_formula(double delta_cost, double Temprature, int scaler);
int annealing_process(instance *inst, int scaler);
double average_delta_cost_between_two_edges(instance *inst);

///////////////////////////////////Genetic Algorithm/////////////////////////////////////////
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
Individual * find_champion_individual(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size);
void calculate_individual_fitness(instance *inst, Individual *individual);
void calculate_population_fitness(instance *inst, Individual *population, int population_size);
void avoid_bad_genes(instance *inst, Individual *children, int children_size);
void print_parents_and_child(instance *inst, Individual *parent1, Individual *parent2, Individual *child);
void crossover(instance *inst, Individual *population, Individual *children, int children_size, int cutting_type );
void mutate_population(instance *inst, Individual *population, int population_size, Individual *mutations, int mutants_size);
void alteration_of_genes(instance *inst, int *genes);
Individual * survival_probabilities_of_generation(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion);
Individual * elitism(Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion);
int compare_individuals(const void *a, const void *b);
void repair_bad_genes(instance *inst, Individual *children, int children_size, int apply_two_opt);
void eliminate_multiple_visits(instance *inst, Individual *indiviual);
void repair_extra_mileage(instance *inst, Individual *indiviual);
