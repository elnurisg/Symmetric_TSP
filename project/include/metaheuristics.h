#include "../include/utilities.h"
#include "../include/heuristics.h"
#include <limits.h>

typedef struct{
	int *genes; // in our case it is tsp solution
	double fitness; // in our case it is cost of tsp and therefore smallest is the fittest
	double probability; // survival probability
						// its fitness divide by fitness of the champion in its generation
} Individual;


///////////////////////////////////Tabu Search/////////////////////////////////////////

/**
 * @brief Tabu Search: A metaheuristic method for solving combinatorial optimization problems by iteratively moving from a current solution to a neighboring solution while avoiding previously visited solutions.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] tenure_mode The tenure mode of the Tabu Search:
 *                         - 0: Reactive step tenure.
 *                         - 1: Reactive line tenure.
 *                         - 2: Random tenure.
 * @param[in] aspiration_criteria Flag indicating whether the Tabu Search should use aspiration criteria or not:
 *                                - 0: Use aspiration criteria.
 *                                - 1: Do not use aspiration criteria.
 * @return Returns 0 if the Tabu Search is successfully applied without encountering any errors.
**/
int tabu_search(instance *inst, int tenure_mode, int aspiration_criteria);

/**
 * @brief Updates the tenure length for Tabu Search.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] current_tenure Current tenure length.
 * @param[in] iteration Current iteration number.
 * @param[in] upper_bound_tenure Upper bound on the tenure length.
 * @param[in] tenure_mode The tenure mode of the Tabu Search.
 * @return Returns the updated tenure length.
**/
int tenure_length_update(instance *inst, int current_tenure, int iteration, int upper_bound_tenure, int tenure_mode);


///////////////////////////////////Variable Neighborhood Search/////////////////////////////////////////	

/**
 * @brief Variable Neighborhood Search: A metaheuristic method for solving combinatorial optimization problems by iteratively exploring different neighborhoods to find better solutions.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] kick_neighborhood The neighborhood type of the kick:
 *                              - 3: 3-opt kick.
 * 								- 4: 4-opt kick.
 *                              - 5: 5-opt kick.
 *                              - etc.
 * @return Returns 0 if the Variable Neighborhood Search is successfully applied without encountering any errors.
 * @note The 1-OPT kick doesn't update the tour and 2-OPT kick is not recommended as 2-OPT refining is applied. Therefore, it is recommended to use an integer greater than or equal to 3 and less than or equal to one-third of the number of nodes in the TSP instance..
**/
int variable_neighborhood_search(instance *inst, int kick_neighborhood);

/**
 * @brief Applies an N-opt kick to the current solution in Variable Neighborhood Search, changing N edges randomly to explore the neighborhood.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] n The size of the neighborhood for the kick.
 * @param[in, out] current_solution The current solution to be improved.
**/
void n_opt_kick(instance *inst, int n, int *current_solution);

/**
 * @brief Generates a new tour from the break positions.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] break_positions Array of break positions indicating segments in the current solution.
 * @param[in] arr_size Size of the break_positions array.
 * @param[in, out] current_solution The current solution from which a new tour will be generated.
**/
void new_tour_from_break_positions(instance *inst, int *break_positions, int arr_size, int *current_solution);

/**
 * @brief Copies a segment of the best solution of the instance (best_sol) into the given current solution.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] current_solution The current solution into which the segment will be copied.
 * @param[in] starting_pos Starting position in the current solution of the segment to be copied.
 * @param[in] ending_pos Ending position in the current solution of the segment to be copied.
 * @param[in] from_pos Starting position in best_sol from which the segment will be copied.
 * @return Returns the new position (from_pos) in best_sol for the next segment.
**/
int copy_segment(instance *inst, int *current_solution, int starting_pos, int ending_pos, int from_pos);

/**
 * @brief Copies a segment of the best solution of the instance (best_sol) in reverse order into the given current solution.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] current_solution The current solution from which the segment will be copied.
 * @param[in] starting_pos Starting position in the current solution of the segment to be copied.
 * @param[in] ending_pos Ending position in the current solution of the segment to be copied.
 * @param[in] from_pos Starting position in best_sol from which the segment will be copied.
 * @return Returns the new position (from_pos) in best_sol for the next segment.
**/
int copy_segment_in_reverse_order(instance *inst, int *current_solution, int starting_pos, int ending_pos, int from_pos);


///////////////////////////////////Simulated Annealing/////////////////////////////////////////

/**
 * @brief Simulated Annealing: A probabilistic metaheuristic algorithm for solving combinatorial optimization problems, inspired by the annealing process in metallurgy.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] annealing_iterations The number of iterations for the simulated annealing process.
 * @return Returns 0 if the Simulated Annealing algorithm is successfully applied without encountering any errors.
**/
int simulated_annealing(instance *inst, int annealing_iterations);

/**
 * @brief Calculates the acceptance probability using the Metropolis formula in Simulated Annealing.
 *
 * @param[in] delta_cost Change in cost between the current and proposed solution.
 * @param[in] temperature Current temperature in the annealing process.
 * @param[in] scaler Scaling factor for the Metropolis formula.
 * @return Returns the acceptance probability.
**/
double metropolis_formula(double delta_cost, double temperature, int scaler);

/**
 * @brief Performs the annealing process in Simulated Annealing.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] scaler Scaling factor for the acceptance probability (calculated by Metropolis formula).
 * @param[in, out] tsp_sol The given solution to which the annealing process will be applied.
**/
void annealing_process(instance *inst, int scaler, int *tsp_sol);

/**
 * @brief Calculates the average change in cost between two edges in the given TSP solution.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] tsp_sol The given TSP solution for which the average is going to be calculated.
 * @return Returns the average change in cost between two edges.
**/
double average_delta_cost_between_two_edges(instance *inst, int *tsp_sol);


///////////////////////////////////Genetic Algorithm/////////////////////////////////////////

/**
 * @brief Genetic Algorithm: A metaheuristic method for solving combinatorial optimization problems by mimicking the process of natural selection.
 *
 * @param[in, out] inst Input instance of the TSP problem.
 * @param[in] repair Switch for the genetic algorithm to use repair or not:
 *                  - 0: Punish their fitness with penalty without repairing.
 *                  - 1: Repair the bad genes.
 *                  - 2: Repair the bad genes and use two-opt refining while repairing the genes.
 * @param[in] cutting_type The cutting type of the crossover:
 *                         - 0: Cutting from the middle.
 *                         - 1: Cutting from a random position.
 *                         - 2: Cutting from one-fourth of the tour.
 * @return Returns 0 if the genetic algorithm is successfully applied without encountering any errors.
 *
 * @note If repair is off (repair=0), then it might not evolve to be a tour if time is not enough. Therefore, it has been decided to repair (with 2-OPT) the last survivor champion. This repair (with 2-OPT) is also applied for the last champion of repair=1 to ensure fairness.
**/
int genetic_algorithm(instance *inst, int repair, int cutting_type);

/**
 * @brief Initializes a population of individuals randomly for the genetic algorithm.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] population_size Size of the population to be initialized.
 * @param[in, out] population Array of individuals representing the population.
**/
void initialize_population_randomly(instance *inst, int population_size, Individual *population);

/**
 * @brief Initializes the genes of the individual randomly for the genetic algorithm.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] individual Individual to be initialized.
**/
void initialize_individual(instance *inst, Individual *individual);

/**
 * @brief Prints the genes of the individuals in the population.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] population_size Size of the population.
 * @param[in] population Array of individuals representing the population.
**/
void print_population(instance *inst, int population_size, Individual *population);

/**
 * @brief Frees the memory allocated for the population.
 *
 * @param[in] population Array of individuals representing the population.
 * @param[in] population_size Size of the population.
**/
void free_population(Individual *population, int population_size);

/**
 * @brief Frees the memory allocated for an individual.
 *
 * @param[in] individual Individual whose memory is to be freed.
**/
void free_Individual(Individual *individual);

/**
 * @brief Calculates the fitness of an individual's genes in the context of the TSP problem, where genes represent the tour and fitness represents the total tour cost.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] individual Individual for which fitness is calculated.
**/
void calculate_individual_fitness(instance *inst, Individual *individual);

/**
 * @brief Calculates the fitness of each individual in a population in the context of the TSP problem, where genes represent the tour and fitness represents the total tour cost.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] population Array of individuals representing the population.
 * @param[in] population_size Size of the population.
**/
void calculate_population_fitness(instance *inst, Individual *population, int population_size);

/**
 * @brief Modifies individuals in a population to avoid bad genes (tours) by detecting defects, such as subtours, missing elements, or incomplete tours where the last element differs from the first, and applying penalties to their fitness accordingly.
 * If defects are detected, the function increases the cost (fitness) of the individual by multiplying it by the number of defects found, discouraging their presence.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] population Array of individuals representing the population.
 * @param[in] population_size Size of the population.
**/
void avoid_bad_genes(instance *inst, Individual *population, int population_size);

/**
 * @brief Prints the genetic details (genes-tours) of the parent individuals (parent1 and parent2) and the resulting child individual.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] parent1 First parent individual.
 * @param[in] parent2 Second parent individual.
 * @param[in] child Resulting child individual.
 * 
 * @note It is designed to assist in visualizing the genetic information during various genetic algorithm operations.
**/
void print_parents_and_child(instance *inst, Individual *parent1, Individual *parent2, Individual *child);

/**
 * @brief Performs crossover operation to generate children individuals from parent individuals.
 * This function applies crossover operation to the parent individuals within the population to produce children individuals.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] population Array of individuals representing the parent population.
 * @param[in] population_size Size of the parent population.
 * @param[in, out] children Array to store the resulting children individuals.
 * @param[in] children_size Size of the children array.
 * @param[in] cutting_type The type of crossover performed:
 *                         - 0: Cutting from the middle. 
 *                         - 1: Cutting from a random position.
 *                         - 2: Cutting from one-fourth of the tour.
**/
void crossover(instance *inst, Individual *population, int population_size, Individual *children, int children_size, int cutting_type );

/**
 * @brief Mutates individuals in a population to introduce diversity and explore new solutions.
 * This function applies mutation operation to introduce random changes in the genetic makeup (tours) of individuals within the population.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] population Array of individuals representing the population.
 * @param[in] population_size Size of the population.
 * @param[in, out] mutations Array to store the mutated individuals.
 * @param[in] mutants_size Size of the mutations array.
**/
void mutate_population(instance *inst, Individual *population, int population_size, Individual *mutations, int mutants_size);

/**
 * @brief Introduces alterations to the genetic makeup (tour) of an individual to promote diversity.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] genes Array representing the genetic makeup (tour) of an individual.
 * 
* @note Alteration is done by swapping random genes (nodes) in the genetic makeup (tours) and the mutation rate is set to 10% internally which determines the upper bound of alterations (maximum number of nodes in tours which can be swapped).
**/
void alteration_of_genes(instance *inst, int *genes);

/**
 * @brief Calculates the survival probabilities of individuals in a generation based on their fitness, sorts the generation in ascending order, and determines the fittest individual (champion) for the next generation.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] Generation Array containing individuals representing the current generation.
 * @param[in] population Array of individuals representing the parent population.
 * @param[in] population_size Size of the parent population.
 * @param[in] children Array of individuals representing the offspring generated through crossover.
 * @param[in] children_size Size of the offspring population.
 * @param[in] mutations Array of individuals representing the mutated individuals.
 * @param[in] mutants_size Size of the mutated population.
 * @param[in, out] champion Fittest individual (champion) in the generation.
**/
void evaluate_generation(instance *inst, Individual *Generation, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion);

/**
 * @brief Calculates the probability of an individual surviving to the next generation relative to the champion individual's fitness.
 *
 * @param[in] individual Pointer to the individual whose survival probability is to be calculated.
 * @param[in] champion Pointer to the fittest individual (champion) in the current generation.
 * @return Returns the probability of the individual surviving to the next generation.
 * 
 * @note In the context of the TSP, where lower fitness (total cost of the tour) indicates better solutions, the survival probability of an individual is inversely proportional to its fitness relative to the champion individual.
**/
double probability_of_individual(Individual *individual, Individual *champion);

/**
 * @brief Calculates the survival probabilities for each individual in a generation relative to the champion individual.
 *
 * @param[in, out] Generation Array containing individuals representing the current generation.
 * @param[in] generation_size Size of the current generation.
 * @param[in] champion Fittest individual (champion) in the generation.
 * 
 * @note In the context of the TSP, where lower fitness (total cost of the tour) indicates better solutions, the survival probability of an individual is inversely proportional to its fitness relative to the champion individual.
**/
void probability_of_population(Individual *Generation, int generation_size, Individual *champion);

/**
 * @brief Implements elitism by giving more probability to the fittest individuals from the current generation to survive for the next generation (mimics natural selection).
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] population Array of individuals representing the current population (adults/parents).
 * @param[in] population_size Size of the current population.
 * @param[in] children Array of individuals representing the offspring generated through crossover.
 * @param[in] children_size Size of the offspring population.
 * @param[in] mutations Array of individuals representing the mutated individuals.
 * @param[in] mutants_size Size of the mutated population.
 * @param[in, out] champion Fittest individual (champion) in the generation.
 * 
 * @note Champion individual is guaranteed to survive to the next generation, and the rest of the population is selected based on their survival probabilities.
**/
void elitism(instance *inst, Individual *population, int population_size, Individual *children, int children_size, Individual *mutations, int mutants_size, Individual *champion);

/**
 * @brief Compares two individuals based on their fitness (total cost of the tour) for sorting purposes.
 *
 * @param[in] a Pointer to the first individual.
 * @param[in] b Pointer to the second individual.
 * @return Returns 1 if the first individual has higher fitness, -1 if the first individual has lower fitness (better tour as cost is lower), and zero if both individuals have the same fitness.
 **/
int compare_individuals(const void *a, const void *b);

/**
 * @brief Repairs bad genes in the individuals of a population to improve the quality of solutions.
 * In the context of TSP, repair mechanisms consist of eliminating defects such as subtours and adding missing nodes.
 * Additionally, if specified, it applies the 2-opt refinement heuristic to further improve the quality of the repaired solutions.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] population Array of individuals representing the population to be repaired.
 * @param[in] population_size Size of the population.
 * @param[in] apply_two_opt Flag indicating whether to apply the 2-opt refinement heuristic (2: Apply).
**/
void repair_bad_genes(instance *inst, Individual *population, int population_size, int apply_two_opt);

/**
 * @brief Ensures each node is visited exactly once by eliminating multiple visits in an individual's tour (genes) representation.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] individual Individual whose tour is processed to remove duplicate visits.
**/
void eliminate_multiple_visits(instance *inst, Individual *indiviual);

/**
 * @brief Repairs the tour representation (genes) of an Individual by inserting any missing nodes, applying the insertion heuristic.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in, out] individual Individual whose tour representation (genes) undergoes the process.
**/
void repair_extra_mileage(instance *inst, Individual *indiviual);

/**
 * @brief Allocates memory for a population of individuals based on the given population size.
 *
 * @param[in] inst Input instance of the TSP problem.
 * @param[in] population_size Size of the population to be allocated.
 * @return Pointer to the allocated array of individuals representing the population.
**/
Individual * allocate_population(instance *inst, int population_size);