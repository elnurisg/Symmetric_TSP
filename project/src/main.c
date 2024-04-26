#include "../include/tsp.h"           
#include "../include/matheuristics.h"

double second();
double random01();     
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst, Config *config); 
void compute_distance(instance *inst);
void generate_random_instances(instance *inst);
void free_instance(instance *inst);


int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	instance inst; 
	Config config;
	inst.tstart = second();

	parse_command_line(argc, argv, &inst, &config);     
	
	if (inst.nnodes == 0) read_input(&inst); 
	else generate_random_instances(&inst);
	compute_distance(&inst);
	// inst.nnodes + 1 is bcz of the closing the tour
	inst.best_sol = (int *) calloc(inst.nnodes+1, sizeof(int)); 

	if(config.greedy_heuristic == 1)
		if ( greedy_heuristic(&inst, config.greedy_starting_mode, config.greedy_grasp) ) print_error(" error within greedy_heuristic()");

	if (config.extra_mileage_heuristic == 1)
		if (extra_mileage_heuristic(&inst, config.extra_mileage_starting_mode)) print_error("Error within extra_mileage_heuristic()");

	if (config.variable_neighborhood_search == 1) 
		if (variable_neighborhood_search(&inst, config.VNS_kick_neighborhood)) print_error("Error within variable_neighborhood_search()");

	if (config.simulated_annealing == 1) 
		if (simulated_annealing(&inst, config.SA_iterations)) print_error("Error within simulated_annealing()");

	if (config.two_opt_refining_heuristic == 1) 
		if (two_opt_refining_heuristic(&inst, inst.best_sol, 0)) print_error("Error within two_opt_refining_heuristic()");

	if (config.tabu_search == 1) 
		if (tabu_search(&inst, config.tabu_tenure_mode, config.tabu_aspiration)) print_error("Error within tabu_search()");

	if (config.genetic_algorithm == 1) 
		if (genetic_algorithm(&inst, config.genetic_repair, config.genetic_cutting)) print_error("Error within genetic_algorithm()");

	if (config.TSPopt == 1) 
		if (TSPopt(&inst, config.TSPopt_model)) print_error("Error within TSPopt()");

	if (config.hard_fixing == 1) 
		if (hard_fixing(&inst, config.hard_fixing_probability)) print_error("Error within hard_fixing()");

	if (config.local_branching == 1) 
		if (local_branching(&inst, config.local_branching_constraint)) print_error("Error within local_branching()");


	double t2 = second();



	calculate_best_val(&inst);
	printf("\n \tbest_val is %f\n", inst.best_val);
	
	printf("\n----------------------------------------------------------------------------------------------");
	printf("\nTook %f seconds \n\n", t2 - inst.tstart);
	if(verify_tour(&inst, inst.best_sol)==0) printf("\tIt is a tour!\n");
	plot_tsp_tour(&inst, 1);
	
	free_instance(&inst);
	return 0; 
}         


void print_help(instance *inst, Config *config)
{
	printf("\n\navailable parameters (vers. 21-Apr-2024) --------------------------------------------------\n");

	printf("-file <input_file> %s\n", inst->input_file);
    printf("-time_limit <timelimit> %lf\n", inst->timelimit);
    printf("-generate_random <nnodes> %d\n", inst->nnodes);
    printf("-seed <random_seed> %d\n", inst->random_seed);
	printf("\nALGORITHMS:\n");
    printf("-greedy <starting_mode> %d <grasp_mode> %d\n", config->greedy_starting_mode, config->greedy_grasp);
	printf("-extra_mileage <starting_mode> %d\n", config->extra_mileage_starting_mode);
	printf("-two_opt\n");
	printf("-VNS <kick_neighborhood> %d\n", config->VNS_kick_neighborhood);
	printf("-SA <iterations> %d\n", config->SA_iterations);
    printf("-tabu_search <tenure_mode> %d <aspiration> %d\n", config->tabu_tenure_mode, config->tabu_aspiration);
    printf("-genetic_algorithm <repair_mode> %d <cutting_type> %d\n", config->genetic_repair, config->genetic_repair);
	printf("-tsp_opt <model_type> %d\n", config->TSPopt_model);
	printf("-hard_fixing <probability> %f\n", config->hard_fixing_probability);
	printf("-local_branching <branching_constraint> %d\n", config->local_branching_constraint);

	printf("\nenter -help or --help for help\n");
	printf("----------------------------------------------------------------------------------------------\n\n");

}

void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->best_sol);
	free(inst->cost);
}

void compute_distance(instance *inst)
{
	inst->cost = (double *) calloc(inst->nnodes*inst->nnodes, sizeof(double));
	for (int i = 0; i < inst->nnodes; i++) // compute distance for every node
	{
		for (int j = 0; j < inst->nnodes; j++) // with all other nodes
		{
			inst->cost[i*inst->nnodes + j] = dist(i, j, inst);
			if (VERBOSE > 1000)
				printf("cost between %d and %d is %f\n", i, j, inst->cost[i*inst->nnodes + j]);
		}
	}
}

void generate_random_instances(instance *inst)
{
	if (inst->random_seed !=0){
		srand(2635623+abs(inst->random_seed));
		for (size_t i = 0; i < 1000; i++) random();
	}

	inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 
	inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double)); 
	for (int i = 0; i < inst->nnodes; i++)
	{
		inst->xcoord[i] = ((double) rand() / RAND_MAX) * 10000;
		inst->ycoord[i] = ((double) rand() / RAND_MAX) * 10000;
	}
	
}

void read_input(instance *inst) // simplified TSP parser, not all SECTIONs detected  
{                      
	FILE *fin = fopen(inst->input_file, "r");  // enter the full path
	if ( fin == NULL ) print_error(" input file not found!");

	inst->nnodes = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int do_print = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{
		if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			token1 = strtok(NULL, "");  
			if ( VERBOSE >= 10 ) printf(" ... solving instance %s \n\n", token1);//with model %d\n\n", token1, inst->model_type);
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "TSP",3) != 0 ) print_error(" format error:  only TYPE == TSP implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}
		

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes); 
			// inst->demand = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));    
			active_section = 0;  
			continue;
		}

		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) 
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "EUC_2D", 6) != 0 ) print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}            
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;   
			continue;
		}
		
		if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}
				
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); 
			continue;
		}    
		  
		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");     
		    
	}                

	fclose(fin);    
	
}

void parse_command_line(int argc, char** argv, instance *inst, Config *config) 
{ 
	
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	inst->nnodes = 0;
	inst->best_val = INFINITY;
	strcpy(inst->input_file, "NULL");
	inst->random_seed = 0; 
	inst->timelimit = CPX_INFBOUND;
	inst->heur_flag = 0;
	
	// can give some optional values if needed
	config->greedy_starting_mode = -1;
    config->greedy_grasp = -1;
    config->extra_mileage_starting_mode = -1;
    config->VNS_kick_neighborhood = -1;
    config->SA_iterations = -1;
    config->tabu_tenure_mode = -1;
	config->tabu_aspiration = -1;
    config->genetic_repair = -1;
    config->genetic_cutting = -1;
    config->TSPopt_model = -1;
    config->hard_fixing_probability = -1;
    config->local_branching_constraint = -1;


    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { 
			if (i + 1 < argc)	inst->timelimit = atof(argv[++i]); // total time limit
			if (inst->timelimit == CPX_INFBOUND) {
				print_help(inst, config);
				print_error("Invalid argument for -time_limit.\n");
			}
			continue; }

		if ( strcmp(argv[i],"-generate_random") == 0 ) { 
			if (i + 1 < argc)	inst->nnodes = atoi(argv[++i]); 
			if (inst->nnodes == 0) {
				print_help(inst, config);
				print_error("Invalid argument for -generate_random.\n");
			}
			continue; } 

		if ( strcmp(argv[i],"-seed") == 0 ) { // random seed
			inst->random_seed = -1;
			if (i + 1 < argc)	inst->random_seed = abs(atoi(argv[++i]));
			if(inst->random_seed == -1) {
				print_help(inst, config);
				print_error("Invalid argument for -seed.\n");
			}
			continue; } 

		if ( strcmp(argv[i],"-greedy") == 0) { 
			config->greedy_heuristic = 1;
			config->greedy_starting_mode = -1;
			config->greedy_grasp = -1;
			if (i + 2 < argc){
				config->greedy_starting_mode = atoi(argv[++i]);
				config->greedy_grasp = atoi(argv[++i]);
			}
			if ((config->greedy_starting_mode!=0 && config->greedy_starting_mode!=1 && config->greedy_starting_mode!=2) || (config->greedy_grasp != 0 && config->greedy_grasp != 1)){
				print_help(inst, config);
				print_error("Wrong parameters for Greedy Heuristic");
			}
			continue;
			}

		if ( strcmp(argv[i],"-extra_mileage") == 0) { 
			config->extra_mileage_heuristic = 1;
			config->extra_mileage_starting_mode = -1;
			if (i + 1 < argc)
				config->extra_mileage_starting_mode = atoi(argv[++i]);
			if (config->extra_mileage_starting_mode!=0 && config->extra_mileage_starting_mode!=1 && config->extra_mileage_starting_mode!=2){
				print_help(inst, config);
				print_error("Wrong starting mode for Extra Mileage Heuristic");
			}
			continue;
			}

		if ( strcmp(argv[i],"-two_opt") == 0) { 
			config->two_opt_refining_heuristic = 1;
			continue;
			}

		if ( strcmp(argv[i],"-VNS") == 0) { 
			config->variable_neighborhood_search = 1;
			config->VNS_kick_neighborhood = -1;
			if (i + 1 < argc)
				config->VNS_kick_neighborhood = atoi(argv[++i]);
			if (config->VNS_kick_neighborhood==-1){
				print_help(inst, config);
				print_error("Wrong kick neighborhood for Variable Neighborhood");
			}
			continue;
			}
		
		if ( strcmp(argv[i],"-SA") == 0) { 
			config->simulated_annealing = 1;
			config->SA_iterations = -1;
			if (i + 1 < argc)
				config->SA_iterations = atoi(argv[++i]);
			if (config->SA_iterations==-1){
				print_help(inst, config);
				print_error("Wrong iteration number for Simulated Annealing");
			}
			continue;
			}

		if ( strcmp(argv[i],"-tabu_search") == 0) { 
			config->tabu_search = 1;
			config->tabu_tenure_mode = -1;
			config->tabu_aspiration = -1;
			if (i + 2 < argc){
				config->tabu_tenure_mode = atoi(argv[++i]);
				config->tabu_aspiration = atoi(argv[++i]);
			}
			if ((config->tabu_tenure_mode!=0 && config->tabu_tenure_mode!=1 && config->tabu_tenure_mode!=2) || (config->tabu_aspiration != 0 && config->tabu_aspiration != 1)){
				print_help(inst, config);
				print_error("Wrong parameters for Tabu Search");
			}
			continue;
			}

		if ( strcmp(argv[i],"-genetic_algorithm") == 0) { 
			config->genetic_algorithm = 1;
			config->genetic_repair = -1;
			config->genetic_cutting = -1;
			if (i + 2 < argc){
				config->genetic_repair = atoi(argv[++i]);
				config->genetic_cutting = atoi(argv[++i]);
			}
			if ((config->genetic_repair!=0 && config->genetic_repair!=1 && config->genetic_repair!=2) || (config->genetic_cutting != 0 && config->genetic_cutting != 1 && config->genetic_cutting != 2)){
				print_help(inst, config);
				print_error("Wrong parameters for Genetic Algorithm");
			}
			continue;
			}

		if ( strcmp(argv[i],"-tsp_opt") == 0) { 
			config->TSPopt = 1;
			config->TSPopt_model = -1;
			if (i + 1 < argc)
				config->TSPopt_model = atoi(argv[++i]);
			if (config->TSPopt_model!=0 && config->TSPopt_model!=1){
				print_help(inst, config);
				print_error("Wrong model type for TSPopt (exact solution)");
			}
			continue;
			}

		if ( strcmp(argv[i],"-hard_fixing") == 0) { 
			config->hard_fixing = 1;
			config->hard_fixing_probability = -1;
			if (i + 1 < argc)
				config->hard_fixing_probability = atof(argv[++i]);
			if (config->hard_fixing_probability < 0 && config->TSPopt_model > 1){
				print_help(inst, config);
				print_error("Wrong probability for Hard-fixing");
			}
			continue;
			}

		if ( strcmp(argv[i],"-local_branching") == 0) { 
			config->local_branching = 1;
			config->local_branching_constraint = -1;
			if (i + 1 < argc)
				config->local_branching_constraint = atoi(argv[++i]);
			if (config->local_branching_constraint==-1){
				print_help(inst, config);
				print_error("Please enter branching constraint for Local Branching");
			}
			continue;
			}


		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		print_help(inst, config);
	}        
	
	if ( help ) exit(1);

}    





