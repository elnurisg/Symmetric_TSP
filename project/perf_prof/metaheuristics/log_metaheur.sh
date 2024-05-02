#!/bin/bash

# Path to your executable
EXECUTABLE="../../tsp"


run_tabu() {

    local NUM_INSTANCES_TABU=$1
    NUM_NODES=2000  
    TIME_LIMIT=600

    directory="./tabu_search/logs"
    mkdir -p $directory

    echo "Tabu Search algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES_TABU; seed++)); do
        for tenure_mode in 0 1 2; do
            for aspiration in 0 1; do
                local log_file="$directory/greedy_1_0_tabu_tenure_mode_${tenure_mode}_aspiration_${aspiration}_seed_${seed}.log"
                echo "Running seed $seed with tenure_mode=$tenure_mode, aspiration=$aspiration which initialized by greedy 1 0..."
                # Run the executable with the specified parameters and save output to log file
                "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -greedy 1 0 "-tabu_search" "$tenure_mode" "$aspiration" -seed "$seed" > "$log_file"
            done
        done
        echo "Done with seed $seed"
        echo
    done
}

run_vns() {

    local NUM_INSTANCES_VNS=$1
    NUM_NODES=2000  
    TIME_LIMIT=600

    directory="./VNS/logs"
    mkdir -p $directory

    echo "VNS algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES_VNS; seed++)); do
        for n_opt in 3 5 7 10; do
            local log_file="$directory/greedy_1_0_VNS_kick_neighborhood_${n_opt}_seed_${seed}.log"
            echo "Running seed $seed with kick_neighborhood=$n_opt which initialized by greedy 1 0..."
            # Run the executable with the specified parameters and save output to log file
            "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -greedy 1 0 "-VNS" "$n_opt" -seed "$seed" > "$log_file"
        done
        echo "Done with seed $seed"
        echo
    done
}

run_simulated_annealing() {

    local NUM_INSTANCES_SA=$1
    NUM_NODES=2000  
    TIME_LIMIT=600

    directory="./simulated_annealing/logs"
    mkdir -p $directory

    echo "Simulated Annealing algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES_SA; seed++)); do
        for iterations in 1 2 3; do
            local log_file="$directory/greedy_1_0_simulated_annealing_iterations_${iterations}_seed_${seed}.log"
            echo "Running seed $seed with iterations=$iterations which initialized by greedy 1 0..."
            # Run the executable with the specified parameters and save output to log file
            "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -greedy 1 0 "-SA" "$iterations" -seed "$seed" > "$log_file"
        done
        echo "Done with seed $seed"
        echo
    done
}

run_genetic_alg() {

    local NUM_INSTANCES_GA=$1
    NUM_NODES=2000  
    TIME_LIMIT=600

    directory="./genetic_algorithm/logs"
    mkdir -p $directory

    echo "Genetic algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES_GA; seed++)); do
        for repair_mode in 0 1 2; do
            for cutting_type in 0 1 2; do
                local log_file="$directory/genetic_algorithm_repair_mode_${repair_mode}_cutting_type_${cutting_type}_seed_${seed}.log"
                echo "Running seed $seed with repair_mode=$repair_mode, cutting_type=$cutting_type..."
                # Run the executable with the specified parameters and save output to log file
                "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -genetic_algorithm "$repair_mode" "$cutting_type" -seed "$seed" > "$log_file"
            done
        done
        echo "Done with seed $seed"
        echo
    done
}

# All algos are initialized by greedy 1 0
seed_num_meta_heuristic=50

run_tabu "$seed_num_meta_heuristic"
# run_vns "$seed_num_meta_heuristic"
# run_simulated_annealing "$seed_num_meta_heuristic"
# run_genetic_alg "$seed_num_meta_heuristic"
