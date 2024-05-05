#!/bin/bash

# Path to executable
EXECUTABLE="../tsp"


find_optimal() {

    local NUM_INSTANCES=$1
    local NUM_NODES=$2  

    directory="./optimal_${NUM_NODES}_${NUM_INSTANCES}/logs"
    mkdir -p $directory

    echo "Branch and Cut algorithm logs are being generated for $NUM_NODES nodes in order to find optimal results for $NUM_INSTANCES "
    for ((seed = 0; seed < NUM_INSTANCES; seed++)); do
        local log_file="$directory/optimal_seed_${seed}.log"
        echo "Running seed $seed ..."
        # Run the executable with the specified parameters and save output to log file
        "$EXECUTABLE" -generate_random "$NUM_NODES" -tsp_opt 1 -seed "$seed" > "$log_file"
        echo "Done with seed $seed"
        echo
    done
}

seed_num=10
nodes_num=500

find_optimal $seed_num $nodes_num


