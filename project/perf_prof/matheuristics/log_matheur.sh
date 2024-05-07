#!/bin/bash

# Path to executable
EXECUTABLE="../../tsp"


run_hard_fixing() {

    local NUM_INSTANCES=$1
    NUM_NODES=1000  
    TIME_LIMIT=60

    directory="./hard_fixing/logs"
    mkdir -p $directory

    echo "Hard fixing algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES; seed++)); do
        for pfix in 0.2 0.3 0.5 0.8; do
            local log_file="$directory/hard_fixing_${pfix}_seed_${seed}.log"
            echo "Running seed $seed with pfix=$pfix..."
            # Run the executable with the specified parameters and save output to log file
            "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -hard_fixing "$pfix" -seed "$seed" > "$log_file"
        done
        echo "Done with seed $seed"
        echo
    done
}


run_local_branching() {

    local NUM_INSTANCES=$1
    NUM_NODES=1000  
    TIME_LIMIT=60

    directory="./local_branching/logs"
    mkdir -p $directory

    echo "Local branching algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES; seed++)); do
        for branching_constraint in 10 20 30; do
            local log_file="$directory/local_branching_${branching_constraint}_seed_${seed}.log"
            echo "Running seed $seed with branching constraint=$branching_constraint..."
            # Run the executable with the specified parameters and save output to log file
            "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -local_branching "$branching_constraint" -seed "$seed" > "$log_file"
        done
        echo "Done with seed $seed"
        echo
    done
}


seed_num_metaheur=10

## running not in parallel bcz it can affect the performance of the system as Cplex is already using a lot of resources
run_hard_fixing "$seed_num_metaheur"
run_local_branching "$seed_num_metaheur"
