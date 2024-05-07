#!/bin/bash

# Path to executable
EXECUTABLE="../../tsp"

run_tsp_opt() {

    local NUM_INSTANCES=$1
    NUM_NODES=500  

    directory="./logs"
    mkdir -p $directory

    echo "TSP-opt logs are being generated for $NUM_NODES nodes."
    for ((seed = 0; seed < NUM_INSTANCES; seed++)); do
        for model_type in 0 1; do
            if [ $model_type -eq 0 ]; then
                local log_file="$directory/exact_method_Benders_loop_seed_${seed}.log"
                echo "Running seed $seed with model_type=$model_type (Benders loop)..."
            else
                local log_file="$directory/exact_method_Branch_and_Cut_seed_${seed}.log"
                echo "Running seed $seed with model_type=$model_type (Branch and Cut)..."
            fi
            # Run the executable with the specified parameters and save output to log file
            "$EXECUTABLE" -generate_random "$NUM_NODES" "-tsp_opt" "$model_type" -seed "$seed" > "$log_file"
        done
        echo "Done with seed $seed"
        echo
    done
}

run_tsp_opt 20