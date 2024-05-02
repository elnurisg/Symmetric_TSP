#!/bin/bash

# Path to your executable
EXECUTABLE="../../tsp"


run_greedy() {

    local NUM_INSTANCES_GREEDY=$1
    local two_opt=$2
    NUM_NODES=2000  
    TIME_LIMIT=300

    if [[ $2 == "two_opt" ]]; then
        directory="./greedy_two_opt/logs"
    else
        directory="./greedy/logs"
    fi

    mkdir -p $directory

    echo "Greedy algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES_GREEDY; seed++)); do
        for starting_mode in 0 1 2; do
            for grasp in 0 1; do
                if [[ $2 == "two_opt" ]]; then
                    local log_file="$directory/greedy_starting_mode_${starting_mode}_grasp_${grasp}_${two_opt}_seed_${seed}.log"
                    echo "Running seed $seed with starting_mode=$starting_mode, grasp=$grasp and 2-OPT..."
                    # Run the executable with the specified parameters and save output to log file
                    "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -greedy "$starting_mode" "$grasp" "-$two_opt" -seed "$seed" > "$log_file"
                else
                    local log_file="$directory/greedy_starting_mode_${starting_mode}_grasp_${grasp}_seed_${seed}.log"
                    echo "Running seed $seed with starting_mode=$starting_mode, grasp=$grasp..."
                    # Run the executable with the specified parameters and save output to log file
                    "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -greedy "$starting_mode" "$grasp" -seed "$seed" > "$log_file"
                fi
            done
        done
        echo "Done with seed $seed"
        echo
    done
}


run_extra_mileage() {

    local NUM_INSTANCES_EXTRA_MILEAGE=$1
    local two_opt=$2
    NUM_NODES=2000  
    TIME_LIMIT=300

    if [[ $2 == "two_opt" ]]; then
        directory="./extra_mileage_two_opt/logs"
    else
        directory="./extra_mileage/logs"
    fi

    mkdir -p $directory

    echo "Extra Mileage algorithm logs are being generated for $NUM_NODES nodes and $TIME_LIMIT seconds time limit."
    for ((seed = 0; seed < NUM_INSTANCES_EXTRA_MILEAGE; seed++)); do
        for starting_mode in 0 1 2; do
            if [[ $2 == "two_opt" ]]; then
                local log_file="$directory/extra_mileage_starting_mode_${starting_mode}_${two_opt}_seed_${seed}.log"
                echo "Running seed $seed with starting_mode=$starting_mode and 2-OPT..."
                # Run the executable with the specified parameters and save output to log file
                "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -extra_mileage "$starting_mode" "-$two_opt" -seed "$seed" > "$log_file"
            else
                local log_file="$directory/extra_mileage_starting_mode_${starting_mode}_seed_${seed}.log"
                echo "Running seed $seed with starting_mode=$starting_mode..."
                # Run the executable with the specified parameters and save output to log file
                "$EXECUTABLE" -generate_random "$NUM_NODES" -time_limit "$TIME_LIMIT" -extra_mileage "$starting_mode" -seed "$seed" > "$log_file"
            fi
        done
        echo "Done with seed $seed"
        echo
    done
}


seed_num_heuristic=50

# run_greedy "$seed_num_heuristic"
# run_greedy "$seed_num_heuristic" "two_opt"
# run_extra_mileage "$seed_num_heuristic"
# run_extra_mileage "$seed_num_heuristic" "two_opt"

