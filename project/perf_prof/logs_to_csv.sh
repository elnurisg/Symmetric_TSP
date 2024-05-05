#!/bin/bash

# Function to create CSV file
create_csv_for_algorithm() {
    local log_dir="$1"/logs
    local csv_file=""$1"/$(basename "$1").csv"

    # Count the number of unique types (excluding the seed part)
    local num_types=$(ls "${log_dir}"/*.log | xargs -n 1 basename | grep -c -E '_seed_0.log')

    # Header for the CSV file (with type names)
    local header=$(ls "${log_dir}"/*.log | xargs -n 1 basename | grep -E '_seed_0.log' | sed 's/_seed_0.log$//' | tr '\n' ', ' | sed 's/,$//')
    echo "$num_types, $header" > "$csv_file"

    # Extract unique types from filenames
    local types=$(ls "${log_dir}"/*.log | xargs -n 1 basename | grep -E '_seed_0.log' | sed 's/_seed_0.log$//' )

    # Iterate over each seed
    local seeds=$(ls "${log_dir}"/*.log | grep -oE 'seed_[0-9]+' | cut -d '_' -f 2 | sort -u | sort -n)

    # Iterate over each seed
    for seed in $seeds; do
        # Initialize the row for this seed
        local row="$seed"

        # Iterate over each type
        for type in $types; do
            # Find the corresponding log file
            local logfile="${log_dir}/${type}_seed_${seed}.log"

            # Extract the number from the log file
            local number=$(grep -oE 'best_val is [0-9]+\.[0-9]+' "$logfile" | cut -d ' ' -f 3)

            # Append the number to the row
            row="$row, $number"
        done
 
        # Append the row to the CSV file
        echo "$row" >> "$csv_file"
    done
}

heuristics() {

    # Call the function with the directory path
    create_csv_for_algorithm "./heuristics/greedy"

    create_csv_for_algorithm "./heuristics/greedy_two_opt"

    create_csv_for_algorithm "./heuristics/extra_mileage"

    create_csv_for_algorithm "./heuristics/extra_mileage_two_opt"
}

metaheuristics() {

    # Call the function with the directory path
    create_csv_for_algorithm "./metaheuristics/genetic_algorithm"

    create_csv_for_algorithm "./metaheuristics/simulated_annealing"

    create_csv_for_algorithm "./metaheuristics/tabu_search"

    create_csv_for_algorithm "./metaheuristics/VNS"
}

optimal() {
    create_csv_for_algorithm "./optimal_500_10"
}


# heuristics
# metaheuristics
# optimal