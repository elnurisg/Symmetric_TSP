#!/bin/bash

# Call the function with the directory path
create_csv_for_algorithm() {
    local log_dir="$1"/logs
    local csv_file=""$1"/$(basename "$1").csv"
    local quality=$2

    if [ "$quality" = "cost" ] || [ -z "$quality" ]; then
        local pattern='best_val is [0-9]+\.[0-9]+'
    elif [ "$quality" = "time" ]; then
        local pattern='Took [0-9]+\.[0-9]+ seconds'
    fi

    # Count the number of unique types (excluding the seed part)
    local num_types=$(ls "${log_dir}"/*.log | xargs -n 1 basename | grep -c -E '_seed_0.log')

    local header=$(ls "${log_dir}"/*.log | xargs -n 1 basename | grep -E '_seed_0.log' | sed 's/_seed_0.log$//' | tr '\n' ', ' | sed 's/,$//')
    echo "$num_types, $header" > "$csv_file"

    # Extract unique types from filenames
    local types=$(ls "${log_dir}"/*.log | xargs -n 1 basename | grep -E '_seed_0.log' | sed 's/_seed_0.log$//' )
    local seeds=$(ls "${log_dir}"/*.log | grep -oE 'seed_[0-9]+' | cut -d '_' -f 2 | sort -u | sort -n)

    for seed in $seeds; do
        # Initialize the row for this seed
        local row="$seed"

        for type in $types; do
            local logfile="${log_dir}/${type}_seed_${seed}.log"
            local number=$(grep -oE "$pattern" "$logfile" | grep -oE '[0-9]+\.[0-9]+')
            row="$row, $number"
        done
 
        echo "$row" >> "$csv_file" # Append the row to the csv file
    done
}

heuristics() {
    create_csv_for_algorithm "./heuristics/greedy"
    create_csv_for_algorithm "./heuristics/greedy_two_opt"
    create_csv_for_algorithm "./heuristics/extra_mileage"
    create_csv_for_algorithm "./heuristics/extra_mileage_two_opt"
}

metaheuristics() {
    create_csv_for_algorithm "./metaheuristics/genetic_algorithm"
    create_csv_for_algorithm "./metaheuristics/simulated_annealing"
    create_csv_for_algorithm "./metaheuristics/tabu_search"
    create_csv_for_algorithm "./metaheuristics/VNS"
}

optimal() {
    create_csv_for_algorithm "./optimal_500_10"
}

exact() {
    create_csv_for_algorithm "./exact_methods" "time"
}

# heuristics
# metaheuristics
# optimal

exact