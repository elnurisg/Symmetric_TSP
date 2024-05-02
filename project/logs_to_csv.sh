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
    create_csv_for_algorithm "./perf_prof/heuristics/greedy"
    python3 ./perf_prof/perfprof.py -D , -T 300 -M 2 -X "Cost Ratio" ./perf_prof/heuristics/greedy/greedy.csv ./perf_prof/heuristics/greedy/greedy.pdf -P "50 instances for Greedy Heuristic"

    create_csv_for_algorithm "./perf_prof/heuristics/greedy_two_opt"
    python3 ./perf_prof/perfprof.py -D , -T 300 -M 2 -X "Cost Ratio" ./perf_prof/heuristics/greedy_two_opt/greedy_two_opt.csv ./perf_prof/heuristics/greedy_two_opt/greedy_two_opt.pdf -P "50 instances for Greedy Heuristic and 2-OPT"

    create_csv_for_algorithm "./perf_prof/heuristics/extra_mileage"
    python3 ./perf_prof/perfprof.py -D , -T 300 -M 2 -X "Cost Ratio" ./perf_prof/heuristics/extra_mileage/extra_mileage.csv ./perf_prof/heuristics/extra_mileage/extra_mileage.pdf -P "50 instances for Extra Mileage Heuristic"

    create_csv_for_algorithm "./perf_prof/heuristics/extra_mileage_two_opt"
    python3 ./perf_prof/perfprof.py -D , -T 300 -M 2 -X "Cost Ratio" ./perf_prof/heuristics/extra_mileage_two_opt/extra_mileage_two_opt.csv ./perf_prof/heuristics/extra_mileage_two_opt/extra_mileage_two_opt.pdf -P "50 instances for Extra Mileage Heuristic and 2-OPT"
}

metaheuristics() {

    # Call the function with the directory path
    create_csv_for_algorithm "./perf_prof/metaheuristics/genetic_algorithm"
    python3 ./perf_prof/perfprof.py -D , -T 600 -M 2 ./perf_prof/metaheuristics/genetic_algorithm/genetic_algorithm.csv ./perf_prof/metaheuristics/genetic_algorithm/genetic_algorithm.pdf -P "50 instances for Genetic algorithm"

    create_csv_for_algorithm "./perf_prof/metaheuristics/simulated_annealing"
    python3 ./perf_prof/perfprof.py -D , -T 600 -M 2 ./perf_prof/metaheuristics/simulated_annealing/simulated_annealing.csv ./perf_prof/metaheuristics/simulated_annealing/simulated_annealing.pdf -P "50 instances for Simulated Annealing"

    create_csv_for_algorithm "./perf_prof/metaheuristics/tabu_search"
    python3 ./perf_prof/perfprof.py -D , -T 600 -M 2 ./perf_prof/metaheuristics/tabu_search/tabu_search.csv ./perf_prof/metaheuristics/tabu_search/tabu_search.pdf -P "50 instances for Tabu Search"

    create_csv_for_algorithm "./perf_prof/metaheuristics/VNS"
    python3 ./perf_prof/perfprof.py -D , -T 600 -M 2 ./perf_prof/metaheuristics/VNS/VNS.csv ./perf_prof/metaheuristics/VNS/VNS.pdf -P "50 instances for Variable Neighborhood Search"
}


# heuristics
# metaheuristics