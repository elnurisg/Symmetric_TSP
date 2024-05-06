#!/bin/bash

heuristics(){
    python3 ./perfprof.py -D , -T 300 -M 1.75 -X "Cost Ratio" ./heuristics/greedy/greedy.csv ./heuristics/greedy/greedy.pdf -P "50 instances for Greedy Heuristic"

    python3 ./perfprof.py -D , -T 300 -M 1.076 -X "Cost Ratio" ./heuristics/greedy_two_opt/greedy_two_opt.csv ./heuristics/greedy_two_opt/greedy_two_opt.pdf -P "50 instances for Greedy Heuristic and 2-OPT"

    python3 ./perfprof.py -D , -T 300 -M 1.03 -X "Cost Ratio" ./heuristics/extra_mileage/extra_mileage.csv ./heuristics/extra_mileage/extra_mileage.pdf -P "50 instances for Extra Mileage Heuristic"

    python3 ./perfprof.py -D , -T 300 -M 1.022 -X "Cost Ratio" ./heuristics/extra_mileage_two_opt/extra_mileage_two_opt.csv ./heuristics/extra_mileage_two_opt/extra_mileage_two_opt.pdf -P "50 instances for Extra Mileage Heuristic and 2-OPT"
}

metaheuristics(){
    python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/genetic_algorithm.csv ./metaheuristics/genetic_algorithm/genetic_algorithm.pdf -P "10 instances for Genetic algorithm"

    python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/simulated_annealing/simulated_annealing.csv ./metaheuristics/simulated_annealing/simulated_annealing.pdf -P "10 instances for Simulated Annealing"

    python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/tabu_search/tabu_search.csv ./metaheuristics/tabu_search/tabu_search.pdf -P "10 instances for Tabu Search"

    python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/VNS/VNS.csv ./metaheuristics/VNS/VNS.pdf -P "10 instances for Variable Neighborhood Search"
}

## updates the csv file of metaheuristics with adding also the optimal values
modify_csv_with_optimal() {
    input_file="$1"
    cut -d ',' -f2- "$input_file" | paste -d ',' ./optimal_500_10/optimal_500_10.csv - > temp.csv
    total_columns=$(head -n 1 "$input_file" | awk -F ',' '{print NF}')
    awk -F',' 'NR==1{$1="'$total_columns'"}1' OFS=',' temp.csv > temp && mv temp "$input_file"
    rm temp.csv
}

metaheuristics_with_optimal() {
    modify_csv_with_optimal ./metaheuristics/genetic_algorithm/genetic_algorithm.csv
    modify_csv_with_optimal ./metaheuristics/simulated_annealing/simulated_annealing.csv
    modify_csv_with_optimal ./metaheuristics/tabu_search/tabu_search.csv
    modify_csv_with_optimal ./metaheuristics/VNS/VNS.csv
}

tabu_more_csv() {

    # Tabu input file path
    input_file="./metaheuristics/tabu_search/tabu_search.csv"

    # output file paths
    output_aspiration_0="./metaheuristics/tabu_search/tabu_with_aspiration_0.csv"
    output_aspiration_1="./metaheuristics/tabu_search/tabu_with_aspiration_1.csv"

    cut -d, -f 1-2,3,5,7 $input_file > ./temp0.csv
    cut -d, -f 1-2,4,6,8 $input_file > ./temp1.csv

    sed '1 s/7/4/' ./temp0.csv > $output_aspiration_0
    sed '1 s/7/4/' ./temp1.csv > $output_aspiration_1

    rm ./temp0.csv
    rm ./temp1.csv
}

# To perform the performance profile for the tabu search with aspiration 0 and 1
further_tabu() {
    tabu_more_csv

    python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/tabu_search/tabu_with_aspiration_0.csv ./metaheuristics/tabu_search/tabu_with_aspiration_0.pdf -P "10 instances for Tabu Search with aspiration criteria"
    python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/tabu_search/tabu_with_aspiration_1.csv ./metaheuristics/tabu_search/tabu_with_aspiration_1.pdf -P "10 instances for Tabu Search without aspiration criteria"
}

# heuristics
# metaheuristics_with_optimal
# metaheuristics

further_tabu