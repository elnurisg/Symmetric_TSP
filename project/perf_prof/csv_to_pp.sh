#!/bin/bash

heuristics(){
    python3 ./perfprof.py -D , -T 300 -M 1.75 -X "Cost Ratio" ./heuristics/greedy/greedy.csv ./heuristics/greedy/greedy.pdf -P "50 instances for Greedy Heuristic"

    python3 ./perfprof.py -D , -T 300 -M 1.076 -X "Cost Ratio" ./heuristics/greedy_two_opt/greedy_two_opt.csv ./heuristics/greedy_two_opt/greedy_two_opt.pdf -P "50 instances for Greedy Heuristic and 2-OPT"

    python3 ./perfprof.py -D , -T 300 -M 1.03 -X "Cost Ratio" ./heuristics/extra_mileage/extra_mileage.csv ./heuristics/extra_mileage/extra_mileage.pdf -P "50 instances for Extra Mileage Heuristic"

    python3 ./perfprof.py -D , -T 300 -M 1.022 -X "Cost Ratio" ./heuristics/extra_mileage_two_opt/extra_mileage_two_opt.csv ./heuristics/extra_mileage_two_opt/extra_mileage_two_opt.pdf -P "50 instances for Extra Mileage Heuristic and 2-OPT"
}

metaheuristics(){
    # do not perfprof genetic algo bcz it has too many parameters so we do deeper analyze
    # python3 ./perfprof.py -D , -T 600 -M 1.1 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/genetic_algorithm.csv ./metaheuristics/genetic_algorithm/genetic_algorithm.pdf -P "10 instances for Genetic algorithm"

    python3 ./perfprof.py -D , -T 600 -M 1.11 -X "Cost Ratio" ./metaheuristics/simulated_annealing/simulated_annealing.csv ./metaheuristics/simulated_annealing/simulated_annealing.pdf -P "10 instances for Simulated Annealing"

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

genetic_analyze_cut() {
    # Genetic input file path
    input_file="./metaheuristics/genetic_algorithm/genetic_algorithm.csv"
    directory=./metaheuristics/genetic_algorithm/cutting_analyze

    mkdir -p $directory 
    # output file paths
    output_cutting_type_0="$directory/GA_cutting_type_0.csv"
    output_cutting_type_1="$directory/GA_cutting_type_1.csv"
    output_cutting_type_2="$directory/GA_cutting_type_2.csv"

    cut -d, -f 1-2,3,6,9 $input_file > ./temp0.csv
    cut -d, -f 1-2,4,7,10 $input_file > ./temp1.csv
    cut -d, -f 1-2,5,8,11 $input_file > ./temp2.csv

    sed '1 s/10/4/' ./temp0.csv > $output_cutting_type_0
    sed '1 s/10/4/' ./temp1.csv > $output_cutting_type_1
    sed '1 s/10/4/' ./temp2.csv > $output_cutting_type_2

    rm ./temp0.csv
    rm ./temp1.csv
    rm ./temp2.csv

}

genetic_analyze_repair() {
    # Genetic input file path
    input_file="./metaheuristics/genetic_algorithm/genetic_algorithm.csv"
    directory=./metaheuristics/genetic_algorithm/repair_analyze

    mkdir -p $directory 
    # output file paths
    output_repair_mode_0="$directory/GA_repair_mode_0.csv"
    output_repair_mode_1="$directory/GA_repair_mode_1.csv"
    output_repair_mode_2="$directory/GA_repair_mode_2.csv"

    cut -d, -f 1-2,3,4,5 $input_file > ./temp0.csv
    cut -d, -f 1-2,6,7,8 $input_file > ./temp1.csv
    cut -d, -f 1-2,9,10,11 $input_file > ./temp2.csv

    sed '1 s/10/4/' ./temp0.csv > $output_repair_mode_0
    sed '1 s/10/4/' ./temp1.csv > $output_repair_mode_1
    sed '1 s/10/4/' ./temp2.csv > $output_repair_mode_2

    rm ./temp0.csv
    rm ./temp1.csv
    rm ./temp2.csv

}

genetic_more_csv() {
    genetic_analyze_cut
    genetic_analyze_repair
}

further_genetic() {
    genetic_more_csv
    
    python3 ./perfprof.py -D , -T 600 -M 1.17 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/cutting_analyze/GA_cutting_type_0.csv ./metaheuristics/genetic_algorithm/cutting_analyze/GA_cuting_type_0.pdf -P "10 instances for GA with cutting from middle"
    python3 ./perfprof.py -D , -T 600 -M 1.2 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/cutting_analyze/GA_cutting_type_1.csv ./metaheuristics/genetic_algorithm/cutting_analyze/GA_cuting_type_1.pdf -P "10 instances for GA with cutting from random position"
    python3 ./perfprof.py -D , -T 600 -M 1.2 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/cutting_analyze/GA_cutting_type_2.csv ./metaheuristics/genetic_algorithm/cutting_analyze/GA_cutting_type_2.pdf -P "10 instances for GA with cutting from one-fourth of the tour"
 
    python3 ./perfprof.py -D , -T 600 -M 1.18 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/repair_analyze/GA_repair_mode_0.csv ./metaheuristics/genetic_algorithm/repair_analyze/GA_repair_mode_0.pdf -P "10 instances for GA without repair"
    python3 ./perfprof.py -D , -T 600 -M 1.17 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/repair_analyze/GA_repair_mode_1.csv ./metaheuristics/genetic_algorithm/repair_analyze/GA_repair_mode_1.pdf -P "10 instances for GA with repairing the bad genes"
    python3 ./perfprof.py -D , -T 600 -M 1.12 -X "Cost Ratio" ./metaheuristics/genetic_algorithm/repair_analyze/GA_repair_mode_2.csv ./metaheuristics/genetic_algorithm/repair_analyze/GA_repair_mode_2.pdf -P "10 instances for GA with repairing the bad genes + two-opt refining"
}

create_csv_best_methods_heuristics() {

    input_file1="./heuristics/extra_mileage/extra_mileage.csv"
    input_file2="./heuristics/extra_mileage_two_opt/extra_mileage_two_opt.csv"
    input_file3="./heuristics/greedy/greedy.csv"
    input_file4="./heuristics/greedy_two_opt/greedy_two_opt.csv"

    directory=./heuristics

    # output file paths
    output="$directory/best_heur_methods.csv"

    paste -d ',' \
        <(cut -d, -f 1,4 "$input_file1") \
        <(cut -d, -f 4 "$input_file2") \
        <(cut -d, -f 6 "$input_file3") \
        <(cut -d, -f 6 "$input_file4") | \
    awk 'NR==1{$1=4}1' FS=',' OFS=',' > "$output"

}

best_methods_heuristics() {
    create_csv_best_methods_heuristics

    python3 ./perfprof.py -D , -T 300 -M 1.25 -X "Cost Ratio" ./heuristics/best_heur_methods.csv ./heuristics/best_heur_methods.pdf -P "50 instances for Best Heuristics"
}

create_csv_best_methods_metaheuristics() {

    input_file1="./metaheuristics/genetic_algorithm/genetic_algorithm.csv"
    input_file2="./metaheuristics/simulated_annealing/simulated_annealing.csv"
    input_file3="./metaheuristics/tabu_search/tabu_search.csv"
    input_file4="./metaheuristics/VNS/VNS.csv"

    directory=./metaheuristics

    # output file paths
    output="$directory/best_metaheur_methods.csv"

    paste -d ',' \
        <(cut -d, -f 1-2,11 "$input_file1") \
        <(cut -d, -f 5 "$input_file2") \
        <(cut -d, -f 5 "$input_file3") \
        <(cut -d, -f 3 "$input_file4") | \
    awk 'NR==1{$1=5}1' FS=',' OFS=',' > "$output"

}

best_methods_metaheuristics() {
    create_csv_best_methods_metaheuristics

    python3 ./perfprof.py -D , -T 600 -M 1.12 -X "Cost Ratio" ./metaheuristics/best_metaheur_methods.csv ./metaheuristics/best_metaheur_methods.pdf -P "10 instances for Best Metaheuristics"

}

exact_method_with_optimal() {
    modify_csv_with_optimal ./exact_methods/exact_methods.csv
}

exact_method() {
    python3 ./perfprof.py -D , -M 10 -X "Time Ratio" ./exact_methods/exact_methods.csv ./exact_methods/exact_methods.pdf -P "20 instances for Exact methods"
}


# heuristics
# metaheuristics_with_optimal
# metaheuristics

# further_tabu
# further_genetic

# best_methods_heuristics
# best_methods_metaheuristics

# exact_method_with_optimal
exact_method