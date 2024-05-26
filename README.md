# Symmetric Traveling Salesman Problem Solver

This C project implements various algorithms for solving the Symmetric Traveling Salesman Problem (TSP), along with performance analysis tools. For more detailed information, please refer to the accompanying [paper.pdf](./paper.pdf) file in repository. The paper provides an in-depth explanation of the project and its analysis.

## Objective

The primary objective of this project is to explore various algorithms and methodologies for efficiently solving the Traveling Salesman Problem (TSP).

## Implemented Algorithms

### Constructive Heuristics:
- Greedy heuristic with GRASP
- Extra-mileage heuristic
- 2-OPT refining heuristic

### Metaheuristics:
- Tabu Search
- VNS (Variable Neighborhood Search)
- Simulated Annealing
- Genetic Algorithm

### Matheuristics:
- Hard fixing
- Local branching

### Exact methods using CPLEX:
- Benders' loop
- Callback method

## Performance Analysis

Performance profiles of the implemented methods are provided in the [perf_prof](./project/perf_prof) folder. These profiles include computational performance metrics, solution quality, and convergence behavior.

## Usage

1. **Installation**: Clone the repository to your local machine.
```sh   
    git clone https://github.com/elnurisg/Symmetric_TSP
```

2. **Compilation**: Use the provided Makefile to compile the project.
```sh
    make
```
3. **Execution**: Run the solver with the desired parameters.
```sh
    ./tsp <parameters>
```

4. **Performance Profiling**: Explore performance profiles in the perf_prof folder for detailed analysis.

## Dependencies

- CPLEX Optimization Studio version 22.1.1
- GNU Plot for visualization

## Documentation

Implemented algorithms are documented with Doxygen for easy reference and understanding.

