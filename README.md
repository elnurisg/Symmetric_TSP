This C project implements various algorithms and their performance analysis for solving the Traveling Salesman Problem (TSP), covering heuristics, metaheuristics, matheuristics, and constructive heuristics for efficient solutions.

**Constructive Heuristics:**
- Greedy heuristic with GRASP
- Extra-mileage heuristic
- 2-OPT refining heuristic

**Metaheuristics:**
- Tabu Search
- VNS (Variable Neighborhood Search)
- Simulated Annealing
- Genetic Algorithm

**Matheuristics:**
- Hard fixing
- Local branching

**Also implemented:**
- Basic TSP model in CPLEX (without SECs)
- Benders' "loop" method with SEC separation for integer solutions
- Patching heuristic
- Branch-and-Cut implementation in CPLEX through callbacks
- Posting heuristic solutions and adding user cuts for fractional solutions
