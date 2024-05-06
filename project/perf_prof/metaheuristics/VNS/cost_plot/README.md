# About
This folder contains a Python script that reads cost values from text files, plots the data, and saves the plots as PNG files. Each text file corresponds to a plot of cost values and is generated using the Variable Neighborhood Search (VNS) algorithm with a time limit of 5 seconds and with different types of kicks, which correspond to the names of the files. VNS is based on a Greedy Heuristic with the starting mode set to 2 (trying all initial nodes and choosing the best one) and without GRASP. Instance used: a280.tsp

# Info
USAGE: python3 plot_cost.py
And after you should add the number of points which you want to visualize. As in our case the smallest file has slightly more than 400 point, internally has been decided to keep the range between 100-400. It will skip some points to equalize the number of points in graph so that comparison among methods can be fair. Having 500 or more points make the plot useless as points become unseperable.