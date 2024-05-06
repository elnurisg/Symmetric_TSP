import os
import matplotlib.pyplot as plt

def plot_costs(file_path, number_of_points=200):
    with open(file_path, 'r') as file:
        costs_refining = []
        costs_kicks = []   

        line_num = 0
        for line in file:
            if line_num % 2 == 0: # 2-opt
                cost1 = float(line.strip().replace('$$', ''))
                line_num += 1
                costs_refining.append(cost1)
            elif line_num % 2 == 1: # kick
                cost2 = float(line.strip().replace('$$', ''))
                line_num += 1
                costs_kicks.append(cost2)

    skip_count = (2*len(costs_refining)) // number_of_points
    costs = []
    for i in range(line_num//2):
        if i % skip_count == 0:
            costs.append(costs_refining[i])
            costs.append(costs_kicks[i])

    plt.figure(figsize=(12, 4)) 
    plt.plot(costs, marker='o', linestyle='-')
    plt.xlabel('Index')
    plt.ylabel('Cost')
    plt.title('Cost Values')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.splitext(file_path)[0] + '_plot.png')
    plt.close()

number_of_points = int(input('Enter the number of points for the plot: '))
if number_of_points < 50 or number_of_points > 400:
    print('Number of points should be between 50 and 400. Because costs_VNS_15-OPT.txt has only slightly more than 400 nodes. Exiting...')
    exit()

script_folder = os.path.dirname(os.path.abspath(__file__))
txt_files = [file for file in os.listdir(script_folder) if file.endswith('.txt')]

for file in txt_files:
    file_path = os.path.join(script_folder, file)
    
    # Automatically calculate skip_count based on the number of lines in the file
    plot_costs(file_path, number_of_points)
    print(f'Plot saved for {file}')
