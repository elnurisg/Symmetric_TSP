import os
import matplotlib.pyplot as plt

def plot_costs(file_path):
    with open(file_path, 'r') as file:
        costs = [float(line.strip().replace('$$', '')) for line in file]

    plt.figure(figsize=(12, 4)) 
    plt.plot(costs, marker='o', linestyle='-')
    plt.xlabel('Index')
    plt.ylabel('Cost')
    plt.title('Cost Values')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.splitext(file_path)[0] + '_plot.png')
    plt.close()


script_folder = os.path.dirname(os.path.abspath(__file__))
txt_files = [file for file in os.listdir(script_folder) if file.endswith('.txt')]

for file in txt_files:
    file_path = os.path.join(script_folder, file)
    plot_costs(file_path)
    print(f'Plot saved for {file}')
