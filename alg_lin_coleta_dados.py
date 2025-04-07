import numpy as np
import matplotlib.pyplot as plt

def make_graph():

    iterations = []
    residuals = []
    with open('resultados_sol_iterativa.txt', 'r') as file:
        for line in file:
            parts = line.strip().split()
            iter = int(parts[3])
            res = float(parts[7])
            iterations.append(iter)
            residuals.append(res)
    
    plt.figure(figsize=(8, 5))
    plt.plot(iterations, residuals, marker='o', linestyle='-', color='blue')
    plt.title('Residuals vs Iterations')
    plt.xlabel('Iterations')
    plt.ylabel('Residuals')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return


make_graph()