import numpy as np
import matplotlib.pyplot as plt

def make_graph():

    loop_matriz = [3,5,7]

    for i in loop_matriz:
        iterations = []
        residuals = []

        with open(f'resultados/resultados_sol_iterativa_dim_{i}.txt', 'r') as file:
            for line in file:
                parts = line.strip().split()
                iter = int(parts[3])
                res = float(parts[7])
                iterations.append(iter)
                residuals.append(res)
        
        plt.figure(figsize=(8, 5))
        plt.plot(iterations, residuals, marker='o', linestyle='-', color='blue')
        plt.title(f'Residuals vs Iterations\nMatrix {i}x{i}')
        plt.xlabel('Iterations')
        plt.ylabel('Residuals')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'graficos/mdim_{i}_wb.png', transparent= False)
        plt.savefig(f'graficos/mdim_{i}_tb.png', transparent= True)

    return


make_graph()