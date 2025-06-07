import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def visualize_results(filepath="results.txt"):
    """
    Reads simulation results from a CSV file and creates heatmap plots
    for pressure and saturation.
    """
    print(f"Reading data from {filepath}...")
    try:
        data = pd.read_csv(filepath)
    except FileNotFoundError:
        print(f"Error: The file {filepath} was not found.")
        print("Please run the C++ simulator first to generate the results file.")
        return

    # Определяем размеры сетки из данных
    width = data['i'].max() + 1
    height = data['j'].max() + 1
    
    print(f"Grid size detected: {width}x{height}")

    # Преобразуем данные из формата [i, j, value] в 2D матрицу
    pressure_grid = data.pivot(index='j', columns='i', values='pressure').values
    saturation_grid = data.pivot(index='j', columns='i', values='saturation').values

    # --- Создание графика для Давления ---
    plt.figure(figsize=(10, 8))
    plt.imshow(pressure_grid, origin='lower', cmap='viridis')
    plt.colorbar(label='Pressure (units)')
    plt.title('Pressure Distribution')
    plt.xlabel('Grid Cell (i)')
    plt.ylabel('Grid Cell (j)')
    pressure_plot_path = "pressure_map.png"
    plt.savefig(pressure_plot_path)
    print(f"Pressure map saved to {pressure_plot_path}")
    plt.close()

    # --- Создание графика для Насыщенности ---
    plt.figure(figsize=(10, 8))
    plt.imshow(saturation_grid, origin='lower', cmap='jet', vmin=0, vmax=1)
    plt.colorbar(label='Water Saturation (fraction)')
    plt.title('Saturation Distribution')
    plt.xlabel('Grid Cell (i)')
    plt.ylabel('Grid Cell (j)')
    saturation_plot_path = "saturation_map.png"
    plt.savefig(saturation_plot_path)
    print(f"Saturation map saved to {saturation_plot_path}")
    plt.close()

if __name__ == "__main__":
    visualize_results() 