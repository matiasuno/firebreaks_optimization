import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def plot_squared_graph(prob_map, highlight_nodes=[]):
    """
    Plots a squared graph colored based on a probability map and highlights
    specific nodes with a circle.

    Parameters:
        prob_map (2D array): Probability map for coloring the squares.
        highlight_nodes (list of integers): List of cell numbers to highlight.
    """
    size = len(prob_map)
    fig, ax = plt.subplots()
    
    # Create grid of squares and color them based on probability map
    for i in range(size):
        for j in range(size):
            square = plt.Rectangle((i, j), 1, 1, fc=(1, 1, 1, prob_map[i][j]), ec='black')
            ax.add_patch(square)
    
    # Highlight specific nodes with a circle
    for node in highlight_nodes:
        x = (node - 1) % size
        y = (node - 1) // size
        circle = Circle((x + 0.5, y + 0.5), 0.4, color='red', fill=False, linewidth=2)
        ax.add_patch(circle)
    
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_aspect('equal', 'box')
    ax.set_xticks(np.arange(0, size, 1))
    ax.set_yticks(np.arange(0, size, 1))
    ax.grid(color='black', linestyle='-', linewidth=1)
    ax.xaxis.tick_top()
    ax.invert_yaxis()
    plt.show()

def adjacency(rows, columns):
    """
    Creates a dictionary representing a grid with the number of rows and columns,
    where each cell contains its number as well as the numbers of its adjacent cells.

    Parameters:
        rows (int): Number of rows in the grid.
        columns (int): Number of columns in the grid.

    Returns:
        dict: A dictionary representing the grid.
    """
    grid_dict = {}
    cell_number = 1

    # Iterate over each cell in the grid
    for row in range(rows):
        for col in range(columns):
            adjacent_cells = []
            if row > 0:
                adjacent_cells.append(cell_number - columns)  # Up
                if col > 0:
                    adjacent_cells.append(cell_number - columns - 1)  # Up Left
                if col < columns - 1:
                    adjacent_cells.append(cell_number - columns + 1)  # Up Right
            if row < rows - 1:
                adjacent_cells.append(cell_number + columns)  # Down
                if col > 0:
                    adjacent_cells.append(cell_number + columns - 1)  # Down Left
                if col < columns - 1:
                    adjacent_cells.append(cell_number + columns + 1)  # Down Right
            if col > 0:
                adjacent_cells.append(cell_number - 1)  # Left
            if col < columns - 1:
                adjacent_cells.append(cell_number + 1)  # Right

            grid_dict[cell_number] = adjacent_cells
            cell_number += 1

    return grid_dict
