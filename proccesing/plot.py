import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from statistics import mean

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


def boxplot(scars):
    # Creating a boxplot       
    # Colors
    BG_WHITE = "#ffffff"
    GREY_LIGHT = "#b4aea9"
    GREY50 = "#7F7F7F"
    BLUE_DARK = "#1B2838"
    BLUE = "#2a475e"
    BLACK = "#282724"
    GREY_DARK = "#747473"
    RED_DARK = "#850e00"

    # Colors taken from Dark2 palette in RColorBrewer R library
    COLOR_SCALE = ["#1B9E77", "#D95F02", "#7570B3"]

    POSITIONS = [0, 1,2]
    # Horizontal lines
    HLINES = [0.2, 0.4, 0.6, 0.8]

    fig, ax = plt.subplots(figsize= (10, 8))

    for h in HLINES:
        ax.axhline(h, color=GREY50, ls=(0, (5, 5)), alpha=0.8, zorder=0)

    fig.patch.set_facecolor(BG_WHITE)
    ax.set_facecolor(BG_WHITE)

    # Add boxplots ---------------------------------------------------
    # Note that properties about the median and the box are passed
    # as dictionaries.

    medianprops = dict(
        linewidth=4, 
        color=GREY_DARK,
        solid_capstyle="butt"
    )
    boxprops = dict(
        linewidth=2, 
        color=GREY_DARK
    )

    ax.boxplot(
        scars,
        positions=POSITIONS, 
        showfliers = False, # Do not show the outliers beyond the caps.
        showcaps = False,   # Do not show the caps
        medianprops = medianprops,
        whiskerprops = boxprops,
        boxprops = boxprops
    )

    means = [mean(y) for y in scars]
    for i, mean in enumerate(means):
        # Add dot representing the mean
        ax.scatter(i, mean, s=100, color=RED_DARK, zorder=3)
        
        # Add line conecting mean value and its label
        ax.plot([i, i + 0.25], [mean, mean], ls="dashdot", color="black", zorder=3)
        
        # Add mean value label.
        ax.text(
            i + 0.25,
            mean,
            r"$\hat{\mu}_{\rm{mean}} = $" + str(round(mean, 2)),
            fontsize=13,
            va="center",
            bbox = dict(
                facecolor="white",
                edgecolor="black",
                boxstyle="round",
                pad=0.15
            ),
            zorder=10 # to make sure the line is on top
        )

    # Displaying the plot
    plt.title('Burn Probability vs. DPV vs. Stochastic MIP  \n (out of sample, n=1000)',size=20)
    ax.set_xticklabels(['BP','DPV','MIP'], size=15, ha="center", ma="center")
    plt.xlabel('Decision models', size=15)
    plt.ylabel('Percentage of burned forest',size=15,labelpad=10)
    plt.show()
