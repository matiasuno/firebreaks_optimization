from read_data import *
from optimization import *
from plot import *


## DATA READING
forest_example = "C:/Users/matia/Documents/source_win/firebreaks_optimization/sample_test/data/CanadianFBP/400cellsC1"
results_example = "C:/Users/matia/Documents/source_win/firebreaks_optimization/sample_test/results/400cellsC1"

forest_path = forest_example
results_path = results_example
nsims = 10

params = read_sims(forest_path,results_path,nsims)
NCells,ignitions_points,avail,scar_graphs = params
adjacency_matrix = adjacency(20,20)
print(adjacency_matrix[1])


## OPTIMIZATION MODEL
i = [0.03]
for intensity in i:
    print("-"*20,"Solving model","-"*20)
    #intensity = 0
    gap = 0
    tlimit = 1800
    w_parameter = 1
    verbose = 0
    fo, fb_list, expected_value, lista_aux = model(intensity,nsims,gap,tlimit,w_parameter,params,verbose)
    print("fo: ",fo)
    print("ev: ",expected_value)
    print("firebreaks: ",fb_list)
    print("total burned forest by scenario: ", lista_aux)

    # Example usage
    prob_map = np.random.rand(20, 20)  # Example probability map
    highlight_nodes = fb_list  # Example highlight node coordinates
    plot_squared_graph(prob_map, highlight_nodes)