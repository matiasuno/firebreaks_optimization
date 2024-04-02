from read_data import *
from optimization import *
from plot import *

## DATA READING
forest_example = "C:/Users/matia/source_win/firebreaks_optimization/sample_test/data/CanadianFBP/Sub40"
results_example = "C:/Users/matia/source_win/firebreaks_optimization/sample_test/Sub40"
harvest_folder = "C:/Users/matia/source_win/firebreaks_optimization/harvest/"

forest_path = forest_example
results_path = results_example
nsims = 50
total_nsims = 1000

params = read_sims(forest_path,results_path,nsims,total_nsims)
NCells,ignitions_points,avail,scar_graphs,bp = params


## OPTIMIZATION MODEL
contador = 1
sol = []
intensities = [0,0.01]
scars = []
l = 1
for i in intensities:
    print("-"*20,"Solving model","-"*20)
    #intensity = 0
    gap = 0
    tlimit = 1800
    w_parameter = 1
    linked = False
    solution = []
    verbose = 0
    
    fo, fb_list, expected_value, lista_aux = model(i,nsims,gap,tlimit,w_parameter,l,params[0:-1],sol,verbose)
    print("alpha=",i,"EV:",expected_value)

    scars.append(lista_aux)
    # Solution visualization
    prob_map = np.random.rand(40, 40)  # Example probability map
    highlight_nodes = fb_list  # Example highlight node coordinates
    plot_squared_graph(prob_map, highlight_nodes)

    #harvested(harvest_folder+'harvest'+str(i)+'.csv',fb_list)
    contador = contador+1