from read_data import *

forest_path = ""
results_path = ""
nsims = 5


params = read_sims(forest_path,results_path,nsims)
NCells,ignitions_points,avail,scar_graphs = params