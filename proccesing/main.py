from read_data import *

forest_example = "C:/Users/matia/source_win/firebreaks_optimization/sample_test/data/CanadianFBP/400cellsC1"
results_example = "C:/Users/matia/source_win/firebreaks_optimization/sample_test/results/400cellsC1"

forest_path = forest_example
results_path = results_example
nsims = 5

params = read_sims(forest_path,results_path,nsims)
NCells,ignitions_points,avail,scar_graphs = params