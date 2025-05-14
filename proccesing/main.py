from read_data import read_asc, read_ignition_points, read_messages, extract_points
from optimization import model
from evaluation import c2f_evaluation
import os

fuels = "/home/matias/Documents/firebreaks_optimization/data/CanadianFBP/Sub40/fuels.asc"
forest = "/home/matias/Documents/firebreaks_optimization/data/CanadianFBP/Sub40"
ev_output = "/home/matias/Documents/firebreaks_optimization/results/results_c/"
ignitions = "/home/matias/Documents/firebreaks_optimization/results/Sub40c/IgnitionsHistory/replication.csv"
msg_path = "/home/matias/Documents/firebreaks_optimization/results/Sub40c/Messages/"

header,data,nodos = read_asc(fuels)
scars_graphs = read_messages(msg_path)
ig_points = extract_points(ignitions)

i = 0.01
nsims = 100
gap = 0.01
tlimit = 60*5
lmbda = 1
solution = None
verbose = 1
ev_nsims = 1000

for i in [0,0.01,0.03,0.05]:
    print("i:", i)

    # Solve the optimisation model
    params = [i, nsims, gap, tlimit, lmbda, solution, verbose]
    fo, fb_list, ev, lista_aux = model(params, nodos, ig_points, scars_graphs)
    #print("ev model:", ev)

    # Carry out the C2F evaluation
    #params = [0, ev_nsims, gap, tlimit, lmbda, solution, verbose]
    ev_output_i = f"{ev_output}i_{i}/"
    c2f_evaluation(forest, ev_output_i, fuels, ev_nsims, fb_list)
    #scars2 = read_messages(f'{ev_output}Messages/')
    #points2 = extract_points(f'{ev_output}IgnitionsHistory/replication.csv')
    #fo,fb_list,ev,lista_aux = model(params, nodos, points2, scars2)
    #print("ev eval:", ev)