from gurobipy import GRB
import gurobipy as gp
import networkx as nx
import csv
import os

def save_results(results, header, output_path="results.csv"):
    file_exists = os.path.isfile(output_path)

    # open file in append mode
    with open(output_path, mode="a", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        # write header only if file does not exist yet
        if not file_exists:
            writer.writerow(header)
        writer.writerows(results)

def get_wc(lista_aux, alpha=0.9):
    """Calculate mean of worst alpha scenarios."""
    lista_aux.sort()
    n_last = int(len(lista_aux)*alpha)
    wc = sum(lista_aux[-n_last:])/(len(lista_aux)*(1-alpha))
    return round(wc,3)


def mip_model(params, n_nodos, scar_graphs, lmbda, n_cfuegos=None, output_flag=False):

    #forest and simulation data
    intensity, nsims, gap, tlimit = params
    sims = list(range(1,nsims+1))
    cortafuegos = int(n_nodos*intensity)
    nodos = list(range(1,n_nodos+1))

    if n_cfuegos:
        cortafuegos = n_cfuegos

    #print("cortafuegos",cortafuegos)

    #optimization parameters
    model = gp.Model()
    model.setParam("OutputFlag", output_flag)

    #model variables
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    eta = model.addVars(sims, vtype=GRB.CONTINUOUS, lb=0)
    phi = model.addVar(vtype=GRB.CONTINUOUS, lb =0)
    
    #model objective function
    f_esperanza = gp.quicksum(x[n,s] for n in nodos for s in sims)/nsims
    f_cvar = phi+(1/(1-0.9))*gp.quicksum(eta[s] for s in sims)/nsims
    model.setObjective(lmbda*f_esperanza+(1-lmbda)*f_cvar, GRB.MINIMIZE)
    
    #firebreak intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    ignitions_points = []
    #fire spread constraint
    for s in sims:
        H = scar_graphs[s-1]
        h_nodes = list(H.nodes())
        if h_nodes:
            ignitions_points.append(h_nodes[0])      
            for n in h_nodes:
                nbrs = list(H.neighbors(n))
                for nbr in nbrs:
                    model.addConstr(x[n,s] <= x[nbr,s]+y[nbr])
        else:
            ignitions_points.append(False)
    
    #starting points constraint            
    for s in sims:
        point = ignitions_points[s-1]
        if point:
            model.addConstr(x[point,s] == 1)

    #CVaR constraint
    for s in sims:
        model.addConstr(gp.quicksum(x[n,s] for n in nodos)-phi <= eta[s])

    #extra options
    model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    
    #model running
    model.optimize()
    
    #results proccesing
    gap = model.MIPGap
    gap = round(gap,3)
    fo = round(model.ObjVal,3)
    time = model.Runtime
    time = round(time/60,3) #minutes
    
    lista_aux=[]
    for s in sims:
        suma = sum(x[n,s].x for n in nodos)
        lista_aux.append(suma/n_nodos)
    ev = sum(lista_aux)/len(lista_aux)
    ev = round(ev,3)

    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
    
    return fo, fb_list, ev, lista_aux, gap, time