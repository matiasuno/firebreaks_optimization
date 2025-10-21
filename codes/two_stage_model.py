from gurobipy import GRB
import gurobipy as gp
from operations_csv import get_all_messages, write_treatment_csv
from operations_raster import read_raster, write_treatment_raster
import networkx as nx

def model(params, n_nodos, scar_graphs, lmbda, n_cfuegos=None):

    #forest and simulation data
    intensity, nsims, gap, tlimit = params
    sims = list(range(1,nsims+1))
    cortafuegos = int(n_nodos*intensity)
    nodos = list(range(1,n_nodos+1))

    if n_cfuegos:
        cortafuegos = n_cfuegos

    print("cortafuegos",cortafuegos)

    #optimization parameters
    model = gp.Model()
    model.setParam("OutputFlag", 1)

        
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
    
    lista_aux=[]
    for s in sims:
        suma = sum(x[n,s].x for n in nodos)
        lista_aux.append(suma/n_nodos)
    ev = sum(lista_aux)/len(lista_aux)

    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
    return fo, fb_list, ev, lista_aux


if __name__ == "__main__":
    
    forest = "sub40"
    ruta_base = "/Users/matiasvilches/Documents/F2A/papers/two_stage"
    msg_path = f"{ruta_base}/sims/{forest}/Messages/"
    fuels = f"{ruta_base}/forest/{forest}/fuels.asc"
    
    data, profile, n_nodos = read_raster(fuels)
    scar_graphs = get_all_messages(msg_path)

    intensity = 0.05
    time_limit = 3600
    params = (intensity, 1000, 0.01, time_limit) #intensity, nsims, gap, tlimit
    lmbda = 1

    treatment_output_raster = f"{ruta_base}/results/{forest}/firebreaks_{forest}_i{intensity}.tif"
    treatment_output_csv = f"{ruta_base}/results/{forest}/firebreaks_{forest}_i{intensity}.csv"

    fo, fb_list, ev, lista_aux = model(params, n_nodos, scar_graphs, lmbda, n_cfuegos=None)

    write_treatment_raster(fuels,fb_list,treatment_output_raster)
    write_treatment_csv(treatment_output_csv,fb_list)