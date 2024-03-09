import gurobipy as gp
from gurobipy import GRB

def model(intensity,nsims,gap,tlimit,w_parameter,params,verbose):
    
    #forest and simulation data
    NCells,ignitions_points,nodos,scar_graphs = params
    sims = list(range(1,nsims+1))
    cortafuegos = int(NCells*intensity)

    #optimization parameters
    model = gp.Model()
    model.setParam("OutputFlag", verbose)
    
    #model variables
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    
    #model objective function
    f_esperanza = gp.quicksum(x[n,s]*w_parameter for n in nodos for s in sims)/nsims
    model.setObjective(f_esperanza, GRB.MINIMIZE)
    
    #firebreak intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    #fire spread constraint
    for s in sims:
        H = scar_graphs[s-1]      
        for n in list(H.nodes):
            nbrs = list(H.neighbors(n))
            for nbr in nbrs:
                model.addConstr(x[n,s] <= x[nbr,s]+y[nbr])
    
    #starting points constraint            
    for s in sims:
        point = ignitions_points[s-1]
        model.addConstr(x[point,s] == 1)
    
    #extra options
    model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    
    #model running
    model.optimize()
    
    #results proccesing
    gap = model.MIPGap
    gap = round(gap,3)
    incumbent = round(model.ObjVal,5)
    bestbd = round(model.ObjBound,5)
    tiempo = round(model.Runtime)
    fo = round(model.ObjVal,3)
    
    lista_aux=[]
    for s in sims:
        suma = sum(x[n,s].x for n in nodos)
        lista_aux.append(suma)
        
    ev = sum(lista_aux)/len(lista_aux)
    #lista_aux.sort(reverse=True)
        
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
      
    
    return fo, fb_list, ev, lista_aux