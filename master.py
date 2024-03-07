#%%
import ReadDataPrometheus
import os
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import networkx as nx
import random as rd
import gurobipy as gp
from gurobipy import GRB
from winsound import Beep
import heapq
from operator import itemgetter
import pandas as pd
import csv
import builtins
import math
import subprocess
import itertools

#%% IMPORTAR DATOS C2F

def read_sims(data,ruta,cantidad,file_list,cmd_bool):

    Folder = os.getcwd()
    Folder = data
    FBPlookup = Folder + '/fbp_lookup_table.csv'
    ForestFile = Folder + '/Forest.asc'
    FBPDict, ColorsDict = ReadDataPrometheus.Dictionary(FBPlookup)
    CellsGrid3, CellsGrid4, Rows, Cols, AdjCells, CoordCells, CellSize = ReadDataPrometheus.ForestGrid(ForestFile,FBPDict)
    
    NCells = Rows * Cols
    m = Rows
    n = Cols
    Colors = []
    
    for i in range(NCells):
        if str(CellsGrid3[i]) not in ColorsDict.keys():
            Colors.append((1.0,1.0,1.0,1.0))
        if str(CellsGrid3[i]) in ColorsDict.keys():
            Colors.append(ColorsDict[str(CellsGrid3[i])])
    
    AvailSet = set()
    NonAvailSet = set()
    
    for i in range(NCells):
        if CellsGrid4[i] != 'NF':
            AvailSet.add(i+1)
        else:
            NonAvailSet.add(i+1)
    
    setDir = ['S', 'SE', 'E', 'NE', 'N', 'NW', 'W', 'SW']
    aux = set([])
    
    Adjacents = {}
    for k in range(NCells):
        aux = set([])
        for i in setDir:
            if AdjCells[k][i] != None:
                if AdjCells[k][i][0] in AvailSet :
                    aux = aux | set(AdjCells[k][i])
        Adjacents[k + 1] = aux & AvailSet
    
    FAG = nx.Graph()
    FAG.add_nodes_from(list(AvailSet))
    
    coord_pos = dict()
    for i in AvailSet:
        coord_pos[i] = CoordCells[i-1]
    
    ColorsAvail = {}
    for i in AvailSet:
        FAG.add_edges_from([(i,j) for j in Adjacents[i]])
        ColorsAvail[i] = Colors[i-1]
    
    Folder = str(ruta)
    
    cmdoutput1 = os.listdir(Folder + '/Messages')
    if ".DS_Store" in cmdoutput1:
        idx = cmdoutput1.index('.DS_Store')
        del cmdoutput1[idx]
    if "desktop.ini" in cmdoutput1:
        idx = cmdoutput1.index('desktop.ini')
        del cmdoutput1[idx]
    pathmessages = Folder + '/Messages/'
    
    MDG = nx.MultiDiGraph()
    
    
    avail = list(AvailSet)
    st_points = []
    sim_edges = []
    burned = dict.fromkeys(AvailSet,0)
    
    contador_sims = 0
    
    histograma = dict.fromkeys(range(0,cantidad))
    
    dpv_count = dict.fromkeys(AvailSet,0)
    if cmd_bool == True:
        file_list = cmdoutput1
    
    for file in file_list:
        
    #for f in range(1, cantidad+1):
    
        #file = str(f).zfill(2) + '.csv'
        #file = os.path.join('MessagesFile' + file)
        
        
        
        if contador_sims < cantidad:
            
            #print(file)
            
            H = nx.read_edgelist(path = pathmessages + file,
                                 delimiter=',',
                                 create_using = nx.DiGraph(),
                                 nodetype = int,
                                 data = [('time', float), ('ros', float)])
            try:
                st_points.append(list(H.edges)[0][0])
            except IndexError:
                st_points.append('na')
                continue
                
        
            MDG.add_weighted_edges_from(H.edges(data='time'), weight='ros')
            contador = dict.fromkeys(list(H.nodes),0)
            for e in H.edges:
                for n in list(H.nodes):
                    if n == e[1] and contador[n] == 0:
                        contador[n]=contador[n]+1
                        burned[e[1]]=burned[e[1]]+1
                        
            area = len(list(H.nodes))
            histograma[contador_sims] = area
            
            dpv_aux = dpv(H)
                  
            for key in dpv_aux:
                dpv_count[key] = dpv_count[key] + dpv_aux[key]
                        
        contador_sims = contador_sims+1
    
    for point in st_points:
        try:
            burned[point] = burned[point]+1
        except:
            continue
    edges=list(MDG.edges)
    
    
    nodos = list(MDG.nodes())
    nodos = list(avail)

    
    params = [NCells,file_list,st_points,nodos,AvailSet,pathmessages,FAG,ColorsAvail,coord_pos,burned,file_list,histograma,CoordCells]
        
    return params,MDG,dpv_count


#%% MODELO 1

def modelo1(intensidad,warmstart,gap,tlimit,xvartype,yvartype,w,params):
    
    NCells,cmdoutput1,st_points,nodos,AvailSet,pathmessages = params
    
    sims = list(range(len(st_points)))
    n_sims = len(sims)
    for i in range(len(sims)):
        sims[i] = sims[i]+1
    
    cortafuegos = int(NCells*intensidad)
    
    model = gp.Model()
    model.setParam("OutputFlag", 0)
    
    x = model.addVars(nodos, sims, vtype=xvartype)
    y = model.addVars(nodos, vtype=yvartype)
    
    model.setObjective(gp.quicksum(x[n,s]*w[n] for n in nodos for s in sims)/n_sims, GRB.MINIMIZE)
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    simulacion = 0
    contador_sims = 0
    for k in cmdoutput1:
        if contador_sims < n_sims:
            simulacion = simulacion +1
            # Leemos la marca
            H = nx.read_edgelist(path = pathmessages + k,
                                 delimiter=',',
                                 create_using = nx.DiGraph(),
                                 nodetype = int,
                                 data = [('time', float), ('ros', float)])        
            for n in list(H.nodes):
                nbrs = list(H.neighbors(n))
                for nbr in nbrs:
                    model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
        contador_sims = contador_sims +1
                
    for s in sims:
        point = st_points[s-1]
        model.addConstr(x[point,s] == 1)
            
    for n in warmstart:
        y[n].Start = 1
            
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()
     
    if yvartype == GRB.BINARY:
        gap = model.MIPGap
        gap = round(gap,3)
        incumbent = round(model.ObjVal,5)
        bestbd = round(model.ObjBound,5)
    else:
        gap = 'n/a'
        incumbent = 'n/a'
        bestbd = 'n/a'
        
    tiempo = round(model.Runtime)
    
    #modelo 1
    fo1 = sum(x[n,s].x*w[n] for n in nodos for s in sims)/n_sims
      
    #modelo 2
    eta_aux=[]
    for n in nodos:
        eta_aux.append(sum(x[n,s].x for s in sims))
    fo2 = max(eta_aux)/n_sims
      
    #modelo 3
    eta_aux2=[]
    for s in sims:
        eta_aux2.append(sum(x[n,s].x for n in nodos))
    fo3 = max(eta_aux2)/NCells
      
    #modelo 4
    fo4 = sum(y[n].x for n in nodos)
    
    b = 'n/a'
    lmbda = 'n/a'
    
    lista_aux=[] 
    for s in sims:
        lista_aux.append(sum(x[n,s].x for n in nodos))
    
    lista_aux = sort(lista_aux)
    
    peores_casos1 = int(round((1-0.9)*n_sims,0))
    peores_casos2 = int(round((1-0.95)*n_sims,0))
    
    peores_casos1 = 1
    peores_casos2 = 1
    
    prom_pcb1= 0
    prom_pcb2= 0
    for i in range(peores_casos1):
        prom_pcb1 = prom_pcb1 + lista_aux[-i-1]
    prom_pcb1 = prom_pcb1/peores_casos1
    
    for i in range(peores_casos2):
        prom_pcb2 = prom_pcb2 + lista_aux[-i-1]
    prom_pcb2 = prom_pcb2/peores_casos2
    
   
    
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
            
    burned_map = dict.fromkeys(AvailSet,0)
    for s in sims:
        for n in nodos:
            if x[n,s].x > 0:
                burned_map[n] = burned_map[n]+1
                
    titulo = 'm1_i'+str(intensidad)+'s'
    results = ['modelo1',n_sims,intensidad,fo1,fo2,fo3,tiempo,gap,lmbda,b,prom_pcb1,prom_pcb2,xvartype,yvartype,incumbent,bestbd] 
    titulo2 = 'm1'+' '+str(contador_cfuegos)+'cf '+str(n_sims)+'sims'+' fo:'+str(round(model.ObjVal,2))
    titulos = [titulo,titulo2]
    
            
    d = dict.fromkeys(nodos,0)
    for n in nodos:
        for s in sims:
            d[n] = d[n]+x[n,s].x
    
    
    return results, fb_list, burned_map,titulos,st_points

#%% MODELO 2       

def modelo2(intensidad,warmstart,gap,tlimit,xvartype,yvartype,w,params):
    
    NCells,cmdoutput1,st_points,nodos,AvailSet,pathmessages = params

    sims = list(range(len(st_points)))
    n_sims = len(sims)
    for i in range(len(sims)):
        sims[i] = sims[i]+1
        
    cortafuegos = int(NCells*intensidad)
    
    model = gp.Model()
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    eta = model.addVar()
    
    model.setObjective(eta, GRB.MINIMIZE)
    
    #firebreaks intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) <= cortafuegos)
    
    #restriccion que genera la cicatriz
    simulacion = 0
    contador_sims = 0
    for k in cmdoutput1:
        if contador_sims < n_sims:
            simulacion = simulacion +1
            # Leemos la marca
            H = nx.read_edgelist(path = pathmessages + k,
                                 delimiter=',',
                                 create_using = nx.DiGraph(),
                                 nodetype = int,
                                 data = [('time', float), ('ros', float)])        
            for n in list(H.nodes):
                nbrs = list(H.neighbors(n))
                for nbr in nbrs:
                    model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
        contador_sims = contador_sims +1
    
    #restriccion starting points            
    for s in sims:
        point = st_points[s-1]
        model.addConstr(x[point,s] == 1)
    
    #restriccion de burn probability
    for n in nodos:
        model.addConstr(gp.quicksum(x[n,s] for s in sims) <= eta)
    
    for n in warmstart:
        y[n].Start = 1
            
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()
     
    if yvartype == GRB.BINARY:
        gap = model.MIPGap
        gap = round(gap,3)
        incumbent = round(model.ObjVal,5)
        bestbd = round(model.ObjBound,5)
    else:
        gap = 'n/a'
        incumbent = 'n/a'
        bestbd = 'n/a'
        
    tiempo = round(model.Runtime)
    
    #modelo 1
    fo1 = sum(x[n,s].x*w[n] for n in nodos for s in sims)/n_sims
      
    #modelo 2
    eta_aux=[]
    for n in nodos:
        eta_aux.append(sum(x[n,s].x for s in sims))
    fo2 = max(eta_aux)/n_sims
      
    #modelo 3
    eta_aux2=[]
    for s in sims:
        eta_aux2.append(sum(x[n,s].x for n in nodos))
    fo3 = max(eta_aux2)/NCells
      
    #modelo 4
    fo4 = sum(y[n].x for n in nodos)
    
    b = 'n/a'
    lmbda = 'n/a'
    
    lista_aux=[] 
    for s in sims:
        lista_aux.append(sum(x[n,s].x for n in nodos))
    
    lista_aux = sort(lista_aux)
    
    peores_casos1 = int(round((1-0.9)*n_sims,0))
    peores_casos2 = int(round((1-0.95)*n_sims,0))
    
    prom_pcb1= 0
    prom_pcb2= 0
    for i in range(peores_casos1):
        prom_pcb1 = prom_pcb1 + lista_aux[-i-1]
    prom_pcb1 = prom_pcb1/peores_casos1
    
    for i in range(peores_casos2):
        prom_pcb2 = prom_pcb2 + lista_aux[-i-1]
    prom_pcb2 = prom_pcb2/peores_casos2
    
    
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
            
    burned_map = dict.fromkeys(AvailSet,0)
    for s in sims:
        for n in nodos:
            if x[n,s].x > 0:
                burned_map[n] = burned_map[n]+1
                
    titulo = 'm2_intensidad'+str(intensidad)
    results = ['modelo2',n_sims,intensidad,fo1,fo2,fo3,tiempo,gap,lmbda,b,prom_pcb1,prom_pcb2,xvartype,yvartype,incumbent,bestbd] 
    titulo2 = 'm2'+' '+str(contador_cfuegos)+'cf '+str(n_sims)+'sims'+' fo:'+str(round(model.ObjVal,2))
    titulos = [titulo,titulo2]
    
    return results, fb_list, burned_map,titulos

#%% MODELO 3 

def modelo3(intensidad,warmstart,gap,tlimit,xvartype,yvartype,w,params):
    
    NCells,cmdoutput1,st_points,nodos,AvailSet,pathmessages = params

    sims = list(range(len(st_points)))
    n_sims = len(sims)
    for i in range(len(sims)):
        sims[i] = sims[i]+1
        
    cortafuegos = int(NCells*intensidad)

    model = gp.Model()
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    eta = model.addVar()
    
    model.setObjective(eta, GRB.MINIMIZE)
    
    model.addConstr(gp.quicksum(y[n] for n in nodos) <= cortafuegos)
    
    simulacion = 0
    contador_sims = 0
    for k in cmdoutput1:
        if contador_sims < n_sims:
            simulacion = simulacion +1
            # Leemos la marca
            H = nx.read_edgelist(path = pathmessages + k,
                                 delimiter=',',
                                 create_using = nx.DiGraph(),
                                 nodetype = int,
                                 data = [('time', float), ('ros', float)])        
            for n in list(H.nodes):
                nbrs = list(H.neighbors(n))
                for nbr in nbrs:
                    model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
        contador_sims = contador_sims +1
    
    for s in sims:
        point = st_points[s-1]
        model.addConstr(x[point,s] == 1)
        model.addConstr(gp.quicksum(x[n,s] for n in nodos) <= eta)
    
    for n in warmstart:
        y[n].Start = 1
            
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()
     
    if yvartype == GRB.BINARY:
        gap = model.MIPGap
        gap = round(gap,3)
        incumbent = round(model.ObjVal,5)
        bestbd = round(model.ObjBound,5)
    else:
        gap = 'n/a'
        incumbent = 'n/a'
        bestbd = 'n/a'
        
    tiempo = round(model.Runtime)
    
    #modelo 1
    fo1 = sum(x[n,s].x*w[n] for n in nodos for s in sims)/n_sims
      
    #modelo 2
    eta_aux=[]
    for n in nodos:
        eta_aux.append(sum(x[n,s].x for s in sims))
    fo2 = max(eta_aux)/n_sims
      
    #modelo 3
    eta_aux2=[]
    for s in sims:
        eta_aux2.append(sum(x[n,s].x for n in nodos))
    fo3 = max(eta_aux2)/NCells
      
    #modelo 4
    fo4 = sum(y[n].x for n in nodos)
    
    b = 'n/a'
    lmbda = 'n/a'
    
    lista_aux=[] 
    for s in sims:
        lista_aux.append(sum(x[n,s].x for n in nodos))
    
    lista_aux = sort(lista_aux)
    
    peores_casos1 = int(round((1-0.9)*n_sims,0))
    peores_casos2 = int(round((1-0.95)*n_sims,0))
    
    prom_pcb1= 0
    prom_pcb2= 0
    for i in range(peores_casos1):
        prom_pcb1 = prom_pcb1 + lista_aux[-i-1]
    prom_pcb1 = prom_pcb1/peores_casos1
    
    for i in range(peores_casos2):
        prom_pcb2 = prom_pcb2 + lista_aux[-i-1]
    prom_pcb2 = prom_pcb2/peores_casos2
        
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
            
    burned_map = dict.fromkeys(AvailSet,0)
    for s in sims:
        for n in nodos:
            if x[n,s].x > 0:
                burned_map[n] = burned_map[n]+1
                
    titulo = 'm3_intensidad'+str(intensidad)
    results = ['modelo3',n_sims,intensidad,fo1,fo2,fo3,tiempo,gap,lmbda,b,prom_pcb1,prom_pcb2,xvartype,yvartype,incumbent,bestbd] 
    titulo2 = 'm3'+' '+str(contador_cfuegos)+'cf '+str(n_sims)+'sims'+' fo:'+str(round(model.ObjVal,2))
    titulos = [titulo,titulo2]
    
    return results, fb_list, burned_map,titulos

#%% MODELO 4

def modelo4(intensidad,n_sims,warmstart,gap,tlimit,xvartype,yvartype,w,params,b,lmbda,fi,forest,cb,firebreaks,dpv_firebreaks):
    
    #forest = nsims
    
    NCells,cmdoutput1,st_points,nodos,AvailSet,pathmessages = params

    sims = list(range(1,n_sims+1))
        
    celdas = len(AvailSet)    
        
    cortafuegos = int(NCells*intensidad)
    

    model = gp.Model()
    model.setParam("OutputFlag", cb)
    
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    eta = model.addVars(sims, vtype=GRB.CONTINUOUS, lb=0)
    phi = model.addVar(vtype=GRB.CONTINUOUS, lb =0)
    
    f_esperanza = gp.quicksum(x[n,s]*w[n] for n in nodos for s in sims)/n_sims
    f_cvar = phi+(1/(1-b))*gp.quicksum(eta[s] for s in sims)/n_sims
    
    model.setObjective(lmbda*f_esperanza+(1-lmbda)*f_cvar, GRB.MINIMIZE)
    
    #firebreak intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    #fire spread constraint
    simulacion = 0
    contador_sims = 0
    for k in cmdoutput1:
        if contador_sims < n_sims:
            simulacion = simulacion +1
            # Leemos la marca
            H = nx.read_edgelist(path = pathmessages + k,
                                 delimiter=',',
                                 create_using = nx.DiGraph(),
                                 nodetype = int,
                                 data = [('time', float), ('ros', float)])        
            for n in list(H.nodes):
                nbrs = list(H.neighbors(n))
                for nbr in nbrs:
                    model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
        contador_sims = contador_sims +1
    
    #starting points constraint            
    for s in sims:
        point = st_points[s-1]
        model.addConstr(x[point,s] == 1)
    
    for s in sims:
        model.addConstr(gp.quicksum(x[n,s] for n in nodos)-phi <= eta[s])
        
    for n in firebreaks:
        model.addConstr(y[n] == 1)
        
    for n in dpv_firebreaks:
        y[n].Start = 1
        #model.addConstr(y[n] == 1)
    
    for n in warmstart:
        y[n].Start = 1
    
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()
     
    if yvartype == GRB.BINARY:
        gap = model.MIPGap
        gap = round(gap,3)
        incumbent = round(model.ObjVal,5)
        bestbd = round(model.ObjBound,5)
    else:
        gap = 'n/a'
        incumbent = 'n/a'
        bestbd = 'n/a'
        
    tiempo = round(model.Runtime)
    
    
    #modelo 1
    fo1 = sum(x[n,s].x*w[n] for n in nodos for s in sims)/n_sims
      
    #modelo 2
    eta_aux=[]
    for n in nodos:
        eta_aux.append(sum(x[n,s].x for s in sims))
    fo2 = max(eta_aux)/n_sims
      
    #modelo 3
    eta_aux2=[]
    for s in sims:
        eta_aux2.append(sum(x[n,s].x for n in nodos))
    fo3 = max(eta_aux2)/NCells
      
    #modelo 4
    fo4 = sum(y[n].x for n in nodos)
    
    lista_aux=[] 
    for s in sims:
        lista_aux.append(sum(x[n,s].x for n in nodos))
    
    lista_aux = sort(lista_aux)
    #print(lista_aux)
    
    # peores_casos1 = int(round((1-0.9)*n_sims,0))
    # peores_casos2 = int(round((1-0.95)*n_sims,0))
    
    prom_pcb1= 0
    prom_pcb2= 0
    
    # for i in range(peores_casos1):
    #     prom_pcb1 = prom_pcb1 + lista_aux[-i-1]
    # prom_pcb1 = prom_pcb1/peores_casos1
    
    # for i in range(peores_casos2):
    #     prom_pcb2 = prom_pcb2 + lista_aux[-i-1]
    # prom_pcb2 = prom_pcb2/peores_casos2
    
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
    
    titulo = 'modelo4'       
    results = [titulo,n_sims,intensidad,cortafuegos,fo1,fo2,fo3,fo4,tiempo,gap,lmbda,b,prom_pcb1,prom_pcb2]
            
    burned_map = dict.fromkeys(AvailSet,0)
    for s in sims:
        for n in nodos:
            if x[n,s].x > 0:
                burned_map[n] = burned_map[n]+1
                
    
    titulo = 'm4_i'+str(fi)+'b'+str(b)+'l'+str(lmbda)+'s'+forest
    results = ['modelo4',forest,fi,fo1,fo2,fo3,tiempo,gap,lmbda,b,prom_pcb1,prom_pcb2,xvartype,yvartype,incumbent,bestbd] 
    titulo2 = 'm4 '+str(contador_cfuegos)+'cf '+str(n_sims)+'sims beta'+str(b)+' lambda'+str(lmbda)
    titulos = [titulo,titulo2]
    
    return results, fb_list, burned_map,titulos,st_points,tiempo,model.ObjVal
    
#%% MODELO 42

def modelo42(intensidad,n_sims,warmstart,gap,tlimit,xvartype,yvartype,w,params,b,lmbda,fi,forest,cb,firebreaks,sim):
    
    #forest = nsims
    
    NCells,cmdoutput1,st_points,nodos,AvailSet,pathmessages = params

    
    sims = list(range(1,n_sims+1))
    
        
    celdas = len(AvailSet)    
        
    #cortafuegos = int(NCells*intensidad)
    cortafuegos = int(intensidad)

    model = gp.Model()
    model.setParam("OutputFlag", cb)
    
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    eta = model.addVars(sims, vtype=GRB.CONTINUOUS, lb=0)
    phi = model.addVar(vtype=GRB.CONTINUOUS, lb =0)
    
    f_esperanza = gp.quicksum(x[n,s]*w[n] for n in nodos for s in sims)/n_sims
    f_cvar = phi+(1/(1-b))*gp.quicksum(eta[s] for s in sims)/n_sims
    
    model.setObjective(lmbda*f_esperanza+(1-lmbda)*f_cvar, GRB.MINIMIZE)
    
    #firebreak intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    #fire spread constraint
    simulacion = 0
    contador_sims = 0
    for f in range(1, n_sims+1):
    
        file = str(f).zfill(2) + '.csv'
        file = os.path.join('MessagesFile' + file)
        #print(file)
        simulacion = simulacion +1
        # Leemos la marca
        H = nx.read_edgelist(path = pathmessages + file,
                             delimiter=',',
                             create_using = nx.DiGraph(),
                             nodetype = int,
                             data = [('time', float), ('ros', float)])        
        for n in list(H.nodes):
            nbrs = list(H.neighbors(n))
            for nbr in nbrs:
                model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
    
    #starting points constraint            
    for s in sims:
        point = st_points[s-1]
        try:
            model.addConstr(x[point,s] == 1)
        except:
            continue
        
    
    for s in sims:
        model.addConstr(gp.quicksum(x[n,s] for n in nodos)-phi <= eta[s])
    
    for n in warmstart:
        y[n].Start = 1
        
        
    for n in firebreaks:
        model.addConstr(y[n] == 1)
    
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()
     
    if yvartype == GRB.BINARY:
        gap = model.MIPGap
        gap = round(gap,3)
        incumbent = round(model.ObjVal,5)
        bestbd = round(model.ObjBound,5)
    else:
        gap = 'n/a'
        incumbent = 'n/a'
        bestbd = 'n/a'
        
    tiempo = round(model.Runtime)
    
    lista_aux=[]
    for s in sims:
        suma = sum(x[n,s].x for n in nodos)/400
        lista_aux.append(suma)
        
    ev = sum(lista_aux)/len(lista_aux)
    print(ev)
    lista_aux.sort(reverse=True)
    wc_scenarios = int(n_sims*0.1)
    if wc_scenarios <1:
        wc_scenarios = n_sims
    wc = sum(lista_aux[i] for i in range(wc_scenarios))
    wc= wc/(wc_scenarios)
        
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
    
    titulo = 'modelo4'       
    results = [titulo,n_sims,fi,cortafuegos,lmbda,b,ev,wc,tiempo,gap]
      
    bp_s = []     
    burned_map = dict.fromkeys(range(1,401),0)
    for s in sims:
        bp_s.append(sum(x[n,s].x for n in nodos)/400)
        for n in nodos:
            if x[n,s].x > 0.9:
                burned_map[n] = burned_map[n]+1
    
    titulo = 'm4_i'+str(fi)+'b'+str(b)+'l'+str(lmbda)+'s'+forest
    results = [sim,forest,fi,lmbda,b,ev,wc,tiempo,gap]
    titulo2 = 'm4 '+str(contador_cfuegos)+'cf '+str(n_sims)+'sims beta'+str(b)+' lambda'+str(lmbda)
    titulos = [titulo,titulo2]
    
    return results, fb_list, burned_map,titulos,st_points,tiempo,model.ObjVal,bp_s

#%% MODELO PATTERN

def modelo_pattern(n_sims,gap,tlimit,xvartype,yvartype,w,warmstart,params,cb,firebreaks,cortes,c,arcos):
    
    #forest = nsims
    NCells,cmdoutput1,st_points,nodos,AvailSet,pathmessages = params

    
    sims = list(range(1,n_sims+1))
    arcos = list(range(1,arcos+1))
    
        
    celdas = len(AvailSet)    
        
    cortafuegos = int(len(arcos)*5)
    
    #nodos = list(range(1,401))

    model = gp.Model()
    model.setParam("OutputFlag", cb)
    
    x = model.addVars(nodos, sims, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    z = model.addVars(arcos,nodos,nodos, vtype=GRB.BINARY)
    h = model.addVars(arcos,nodos,vtype=GRB.BINARY)

    f_esperanza = gp.quicksum(x[n,s]*w[n]*c[n] for n in nodos for s in sims)
    
    model.setObjective(f_esperanza, GRB.MINIMIZE)
    
    #firebreak intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    #fire spread constraint
    simulacion = 0
    contador_sims = 0
    for f in range(1, n_sims+1):
    
        file = str(f).zfill(2) + '.csv'
        file = os.path.join('MessagesFile' + file)
        #print(file)
        simulacion = simulacion +1
        # Leemos la marca
        H = nx.read_edgelist(path = pathmessages + file,
                             delimiter=',',
                             create_using = nx.DiGraph(),
                             nodetype = int,
                             data = [('time', float), ('ros', float)])        
        for n in list(H.nodes):
            nbrs = list(H.neighbors(n))
            for nbr in nbrs:
                model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
    
    #starting points constraint            
    for s in sims:
        point = st_points[s-1]
        try:
            model.addConstr(x[point,s] == 1)
        except:
            continue


    for a in arcos:
        model.addConstr(gp.quicksum(h[a,i] for i in nodos) == 5)
        
    dmax = math.sqrt(5)
    for i in nodos:
        #model.addConstr(gp.quicksum(z[a,i,j] for a in arcos for j in nodos  if j != i) >= y[i])
        model.addConstr(gp.quicksum(h[a,i] for a in arcos) >= y[i])
        model.addConstr(gp.quicksum(h[a,i] for a in arcos) <= 1)
        for j in nodos:
            if i != j:
                model.addConstr(gp.quicksum(z[a ,i, j] for a in arcos) <= 1)
                #model.addConstr(gp.quicksum(z[a, i, j] for a in arcos) >= y[i] + y[j] - 1)
                for a in arcos:
                    model.addConstr(z[a, i, j] <= y[i])
                    model.addConstr(z[a, i, j] <= y[j])
                    model.addConstr(z[a,i,j] >= h[a,i]+h[a,j]-1)
                    model.addConstr(z[a, i, j]*distance(400,i,j)[0] <= dmax)
                
    
        
    for f in cortes:
        model.addConstr(gp.quicksum(y[n] for n in f) <= 4)
        
    for n in firebreaks:
        model.addConstr(y[n] == 1)
        
    for n in warmstart:
        y[n].Start = 1
    
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()
     
    if yvartype == GRB.BINARY:
        gap = model.MIPGap
        gap = round(gap,3)
        incumbent = round(model.ObjVal,5)
        bestbd = round(model.ObjBound,5)
    else:
        gap = 'n/a'
        incumbent = 'n/a'
        bestbd = 'n/a'
        
    tiempo = round(model.Runtime)
    
    lista_aux=[]
    for s in sims:
        suma = sum(x[n,s].x for n in nodos)/400
        lista_aux.append(suma)
        
    ev = sum(lista_aux)/len(lista_aux)
    #print(ev)
    lista_aux.sort(reverse=True)
    wc_scenarios = int(n_sims*0.1)
    if wc_scenarios <1:
        wc_scenarios = n_sims
    wc = sum(lista_aux[i] for i in range(wc_scenarios))
    wc= wc/(wc_scenarios)
        
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
    
      
    bp_s = []     
    burned_map = dict.fromkeys(range(1,401),0)
    for s in sims:
        bp_s.append(sum(x[n,s].x for n in nodos)/400)
        for n in nodos:
            if x[n,s].x > 0.9:
                burned_map[n] = burned_map[n]+1
    
    
    resulting_arcs = []
    
    contador=0
    for a in arcos:
        resulting_arcs.append([])
        for i in nodos:
            if h[a,i].x == 1:
                #print(a,i)
                resulting_arcs[contador].append(i)
        contador=+1
            
    
    return fb_list, burned_map,model.ObjVal,bp_s,resulting_arcs

#%% MODELO DINAMICO

def modelo_dinamico(intensidad,n_sims,gap,tlimit,params,cb,firebreaks,configs,path):
       
    NCells,_,st_points,nodos,AvailSet,_ = params

    sims = list(range(1,len(n_sims)+1))
    print(sims)
    celdas = len(AvailSet)    
    cortafuegos = int(NCells*intensidad)
    w = dict.fromkeys(nodos,1/NCells)
    
    nodos = list(range(1,10))

    model = gp.Model()
    model.setParam("OutputFlag", cb)
    
    x = model.addVars(nodos, sims, configs, vtype=GRB.BINARY)
    y = model.addVars(nodos, vtype=GRB.BINARY)
    c = model.addVars(configs, vtype=GRB.BINARY)
    z = model.addVars(configs,sims, vtype=GRB.BINARY)
    
    f_esperanza = gp.quicksum(x[n,s]*w[n] for n in nodos for s in sims)
    model.setObjective(f_esperanza, GRB.MINIMIZE)
    
    #firebreak intensity constraint
    model.addConstr(gp.quicksum(y[n] for n in nodos) == cortafuegos)
    
    
    
    df = pd.DataFrame(columns = ['scenario','firebreaks'])
    contador = 1
    for conf_folder in configs:
        
        print('conf:',conf_folder)
        
        df.loc[len(df)] = [contador,list(conf_folder)]
        contador = contador+1
            
        pathmessages = path+'/'+conf_folder+'/Messages'
    
        #fire spread constraint
        simulacion = 0
        contador_sims = 0
        for f in range(1, n_sims+1):
        
            file = str(f).zfill(2) + '.csv'
            file = os.path.join('MessagesFile' + file)
            #print(file)
            simulacion = simulacion +1
            # Leemos la marca
            H = nx.read_edgelist(path = pathmessages +'/'+ file,
                                 delimiter=',',
                                 create_using = nx.DiGraph(),
                                 nodetype = int,
                                 data = [('time', float), ('ros', float)])        
            for n in list(H.nodes):
                nbrs = list(H.neighbors(n))
                for nbr in nbrs:
                    model.addConstr(x[n,simulacion] <= x[nbr,simulacion]+y[nbr])
    
    #starting points constraint            
    for s in sims:
        point = st_points[s-1]
        try:
            model.addConstr(x[point,s] == c[s])
        except:
            continue
        
    for s in sims:
        for n in nodos:
            x[n,s] <= c[s]
        
    for n in firebreaks:
        model.addConstr(y[n] == 1)
        
    for cfg in range(0,len(configs)+1):
        config = df.loc[cfg,'firebreaks']
        for n in config:
            model.addConstr(y[n] >= c[s])
            
    
    if yvartype == GRB.BINARY:
        model.Params.MIPGap = gap
    model.Params.TimeLimit = tlimit
    model.optimize()

        
    tiempo = round(model.Runtime)
    
    lista_aux=[]
    for s in sims:
        suma = sum(x[n,s].x for n in nodos)/400
        lista_aux.append(suma)
        
    ev = sum(lista_aux)/len(lista_aux)
    print(ev)
    lista_aux.sort(reverse=True)
    wc_scenarios = int(n_sims*0.1)
    if wc_scenarios <1:
        wc_scenarios = n_sims
    wc = sum(lista_aux[i] for i in range(wc_scenarios))
    wc= wc/(wc_scenarios)
        
    contador_cfuegos=0
    fb_list = []
    for n in nodos:
        if y[n].x > 0.9:
            contador_cfuegos = contador_cfuegos+1
            fb_list.append(n)
    
    titulo = 'modelo4'       
    results = [titulo,n_sims,fi,cortafuegos,lmbda,b,ev,wc,tiempo,gap]
      
    bp_s = []     
    burned_map = dict.fromkeys(range(1,401),0)
    for s in sims:
        bp_s.append(sum(x[n,s].x for n in nodos)/400)
        for n in nodos:
            if x[n,s].x > 0.9:
                burned_map[n] = burned_map[n]+1
    
    results = [n_sims,intensidad,ev,wc,tiempo,gap,model.ObjVal]
    
    return results, fb_list, burned_map,st_points,bp_s


#%% save results

def save_results(lista,path,f):
    
    df_results = pd.DataFrame(index=None)
    cols = ['seed','modelo','nsims','fsize','intensidad','fo1','fo2','fo3','ptime','gap','lambda','beta','10 wc avg','5 wc avg','dom x','dom y','incumbent','bestbd']
    dic = dict.fromkeys(cols,0)
    lista.insert(3,f)
    i = 0
    for key in dic:
        dic[key] = lista[i]
        i = i+1
    df_results = df_results.append(dic,ignore_index = True)
    
    with pd.ExcelWriter(path,engine = 'openpyxl',mode='a',if_sheet_exists='overlay') as writer:
        df_results.to_excel(writer, sheet_name='Sheet1',startrow=writer.sheets['Sheet1'].max_row, header=None,index=False)
        
def save_results2(lista,path,f):
    
    df_results = pd.DataFrame(index=None)
    cols = ['seed','nsims','forest','intensidad','lambda','b','ev','wc','ptime','gap']
    dic = dict.fromkeys(cols,0)
    lista.insert(0,f)
    i = 0
    for key in dic:
        dic[key] = lista[i]
        i = i+1
    df_results = df_results.append(dic,ignore_index = True)
    
    with pd.ExcelWriter(path,engine = 'openpyxl',mode='a',if_sheet_exists='overlay') as writer:
        df_results.to_excel(writer, sheet_name='Sheet1',startrow=writer.sheets['Sheet1'].max_row, header=None,index=False)
        
def harvested(folder,l,titulos,i): #funcion que pasa una lista de elementos a un archivo .csv que contiene a los cortafuegos
    datos=[np.insert(l,0,1)] #inserto el elemento 1 que corresponde al ano que necesita el archivo 
    if len(l)==0: #si no hay cortafuegos
        cols=['Year Number'] #creo solamente una columna correspondiente al ano
    else: #si hay cortafuegos
        colu=['Year Number',"Cell Numbers"] #creo 2 columnas
        col2=[""]*(len(l)-1) #creo el resto de columnas correspondientes a los otros nodos
        cols=colu+col2 #junto ambas columnas
    titulo = titulos[0]
    titulo2 = titulos[1]
    
    #folder = 'eval2/'+titulo[0:2]+'/'
    df = pd.DataFrame(datos,columns=cols) #creo el dataframe
    df.to_csv(folder,index=False) #guardo el dataframe
    
def harvested2(folder,l,titulos): #funcion que pasa una lista de elementos a un archivo .csv que contiene a los cortafuegos
    datos=[np.insert(l,0,1)] #inserto el elemento 1 que corresponde al ano que necesita el archivo 
    if len(l)==0: #si no hay cortafuegos
        cols=['Year Number'] #creo solamente una columna correspondiente al ano
    else: #si hay cortafuegos
        colu=['Year Number',"Cell Numbers"] #creo 2 columnas
        col2=[""]*(len(l)-1) #creo el resto de columnas correspondientes a los otros nodos
        cols=colu+col2 #junto ambas columnas
    titulo = titulos[0]
    titulo2 = titulos[1]
    
    #folder = 'eval2/'+titulo[0:2]+'/'
    df = pd.DataFrame(datos,columns=cols) #creo el dataframe
    df.to_csv(folder+"harvested_"+titulo+".csv",index=False) #guardo el dataframe
    
def harvested3(folder,l): #funcion que pasa una lista de elementos a un archivo .csv que contiene a los cortafuegos
    datos=[np.insert(l,0,1)] #inserto el elemento 1 que corresponde al ano que necesita el archivo 
    if len(l)==0: #si no hay cortafuegos
        cols=['Year'] #creo solamente una columna correspondiente al ano
    else: #si hay cortafuegos
        colu=['Year',"Ncell"] #creo 2 columnas
        col2=[""]*(len(l)-1) #creo el resto de columnas correspondientes a los otros nodos
        cols=colu+col2 #junto ambas columnas
    
    #folder = 'eval2/'+titulo[0:2]+'/'
    df = pd.DataFrame(datos,columns=cols) #creo el dataframe
    df.to_csv(folder,index=False)
    
def model_map(fb_list,burned_results,burned,FAG,titulos,b,lmbda,colors,folder,coord_pos,st_points,show,maxi):

    min_node = builtins.min(burned_results, key=burned_results.get)
    burned_results[min_node] = max(burned.values())
    maximo = builtins.max(burned.values())
    #burned_results[min_node] = maxi
    
    figure(1, figsize=(11*2.5,9*2.5))
    nx.draw_networkx(FAG, pos = coord_pos,
                             node_size = 50,
                             with_labels = False,
                             node_shape = 's',
                             node_color = 'k')
    
    nx.draw_networkx_nodes(FAG, pos = coord_pos,
                       node_size = 50,
                       nodelist = list(FAG.nodes),
                       node_shape='s',
                       node_color = list(colors.values()))
    
    
    if show == True:
        nc = nx.draw_networkx_nodes(FAG, pos = coord_pos,
                                        # node_size = nodesize,  # array
                                        # nodelist=nodelist,  # list
                                    linewidths = 2.0,
                                    node_size = 15,
                                    cmap = 'Reds',
                                    node_shape = 's',
                                    node_color = list(burned_results.values())
                                    )
    
    cf = nx.Graph()
    nx.draw_networkx_nodes(cf, pos = coord_pos,
                       node_size = 10,
                       nodelist = fb_list,
                       node_shape='s',
                       node_color = 'b')
    
    
    # sp = nx.Graph()
    # nx.draw_networkx_nodes(cf, pos = coord_pos,
    #                     node_size = 15,
    #                     nodelist = st_points,
    #                     node_shape='s',
    #                     node_color = 'cyan')
    
    titulo = titulos[0]
    folder = folder+'/'
    
    #plt.title(titulos[1])
    if show == True:
        plt.colorbar(nc)
     
    
    plt.savefig(folder+titulo+'.png', bbox_inches='tight')
    plt.clf()
    plt.show()

    return maximo

def random_solution(NCells,nodos,folder,intensidad):

    intensidades = intensidad
    
    fb_list = []
    for i in intensidades:
        harvested2(folder,fb_list,['rd_solutions'+str(i),'a'])
        elements = rd.sample(nodos,int(0.01*NCells))
        for e in elements:
            nodos.remove(e)

        fb_list = fb_list + elements
        
    
#%% find most burned

def find_scenarios(folder_sc):
    
    # Ruta de la carpeta que contiene los archivos Excel
    folder_path = folder_sc+'/Messages/'

    # Obtener la lista de archivos en la carpeta
    find_list = os.listdir(folder_path)

    # Filtrar la lista para quedarnos solo con los archivos Excel
    excel_files = [file for file in find_list if file.endswith('.csv')]

    # Ordenar la lista de archivos por tamaño en orden descendente
    sorted_files = sorted(excel_files, key=lambda file: os.path.getsize(os.path.join(folder_path, file)), reverse=True)

    #sorted_files = sorted_files[0:n]
    # print("Archivos ordenados por tamaño:")
    for file in range(len(sorted_files)):
         sorted_files[file] = ''+sorted_files[file]
    
        
    return sorted_files

#%% succesors

def dpv_f(G):
    nodes = list(G.nodes())
    dic = dict.fromkeys(nodes,0)
    for n in nodes:
        subgraph = nx.ego_graph(G,n,radius=10000)
        neighbors= list(subgraph.nodes())
        dpv = len(neighbors)
        dic[n] = dic[n]+dpv
    return dic

def get_fb_from_dpv(dictionary,k):
    sorted_items = sorted(dictionary.items(), key=lambda x: x[1], reverse=True)
    top_k_items = sorted_items[:k]
    top_k_keys = [item[0] for item in top_k_items]
    return top_k_keys

#READ SOLUTION
# path = 'C:/Users/matia/Documents/proyecto_magister/resultados/expfinal2/sols1002/harvested_0.04.csv'
# df = pd.read_csv(path,header=1)
# firebreaks =[]
# for i in df:
#     print(i)
#     firebreaks.append(int(i))

def distance(n,i,j):
    dim = n**(1/2)
    p1 = [(i-1)//(dim),i-dim*((i-1)//(dim))-1];
    p2 = [(j-1)//(dim),j-dim*((j-1)//(dim))-1];

    c1 = (i-1)//(dim)
    c2 = i-dim*((i-1)//(dim))-1

    distancia = math.dist(p1,p2)
    return distancia,c1,c2

def arc_pattern(lista):
    d1 = 4+math.sqrt(5)
    d2 = 3+math.sqrt(2)+math.sqrt(5)
    d3 = 2+2*math.sqrt(2)
    lista_aux = [d1,d2,d3,d1,d2]
    distancias = []
    for i in lista:
        suma = sum(distance(400,i,j)[0] for j in lista)
        distancias.append(suma)
    prueba = False
    for i in distancias:
        if i in lista_aux:
            lista_aux.remove(i)
    if len(lista_aux) == 0:
        prueba = True
    return prueba

def c2f(results_folder,ev_sims,harvest_file,seed,forest):
    command = "wsl cd ../C2F-W;"
    c = 'python3 main.py --input-instance-folder data/'+forest+'/ --output-folder '+results_folder +' --sim-years 1 --nsims '+str(ev_sims)+' --finalGrid --weather random --nweathers 2 --Fire-Period-Length 1.0 --output-messages --ROS-CV 0.1 --seed '+str(seed)#+' --HarvestedCells ../../../'+harvest_file
    command = command + c
    result = subprocess.check_output(command, shell=True)
    #print(result.decode("utf-8"))
    
def c2fw(results_folder,ev_sims,harvest_file,seed,forest):
    command = "wsl cd ../C2F-W;"
    c = 'python3 main.py --input-instance-folder data/'+forest+'/ --output-folder '+results_folder +' --sim-years 1 --nsims '+str(ev_sims)+' --weather rows --Fire-Period-Length 1.0 --output-messages --ROS-CV 0.0 --seed '+str(seed)+' --cros --sim C --out-cfb --out-sfb --FirebreakCells '+harvest_file
    command = command + c
    result = subprocess.check_output(command, shell=True)
    #print(result.decode("utf-8"))
    
def process_file(file_path,vacios,dim):
    
    try:
        df = pd.read_csv(file_path, header=None)
    except pd.errors.EmptyDataError:
        vacios = vacios+1
        # El archivo CSV está vacío, puedes manejar esta situación como desees
        return 0.0,vacios,[]
    
    df = pd.read_csv(file_path, header=None)
    w = {key: 0 for key in range(1, dim+1)}
    burned = set()
    for i in range(len(df)):
        q = df.loc[i][1]
        burned.add(q)
    q = df.loc[0][0]
    burned.add(q)
    return len(burned) / dim,vacios,burned


def read_firescar(path,folder,sims,cols,k,fi,l):

    dim = cols*cols

    lista = []
    vacios = 0
    for f in range(1, sims+1):
        #print(f)
        file = str(f).zfill(2) + '.csv'
        ruta = os.path.join(folder,'Messages', 'MessagesFile' + file)
        p1,vacios,burned = process_file(ruta,vacios,dim)
        lista.append(p1)
        #print(burned)
    
    #print(lista)
    suma = sum(lista) / (sims)
    #print(suma)
    lista.sort(reverse=True)
    wc_scenarios = int(sims*0.1)
    if wc_scenarios <1:
        wc_scenarios = sims
    wc = sum(lista[i] for i in range(wc_scenarios))
    wc= wc/(wc_scenarios)
    
    lista = [k,'100',100,'20x20',l,fi,suma,wc]
    cols = ['seed','modelo','nsims','fsize','lambda','intensidad','fo1','10wc']
    
    df_results = pd.DataFrame([lista], columns=cols)
    print(wc)
    
    with pd.ExcelWriter(path,engine = 'openpyxl',mode='a',if_sheet_exists='overlay') as writer:
        df_results.to_excel(writer, sheet_name='Sheet1',startrow=writer.sheets['Sheet1'].max_row, header=None,index=False)
    
    return suma,wc