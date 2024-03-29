import networkx as nx
import pandas as pd
import ReadDataPrometheus
import os
import numpy as np

def read_sims(forest_path,results_path,nsims,total_nsims):

    Folder = forest_path
    FBPlookup = Folder + '/fbp_lookup_table.csv'
    ForestFile = Folder + '/fuels.asc'
    IgnitionsFolder = results_path + '/IgnitionsHistory/ignitions_log.csv'
    msg_path = results_path + '/Messages/'
    FBPDict, ColorsDict = ReadDataPrometheus.Dictionary(FBPlookup)
    CellsGrid3, CellsGrid4, Rows, Cols, AdjCells, CoordCells, CellSize = ReadDataPrometheus.ForestGrid(ForestFile,FBPDict)
    
    NCells = Rows * Cols
    m = Rows
    n = Cols
    
    AvailSet = set()
    NonAvailSet = set()
    
    for i in range(NCells):
        if CellsGrid4[i] != 'NF':
            AvailSet.add(i+1)
        else:
            NonAvailSet.add(i+1)
    
    # Read ignitions points folder
    df = pd.read_csv(IgnitionsFolder, header=None, names=['points'])
    # Now, you can access the single column DataFrame
    single_column_df = df['points']
    # Convert the DataFrame column to a list
    ignitions_points_aux = single_column_df.tolist()
    ignitions_points = []
    
    cmdoutput1 = os.listdir(msg_path)
    if ".DS_Store" in cmdoutput1:
        idx = cmdoutput1.index('.DS_Store')
        del cmdoutput1[idx]
    if "desktop.ini" in cmdoutput1:
        idx = cmdoutput1.index('desktop.ini')
        del cmdoutput1[idx]
    
    avail = list(AvailSet)
    burn_probability = dict.fromkeys(list(range(1,Rows*Cols+1)),0)

    scar_graphs = []
    contador_sims = 1
    for f in range(1, total_nsims+1):
        
        ignitions_points.append(ignitions_points_aux[f-1])
        file = 'MessagesFile' + str(f).zfill(4) + '.csv'
        H = nx.read_edgelist(path = msg_path + file,
                                delimiter=',',
                                create_using = nx.DiGraph(),
                                nodetype = int,
                                data = [('time', float), ('ros', float)])

        scar_graphs.append(H)
        contador_sims = contador_sims+1

        for nodo in list(H.nodes()):
            burn_probability[nodo] = burn_probability[nodo]+1        

        if f == nsims:
            break

    for key in burn_probability:
        burn_probability[key] = burn_probability[key]/nsims

    
    params = [NCells,ignitions_points,avail,scar_graphs,burn_probability]
        
    return params

def harvested(folder,l): #funcion que pasa una lista de elementos a un archivo .csv que contiene a los cortafuegos
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

def dpv(graphs,ncells,bp):
    dpv_dic = dict.fromkeys(list(range(1,ncells+1)),0)

    for g in graphs:
        g_nodes = g.nodes()
        for n in g_nodes:
                ego_g = nx.ego_graph(g, n, 1000000)
                ego_nodes = list(ego_g.nodes)
                dpv_n = 0
                #print(ego_nodes)
                for j in ego_nodes:
                    dpv_n = dpv_n+bp[j]
                
                
        dpv_dic[n] = dpv_dic[n]+dpv_n
    return dpv_dic