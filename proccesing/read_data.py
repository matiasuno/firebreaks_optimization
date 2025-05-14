import networkx as nx
import pandas as pd
import os
import numpy as np
from heapq import nlargest
import csv

def get_graph(msg_path):
    H = nx.read_edgelist(path = msg_path,
                            delimiter=',',
                            create_using = nx.DiGraph(),
                            nodetype = int,
                            data = [('time', float), ('ros', float)])

    return H

def read_asc(file_path):
    with open(file_path, 'r') as file:
        # Leer las primeras seis líneas y almacenarlas
        header = [next(file) for _ in range(6)]
        # Leer los datos numéricos
        data = np.loadtxt(file, dtype=np.float32)
        if data.ndim != 2:
            raise ValueError(f"Expected a 2D array, but got {data.ndim}D array.")
        xdim, ydim = np.shape(data)
        nodos = int(xdim*ydim)
    return header, data, nodos

def get_dpv(graphs,ncells,risk_values,intensity):
    dpv_dic = dict.fromkeys(list(range(1,ncells+1)),0)
    for g in graphs:
        g_nodes = g.nodes()
        for n in g_nodes:
                ego_g = nx.ego_graph(g, n, 1000000)
                ego_nodes = list(ego_g.nodes)
                dpv_n = 0
                #print(ego_nodes)
                for j in ego_nodes:
                    dpv_n = dpv_n+risk_values[j]
        dpv_dic[n] = dpv_dic[n]+dpv_n
    dpv_sol = nlargest(int(ncells*intensity), dpv_dic, key=dpv_dic.get)  
    return dpv_sol

def get_bp(bp,ncells,intensity):
    bp_sol = nlargest(int(ncells*intensity), bp, key=bp.get)
    return bp_sol

def read_ignition_points(ignitions):
    # Read the ignitions data
    ignitions_points = []
    with open(ignitions, 'r') as file:
        for line in file:
            ignitions_points.append(int(line.strip()))
    return ignitions_points

def read_messages(msg_path):
    # Read the messages data
    messages = []
    for file in os.listdir(msg_path):
        if file.endswith('.csv'):
            graph = get_graph(f"{msg_path}/{file}")
            messages.append(graph)
    return messages

def write_treatment(header_file,firebreaks,output_path):
    header = {}
    for line in header_file:
        key, value = line.strip().split(maxsplit=1)
        header[key] = float(value) if '.' in value or 'e' in value.lower() else str(value)
    ncols = int(header['ncols'])
    nrows = int(header['nrows'])
    data = np.zeros((nrows, ncols))
    for n in firebreaks:
        # Ajusta 'n' para que sea base 0 en lugar de base 1
        n -= 1  
        row = n // ncols  # Índice de fila
        col = n % ncols   # Índice de columna

        # Solo actualiza si el índice está dentro del rango
        if 0 <= row < nrows and 0 <= col < ncols:
            data[row][col] = 1
            
    with open(output_path, 'w') as file:
        # Escribir el encabezado
        file.writelines(header_file)
        # Escribir los datos numéricos
        np.savetxt(file, data, fmt='%.6f')

def harvested(output,fbs): #funcion que pasa una lista de elementos a un archivo .csv que contiene a los cortafuegos
    datos=[np.insert(fbs,0,1)] #inserto el elemento 1 que corresponde al ano que necesita el archivo 
    if len(fbs)==0: #si no hay cortafuegos
        cols=['Year'] #creo solamente una columna correspondiente al ano
    else: #si hay cortafuegos
        colu=['Year',"Ncell"] #creo 2 columnas
        col2=[""]*(len(fbs)-1) #creo el resto de columnas correspondientes a los otros nodos
        cols=colu+col2 #junto ambas columnas
    df = pd.DataFrame(datos,columns=cols) #creo el dataframe
    df.to_csv(output,index=False)
    
def extract_points(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        points = [int(row['point']) for row in reader]
    return points