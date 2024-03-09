import networkx as nx
import pandas as pd
import ReadDataPrometheus
import os

def read_sims(forest_path,results_path,nsims):

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

    scar_graphs = []
    contador_sims = 1
    for f in range(1, nsims+1):
        
        ignitions_points.append(ignitions_points_aux[f-1])
        file = 'MessagesFile' + str(f).zfill(2) + '.csv'
        H = nx.read_edgelist(path = msg_path + file,
                                delimiter=',',
                                create_using = nx.DiGraph(),
                                nodetype = int,
                                data = [('time', float), ('ros', float)])

        scar_graphs.append(H)
        contador_sims = contador_sims+1
    
    params = [NCells,ignitions_points,avail,scar_graphs]
        
    return params

