from master import *
import pandas as pd
import copy



folder2 = 'C:/Users/matia/Documents/proyecto_magister/resultados/Sub20mod/0'
data = 'C:/Users/matia/Documents/proyecto_magister/simulador/Cell2fire-main/data/20x20'

rd.seed(123)
archivos=[]
#archivos = find_scenarios(folder2)

xvar = [GRB.BINARY,GRB.CONTINUOUS]
yvar = [GRB.BINARY,GRB.CONTINUOUS]

x = GRB.BINARY
y = GRB.BINARY

sims = 1000

intensities = [0.01,0.02,0.03,0.04,0.05]
#intensities = [0]
lbda = [1]#,0.5,0]
seeds = ['seed1']#,'seed2','seed3','seed4','seed5']
nsims = [20]#,100,180]
 
firebreaks=[]  
contador = 1
instancias = len(intensities)*len(nsims)*len(seeds)*len(lbda)
gap = 0
tlimit = 3600*(1/2)
b = 0.9

ruta = folder2
params,mdg,dpv_dic = read_sims(data,ruta,sims,archivos,True)
avail = params[4]
nodos = params[3]
cells = params[0]
FAG = params[6]
colors = params[7]
coordpos = params[8]
burned = params[9]
cc = params[12]
st_points = params[2]
fb_list = []
for s in seeds: 
    for l in lbda:
        for sim in nsims:
            for i in intensities:
                
                # path_sol = 'C:/Users/matia/Documents/proyecto_magister/resultados/Sub20mod/sols/harvested_'
                # path_sol2 = path_sol+str(i)+'_l'+str(l)+'_'+s+'_s'+str(sim)+".csv"
                # df = pd.read_csv(path_sol2,header=1)
                # firebreaks =[]
                # for fb in df:
                #     #print(fb)
                #     firebreaks.append(int(fb))
                # del firebreaks[0]
                
                # print('-'*30,s,'-'*30)
                # print('-'*30,'lambda:',l,'-'*30)
                # print('-'*30,'nsims:',sim,'-'*30)
                # print('-'*30,'alpha:',i,'-'*30)
                
                if i != 0:
                    nodos2 = copy.deepcopy(nodos)
                    elements = rd.sample(nodos2,int(0.01*400))
                    for e in elements:
                        nodos2.remove(e)
                    fb_list = fb_list + elements
                
                print(contador,'de: ',instancias)
                
                cortafuegos = int(cells*i)                

                w = dict.fromkeys(nodos,1/cells)
                
                # for sp in st_points:
                #     if sp != 'na':
                #         dpv_dic[sp] = 0   
                # ws = get_fb_from_dpv(dpv_dic, cortafuegos)
                    
                    
                r_folder = 'C:/Users/matia/Documents/proyecto_magister/resultados/Sub20mod'
                
                
                # #modelo4(intensidad, warmstart, gap, tlimit, xvartype, yvartype, w, params, b, lmbda, fi, forest)
                # #results = modelo4(i,warmstart,gap,tlimit,x,y,w,params[0:6], b, l,i,str(sim),1,firebreaks,[])
                results = modelo42(len(fb_list), sims,[], gap, tlimit, x, y, w, params[0:6], b, l, i, '20x20', 0,fb_list,sim)
                #return results, fb_list, burned_map,titulos,st_points,tiempo,model.ObjVal,bp_s
                # #results[0].insert(0,s)
                
                save_results2(results[0],r_folder+'/ev_static1000.xlsx',s)
                
                #print(results[1])
                #print(results[-2])
                
                #name = "harvested_"+str(i)+'_l'+str(l)+'_'+s+'_s'+str(sim)+".csv"
                #harvested(r_folder+'/sols/'+name, results[1], results[3],i)
                
                # #tlimit = tlimit-results[5]
                
                # #print(results[5],results[0][4])
                
                # #firebreaks = results[1]
                contador = contador +1
    

