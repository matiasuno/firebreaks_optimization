from read_data import *
from optimization import *
from plot import *
from heapq import nlargest
from statistics import mean


## DATA READING
scars = []
results_folders = ['bp','dpv','mip']
for r in results_folders:
    forest_example = "C:/Users/matia/source_win/firebreaks_optimization/sample_test/data/CanadianFBP/Sub20"
    results_example = "C:/Users/matia/source_win/firebreaks_optimization/sample_test/Sub20_"
    harvest_folder = "C:/Users/matia/source_win/firebreaks_optimization/harvest/"

    forest_path = forest_example
    results_path = results_example+r
    nsims = 1000
    total_nsims = 1000

    params = read_sims(forest_path,results_path,nsims,total_nsims)
    NCells,ignitions_points,avail,scar_graphs,bp = params
    #adjacency_matrix = adjacency(20,20)
    #dpv_dic = dpv(scar_graphs,NCells,bp)


    contador = 1
    #'''
    ## OPTIMIZATION MODEL
    i = [0]
    for intensity in i:
        #bp_sol = nlargest(int(NCells*intensity), bp, key=bp.get)
        #dpv_sol = nlargest(int(NCells*intensity), dpv_dic, key=dpv_dic.get)
        #mip_sol = []
        #sols = [bp_sol,dpv_sol,mip_sol]
        sols = [[]]
        for s in sols:
            print("-"*20,"Solving model","-"*20)
            #intensity = 0
            gap = 0
            tlimit = 1800
            w_parameter = 1
            linked = False
            solution = []
            verbose = 0
            
            fo, fb_list, expected_value, lista_aux = model(intensity,nsims,gap,tlimit,w_parameter,params[0:-1],linked,s,verbose)
            scars.append(lista_aux)
            #print("fo: ",fo)
            print("alpha: ", intensity,"ev: ",expected_value)
            #print("firebreaks: ",fb_list)
            #print("total burned forest by scenario: ", lista_aux)

            # Example usage
            prob_map = np.random.rand(40, 40)  # Example probability map
            highlight_nodes = fb_list  # Example highlight node coordinates
            #plot_squared_graph(prob_map, highlight_nodes)

            #harvested(harvest_folder+'harvest'+str(contador)+'.csv',fb_list)
            contador = contador+1
        #'''


# Creating a boxplot
        
# Colors
BG_WHITE = "#ffffff"
GREY_LIGHT = "#b4aea9"
GREY50 = "#7F7F7F"
BLUE_DARK = "#1B2838"
BLUE = "#2a475e"
BLACK = "#282724"
GREY_DARK = "#747473"
RED_DARK = "#850e00"

# Colors taken from Dark2 palette in RColorBrewer R library
COLOR_SCALE = ["#1B9E77", "#D95F02", "#7570B3"]

POSITIONS = [0, 1,2]
# Horizontal lines
HLINES = [0.2, 0.4, 0.6, 0.8]

fig, ax = plt.subplots(figsize= (10, 8))

for h in HLINES:
    ax.axhline(h, color=GREY50, ls=(0, (5, 5)), alpha=0.8, zorder=0)

fig.patch.set_facecolor(BG_WHITE)
ax.set_facecolor(BG_WHITE)

# Add boxplots ---------------------------------------------------
# Note that properties about the median and the box are passed
# as dictionaries.

medianprops = dict(
    linewidth=4, 
    color=GREY_DARK,
    solid_capstyle="butt"
)
boxprops = dict(
    linewidth=2, 
    color=GREY_DARK
)

ax.boxplot(
    scars,
    positions=POSITIONS, 
    showfliers = False, # Do not show the outliers beyond the caps.
    showcaps = False,   # Do not show the caps
    medianprops = medianprops,
    whiskerprops = boxprops,
    boxprops = boxprops
)

means = [mean(y) for y in scars]
for i, mean in enumerate(means):
    # Add dot representing the mean
    ax.scatter(i, mean, s=100, color=RED_DARK, zorder=3)
    
    # Add line conecting mean value and its label
    ax.plot([i, i + 0.25], [mean, mean], ls="dashdot", color="black", zorder=3)
    
    # Add mean value label.
    ax.text(
        i + 0.25,
        mean,
        r"$\hat{\mu}_{\rm{mean}} = $" + str(round(mean, 2)),
        fontsize=13,
        va="center",
        bbox = dict(
            facecolor="white",
            edgecolor="black",
            boxstyle="round",
            pad=0.15
        ),
        zorder=10 # to make sure the line is on top
    )

# Displaying the plot
plt.title('Burn Probability vs. DPV vs. Stochastic MIP  \n (out of sample, n=1000)',size=20)
ax.set_xticklabels(['BP','DPV','MIP'], size=15, ha="center", ma="center")
plt.xlabel('Decision models', size=15)
plt.ylabel('Percentage of burned forest',size=15,labelpad=10)
plt.show()