import subprocess as subp
import os
from read_data import read_asc, write_treatment,harvested

def c2f_evaluation(forest,output,fuels,nsims,fb_sol):

    print("simulating")

    #INITIALIZATION
    
    #PARAMETERS
    input_folder = f'--input-instance-folder {forest}' 
    output_folder = f'--output-folder {output}'
    nsims = f'--nsims {nsims}'
    nthreads = '--nthreads 28'
    weather_opt = '--nweathers 86 --weather random'
    seed = '--seed 333'
    extra = '--ignitionsLog'
    outputs = '--output-messages'
    ignitions = '--ignitions-random'


    #SAVE FIREBREAK SOLUTIONS
    try:
        os.makedirs(output)
    except FileExistsError:
        print(f"Directory {output} already exists. Skipping creation.")
    harvest_file_asc = f'{output}/harvested.asc'
    harvest_file_csv = f'{output}/harvested.csv'
    header,data,nodos = read_asc(fuels)
    write_treatment(header,fb_sol,harvest_file_asc)
    harvest_csv = harvested(harvest_file_csv,fb_sol)
    
    #CALL C2F
    firebreaks_opt = f'--FirebreakCells {harvest_file_csv}'
    options = " ".join([input_folder,output_folder,nsims,nthreads,seed,extra,outputs,weather_opt,firebreaks_opt,ignitions])
    c2f_call = '/home/matias/source/C2F-W/Cell2Fire/Cell2Fire --sim C  '+options
    subp.call(c2f_call, shell=True, stdout=subp.DEVNULL)