from operations_csv import get_all_messages, write_treatment_csv
from operations_raster import read_raster, write_treatment_raster
from two_stage_model import mip_model, save_results, get_wc

if __name__ == "__main__":
    
    #ruta_base = "/Users/matiasvilches/Documents/F2A/source/two_stage"
    ruta_base = "/home/matias/Documents/source/two_stage"
    forest_list = ["sub20", "sub40", "sub100"]
    header = ["forest", "nsims","intensity", "lambda", "instance", "time","gap", "fo", "ev","ws",]
    results_output_csv = f"{ruta_base}/results/summary_results.csv"

    time_limit = 3600
    time_limit = 300
    min_gap = 0.01
    lambda_list = [1,0.5,0]
    intensity_list = [0.01,0.02,0.03,0.04,0.05]
    nsims_list = [20,60,100,140,180]
    instances = [i for i in range(1,6)]

    """
    forest_list = ["sub40"]
    lambda_list = [1]
    intensity_list = [0.02]
    nsims_list = [100]
    n_instances = 1
    instances = [i for i in range(1,1+n_instances)]
    """

    total_experiments = len(forest_list)*len(nsims_list)*len(lambda_list)*len(intensity_list)*len(instances)

    experiment_counter = 0
    for forest in forest_list:
        for nsims in nsims_list:
            for lmda in lambda_list:
                for intensity in intensity_list:
                    for i in instances:
                        experiment_counter += 1
                        print(f"Running",f"{experiment_counter}/{total_experiments}", f"forest: {forest}, lambda: {lmda}, intensity: {intensity}, nsims: {nsims}, instance: {i}")
                    
                        msg_path = f"{ruta_base}/sims/{forest}/msg{i}/"
                        fuels = f"{ruta_base}/forest/{forest}/fuels.asc"

                        data, profile, n_nodos = read_raster(fuels)
                        scar_graphs = get_all_messages(msg_path)

                        params = (intensity, nsims, min_gap, time_limit) #intensity, nsims, gap, tlimit

                        treatment_output_raster = f"{ruta_base}/results/{forest}/firebreaks_{forest}_i{intensity}.tif"
                        treatment_output_csv = f"{ruta_base}/results/{forest}/{forest}_i{intensity}_l{lmda}_inst{i}_n{nsims}.csv"

                        fo, fb_list, ev, lista_aux, gap, time = mip_model(params, n_nodos, scar_graphs, lmda, n_cfuegos=None)
                        wc = get_wc(lista_aux, alpha=0.1)

                        print(ev)

                        write_treatment_csv(treatment_output_csv,fb_list)

                        results = [forest, nsims, intensity, lmda, i, time, gap, fo, ev, wc]
                        save_results([results], header=header, output_path=results_output_csv)