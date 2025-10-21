import networkx as nx
import pandas as pd
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from operations_csv import convert_csv_to_pickle
from collections import defaultdict
import rasterio
import numpy as np

def process_file(file_path, nsims):
    """
    Procesa un archivo Pickle y devuelve la contribución a bp_dic.
    """

    try:
        df = pd.read_pickle(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return {}

    df.columns = ['source', 'target', 'time', 'ros']
    H = nx.from_pandas_edgelist(df, source='source', target='target', edge_attr=True, create_using=nx.DiGraph())
    nodes = list(H.nodes())
    r = {n: 1 / nsims for n in nodes}
    #print("Number of nodes:", len(r))
    return r

def process_pickle_parallel(pickle_path, nsims, n_cores):
    """
    Procesa múltiples archivos Pickle en paralelo y combina los resultados.
    """
    # Lista de archivos Pickle con rutas completas
    # only save number of files < nsims
    files = [os.path.join(pickle_path, file) for file in os.listdir(pickle_path) if file.endswith('.pkl')][:nsims]
    bp_dic = defaultdict(float)
    # Paralelismo con número de núcleos configurables
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = executor.map(process_file, files, [nsims] * len(files))
        for i,bp_contrib in enumerate(results,1):
            print(f"Processed file {i}/{len(files)}")
            for n, contrib in bp_contrib.items():
                bp_dic[n] = bp_dic.get(n, 0) + contrib
    
    return bp_dic


def write_raster_from_dict_ids(base_raster, file_path, raster_dict, default_nodata=-9999):
    """
    Crear un raster a partir de un diccionario {cell_id: valor}.
    - base_raster: ruta a raster base (se usa para crs, transform, shape, nodata por defecto)
    - file_path: ruta de salida (.tif)
    - raster_dict: dict {cell_id (1-based): value}
    - default_nodata: valor nodata si el raster base no lo define
    """
    # --- abrir base raster y leer metadatos ---
    with rasterio.open(base_raster) as src:
        nrows, ncols = src.height, src.width
        transform = src.transform
        crs = src.crs
        src_nodata = src.nodata
        profile = src.profile.copy()

    # nodata a usar
    nodata = src_nodata if src_nodata is not None else default_nodata

    # crear array con nodata
    # usar float32 por compatibilidad general; si quieres mantener enteros, modifica dtype
    out_dtype = "float32"
    array = np.full((nrows, ncols), nodata, dtype=out_dtype)

    # rellenar valores desde raster_dict (verificando rangos)
    total = nrows * ncols
    for cell_id, value in raster_dict.items():
        try:
            cid = int(cell_id)
        except Exception:
            print(f"Warning: cell_id {cell_id} no es convertible a int -> se ignora.")
            continue

        if cid < 1 or cid > total:
            print(f"Warning: cell_id {cid} fuera de rango [1, {total}] -> se ignora.")
            continue

        row = (cid - 1) // ncols   # 0-based
        col = (cid - 1) % ncols
        array[row, col] = value

    # asegurar carpeta de salida
    out_dir = os.path.dirname(file_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # preparar perfil de salida (basado en profile del input, pero con cambios)
    profile_out = profile
    profile_out.update({
        "driver": "GTiff",
        "dtype": out_dtype,
        "count": 1,
        "nodata": nodata,
        "height": nrows,
        "width": ncols,
        "crs": crs,
        "transform": transform
    })

    # escribir raster
    with rasterio.open(file_path, "w", **profile_out) as dst:
        dst.write(array.astype(out_dtype), 1)

    print(f"Raster escrito en: {file_path}")
    print(f"shape: {array.shape}, dtype: {array.dtype}, nodata: {nodata}")


def bp_calculation(fuels,csv_path,pickle_path,output_path,nsims,ncores):

    # Paso 1: Convertir CSV a Pickle
    print("Convirtiendo CSV a Pickle...")
    convert_csv_to_pickle(csv_path, pickle_path)

    # Paso 2: Procesar Pickles en paralelo
    print("Procesando Pickles en paralelo...")
    bp = process_pickle_parallel(pickle_path, nsims, ncores)

    # Paso 3: Escribir resultados
    print("writing results...")
    write_raster_from_dict_ids(fuels, output_path, bp,default_nodata=-9999)
    #print(f"Cálculo completado. Resultado guardado en: {output_path}")

    # Limpiar Pickle (opcional)
    print("cleaning pickle...")
    shutil.rmtree(pickle_path)

# Ejecución principal
if __name__ == "__main__":

    # Ejemplo de uso:
    ruta_base = "/Users/matiasvilches/Documents/F2A/papers/two_stage"
    #ruta_base = "/home/matias/Documents/ITREND"
    
    forest = "sub40"
    fuels_path = f"{ruta_base}/forest/{forest}/fuels.asc"
    csv_path = f"{ruta_base}/sims/{forest}_i01/Messages/"
    pickle_path = f"{ruta_base}/sims/{forest}/Pickles/"
    nsims = 1000
    ncores = 1
    output_path = f"{ruta_base}/results/{forest}/bp_{forest}_i01.tif"

    bp_calculation(fuels_path, csv_path, pickle_path, output_path, nsims, ncores)