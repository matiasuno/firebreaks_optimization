import os
import numpy as np
import rasterio
import pandas as pd
import numpy as np
from rasterio import features
import fiona

# read raster using rasterio
def read_raster(file_path):
    with rasterio.open(file_path) as src:
        data = src.read(1)  # Read the first band
        profile = src.profile
        if data.ndim != 2:
            raise ValueError(f"Expected a 2D array, but got {data.ndim}D array.")
        xdim, ydim = data.shape
        nodos = int(xdim*ydim)
        
    return data, profile, nodos

# write raster from array using rasterio
def write_raster(file_path, array, reference_raster):
    with rasterio.open(reference_raster) as src:
        profile = src.profile
        profile.update(
            driver='AAIGrid',
            dtype=rasterio.float32,
            nodata=-9999,
            count=1,
            compress=None
        )

        with rasterio.open(file_path, 'w', **profile) as dst:
            dst.write(array.astype(rasterio.float32), 1)



def write_raster_from_id_list(base_raster, output_path, id_list, value=1.0, nodata=0.0):
    """
    Crea un raster a partir de una lista de IDs de celdas, usando un raster base
    para la geometría, transform, CRS y nodata.

    Args:
        base_raster (str): Ruta al raster base (.tif o .asc).
        output_path (str): Ruta del archivo de salida (.tif o .asc).
        id_list (list[int]): Lista de IDs de celdas (1 = (0,0)).
        value (float): Valor a asignar a las celdas de la lista. Default 1.0.
        nodata (float): Valor para las celdas vacías. Default 0.0.
    """
    with rasterio.open(base_raster) as src:
        nrows, ncols = src.height, src.width
        transform = src.transform
        crs = src.crs

    # Crear una matriz llena con nodata
    array = np.full((nrows, ncols), nodata, dtype=np.float32)

    # Asignar el valor a las celdas indicadas
    for cell_id in id_list:
        if cell_id <= 0 or cell_id > nrows * ncols:
            continue  # ignorar IDs fuera de rango
        row = (cell_id - 1) // ncols
        col = (cell_id - 1) % ncols
        array[row, col] = value

    # Elegir driver según extensión
    driver = "AAIGrid" if output_path.lower().endswith(".asc") else "GTiff"

    profile = {
        "driver": driver,
        "dtype": "float32",
        "nodata": nodata,
        "width": ncols,
        "height": nrows,
        "count": 1,
        "crs": crs,
        "transform": transform
    }

    # Escribir raster
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(array, 1)

    print(f"Raster escrito en {output_path} con {len(id_list)} celdas activas.")


def write_raster_from_dict_ids(base_raster, file_path, raster_dict=None, csv_file=None):
    """
    Crear un raster a partir de un diccionario o un archivo CSV con IDs de celdas -> valores,
    usando un raster base para la geometría, transform, CRS y nodata.
    
    Args:
        base_raster (str): Ruta del raster base (.tif o .asc).
        file_path (str): Ruta de salida del nuevo raster.
        raster_dict (dict, optional): Diccionario {cell_id: valor}.
        csv_file (str, optional): CSV con columnas [cell_id, value].
    """
    if raster_dict is None and csv_file is None:
        raise ValueError("Debe especificar raster_dict o csv_file")

    # Cargar datos desde CSV si se da un archivo
    if csv_file:
        df = pd.read_csv(csv_file)
        raster_dict = dict(zip(df["cell_id"], df["value"]))

    with rasterio.open(base_raster) as src:
        nrows, ncols = src.height, src.width
        transform = src.transform
        crs = src.crs
        nodata = src.nodata if src.nodata is not None else -9999
        dtype = src.dtypes[0]

    # inicializar matriz con NODATA
    array = np.full((nrows, ncols), nodata, dtype=dtype)

    # asignar valores
    for cell_id, value in raster_dict.items():
        row = (cell_id - 1) // ncols
        col = (cell_id - 1) % ncols
        array[row, col] = value

    # escribir nuevo raster
    with rasterio.open(
        file_path,
        "w",
        driver="GTiff",
        height=nrows,
        width=ncols,
        count=1,
        dtype=dtype,
        crs=crs,
        transform=transform,
        nodata=nodata
    ) as dst:
        dst.write(array, 1)

def average_asc_files(asc_files, output_path):

    # Obtener lista de archivos .asc
    #asc_files = sorted(glob.glob(os.path.join(input_folder, '*.asc')))

    # Lista para guardar los datos de cada raster
    arrays = []

    # Leer y apilar los datos
    for f in asc_files:
        print(f"Procesando archivo: {f}")
        with rasterio.open(f) as src:
            data = src.read(1).astype(float)
            arrays.append(data)
            profile = src.profile  # Guardamos el profile del primero

    # Stack y promedio
    stack = np.stack(arrays)
    mean_array = np.max(stack, axis=0)

    # Actualizamos el profile para guardar como ASCII
    profile.update(
        driver='AAIGrid',
        dtype='float32',
        count=1
    )

    # Guardar el raster promedio
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(mean_array.astype(np.float32), 1)

    print(f'Promedio guardado en: {output_path}')


def sum_raster_values(raster_path):
    with rasterio.open(raster_path) as src:
        data = src.read(1)  # Read the first band (assuming single-band raster)
        total_sum = np.nansum(data)  # Sum all values, ignoring NaNs (if any)
    return float(total_sum)

def raster_to_dict(raster_path):
    with rasterio.open(raster_path) as src:
        raster_data = src.read(1)  # Leer la primera banda
        rows, cols = raster_data.shape

        # Crear el diccionario con el ID de la celda como clave
        raster_dict = {}
        cell_id = 1

        for row in range(rows):
            for col in range(cols):
                raster_dict[cell_id] = raster_data[row, col]
                cell_id += 1
    return raster_dict

def remove_last_two_rows(directory):
    """Remove the last two rows from all CSV files in a directory."""
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):  # Only process CSV files
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath)
            
            if len(df) > 2:
                df = df.iloc[:-2]  # Remove last two rows
            else:
                print(f"Skipping {filename}, not enough rows.")
                continue
            
            df.to_csv(filepath, index=False)  # Save without index
            print(f"Updated: {filename}")

def create_ascii_with_ones(input_ascii, output_ascii):
    """Crea un archivo ASCII con todos los valores establecidos en 1, basado en un archivo ASCII de entrada."""
    # === Leer el archivo ASCII original ===
    with rasterio.open(input_ascii) as src:
        data = src.read(1)  # Leer los datos originales (solo para obtener la forma)
        profile = src.profile.copy()  # Copiar metadatos

        # === Crear un array del mismo tamaño lleno de 1s ===
        new_data = np.ones_like(data, dtype=np.int32)

        # === Actualizar el perfil para salida ASCII ===
        profile.update(
            driver='AAIGrid',
            dtype=rasterio.int32,
            nodata=-9999,
            count=1,
            compress=None
        )

        # === Escribir el nuevo archivo con todos 1s ===
        with rasterio.open(output_ascii, 'w', **profile) as dst:
            dst.write(new_data, 1)

    print(f'Archivo ASCII generado con todos los valores = 1: {output_ascii}')


def multiplicar_raster_por_valor(input_path, factor):
    with rasterio.open(input_path) as src:
        data = src.read(1)
        nodata = src.nodata

        # Proteger nodata si está definido
        if nodata is not None:
            resultado = np.where(data != nodata, data * factor, nodata)
        else:
            resultado = data * factor

    return resultado

def reciproco_raster_por_valor(input_path: str,
                               output_path: str,
                               zero_replace: float = 100000.0):
    """
    Lee un raster, calcula 1/value donde value != 0 y value != nodata.
    Donde value == 0 o value == nodata deja el valor `zero_replace` (o nodata
    si la celda era nodata) y guarda un raster de salida.
    
    Args:
        input_path: ruta al raster de entrada.
        output_path: ruta al raster de salida (.tif o .asc compatible).
        zero_replace: valor a asignar cuando value == 0 (por defecto 100000.0).
    """
    # Abrir como masked array para que nodata se convierta en máscara
    with rasterio.open(input_path) as src:
        # Leer banda 1 como masked array: las celdas nodata estarán en .mask == True
        band = src.read(1, masked=True)
        src_crs = src.crs
        src_transform = src.transform
        src_width = src.width
        src_height = src.height

    # band es un numpy MaskedArray
    data = band.data.astype(np.float32)      # los datos (con posible contenido nodata)
    mask = band.mask                         # True donde era nodata

    # Inicializar resultado con el valor por defecto (zero_replace)
    resultado = np.full(data.shape, fill_value=float(zero_replace), dtype=np.float32)

    # Condición para calcular reciprocidad: no nodata y no cero
    valid = (~mask) & (data != 0)

    # Evitar división por cero calculando solo donde valid == True
    if np.any(valid):
        resultado[valid] = 1.0 / data[valid]

    # Para las celdas que eran nodata, queremos mantener nodata o asignar
    # un valor de nodata en el archivo de salida. Vamos a definir nodata_out.
    nodata_out = -9999.0  # valor por defecto de nodata de salida

    # Preparar perfil limpio para evitar campos incompatibles con drivers como AAIGrid
    # Elegir driver según extensión de salida
    if output_path.lower().endswith(".asc"):
        driver = "AAIGrid"
    else:
        driver = "GTiff"

    profile = {
        "driver": driver,
        "dtype": "float32",
        "nodata": nodata_out,
        "width": src_width,
        "height": src_height,
        "count": 1,
        "crs": src_crs,
        "transform": src_transform
    }

    # Poner nodata_out en las posiciones originalmente nodata para que se escriban como nodata
    resultado_masked = resultado.copy()
    resultado_masked[mask] = nodata_out

    # Escribir el raster
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(resultado_masked.astype(np.float32), 1)

    print(f"Output written to {output_path}. nodata={nodata_out}, zero replaced by {zero_replace}")

def reciproco_raster_multibanda(input_path: str,
                                output_path: str,
                                zero_replace: float = 100000.0):
    """
    Lee un raster multibanda, calcula 1/value donde value != 0 y value != nodata
    en cada banda. Donde value == 0 o value == nodata deja el valor `zero_replace`
    (o nodata si la celda era nodata) y guarda un raster multibanda de salida.
    
    Args:
        input_path: ruta al raster de entrada.
        output_path: ruta al raster de salida (.tif o .asc compatible).
        zero_replace: valor a asignar cuando value == 0 (por defecto 100000.0).
    """
    with rasterio.open(input_path) as src:
        src_crs = src.crs
        src_transform = src.transform
        src_width = src.width
        src_height = src.height
        src_count = src.count
        src_nodata = src.nodata

        # definir un nodata para salida (puede ser el mismo que input si existe)
        nodata_out = src_nodata if src_nodata is not None else -9999.0

        if output_path.lower().endswith(".asc"):
            driver = "AAIGrid"
            if src_count > 1:
                raise ValueError("El formato .asc solo admite una banda, no multibanda")
        else:
            driver = "GTiff"

        profile = {
            "driver": driver,
            "dtype": "float32",
            "nodata": nodata_out,
            "width": src_width,
            "height": src_height,
            "count": src_count,
            "crs": src_crs,
            "transform": src_transform
        }

        with rasterio.open(output_path, "w", **profile) as dst:
            for i in range(1, src_count + 1):
                band = src.read(i, masked=True)   # leer banda i como masked array

                data = band.data.astype(np.float32)
                mask = band.mask

                # inicializar resultado con zero_replace
                resultado = np.full(data.shape, float(zero_replace), dtype=np.float32)

                # condición válida: no nodata y no cero
                valid = (~mask) & (data != 0)
                if np.any(valid):
                    resultado[valid] = 1.0 / data[valid]

                # poner nodata donde corresponde
                resultado[mask] = nodata_out

                dst.write(resultado, i)

    print(f"Output written to {output_path}. nodata={nodata_out}, zero replaced by {zero_replace}")


def calculate_raster_mean(raster_path):
    """
    Calculate the mean of values in a raster file, ignoring NaN values.
    
    Parameters:
    raster_path (str): Path to the raster file
    
    Returns:
    float: Mean value of the raster (excluding NaN values)
    """
    try:
        # Open the raster file
        with rasterio.open(raster_path) as src:
            # Read all bands into a numpy array
            data = src.read()
            
            # Convert to float to handle NaN values properly
            data = data.astype(float)
            
            # Replace NoData values with NaN (if specified in metadata)
            if src.nodata is not None:
                data[data == src.nodata] = np.nan
            
            # Calculate mean ignoring NaN values
            mean_value = np.nanmean(data)

            # calculate max
            max_value = np.nanmax(data)

            forest = os.path.basename(raster_path)
            print(f"Mean, max value of {forest} : {mean_value:.4f}, {max_value:.4f}")
            return mean_value,max_value
            
    except Exception as e:
        print(f"Error processing raster file: {e}")
        return None    

def get_ids_in_and_out_polygon(raster_path: str, polygon_path: str, layer: str = None):
    """
    Retorna dos listas de IDs (basados en 1):
    - ids_inside: celdas cubiertas por el polígono
    - ids_outside: celdas fuera del polígono

    Args:
        raster_path (str): ruta al raster base.
        polygon_path (str): ruta al polígono (puede ser .shp o .gpkg).
        layer (str): nombre del layer en caso de .gpkg con múltiples capas (opcional).
    """
    # Abrir raster
    with rasterio.open(raster_path) as src:
        nrows, ncols = src.height, src.width
        transform = src.transform
        crs = src.crs

    # Leer geometrías del polígono
    with fiona.open(polygon_path, layer=layer) as shp:
        geometries = [feature["geometry"] for feature in shp]

    # Crear máscara booleana (True dentro del polígono)
    mask = features.geometry_mask(
        geometries,
        transform=transform,
        invert=True,  # True = dentro del polígono
        out_shape=(nrows, ncols)
    )

    # Crear arrays de índices
    rows, cols = np.indices(mask.shape)
    ids = rows * ncols + cols + 1  # IDs desde 1 (0,0) = 1

    # Extraer listas según la máscara
    ids_inside = ids[mask].tolist()
    ids_outside = ids[~mask].tolist()

    return ids_inside, ids_outside

def write_treatment_raster(fuels,firebreaks,output_path):

    #fuels: .asc file of forest
    #firebreaks: list of firebreaks
    #output_path = path in which will be stored dpv.asc file

    data, profile, nodos = read_raster(fuels)

    nrows = profile["height"]
    ncols = profile["width"]
    data = np.zeros((nrows, ncols))

    for n in firebreaks:

        row = n // ncols  # Calculate row index
        col = n % ncols   # Calculate column index

        data[row][col-1] = 1

    write_raster(output_path,data,fuels)