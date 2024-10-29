"""
This file contains the functions to collect all the data needed to create the CHM and DSM: downloading the LAZ subtiles,
downloading the AHN tiles (DTM and DSM), downloading the building geometries from 3DBAG.
"""

import rasterio
from rasterio.merge import merge
import glob
import os
import requests
import zipfile
import pandas as pd
import geopandas as gpd
import gzip
from io import BytesIO

def download_las_tiles(tile_list_file, output_folder):
    """
    Download LAZ files for each subtile specified in a text file.
    ------
    Input:
    - tile_list_file (str): Path to the text file containing the list of subtiles to download.
    - output_folder (str): Path to the folder where the downloaded files should be saved.
    Output:
    - none:  The function writes output files directly to the specified `output_folder{subtile}`.
    """

    base_url = "https://geotiles.citg.tudelft.nl/AHN4_T"
    os.makedirs(output_folder, exist_ok=True)

    # Read tile names from the text file
    with open(tile_list_file, 'r') as f:
        tile_names = [line.strip() for line in f.readlines() if line.strip()]

    for full_tile_name in tile_names:
        # Extract tile name and sub-tile number
        if '_' in full_tile_name:
            tile_name, sub_tile = full_tile_name.split('_')
        else:
            print(f"Skipping invalid tile entry: {full_tile_name}")
            continue

        # Define folder structure within the specified output folder
        tile_folder = os.path.join(output_folder, tile_name)
        os.makedirs(tile_folder, exist_ok=True)

        # Construct the URL and file path
        sub_tile_str = f"_{int(sub_tile):02}"  # Convert sub_tile to zero-padded two-digit number
        url = f"{base_url}/{tile_name}{sub_tile_str}.LAZ"
        filename = f"{tile_name}{sub_tile_str}.LAZ"
        file_path = os.path.join(tile_folder, filename)

        # Check if the file already exists to avoid re-downloading
        if os.path.exists(file_path):
            print(f"File {file_path} already exists, skipping download.")
            continue

        # Attempt to download the file
        try:
            response = requests.get(url)
            response.raise_for_status()  # Ensure the request was successful
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded and saved {file_path}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {url}: {e}")


def download_and_extract(url, file_path, output_folder):
    """Download a ZIP file from a URL, extract its contents, and delete the ZIP file."""
    if not os.path.exists(file_path):
        try:
            print(f"Downloading from {url}...")
            response = requests.get(url, verify=False)  # Bypass SSL verification
            response.raise_for_status()  # Raise an error for bad responses
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded and saved to {file_path}")

            # Extract the ZIP file
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(output_folder)  # Extract to the output folder
            print(f"Extracted {file_path}")

            # Remove the ZIP file after extraction
            os.remove(file_path)
            print(f"Deleted the ZIP file: {file_path}")

        except requests.exceptions.RequestException as e:
            print(f"Failed to download {url}: {e}")
        except zipfile.BadZipFile as e:
            print(f"Failed to extract {file_path}: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")


def merge_tif_files(input_folder, output_file, file_prefix, nodata_value=-9999):
    """Merge TIF files in a folder with a specific prefix into a single raster file.
    ------
    Input:
    - input_folder (str): Path to the folder containing the TIF files.
    - output_file (str): Path to the output raster file.
    - file_prefix (str): Prefix of files to merge ("M_" for DTM or "R_" for DSM).
    - nodata_value (int, optional): Value to replace the nodata value. Defaults to -9999.

    Output:
    - tile_bounds (list): List containing the name of the tile and the extent of that tile, for downloading the buidling data.
    - The function writes the merged file directly to the specified `output_file`.
    """
    # Find TIF files that match the given prefix
    tif_files = glob.glob(os.path.join(input_folder, f'{file_prefix}*.TIF'))

    tile_bounds = []

    if not tif_files:
        print(f"No TIF files found for prefix {file_prefix}.")
        return

    src_files_to_mosaic = []

    for tif_file in tif_files:
        src = rasterio.open(tif_file)

        # Extract the tile name &
        base_name = os.path.basename(tif_file)
        tile_name = base_name.split('_', 1)[-1].split('.')[0]
        tile_bounds.append([tile_name, src.bounds])

        src_files_to_mosaic.append(src)

    mosaic, out_transform = merge(src_files_to_mosaic)

    # Replace existing nodata values with new nodata value
    for src in src_files_to_mosaic:
        if 'nodata' in src.meta:
            mosaic[mosaic == src.nodata] = nodata_value

    out_meta = src_files_to_mosaic[0].meta.copy()

    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform,
        "count": mosaic.shape[0],
        "nodata": nodata_value
    })

    with rasterio.open(output_file, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Merged {len(tif_files)} TIF files with prefix {file_prefix} into {output_file}")

    # close & delete original TIF files & save bounds
    for src in src_files_to_mosaic:
        src.close()
    for tif_file in tif_files:
        os.remove(tif_file)
        print(f"Deleted original TIF file: {tif_file}")

    return tile_bounds


def download_raster_tiles(tile_list_file, output_folder, name):
    """
    Download DSM and DTM files for each tile specified in a text file.
    ------
    Input:
    - tile_list_file (str): Path to the text file containing the list of tiles to download.
    - output_folder (str): Directory where the downloaded and unzipped files will be saved.
    - name (str): Name of the output raster file.
    Output:
    - tile_bounds (list): List containing the name of the tile and the extent of that tile, for downloading the buidling data.
    - The function writes output files directly to the specified `output_folder`.
    """

    # Base URLs for downloading DTM and DSM tiles
    base_url_dtm = "https://ns_hwh.fundaments.nl/hwh-ahn/ahn4/02a_DTM_0.5m"
    base_url_dsm = "https://ns_hwh.fundaments.nl/hwh-ahn/ahn4/03a_DSM_0.5m"

    # Read tile names from the text file
    with open(tile_list_file, 'r') as f:
        tile_names = [line.strip() for line in f.readlines() if line.strip()]

    os.makedirs(output_folder, exist_ok=True)

    for tile_name in tile_names:
        # Construct the URLs for DTM and DSM
        dtm_url = f"{base_url_dtm}/M_{tile_name}.zip"
        dsm_url = f"{base_url_dsm}/R_{tile_name}.zip"

        # Define file paths for downloaded zip files
        dtm_file_path = os.path.join(output_folder, f"M_{tile_name}.zip")
        dsm_file_path = os.path.join(output_folder, f"R_{tile_name}.zip")

        # Download and extract DTM and DSM files
        download_and_extract(dtm_url, dtm_file_path, output_folder)
        download_and_extract(dsm_url, dsm_file_path, output_folder)

    # Merge all DTM and DSM files separately after all downloads are done
    merged_dtm_output_file = os.path.join(output_folder, f"{name}_DTM.TIF")
    merged_dsm_output_file = os.path.join(output_folder, f"{name}_DSM.TIF")

    # Merge files starting with "M_" for DTM and "R_" for DSM
    merge_tif_files(output_folder, merged_dtm_output_file, "M_")
    tile_bounds = merge_tif_files(output_folder, merged_dsm_output_file, "R_")

    print(f"All DTM and DSM files processed and merged.")

    return tile_bounds


def download_wfs_data(wfs_url, layer_name, bbox, gpkg_name, output_folder, tile_name):
    """
    Download data from a WFS server in batches and save it to a GeoPackage.
    -----------------------------------------------------
    Input:
    -   wfs_url (str): URL of the WFS service.
    -   layer_name (str): The layer name to download.
    -   bbox (tuple): Bounding box as (minx, miny, maxx, maxy).
    -   gpkg_name (str): Name for the output GeoPackage file.
    -   tile_name (str): Layer name for saving in the GeoPackage.
    Output:
    -   None: saves a GeoPackage file to the given {output_gpkg} at layer {tile_name}.
    """
    # Initialize variables for feature collection, max requestable amount from server is 10000
    all_features = []
    start_index = 0
    count = 10000

    while True:
        params = {
            "SERVICE": "WFS",
            "REQUEST": "GetFeature",
            "VERSION": "2.0.0",
            "TYPENAMES": layer_name,
            "SRSNAME": "urn:ogc:def:crs:EPSG::28992",
            "BBOX": f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]},urn:ogc:def:crs:EPSG::28992",
            "COUNT": count,
            "STARTINDEX": start_index
        }

        # Mimicking a QGIS request
        headers = {
            "User-Agent": "Mozilla/5.0 QGIS/33411/Windows 11 Version 2009"
        }

        response = requests.get(wfs_url, params=params, headers=headers)

        # Check if the request was successful & donwload data
        if response.status_code == 200:
            if response.headers.get('Content-Encoding', '').lower() == 'gzip' and response.content[:2] == b'\x1f\x8b':
                data = gzip.decompress(response.content)
            else:
                data = response.content

            with BytesIO(data) as f:
                gdf = gpd.read_file(f)

            all_features.append(gdf)

            # Check if the number of features retrieved is less than the requested count: then we can stop
            if len(gdf) < count:
                break

                # Start index for next request
            start_index += count

        else:
            print(f"Failed to download WFS data. Status code: {response.status_code}")
            print(f"Error message: {response.text}")
            break

            # Concatenate all features into a single GeoDataFrame
    if all_features:
        full_gdf = gpd.GeoDataFrame(pd.concat(all_features, ignore_index=True))

        os.makedirs(output_folder, exist_ok=True)
        output_gpkg = os.path.join(output_folder, f"{gpkg_name}.gpkg")

        # Save the GeoDataFrame to a GeoPackage with the specified tile name
        full_gdf.to_file(output_gpkg, layer=tile_name, driver="GPKG")
        print(f"Downloaded and saved layer '{tile_name}' to {output_gpkg}")

        # Saving to the GeoPackage with tile name
        os.path.join(output_folder, )
        full_gdf.to_file(output_gpkg, layer=tile_name, driver="GPKG")
        print(f"Downloaded and saved layer '{tile_name}' to {output_gpkg}")
    else:
        print("No features were downloaded.")

def setup_WFS_download(gpkg_name, tile_bounds, output_folder):
    """
    Collecting the needed information for downloading the WFS data.
    --------------
    input:
    -   gpkg_name (str): Name for the output GeoPackage file.
    -   tile_bounds (list): list containing nested lists where at [0] the tile name is saved and at [1] the bounds of the tile
    -   output_folder (str): Path to the output GeoPackage file.
    Output:
    -   none: saves a GeoPackage file to the given {gpkg_name} at layer {tile_name}.
    """
    wfs_url = "https://data.3dbag.nl/api/BAG3D/wfs"
    layer_name = "BAG3D:lod13"
    for i in range(len(tile_bounds)):
        tile_name = tile_bounds[i][0]
        bbox = tile_bounds[i][1]
        bbox_tuple = (bbox.left - 40, bbox.bottom - 40, bbox.right + 40, bbox.top + 40)
        print(bbox_tuple)

        download_wfs_data(
            wfs_url=wfs_url,
            layer_name=layer_name,
            bbox=bbox_tuple,
            gpkg_name=gpkg_name,
            output_folder=output_folder,
            tile_name=tile_name
        )