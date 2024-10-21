"""
This file contains the functions to collect all the data needed to create the CHM and DSM.
"""

import rasterio
from rasterio.merge import merge
import glob
import os
import requests
import zipfile

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
    - none: The function writes the merged file directly to the specified `output_file`.
    """
    # Find TIF files that match the given prefix
    tif_files = glob.glob(os.path.join(input_folder, f'{file_prefix}*.TIF'))
    bounds = []

    if not tif_files:
        print(f"No TIF files found for prefix {file_prefix}.")
        return

    src_files_to_mosaic = []

    for tif_file in tif_files:
        src = rasterio.open(tif_file)
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


def download_raster_tiles(tile_list_file, output_folder, name):
    """
    Download DSM and DTM files for each tile specified in a text file.
    ------
    Input:
    - tile_list_file (str): Path to the text file containing the list of tiles to download.
    - output_folder (str): Directory where the downloaded and unzipped files will be saved.
    - name (str): Name of the output raster file.
    Output:
    - none: The function writes output files directly to the specified `output_folder`.
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
    merge_tif_files(output_folder, merged_dsm_output_file, "R_")

    print(f"All DTM and DSM files processed and merged.")
