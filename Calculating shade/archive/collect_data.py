"""
This file contains the functions to collect all the data needed to create the CHM and DSM.
"""

import rasterio
from rasterio.merge import merge
import glob
import cProfile
import io
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
        # Extract tile name and sub-tile number (e.g., "25BZ1_01")
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


def download_raster_tiles(tile_list_file, output_folder):
    """
    Download DSM and DTM files for each tile specified in a text file.
    ------
    Input:
    - tile_list_file (str): Path to the text file containing the list of tiles to download.
    - output_folder (str): Directory where the downloaded and unzipped files will be saved.
    Output:
    - none:  The function writes output files directly to the specified `output_folder`.
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