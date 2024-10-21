from src.collect_data import download_raster_tiles, download_las_tiles, setup_WFS_download
from src.create_first_chm_parallel import process_laz_files as process_laz_files_parallel
from src.create_dsm_and_chm_parallel import process_folders as process_dsm_chm_parallel
from src.shade_parallel import process_folders as process_shade
from src.merge_shademaps import merge_tif_files_by_time

import argparse
import os
import json
from datetime import datetime


def flatten_dict(nested_dict, parent_key='', sep='_'):
    """
    Flatten a nested dictionary ignoring the first layer.
    ------------------------------------------------------------------
    Input:
    - nested_dict (dict): The dictionary to flatten.
    - parent_key (str): The base key string for the flattened keys.
    - sep (str): Separator to use between concatenated keys.

    Output:
    - items: A flat dictionary.
    """
    items = {}

    # flattening the nested dictionary without the first layers to get the required parameters
    for key, value in nested_dict.items():
        if isinstance(value, dict):
            for sub_key, sub_value in flatten_dict(value, parent_key=key, sep=sep).items():
                items[sub_key] = sub_value
        else:
            items[key] = value

    return items

def read_config(file_path):
    """
    Read the contents of the JSON configuration file and save them to a dictionary
    -------------------------------------
    Input:
    - file_path (str): The path to the configuration file.
    output:
    - params (dict): The configuration dictionary with all the parameters needed.
    """
    # Load parameters from the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

        # Flatten the nested structure while ignoring the first layer
        params = {}
        for key in data:
            params.update(flatten_dict(data[key]))

    # Default values to use if not provided in the JSON file
    default_values = {
        'ndvi_threshold': 0.0,
        'resolution': 0.5,
        'filter_size': 3,
        'chm_max_workers': 4,
        'dsm_max_workers': 4,
        'min_vegetation_height': 2,
        'max_vegetation_height': 40,
        'max_shade_workers': 4,
        'start_time': 9,
        'end_time': 20,
        'interval': 30,
        'trans': 10,
        'trunkheight': 25,
        'files_start_time': 900,
        'files_end_time': 2000,
        'remove_las': False,
        'smooth_chm': False,
        'pre_filter': False,
        'speed_up': False,
        'use_chm' : True,
        'delete_input_shade': False
    }

    # Set default values for missing parameters
    for key, default_value in default_values.items():
        params.setdefault(key, default_value)

    # Convert the 'date' parameter to a datetime object if it's provided as a string
    if 'date' in params and isinstance(params['date'], str):
        try:
            params['date'] = datetime.strptime(params['date'], '%Y-%m-%d')
        except ValueError:
            raise ValueError(f"Invalid date format for 'date': {params['date']}")

    return params


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creation of Canopy Height Model (CHM), Digital Surface Model '
                                                 'with ground and buildings (DSM) and shade maps based on AHN '
                                                 'data through a configuration JSON file.')
    parser.add_argument('config_file', type=str, help='Path to the configuration file')

    args = parser.parse_args()

    config_file = args.config_file

    # Check if the config file exists
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")

    # Read parameters from the config file
    params = read_config(config_file)
    for param in params:
        print(param + ' ' + str(params[param]))


    if params.get('download_las'):
        print("Downloading LAS tiles...")
        download_las_tiles(params['subtile_list_file'], params['las_output_folder'])

    if params.get('download_dsm_dtm'):
        print("Downloading and merging DSM and DTM tiles...")
        tile_bounds = download_raster_tiles(params['tile_list_file'], params['ahn_output_folder'], params['merged_name'])
        if params.get('download_buildings'):
            print("Downloading buildings geometries...")
            setup_WFS_download( params["buildings_name"], tile_bounds, params['buildings_output_folder'])

    if params.get('create_chm'):
        print("Creating first CHM files...")
        process_laz_files_parallel(params['las_output_folder'], params['chm_output_folder'], params['ndvi_threshold'],
                                   params['resolution'], params['remove_las'], params['smooth_chm'],
                                   params['filter_size'],  params['pre_filter'], params['chm_max_workers'])

    if params.get('create_final_dsm_chm'):
        # First check if a merged dtm and dsm is already provided, and a buildings geopackage
        output_folder = params['ahn_output_folder']
        name = params['merged_name']

        buildings_output_folder = params['buildings_output_folder']
        buildings_name = params['buildings_name']

        merged_dtm_output_file = os.path.join(output_folder, f"{name}_DTM.TIF")
        merged_dsm_output_file = os.path.join(output_folder, f"{name}_DSM.TIF")
        buildings_file = os.path.join(buildings_output_folder, f"{buildings_name}.gpkg")

        if params.get('merged_dtm') is None or "path":
            params['merged_dtm'] = merged_dtm_output_file
        if params.get('merged_dsm') is None or "path":
            params['merged_dsm'] = merged_dsm_output_file
        if params.get('buildings_path') is None or "path":
            params['buildings_path'] = buildings_file

        # Then we can start processing
        print("Creating final DSM and DTM files...")
        process_dsm_chm_parallel(params['chm_output_folder'], params['merged_dtm'], params['merged_dsm'],
                                 params['buildings_path'], params['output_dsm_chm'], params['dsm_max_workers'],
                                 params['speed_up'], params['min_vegetation_height'], params['max_vegetation_height'])

    if params.get('create_shade'):
        print("Creating shade maps...")
        process_shade(params['output_dsm_chm'], params['output_base_shademap'], params['date'], params['start_time'],
                      params['end_time'], params['interval'], params['use_chm'], params['trans'], params['trunkheight'],
                      params['max_shade_workers'])

    if params.get('merge_shademaps'):
        print("Merging shademaps...")
        print(params["delete_input_shade"])
        merge_tif_files_by_time(params['output_base_shademap'], params['output_folder_merged_shademaps'],
                                params['merged_name'], params['files_start_time'], params['files_end_time'],
                                params['delete_input_shade'])
