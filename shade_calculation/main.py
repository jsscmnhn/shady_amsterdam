from gdown import download

from src.collect_data import download_raster_tiles, download_las_tiles
from src.create_first_chm_parallel import process_laz_files as process_laz_files_parallel
from src.create_dsm_and_chm_parallel import process_folders as process_dsm_chm_parallel
from src.shade_parallel import process_folders as process_shade
from src.merge_shademaps import merge_tif_files_by_time

import configparser
import os
from datetime import datetime

def read_config(file_path):
    params = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            key_value = line.split('=', 1)
            if len(key_value) == 2:
                key = key_value[0].strip()
                value = key_value[1].strip()

                # Store the parameters in the dictionary
                params[key] = value

    # Convert types where necessary
    if 'ndvi_threshold' in params:
        params['ndvi_threshold'] = float(params['ndvi_threshold'])
    if 'resolution' in params:
        params['resolution'] = float(params['resolution'])
    if 'filter_size' in params:
        params['filter_size'] = int(params['filter_size'])
    if 'chm_max_workers' in params:
        params['chm_max_workers'] = int(params['chm_max_workers'])
    if 'dsm_max_workers' in params:
        params['dsm_max_workers'] = int(params['dsm_max_workers'])
    if 'start_time' in params:
        params['start_time'] = int(params['start_time'])
    if 'end_time' in params:
        params['end_time'] = int(params['end_time'])
    if 'interval' in params:
        params['interval'] = int(params['interval'])
    if 'trans' in params:
        params['trans'] = int(params['trans'])
    if 'trunkheight' in params:
        params['trunkheight'] = int(params['trunkheight'])
    if 'files_start_time' in params:
        params['files_start_time'] = int(params['files_start_time'])
    if 'files_end_time' in params:
        params['files_end_time'] = int(params['files_end_time'])
    if 'date' in params:
        try:
            params['date'] = datetime.strptime(params[key], '%Y-%m-%d')
        except ValueError:
            raise ValueError(f"Invalid date format for {'date'}: {params['date']}")

    # Convert boolean strings to actual boolean values
    boolean_keys = ['download_las', 'download_dsm_dtm', 'create_chm',
                    'create_final_dsm_chm', 'create_shade', 'merge_shademaps',
                    'remove_las', 'smooth_chm', 'speed_up', 'delete_input']

    for key in boolean_keys:
        if key in params:
            params[key] = params[key].lower() == 'yes'

    return params


if __name__ == '__main__':

    config_file = 'Shade_config.txt'
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")

    params = read_config(config_file)

    if params.get('download_las'):
        print("Downloading LAS tiles...")
        download_las_tiles(params['subtile_list_file'], params['las_output_folder'])

    if params.get('download_dsm_dtm'):
        print("Downloading and merging DSM and DTM tiles...")
        download_raster_tiles(params['tile_list_file'], params['ahn_output_folder'], params['merged_name'])

    #TODO MAKE SURE THAT OPTIONAL CONFIG ARE MADE EVEN IF NOT INPUT BY USER TO MAKE SURE NO VALUE ERROR
    if params.get('create_chm'):
        print("Creating first CHM files...")
        process_laz_files_parallel(params['las_output_folder'], params['chm_output_folder'], params['ndvi_threshold'],
                                   params['resolution'], params['remove_las'], params['smooth_chm'],
                                   params['filter_size'], params['chm_max_workers'])

    if params.get('create_final_dsm_chm'):
        # First check if a merged dtm and dsm is already provided
        output_folder = params['ahn_output_folder']
        name = params['merged_name']

        merged_dtm_output_file = os.path.join(output_folder, f"{name}_DTM.TIF")
        merged_dsm_output_file = os.path.join(output_folder, f"{name}_DSM.TIF")

        if params.get('merged_dtm') is None:
            params['merged_dtm'] = merged_dtm_output_file
        if params.get('merged_dsm') is None:
            params['merged_dsm'] = merged_dsm_output_file

        # Then we can start processing
        print("Creating final DSM and DTM files...")
        process_dsm_chm_parallel(params['chm_output_folder'], params['merged_dtm'], params['merged_dsm'],
                                 params['buildings_path'], params['output_dsm_chm'], params['dsm_max_workers'],
                                 params['speed_up'], params['min_vegetation_height'], params['max_vegetation_height'])

    if params.get('create_shade'):
        print("Creating shade maps...")
        process_shade(params['output_dsm_chm'], params['output_base_shademap'], params['date'], params['start_time'],
                      params['end_time'], params['interval'], params['trans'], params['trunkheight'],
                      params['max_shade_workers'])

    if params.get('merge_shademaps'):
        print("Merging shademaps...")
        merge_tif_files_by_time(params['output_base_shademap'], params['output_folder_merged_shademaps'],
                                params['merged_name'], params['files_start_time'], params['files_end_time'],
                                params['delete_input_shade'])
