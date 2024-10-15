from src.collect_data import download_raster_tiles, download_las_tiles
from src.create_first_chm import process_laz_files
from src.create_first_chm_parallel import process_laz_files as process_laz_files_parallel
from src.create_dsm_and_chm_parallel import process_folders as process_dsm_chm_parallel
from src.shade_parallel import process_folders as process_shade
from src.merge_shademaps import merge_tif_files_by_time

import argparse
from datetime import datetime


def optional_params(input_str):
    """Parse the input string for optional parameters."""
    params = {}
    for param in input_str.split(","):
        try:
            key, value = param.split("=")
            params[key.strip()] = float(value.strip()) if '.' in value else int(value.strip())
        except ValueError:
            print(f"Invalid parameter format for '{param}'. Expected format: key=value.")
    return params

def main():
    parser = argparse.ArgumentParser(description='Creating a shade map from a CHM and DSM file.')
    args = parser.parse_args()  # No arguments needed now

    user_choice = input(
        "Do you want to perform vegetation extraction, ground extraction, or CHM creation? Enter 'vegetation', 'ground', or 'chm': ").strip().lower()

    # Check for existing CHM and DSM folder
    existing_chm_dsm_path = input(
        "Is there already a main folder with subfolders containing the CHM and DSM for each subtile? (yes/no): ").strip().lower()

    if existing_chm_dsm_path == 'yes':
        chm_dsm_folder = input("Provide path to this folder: ")
        output_base_folder = input("Provide the output base folder path: ")

        while True:
            date_input = input("Provide the date for processing (format YYYY-MM-DD): ")
            try:
                date = datetime.strptime(date_input, "%Y-%m-%d")  # Convert string to datetime object
                break  # Exit the loop if the date is valid
            except ValueError:
                print("Invalid date format. Please enter the date in YYYY-MM-DD format.")

        # Ask user if they want to input optional values
        optional_input = input("Do you want to input optional parameters? (yes/no): ").strip().lower()

        # Set default values for optional parameters
        params = {
            'start_time': 10,
            'end_time': 21,
            'interval': 30,
            'trans': 10,
            'trunkheight': 25,
            'max_workers': 4
        }

        if optional_input == 'yes':
            input_params = input(
                "Enter optional parameters in the format: start_time={int}, end_time={int}, "
                "interval={int}, trans={int}, trunkheight={int}, max_workers={int} ")
            params.update(optional_params(input_params))

        # Call shade processing function
        process_shade(
            chm_dsm_folder,
            output_base_folder,
            date,
            start_time=params['start_time'],
            end_time=params['end_time'],
            interval=params['interval'],
            trans=params['trans'],
            trunkheight=params['trunkheight'],
            max_workers=params['max_workers']
        )

    else:
        # Check if first CHM was created from LAZ data
        direct_chm_creation = input("Has there already been a CHM created from LAZ data? (yes/no): ").strip().lower()
        if direct_chm_creation == 'yes':
            merged_dsm_dtm = input("Is there already a merged DSM and DTM tif file? (yes/no): ").strip().lower()
            if merged_dsm_dtm == 'yes':
                dsm_path = input("Provide the path to the merged DSM file: ")
                dtm_path = input("Provide the path to the merged DTM file: ")
            else:
                # TODO: ADD FUNCTIONALITY TO RUN DSM & DTM DOWNLOAD AND MERGE
                dsm_path = input("Provide the path to the merged DSM file: ")
                dtm_path = input("Provide the path to the merged DTM file: ")

            chm_folder = input("Provide the path to the folder containing the subfolders for each tile: ")
            buildings_path = input("Provide the path to the geopackage containing the buildings: ")
            output_base_folder = input("Provide the output base folder path: ")
            process_dsm_chm_parallel(chm_folder, dtm_path, dsm_path, buildings_path, output_base_folder)

        else:
            # Check if user has downloaded LAZ files and the AHN DSM and DTM
            laz_downloaded = input("Have you already downloaded the LAZ files and do you have a merged DSM and DTM"
                                   "from AHN? (yes/no): ").strip().lower()
            if laz_downloaded == 'yes':
                laz_folder = input("Provide the path to the folder containing folders for each LAZ subtile file: ")
                dsm_path = input("Provide the path to the merged DSM file: ")
                dtm_path = input("Provide the path to the merged DTM file: ")
                buildings_path = input("Provide the path to the geopackage containing the buildings: ")

                output_base_folder = input("Provide the output base folder path: ")
                # Set default values for optional parameters
                params = {
                    'ndvi_threshold': 0.0,
                    'resolution': 0.5,
                    'remove': False,
                    'smooth_chm': False,
                    'filter_size': 3,
                    'max_workers': 4
                }

                optional_input = input(
                    "Do you want to input optional parameters for CHM processing? (yes/no): ").strip().lower()
                if optional_input == 'yes':
                    input_params = input(
                        "Enter optional parameters in the format: ndvi_threshold={float}, resolution={float}, "
                        "remove={bool}, smooth_chm={bool}, filter_size={int}, max_workers={int}: ")
                    params.update(optional_params(input_params))

                process_laz_files_parallel(
                    laz_folder,
                    output_base_folder,
                    ndvi_threshold=params['ndvi_threshold'],
                    resolution=params['resolution'],
                    remove=params['remove'],
                    smooth_chm=params['smooth_chm'],
                    filter_size=params['filter_size'],
                    max_workers=params['max_workers']
                )
                process_dsm_chm_parallel(output_base_folder, dtm_path, dsm_path, buildings_path, output_base_folder)


            else:
                # Ask about missing files
                missing_dsm_dtm = input("Is the DSM and DTM missing? (yes/no): ").strip().lower()
                missing_chm = input("Is the CHM converted from LAZ missing? (yes/no): ").strip().lower()
                output_folder = input("Provide the output folder path: ")

                # Initialize a list to hold download functions and their arguments
                download_tasks = []

                if missing_dsm_dtm == 'yes':
                    # Add DSM/DTM download task
                    dsm_dtm_tiles_file = input(
                        "Provide the path to the text file containing the needed DSM/DTM tiles: ")
                    name = input("Provide the common name for the merged TIF files: ")
                    download_tasks.append((download_raster_tiles, (dsm_dtm_tiles_file, output_folder, name)))

                if missing_chm == 'yes':
                    # Add CHM download task
                    chm_subtiles_file = input("Provide the path to the text file containing the needed CHM subtiles: ")
                    download_tasks.append((download_las_tiles, (chm_subtiles_file, output_folder)))

                    # Set default values for optional parameters
                    params = {
                        'ndvi_threshold': 0.0,
                        'resolution': 0.5,
                        'remove': False,
                        'smooth_chm': False,
                        'filter_size': 3,
                        'max_workers': 4
                    }

                    # Ask if the user wants to input optional parameters
                    optional_input = input(
                        "Do you want to input optional parameters for CHM processing? (yes/no): ").strip().lower()

                    if optional_input == 'yes':
                        input_params = input(
                            "Enter optional parameters in the format: ndvi_threshold={float}, resolution={float}, "
                            "remove={bool}, smooth_chm={bool}, filter_size={int}, max_workers={int}: ")

                        params.update(optional_params(input_params))

                    download_tasks.append((
                        process_laz_files_parallel,
                        (
                            output_folder,
                            output_folder,
                            params['ndvi_threshold'],
                            params['resolution'],
                            params['remove'],
                            params['smooth_chm'],
                            params['filter_size'],
                            params['max_workers']
                        )
                    ))

                # Run all download tasks
                for download_func, download_args in download_tasks:
                    download_func(*download_args)


if __name__ == '__main__':
    main()