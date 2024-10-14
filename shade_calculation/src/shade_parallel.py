""" This file contains the functions that can take a DSM and CHM of an AHN subtile, and calculate shade maps for
them. The process can run in parallel."""

import tqdm
import os
import datetime as dt
import shade_setup as shade

import os
import concurrent.futures
import tqdm

def process_chm_dsm(chm_filename, dsm_filename, folder_path, output_dir, date, start_time=10, end_time=21,
                    interval=30, trans=10, trunkheight=25):
    """
     Function to process a single DSM and CHM file pair.
     ----
     Input:
     - chm_filename (str): Name of the CHM file to be processed.
     - dsm_filename (str): Name of the DSM file to be processed.
     - folder_path (str): Path to the folder containing the CHM and DSM files.
     - output_dir (str): Directory where the output results will be saved.
     - date (str): Date parameter used in the shade calculation.
     - start_time (int): Starting hour for shade calculation (default: 10).
    - end_time (int): Ending hour for shade calculation (default: 21).
    - interval (int): Time interval for the calculations (minutes) (default: 30).
    - trans (int): Transmissivity value for the calculations, percentage for tree shade (default: 10).
    - trunkheight (int): Trunk height for vegetation calculations, percentage of CHM height representing trunk
                         height (default: 25).

     Output:
     - None: This function prints a message upon completion of shade calculation.
     """
    chm_path = os.path.join(folder_path, chm_filename)
    dsm_path = os.path.join(folder_path, dsm_filename)
    title = "_".join(chm_filename.split('_')[1:3]).replace('.tif', '')

    # Call the shade calculation setup function
    shade.shadecalculation_setup(
        filepath_dsm=dsm_path,
        filepath_veg=chm_path,
        tile_no=title,
        date=date,
        intervalTime=interval,
        onetime=0,
        filepath_save=output_dir,
        UTC=2,
        dst=1,
        useveg=1,
        trunkheight=trunkheight,
        transmissivity=trans,
        start_time=start_time,
        end_time=end_time
        )

    print(f"Completed shade calculation for: {chm_filename} and {dsm_filename}.")

def run_shade_calculation(folder_path, output_base_folder, date, start_time, end_time, interval,
                          trans, trunkheight, max_workers=4):
    """
      Function to run shade calculations for pairs of CHM and DSM files in a specified folder.
      ----
      Input:
      - folder_path (str): Path to the folder containing CHM and DSM files.
      - output_base_folder (str): Base path for saving the output files.
      - date (str): Date parameter used in the shade calculation.
      - max_workers (int): Maximum number of worker threads to use for concurrent processing (default: 4).
      - start_time (int): Starting hour for shade calculation (default: 10).
      - end_time (int): Ending hour for shade calculation (default: 21).
      - interval (int): Time interval for the calculations (minutes) (default: 30).
      - trans (int): Transmissivity value for the calculations, percentage for tree shade (default: 10).
      - trunkheight (int): Trunk height for vegetation calculations, percentage of CHM height representing trunk
                         height (default: 25).

      Output:
      - None: This function prints messages regarding the processing status and creates output directories.
      """

    chm_files = [f for f in os.listdir(folder_path) if f.startswith('CHM')]
    dsm_files = [f for f in os.listdir(folder_path) if f.startswith('DSM')]

    if not chm_files or not dsm_files:
        print("No CHM or DSM files found in the specified folder.")
        return

    # Sort files
    chm_files.sort()
    dsm_files.sort()

    if len(chm_files) != len(dsm_files):
        print("Mismatch between the number of CHM and DSM files.")
        return

    first_chm_filename = chm_files[0]
    tile = "_".join(first_chm_filename.split('_')[1:2]).replace('.TIF', '')

    output_dir = os.path.join(output_base_folder, f'{tile}/')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    else:
        print(f"Output directory already exists: {output_dir}")

    # Process files concurrently within the folder
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for chm_filename, dsm_filename in zip(chm_files, dsm_files):
            futures.append(
                executor.submit(process_chm_dsm, chm_filename, dsm_filename, folder_path, output_dir,
                                date, start_time, end_time, interval, trans, trunkheight))

        for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures),
                                desc=f"Processing files in {folder_path}"):
            try:
                future.result()
            except Exception as e:
                print(f"Error occurred while processing: {e}")


def process_folders(base_folder, output_base_folder, date, start_time=10, end_time=21, interval=30, trans=10, trunkheight=25, max_workers=4):
    """
    Function to process all subfolders in a base folder for shade calculations.
    ----
    Input:
    - base_folder (str): Path to the base folder containing subfolders with CHM and DSM files.
    - output_base_folder (str): Base path for saving the output files.
    - date (str): Date parameter used in the shade calculation.
    - max_workers (int): Maximum number of worker threads to use for concurrent processing (default: 4).
    - start_time (int): Starting hour for shade calculation (default: 10).
    - end_time (int): Ending hour for shade calculation (default: 21).
    - interval (int): Time interval for the calculations (minutes) (default: 30).
    - trans (int): Transmissivity value for the calculations, percentage for tree shade (default: 10).
    - trunkheight (int): Trunk height for vegetation calculations, percentage of CHM height representing trunk
                         height (default: 25).

    Output:
    - None: This function prints messages regarding the processing status for each folder.
    """
    # Get a list of all folders in the base folder
    folder_paths = [os.path.join(base_folder, f) for f in os.listdir(base_folder) if
                    os.path.isdir(os.path.join(base_folder, f))]

    for folder_path in tqdm.tqdm(folder_paths, desc="Processing folders"):
        print(f"Processing folder: {folder_path}")
        run_shade_calculation(folder_path, output_base_folder, date, start_time, end_time, interval, max_workers)


# base_folder = "G:/Geomatics/final"
# output_base_folder = "G:/Geomatics/output"
# date = dt.datetime(2015, 7, 1)
#
# process_folders(base_folder, output_base_folder, date, max_workers=50)
