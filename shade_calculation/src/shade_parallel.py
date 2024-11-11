""" This file contains the functions that can take a DSM and CHM of an AHN subtile, and calculate shade maps for
them. The process can run in parallel."""
from turtledemo.penrose import start

from shade_calculation.extra import shade_setup as shade
import datetime
import os
import concurrent.futures
import tqdm

def process_chm_dsm(chm_filename, dsm_filename, folder_path, output_dir, date, start_time=9, end_time=20,
                    interval=30, use_chm=True, trans=10, trunkheight=25):
    """
     Function to process a single DSM and CHM file pair.
     ----
     Input:
     - chm_filename (str): Name of the CHM file to be processed.
     - dsm_filename (str): Name of the DSM file to be processed.
     - folder_path (str): Path to the folder containing the CHM and DSM files.
     - output_dir (str): Directory where the output results will be saved.
     - date (str): Date parameter used in the shade calculation.
     - start_time (int): Starting hour for shade calculation (default: 9).
    - end_time (int): Ending hour for shade calculation (default: 20).
    - interval (int): Time interval for the calculations (minutes) (default: 30).
    - use_chm (bool): Use the CHM file in the shade calculation (default: True).
    - trans (int): Transmissivity value for the calculations, percentage for tree shade (default: 10).
    - trunkheight (int): Trunk height for vegetation calculations, percentage of CHM height representing trunk
                         height (default: 25).

     Output:
     - None: This function prints a message upon completion of shade calculation.
     """

    # update start time and endtime for UMEP tool
    start_time += 1
    end_time += 1

    chm_path = os.path.join(folder_path, chm_filename)
    dsm_path = os.path.join(folder_path, dsm_filename)
    print(chm_path)
    print(dsm_path)
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
        useveg=use_chm,
        trunkheight=trunkheight,
        transmissivity=trans,
        start_time=start_time,
        end_time=end_time
        )

    print(f"Completed shade calculation for: {chm_filename} and {dsm_filename}.")

def run_shade_calculation(file_pairs, output_base_folder, date, start_time, end_time, interval,
                          use_chm=True, trans=10, trunkheight=25, max_workers=4):
    """
     Process CHM and DSM file pairs in parallel.
      ----
      Input:
     - file_pairs (list): List of tuples containing (CHM file, DSM file, folder_path, tile).
      - output_base_folder (str): Base path for saving the output files.
      - date (str): Date parameter used in the shade calculation.
      - max_workers (int): Maximum number of worker threads to use for concurrent processing (default: 4).
      - start_time (int): Starting hour for shade calculation.
      - end_time (int): Ending hour for shade calculation .
      - interval (int): Time interval for the calculations (minutes).
      - trans (int): Transmissivity value for the calculations, percentage for tree shade (default: 10).
      - trunkheight (int): Trunk height for vegetation calculations, percentage of CHM height representing trunk
                         height (default: 25).

      Output:
      - None: This function prints messages regarding the processing status and creates output directories.
      """

    if not file_pairs:
        print("No file pairs.")
        return

    # Create output directories for each tile
    for _, _, folder_path, tile in file_pairs:
        output_dir = os.path.join(output_base_folder, f'{tile}/')

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    # Process files
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        list(tqdm.tqdm(executor.map(
            process_chm_dsm,
            [chm for chm, _, _, _ in file_pairs],
            [dsm for _, dsm, _, _ in file_pairs],
            [folder_path for _, _, folder_path, _ in file_pairs],
            [os.path.join(output_base_folder, f'{tile}/') for _, _, _, tile in file_pairs],
            [date] * len(file_pairs),
            [start_time] * len(file_pairs),
            [end_time] * len(file_pairs),
            [interval] * len(file_pairs),
            [use_chm] * len(file_pairs),
            [trans] * len(file_pairs),
            [trunkheight] * len(file_pairs)
        ), total=len(file_pairs), desc="Processing files", unit="file"))


def process_folders(base_folder, output_base_folder, date, start_time=9, end_time=20, interval=30, use_chm=True, trans=10,
                    trunkheight=25, max_workers=4):
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
    file_pairs = []

    # Collect CHM & DSM file pairs from folders
    for folder_path, _, files in os.walk(base_folder):
        chm_files = [f for f in files if f.startswith('CHM') and (f.endswith('.tif') or f.endswith('.tiff'))]
        dsm_files = [f for f in files if f.startswith('DSM') and (f.endswith('.tif') or f.endswith('.tiff'))]

        # Sort to pair by file order
        chm_files.sort()
        dsm_files.sort()

        if len(chm_files) != len(dsm_files):
            print(f"Warning: Mismatched CHM and DSM file counts in folder {folder_path}")
            continue

        # Add each pair of CHM and DSM files along with the tile name and folder path to the file_pairs list
        for chm, dsm in zip(chm_files, dsm_files):
            # Extract utile name from the CHM filename
            tile = "_".join(chm.split('_')[1:2]).replace('.TIF', '').replace('.tif', '')

            file_pairs.append((chm, dsm, folder_path, tile))

    if file_pairs:
        print(f"Found {len(file_pairs)} CHM/DSM file pairs to process.")
        run_shade_calculation(file_pairs, output_base_folder, date, start_time, end_time, interval,
                                       use_chm, trans, trunkheight, max_workers)
    else:
        print("No CHM/DSM file pairs found.")
