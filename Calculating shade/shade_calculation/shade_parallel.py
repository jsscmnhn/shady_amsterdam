import tqdm
import os
import datetime as dt
import shade_setup as shade

import os
import concurrent.futures
import tqdm

def process_file_pair(chm_filename, dsm_filename, folder_path, output_dir, date):
    """Function to process a single DSM and CHM file pair."""
    chm_path = os.path.join(folder_path, chm_filename)
    dsm_path = os.path.join(folder_path, dsm_filename)
    title = "_".join(chm_filename.split('_')[1:3]).replace('.tif', '')

    # Call the shade calculation setup function
    shade.shadecalculation_setup(
        filepath_dsm=dsm_path,
        filepath_veg=chm_path,
        tile_no=title,
        date=date,
        intervalTime=30,
        onetime=0,
        filepath_save=output_dir,
        UTC=2,
        dst=1,
        useveg=1,
        trunkheight=25,
        transmissivity=10,
        start_time=10,
        end_time=21
        )

    print(f"Completed shade calculation for: {chm_filename} and {dsm_filename}.")


def run_shade_calculation(folder_path, output_base_folder, date, max_workers=4):
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
                executor.submit(process_file_pair, chm_filename, dsm_filename, folder_path, output_dir, date))

        for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures),
                                desc=f"Processing files in {folder_path}"):
            try:
                future.result()
            except Exception as e:
                print(f"Error occurred while processing: {e}")


def process_folders_sequentially(base_folder, output_base_folder, date, max_workers=4):
    # Get a list of all folders in the base folder
    folder_paths = [os.path.join(base_folder, f) for f in os.listdir(base_folder) if
                    os.path.isdir(os.path.join(base_folder, f))]

    for folder_path in tqdm.tqdm(folder_paths, desc="Processing folders"):
        print(f"Processing folder: {folder_path}")
        run_shade_calculation(folder_path, output_base_folder, date, max_workers)


base_folder = "C:/Users/Admin/Documents/Jessica/test"
output_base_folder = "C:/Users/Admin/Documents/Jessica/output"
date = dt.datetime(2015, 7, 1)

process_folders_sequentially(base_folder, output_base_folder, date, max_workers=24)
