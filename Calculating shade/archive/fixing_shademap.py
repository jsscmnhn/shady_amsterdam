import os
import numpy as np
import rasterio
from tqdm import tqdm

# Define the path to the folder containing the TIFF files

folder_path = 'G:\Geomatics\merged_correct_value'

# Get a list of all TIFF files in the folder
tif_files = [f for f in os.listdir(folder_path) if f.endswith('.TIF')]

# Iterate through the TIFF files with a progress bar
for filename in tqdm(tif_files, desc="Processing TIFF files"):
    file_path = os.path.join(folder_path, filename)

    # Open the TIFF file with rasterio
    with rasterio.open(file_path, 'r+') as src:
        # Read the data
        data = src.read(1)  # Read the first band

        # Get the no-data value
        no_data_value = src.nodata

        # Perform the transformation: 1 - value for non-no-data values
        # If no_data_value is set, we exclude it from the operation
        if no_data_value is not None:
            data = np.where(data != no_data_value, 1 - data, no_data_value)
        else:
            # If no_data_value is not set, just perform the operation on all values
            data = 1 - data

        # Write the modified data back to the same file
        src.write(data, 1)

    # Optional: Print progress for each file
    tqdm.write(f"Processed {filename}")

print("All files processed.")

