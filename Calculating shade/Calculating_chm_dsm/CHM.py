import functions
import os
import laspy
import time  # Import time module to track time
from tqdm import tqdm
import numpy as np

def extract_vegetation_points(LasData, ndvi_threshold=0.1):
    """
    Extract vegetation points based on classification and NDVI threshold.

    Parameters:
    - LasData (laspy.LasData): Input point cloud data.
    - ndvi_threshold (float): The NDVI threshold to classify vegetation.

    Returns:
    - veg_points (laspy.LasData): Filtered vegetation points as a new LasData object.
    """

    # Filter points based on classification (vegetation-related classes)
    possible_vegetation_points = LasData[(LasData.classification == 1) |  # Unclassified
                                         (LasData.classification == 3) |  # Low vegetation
                                         (LasData.classification == 4) |  # Medium vegetation
                                         (LasData.classification == 5)]  # High vegetation

    # Calculate NDVI
    red = possible_vegetation_points.red
    nir = possible_vegetation_points.nir
    ndvi = (nir.astype(float) - red) / (nir + red)

    # Filter the points whose NDVI is greater than the threshold
    veg_points = possible_vegetation_points[ndvi > ndvi_threshold]

    return veg_points

def chm_creation(LasData, veg_points, output_filename, resolution=0.5, search_radius=1.0, idw_power=2.0):

    veg_points = extract_vegetation_points(LasData)
    vegetation_data = vegetation_raster_idw(LasData, veg_points, resolution,
                                                      search_radius, idw_power)
    veg_raster = vegetation_data[0]
    grid_centers = vegetation_data[1]
    top_left_x = grid_centers[0][0, 0] - resolution / 2
    top_left_y = grid_centers[1][0, 0] + resolution / 2

    transform = functions.create_affine_transform(top_left_x, top_left_y, resolution)

    functions.write_output(veg_points, veg_raster, transform, output_filename)


def vegetation_raster_idw(LasData, points, resolution, search_radius=1, power=2):
    """
    Create a vegetation raster using Inverse Distance Weighting (IDW) interpolation.

    Parameters:
    - LasData (laspy.LasData): Input LAS data.
    - points (laspy.LasData): Vegetation points to be interpolated.
    - center_point (tuple): Center coordinates (x, y) of the region of interest.
    - resolution (float): Resolution of the raster.
    - half_size (float): Half of the size of the tile (defines the extent of the raster).
    - search_radius (float): Radius for the IDW search. Default is 1.
    - power (float): Power parameter for the IDW algorithm. Default is 2.

    Returns:
    - vege_raster (np.ndarray): Generated raster for vegetation.
    - grid_center_xy (tuple): Grid of x, y center coordinates for each raster cell.
    """

    # Extents of the pc
    min_x, max_x = round(LasData.x.min()), round(LasData.x.max())
    min_y, max_y = round(LasData.y.min()), round(LasData.y.max())

    # Define size of the region
    x_length = max_x - min_x
    y_length = max_y - min_y

    # Number of rows and columns
    cols = round(x_length / resolution)
    rows = round(y_length / resolution)

    # Initialize raster grid
    vege_raster = np.full((rows, cols), np.nan, dtype=np.float32)

    # Calculate center coords for each grid cell
    grid_center_xy = functions.raster_center_coords(min_x, max_x, min_y, max_y, resolution)

    if len(points) == 0:
        print("step4, vegetation_raster_idw: There are no vegetation points in the current area.")
        return vege_raster, grid_center_xy

    # Convert LAS points to NumPy array for processing
    points_list = np.array(points.xyz)

    # Interpolation
    functions.interpolation_idw_kdtree(points_list, vege_raster, grid_center_xy, search_radius, power)

    return vege_raster, grid_center_xy

def process_laz_files(input_folder, output_folder, ndvi_threshold=0.1, resolution=0.5):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    laz_files = [f for f in os.listdir(input_folder) if f.endswith(".LAZ")]

    total_start_time = time.time()

    # Iterate over files
    for file_name in tqdm(laz_files, desc="Processing files", unit="file"):
        file_path = os.path.join(input_folder, file_name)
        tile_name = os.path.splitext(file_name)[0]
        output_filename = os.path.join(output_folder, f"CHM_{tile_name}.TIF")

        start_time = time.time()

        # Load LAS data
        with laspy.open(file_path) as las:
            LasData = las.read()

        # Extract vegetation points
        veg_points = extract_vegetation_points(LasData, ndvi_threshold=ndvi_threshold)

        # Create the CHM and save it
        chm_creation(LasData, veg_points, output_filename, resolution=resolution)

        # Calculate time taken for current file
        elapsed_time = time.time() - start_time
        print(f"\nProcessed {file_name} in {elapsed_time:.2f} seconds and saved output to {output_filename}")

        os.remove(file_path)
        print(f"Deleted file: {file_path}")

    # Calculate and print total processing time
    total_elapsed_time = time.time() - total_start_time
    print(f"\nAll files processed in {total_elapsed_time:.2f} seconds.")

input_folder = "data/25DN2"
output_folder = "output/25DN2"
process_laz_files(input_folder, output_folder)