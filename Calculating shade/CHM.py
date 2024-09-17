import cProfile
import io
import pstats

import numpy as np
from matplotlib import pyplot as plt

import functions

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

def chm_creation(LasData, veg_points, resolution=0.5, search_radius=1.0, idw_power=2.0):

    veg_points = extract_vegetation_points(LasData)
    vegetation_data = functions.vegetation_raster_idw(LasData, veg_points, resolution,
                                             search_radius, idw_power)
    veg_raster = vegetation_data[0]
    grid_centers = vegetation_data[1]
    top_left_x = grid_centers[0][0, 0] - resolution / 2
    top_left_y = grid_centers[1][0, 0] + resolution / 2

    transform = functions.create_affine_transform(top_left_x, top_left_y, resolution)
    functions.write_output(veg_points, veg_raster, transform, "../data/chm.tiff")