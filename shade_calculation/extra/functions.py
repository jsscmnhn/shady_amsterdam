"""
Contains relevant functions for completing tasks or visualizations.

----
Script containing some functions necessary for creating the rasters.
- raster_center_coords:
  To compute the center xy coordinates of a grid

- create_affine_transform:
  To create an Affine object which is a transformer necessary for the write_output function.

- write_output:
  For writing the grid to an output .tiff file, original data can be either from a LAZ file or from a TIF.
"""

import laspy
import numpy as np
import rasterio
from matplotlib import pyplot as plt
from rasterio.transform import Affine
from scipy.spatial import cKDTree


def raster_center_coords(min_x, max_x, min_y, max_y, resolution):
    """
    Compute the center xy coordinates of a grid.
    ----
    Input:
    - min_x, max_x, min_y, max_y(float): Minimum and maximum x and y coordinates of the grid.
    - resolution (float): The length of each cell, function can only be used for square cells.

    Output:
    - grid_center_x: a grid where each cell contains the value of its center point's x coordinates.
    - grid_center_y: a grid where each cell contains the value of its center point's y coordinates.
    """
    # create coordinates for the x and y border of every cell.
    x_coords = np.arange(min_x, max_x, resolution)  # x coordinates expand from left to right.
    y_coords = np.arange(max_y, min_y, -resolution)  # y coordinates reduce from top to bottom.

    # create center point coordinates for evey cell.
    grid_x, grid_y = np.meshgrid(x_coords, y_coords)
    grid_center_x = grid_x + resolution / 2
    grid_center_y = grid_y - resolution / 2
    return grid_center_x, grid_center_y

def create_affine_transform(top_left_x, top_left_y, res):
    """
    Create Affine transform for the write_output function.
    """
    transform = Affine.translation(top_left_x, top_left_y) * Affine.scale(res, -res)
    return transform

def write_output(dataset, output, transform, name, change_nodata=False):
    """
    Write grid to .tiff file.
    ----
    Input:
    - dataset: Can be either a rasterio dataset (for rasters) or laspy dataset (for point clouds)
    - output (Array): the output grid, a numpy grid.
    - name (String): the name of the output file.
    - transform:
      a user defined rasterio Affine object, used for the transforming the pixel coordinates
      to spatial coordinates.
    - change_nodata (Boolean): true: use a no data value of -9999, false: use the datasets no data value
    """
    output_file = name

    #  Determine the CRS: If it's a laspy dataset with a header, use parse_crs; otherwise, use the dataset's crs.
    if hasattr(dataset, 'header') and hasattr(dataset.header, 'parse_crs'):
        crs = dataset.header.parse_crs()
    else:
        crs = dataset.crs

    output = np.squeeze(output)

    # Set the nodata value: use -9999 if nodata_value is True or dataset does not have nodata.
    if change_nodata:
        nodata_value = -9999
    else:
        try:
            nodata_value = dataset.nodata
            if nodata_value is None:
                raise AttributeError("No no data value found in dataset.")
        except AttributeError as e:
            print(f"Warning: {e}. Defaulting to -9999.")
            nodata_value = -9999

    # output the dataset
    with rasterio.open(output_file, 'w',
                       driver='GTiff',
                       height=output.shape[0],  # Assuming output is (rows, cols)
                       width=output.shape[1],
                       count=1,
                       dtype=output.dtype,
                       crs=crs,
                       nodata=nodata_value,
                       transform=transform) as dst:
        dst.write(output, 1)
    print("File written to '%s'" % output_file)

