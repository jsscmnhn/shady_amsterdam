"""
Contains relevant functions for completing tasks or visualizations.

----

- crop_laz:

  Necessary function for crop the LAZ/LAS file to a certain bounding box.

- plot_pointcloud:

  Used for the 3D visualization of point cloud.

- plot_grid:

  Used for the visualization of grid.

- create_affine_transform:

  To create an Affine object which is a transformer necessary for the write_output function.

- write_output:

  For writing the grid to an output .tiff file.
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

    Parameters:

    - resolution: The length of each cell.

    Return:

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


def interpolation_idw_kdtree(points_list, grid, grid_center_xy, search_radius=1, power=2):
    """
    Perform points to raster manipulation using IDW interpolation where the default power is 2.

    *This is the improved version by using KD-Tree method to search points, which takes about
    7 seconds to perform a 500m * 500m raster (resolution=0.5m).*

    ----

    Parameters:

    - points: a numpy 2D array (n,3), it has n rows (n points) and each row contains each point's xyz values.
    - grid: a numpy grid for interpolation, which will be manipulated directly.
    - grid_center_xy:

      a tuple containing two numpy grid (grid_center_x, grid_center_y), both of them has the
      same shape as the grid and contain the xy coordinates of the center points of grid cells.
      (e.g.  the center point of grid[i][j] --> grid_center_xy[0][i][j], grid_center_xy[1][i][j])

    - search_radius:

      the search radius of each cell. Default is 1m. Because the IDW interpolation is used, even though there
      may be some points that shouldn't be included, they shouldn't have significant effect to the assigned z
      value because their weight is low.

    - power: the power for IDW interpolation equation.

    Return:

    - grid: the output is the manipulated/changed input grid.
    """
    rows, cols = grid.shape
    grid_center_x, grid_center_y = grid_center_xy
    points_x = points_list[:, 0]
    points_y = points_list[:, 1]
    points_z = points_list[:, 2]
    tree = cKDTree(points_list[:, :2])  # Use xy coordinates to build KD-Tree.
    for i in range(rows):
        for j in range(cols):
            cell_center = [grid_center_x[i, j], grid_center_y[i, j]]
            # Use query_ball_point method to search points that are inside the range. The
            # idx is a list contains row numbers, indicating the corresponding points.
            idx = tree.query_ball_point(cell_center, search_radius)

            if idx:
                epsilon = 1e-10  # create a very small number to prevent divison by zero.
                distances = np.sqrt((points_x[idx] - cell_center[0]) ** 2 + (points_y[idx] - cell_center[1]) ** 2)
                weights = 1 / (distances ** power + epsilon)  # Add the epsilon here.
                if np.all(weights > 0):  # check there is no zero inside the weights list.
                    grid[i, j] = np.sum(weights * points_z[idx]) / np.sum(weights)
                else:
                    raise ValueError("weights == 0, cannot perform division.")
    return grid


def create_affine_transform(top_left_x, top_left_y, res):
    """
    Create Affine transform for the write_output function.
    """
    transform = Affine.translation(top_left_x, top_left_y) * Affine.scale(res, -res)
    return transform

def write_output(dataset, output, transform, name, nodata_value=False):
    """
    Write grid to .tiff file.

    ----

    Parameters:

    - dataset: Can be either a rasterio dataset (for rasters) or laspy dataset (for point clouds)
    - output: the output grid, a numpy grid.
    - name: the name of the output file.
    - transform:
      a user defined rasterio Affine object, used for the transforming the pixel coordinates
      to spatial coordinates.
    """
    output_file = name

    # Check if dataset is laspy pc
    if hasattr(dataset, 'header') and hasattr(dataset.header, 'parse_crs'):
        crs = dataset.header.parse_crs()
    else:
        crs = dataset.crs

    output = np.squeeze(output)

    if nodata_value == True:
        nodata_value = -9999
    else: nodata_value = dataset.nodata

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

