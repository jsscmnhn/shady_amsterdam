import laspy
import numpy as np


def extract_ground_points(LasData):
    ground_points = LasData[(LasData.classification == 2)]

    return ground_points

def extract_building_points(LasData):
    building_points = LasData[(LasData.classification == 6)]

    return building_points


file_path = 'to_laz_file'
with laspy.open(file_path) as las:
    LasData = las.read()


def save_points_as_las(LasData, output_file, new_points):
    """
    Extract vegetation points and save them as a new LAS file.

    Parameters:
    - LasData (laspy.LasData): Input point cloud data.
    - output_file (str): Path to the output LAS file.
    - ndvi_threshold (float): NDVI threshold for vegetation points.

    Returns:
    - None
    """
    # Extract vegetation points

    # Check if any points were filtered
    if len(new_points) == 0:
        print("No  points found .")
        return

    # Create a new LasData object for the vegetation points
    veg_las_data = laspy.LasData(LasData.header)
    veg_las_data.points = new_points.points  # Assign the filtered vegetation points

    # Save the vegetation points to a new LAS file
    with laspy.open(output_file, mode="w", header=veg_las_data.header) as writer:
        writer.write_points(veg_las_data.points)

    print(f"Vegetation points saved to {output_file}")

