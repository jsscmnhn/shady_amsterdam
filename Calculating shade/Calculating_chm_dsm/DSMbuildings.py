import os
import glob
import re
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.fill import fillnodata
from shapely.geometry import mapping
from rasterio.features import geometry_mask
from functions import write_output
import tqdm

def get_bbox(raster_paths):
    bboxes = []
    for raster_path in raster_paths:
        with rasterio.open(raster_path) as src:
            bbox = src.bounds
            bboxes.append(bbox)
    left = max([bbox.left for bbox in bboxes])
    bottom = max([bbox.bottom for bbox in bboxes])
    right = min([bbox.right for bbox in bboxes])
    top = min([bbox.top for bbox in bboxes])
    return rasterio.coords.BoundingBox(left, bottom, right, top)

def crop_raster(raster_path, bbox, no_data=-9999):
    with rasterio.open(raster_path) as src:
        window = src.window(bbox.left, bbox.bottom, bbox.right, bbox.top)

        # Read the data for the specified window
        cropped_data = src.read(window=window)

        # Read the full mask for the entire dataset
        full_mask = src.read_masks(1)

        # Ensure window attributes are integers
        row_off = int(window.row_off)
        col_off = int(window.col_off)
        height = int(window.height)
        width = int(window.width)

        # Slice the mask to match the window dimensions
        cropped_mask = full_mask[row_off:row_off + height, col_off:col_off + width]

        # Replace no-data values in the data
        if src.nodata is not None:
            cropped_data[cropped_data == src.nodata] = no_data
        return cropped_data, src.window_transform(window), src.crs, cropped_mask


def fill_dtm(cropped_data, transform, crs, src_mask, max_dis, smooth_it, nodata=-9999):
    # Replace any specified no_data values with np.nan for filling
    cropped_data = np.where(cropped_data == nodata, np.nan, cropped_data)

    filled_array = fillnodata(cropped_data, mask=src_mask, max_search_distance=max_dis, smoothing_iterations=smooth_it)

    # Optionally, set filled values back to the nodata value
    # filled_array[np.isnan(filled_array)] = nodata

    return filled_array[0]

def chm_finish(chm_array, dtm_array):
    result_array = chm_array - dtm_array
    result_array[result_array < 0] = 0
    return result_array

def replace_buildings(filled_dtm, dsm_buildings, buildings_geometries, transform, nodata_value=-9999):
    print("filled_dtm shape:", filled_dtm.shape)

    building_mask = geometry_mask(buildings_geometries, transform=transform, invert=True, out_shape=filled_dtm.shape)
    print("Building mask shape:", building_mask.shape)

    # Apply the mask to the filled DTM
    final_dtm = np.where(building_mask, filled_dtm, dsm_buildings)

    return final_dtm

def load_buildings(buildings_path):
    buildings_gdf = gpd.read_file(buildings_path)
    return [mapping(geom) for geom in buildings_gdf.geometry]

def extract_tilename(filename):
    match = re.match(r'CHM_(\w+)_\d+\.TIF', filename)
    if match:
        return match.group(1)
    return None


def process_files(chm_files, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value=-9999,
                  max_search_distance=80.0, smoothing_iterations=2):
    building_geometries = load_buildings(buildings_path)

    # Use tqdm to wrap the file list for progress tracking
    for chm_path in tqdm.tqdm(chm_files, desc="Processing Files", unit="file"):
        chm_filename = os.path.basename(chm_path)

        # Extracting the tile name
        tile = extract_tilename(chm_filename)
        if not tile:
            print(f"Skipping {chm_filename}, couldn't extract common part.")
            continue

        # creating output folder for tile name
        output_folder = os.path.join(output_base_folder, tile)
        os.makedirs(output_folder, exist_ok=True)

        # overlap BBOX dtm, dsm & chm
        raster_paths = [dtm_path, chm_path, dsm_path]
        overlapping_bbox = get_bbox(raster_paths)

        # Cropping rasters to bbox
        dtm_cropped, dtm_transform, dtm_crs, dtm_mask = crop_raster(dtm_path, overlapping_bbox, no_data=nodata_value)
        dsm_cropped, _, _, _ = crop_raster(dsm_path, overlapping_bbox, no_data=nodata_value)
        chm_cropped, _, _, _ = crop_raster(chm_path, overlapping_bbox, no_data=nodata_value)

        # no data values fill
        filled_dtm = fill_dtm(dtm_cropped, dtm_transform, dtm_crs, dtm_mask, max_search_distance, smoothing_iterations,
                              nodata=nodata_value)

        print(filled_dtm.shape)

        # norm CHM calculation
        chm_result = chm_finish(chm_cropped, filled_dtm)

        # Insert buildings from DSM in DTM
        final_dtm_with_buildings = replace_buildings(filled_dtm, dsm_cropped, building_geometries, dtm_transform,
                                                     nodata_value)

        # get number subtile for output
        file_number = re.search(r'_(\d+)\.TIF', chm_filename).group(1)

        # output file names
        output_dtm_filename = f"DTM_{tile}_{file_number}.tif"
        output_chm_filename = f"CHM_{tile}_{file_number}.tif"

        # Saving final DTM + buildings
        output_dtm_path = os.path.join(output_folder, output_dtm_filename)
        write_output(rasterio.open(dtm_path), final_dtm_with_buildings, dtm_transform, output_dtm_path)

        # FOR NOW SAVING TO NEW PATH FOR TESTING -> CHANGE THIS WHEN IT WORKS
        output_chm_path = os.path.join(output_folder, output_chm_filename)
        write_output(rasterio.open(chm_path), chm_result, dtm_transform, output_chm_path)

        print(f"Processed {chm_path} -> {output_dtm_path}")

#%% Run the processing function
geotiff_dtm = "data/DTM_ams.tif"
geotiff_dsm = "data/DSM_ams.tif"
buildings = "data/ams_buildings.gpkg"
chm_folder = "output/25DN2_TEST"
output_base_folder = "final"

chm_files = glob.glob(os.path.join(chm_folder, "*.TIF"))
process_files(chm_files, geotiff_dtm, geotiff_dsm, buildings, output_base_folder)