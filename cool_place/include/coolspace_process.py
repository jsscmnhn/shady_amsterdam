from include.road_process import Road
from include.identification import CoolSpace
from include.building import Building
from include.evaluation import CoolEval
from datetime import datetime, timedelta
from shapely.geometry import base
from shapely import wkt
from rich.progress import Progress
import geopandas as gpd
import pandas as pd
import fiona
import rasterio
import glob
import os
import configparser
import ast
import time
import json

def read_config(filename):
    with open(filename, 'r') as f:
        config = json.load(f)
    return config

def convert_to_datetime(time: int) -> datetime:
    hours, minutes = divmod(time, 100)
    return datetime(2024, 10, 16, hours, minutes)  # the [year, month, day] doesn't matter


def compute_search_range(search_start_time: int,
                         search_end_time: int,
                         start_time: int,
                         end_time: int,
                         time_interval: int) -> list:

    if search_start_time < start_time or search_end_time > end_time:
        try:
            raise ValueError(f"The search range is: {search_start_time} to {search_end_time}, "
                             f"which exceeds the time range of"
                             f" shade maps: {start_time} to {end_time}")
        except ValueError as e:
            print(f"Error: {e}, empty time range will be returned.")
            return [None, None]

    # Convert integers to date-time format
    start_time_dt = convert_to_datetime(start_time)
    search_start_time_dt = convert_to_datetime(search_start_time)
    search_end_time_dt = convert_to_datetime(search_end_time)
    # print(start_time_dt)
    # print(search_start_time_dt)
    # print(search_end_time_dt)

    # Calculate the time interval as timedelta
    interval = timedelta(minutes=time_interval)

    # Compute the index range for shadows based on the input time range
    start_idx = int((search_start_time_dt - start_time_dt) / interval)
    end_idx = int((search_end_time_dt - start_time_dt) / interval)

    return [start_idx, end_idx]


def identification(coolspace_file: gpd.geodataframe,
                   road_file: gpd.geodataframe,
                   building_file: gpd.geodataframe,
                   shademaps_path: str,
                   road_buffer_attri: str,
                   building_buffer_num: int = 4,
                   mode: str = 'single-day',
                   num_days: int = None,
                   single_day_time_range: list = None,
                   time_interval: int = None,
                   search_range: list = None,
                   morning_range: list = None,
                   afternoon_range: list = None,
                   late_afternoon_range: list = None,
                   output_coolspace_type: str = 'land-use') -> gpd.geodataframe:

    coolSpace = CoolSpace(coolspace_file)
    road = Road(road_file)
    building = Building(building_file)

    # Read datasets
    shadows = []
    shadow_files = glob.glob(os.path.join(shademaps_path, '*.TIF'))
    start_time = single_day_time_range[0]
    end_time = single_day_time_range[1]
    daytime = compute_search_range(start_time, end_time, start_time, end_time, time_interval)
    print(f"Process range: {daytime}")
    with Progress() as progress:
        task = progress.add_task("Loading shade maps...", total=daytime[1] + 1)
        for shadow_file in shadow_files:
            shadow_map = rasterio.open(shadow_file, crs=coolSpace.data.crs)
            shadows.append(shadow_map)
            progress.advance(task)
    shadows = shadows[daytime[0]:daytime[1]]
    print(f"Total shade maps: {len(shadows)}")

    # Set the calculation mode.
    # For single-day, it allows query for different time-range after the calculation for the whole day
    # For multi-days, it only calculates for each day and doesn't allow ANY query
    if mode not in ['single-day', 'multi-days']:
        raise ValueError(f"Invalid mode: {mode}, expect 'single-day' or 'multi-days'")

    # Create road buffer based on the given buffer attribute name
    road.create_attribute('typeweg', road_buffer_attri)  # This line should be included in final program
    road.create_buffer(road_buffer_attri)
    road.data.set_geometry("buffered", inplace=True)

    # Create building buffer
    building.create_buffer(building_buffer_num)
    building.data.set_geometry("buffered", inplace=True)

    # Clip the input land-use data to create public space data
    coolSpace.clip(road.data, use_clip=False)
    coolSpace.clip(building.data, use_clip=True, filter_thin=True)

    # Perform shade calculation for ALL shade maps
    coolSpace.calculate_shade(shadows, use_clip=True)

    # Evaluate the shade coverage based on different modes
    if mode == 'single-day':
        coolSpace.evaluate_shade_coverage(attri_name="Day", start=daytime[0], end=daytime[1] - 1)

        if search_range is not None:
            search_start_time = search_range[0]
            search_end_time = search_range[1]
            search = compute_search_range(search_start_time, search_end_time, start_time, end_time, time_interval)
            start = search[0]
            end = search[1] - 1 if search[1] is not None else None
            coolSpace.evaluate_shade_coverage(attri_name="Query", start=start, end=end)

        if morning_range is not None:
            m_start = morning_range[0]
            m_end = morning_range[1]
            morning = compute_search_range(m_start, m_end, start_time, end_time, time_interval)
            start = morning[0]
            end = morning[1] - 1 if morning[1] is not None else None
            coolSpace.evaluate_shade_coverage(attri_name="Morn", start=start, end=end)

        if afternoon_range is not None:
            a_start = afternoon_range[0]
            a_end = afternoon_range[1]
            afternoon = compute_search_range(a_start, a_end, start_time, end_time, time_interval)
            start = afternoon[0]
            end = afternoon[1] - 1 if afternoon[1] is not None else None
            coolSpace.evaluate_shade_coverage(attri_name="Aftrn", start=start, end=end)

        if late_afternoon_range is not None:
            la_start = late_afternoon_range[0]
            la_end = late_afternoon_range[1]
            late_aftrn = compute_search_range(la_start, la_end, start_time, end_time, time_interval)
            start = late_aftrn[0]
            end = late_aftrn[1] - 1 if late_aftrn[1] is not None else None
            coolSpace.evaluate_shade_coverage(attri_name="LtAftrn", start=start, end=end)

    elif mode == 'multi-days':
        # Shade maps number for one day, it assumes that every day has the same number of shade maps
        n = len(shadows) / num_days
        for day in range(num_days):
            search_range = [int(day * n), int(day * n + n - 1)]
            coolSpace.evaluate_shade_coverage(attri_name=f"{day}", start=search_range[0], end=search_range[1])
        coolSpace.evaluate_shade_coverage(attri_name='Alldays', start=0, end=len(shadows) - 1)

    # Get cool space polygons of all-daytime search range. Based on the settings, the output
    # polygons will be either land-use polygons or public space polygons. All the other geometries
    # will be transformed into WKT and stored in columns.
    with Progress() as progress:
        output_task = progress.add_task("Processing output data...", total=1)
        if output_coolspace_type == 'land-use':
            output_gdf = coolSpace.get_cool_spaces(geom_type='geometry')
        elif output_coolspace_type == 'public-space':
            output_gdf = coolSpace.get_cool_spaces(geom_type='clipped')
        else:
            print(f"Invalid type {output_coolspace_type}. Expected types are 'land-use' or 'public-space', "
                  f"to continue the process, 'land-use' will be used for output.")
            output_gdf = coolSpace.get_cool_spaces(geom_type='geometry')
        progress.advance(output_task)

    return output_gdf


def evaluation(coolspace: gpd.geodataframe,
               building_population_file: gpd.geodataframe,
               bench_file: gpd.geodataframe,
               heatrisk_file: gpd.geodataframe,
               pet_file: str,
               gpkg_file: str,
               output_layer: str,
               search_buffer: int = 700) -> gpd.geodataframe:

    # coolspace = gpd.read_file(output_file)
    # building_population = gpd.read_file(building_population_file)
    # bench = gpd.read_file(bench_file)
    # heatrisk = gpd.read_file(heatrisk_file)



    # Prompt the user to specify the start and end indices for sdgeom columns
    start_layer = int(input("Enter the starting layer index (e.g., 1 for sdgeom1): "))
    end_layer = int(input("Enter the ending layer index (e.g., 3 for sdgeom3): "))

    # Generate the WKT column names based on user input
    wkt_columns = [f"sdGeom{i}" for i in range(start_layer, end_layer + 1)]
    cool_eval = CoolEval(cool_places=coolspace,
                         buildings=building_population_file,
                         bench=bench_file,
                         heatrisk=heatrisk_file,
                         pet=pet_file,
                         search_buffer=search_buffer)

    # Iterate over each WKT column provided by the user input
    for col in wkt_columns:
        if col in coolspace.columns:
            print(f"Processing {col}...")


            # Function to safely load WKT geometries and log errors
            def safe_load_wkt(x):
                if pd.notnull(x) and isinstance(x, str) and x.strip():
                    try:
                        return wkt.loads(x)
                    except Exception as e:
                        print(f"Error parsing WKT: '{x}' | Error: {e}")
                return None


            # Convert WKT column to geometries, handling invalid or empty geometries
            coolspace[col] = coolspace[col].apply(safe_load_wkt)

            # Remove rows where geometry is None
            valid_shade = coolspace[coolspace[col].notnull()]

            shade = gpd.GeoDataFrame(valid_shade[['id', col]], geometry=col, crs=coolspace.crs)
            shade.rename(columns={col: 'geometry'}, inplace=True)
            shade.set_geometry('geometry', inplace=True)

            # Perform the evaluations on the current shade geometry
            shade = cool_eval.evaluate_capacity(shade, col)
            shade = cool_eval.evaluate_sfurniture(shade, col)
            shade = cool_eval.evaluate_heatrisk(shade, col)
            shade = cool_eval.eval_pet(shade, col)

            # Store the evaluated shade geometries for aggregation
            cool_eval.eval_shades.append(shade)

    # Aggregate results to cool places
    cool_eval.aggregate_to_coolspaces()

    # Export the final results
    cool_eval.export_eval_gpkg(gpkg_file, layer_name=output_layer)

    print("Processing complete!")

    return cool_eval


def list_to_string(gdf: gpd.geodataframe) -> None:
    for col in gdf.columns:
        if gdf[col].apply(lambda x: isinstance(x, list)).any():
            gdf[col] = gdf[col].apply(lambda x: str(x) if isinstance(x, list) else x)


def drop_or_wkt(gdf: gpd.geodataframe, mode='to_wkt') -> None:
    geometry_columns = [col for col in gdf.columns if col != gdf.geometry.name]
    for col in geometry_columns:
        if isinstance(gdf[col].iloc[0], base.BaseGeometry):
            if mode == 'to_wkt':
                gdf[col] = gpd.GeoSeries(gdf[col]).to_wkt()
                if col == 'geometry':
                    gdf.rename(columns={col: 'orig_geom'}, inplace=True)
            else:
                gdf.drop(columns=col, inplace=True)


def output_all_shade_geoms(gdf: gpd.geodataframe, folder_path: str) -> None:
    """
    This method is used for output all shade goemetries to a GeoPackage for routing process.
    It filters out the 'overig' type and 'bebouwd gebied' type which are not public area, thus
    it is a hard-coded method. Be careful when using this method.

    :param gdf: Input GeoDataFrame.
    :param folder_path: Output folder path for the GeoPackage (shadeGeoms.gpkg)
    """
    sd_geom_cols = [col for col in gdf.columns if col.startswith('sdGeom')]
    sd_area_cols = [col for col in gdf.columns if col.startswith('sdArea')]
    sd_avg_cols = [col for col in gdf.columns if col.startswith('sdAvg')]

    with Progress() as progress:
        task = progress.add_task("Outputing shade geometries...", total=len(sd_geom_cols))
        filtered_gdf = gdf[~gdf['typelandge'].isin(['overig', 'bebouwd gebied'])]
        for col1, col2, col3 in zip(sd_geom_cols, sd_area_cols, sd_avg_cols):
            single_layer_df = filtered_gdf[['id', col3, col2, col1]].copy()
            single_layer_df[col1] = single_layer_df[col1].apply(wkt.loads)
            single_layer_gdf = gpd.GeoDataFrame(single_layer_df, geometry=col1, crs=gdf.crs)

            if single_layer_gdf[col1].notna().any():
                output_gpkg = folder_path + "shadeGeoms.gpkg"
                single_layer_gdf.to_file(output_gpkg, layer=col1, driver="GPKG")

            progress.advance(task)