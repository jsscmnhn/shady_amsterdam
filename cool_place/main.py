from include.road_process import Road
from include.identification import CoolSpace
from include.building import Building
from include.evaluation import CoolEval
from datetime import datetime, timedelta
from shapely.geometry import base
from rich.progress import Progress
import geopandas as gpd
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
        raise ValueError(f"The search range is: {search_start_time} to {search_end_time}, "
                         f"which exceeds the time range of"
                         f"shade maps: {start_time} to {end_time}")
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
    with Progress() as progress:
        task = progress.add_task("Loading shade maps...", total=daytime[1] + 1)
        for shadow_file in shadow_files:
            shadow_map = rasterio.open(shadow_file, crs=coolSpace.data.crs)
            shadows.append(shadow_map)
            progress.advance(task)
    shadows = shadows[daytime[0]:daytime[1]]

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
    coolSpace.clip(building.data, use_clip=True)

    # Perform shade calculation for ALL shade maps
    coolSpace.calculate_shade(shadows, use_clip=True)

    # Evaluate the shade coverage based on different modes
    if mode == 'single-day':
        coolSpace.evaluate_shade_coverage(attri_name="Day", start=daytime[0], end=daytime[1])

        if search_range is not None:
            search_start_time = search_range[0]
            search_end_time = search_range[1]
            search = compute_search_range(search_start_time, search_end_time, start_time, end_time, time_interval)
            coolSpace.evaluate_shade_coverage(attri_name="Query", start=search[0], end=search[1])

        if morning_range is not None:
            m_start = morning_range[0]
            m_end = morning_range[1]
            morning = compute_search_range(m_start, m_end, start_time, end_time, time_interval)
            coolSpace.evaluate_shade_coverage(attri_name="Morn", start=morning[0], end=morning[1])

        if afternoon_range is not None:
            a_start = afternoon_range[0]
            a_end = afternoon_range[1]
            afternoon = compute_search_range(a_start, a_end, start_time, end_time, time_interval)
            coolSpace.evaluate_shade_coverage(attri_name="Aftrn", start=afternoon[0], end=afternoon[1])

        if late_afternoon_range is not None:
            la_start = late_afternoon_range[0]
            la_end = late_afternoon_range[1]
            late_aftrn = compute_search_range(la_start, la_end, start_time, end_time, time_interval)
            coolSpace.evaluate_shade_coverage(attri_name="LtAftrn", start=late_aftrn[0], end=late_aftrn[1])

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


def evaluation(coolspace_file: gpd.geodataframe,
               building_population_file: str,
               bench_file: str,
               heatrisk_file: str,
               pet_file: str,
               search_buffer: int = 700) -> gpd.geodataframe:

    coolspace = coolspace_file
    building_population = gpd.read_file(building_population_file)
    bench = gpd.read_file(bench_file)
    heatrisk = gpd.read_file(heatrisk_file)

    cool_eval = CoolEval(cool_places=coolspace,
                         buildings=building_population,
                         bench=bench,
                         heatrisk=heatrisk,
                         pet=pet_file,
                         search_buffer=search_buffer)

    cool_eval.evaluate_capacity()
    cool_eval.evaluate_sfurniture()
    cool_eval.evaluate_heatrisk()
    cool_eval.eval_pet()

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
            else:
                gdf.drop(columns=col, inplace=True)


# entry
if __name__ == '__main__':

    config = read_config("coolspaceConfig.json")

    # Read config
    with Progress() as progress:
        task = progress.add_task("Reading config parameters...", total=1)
        # read GeoPackage and layers
        gpkg_file                 = config['files']['gpkg_file']
        landuse_layer             = config['files']['landuse_file']
        road_layer                = config['files']['road_file']
        building_layer            = config['files']['building_file']
        building_population_layer = config['files']['building_population_file']
        street_furniture_layer    = config['files']['street_furniture_file']
        heatrisk_layer            = config['files']['heatrisk_file']
        output_layer              = config['files']['output_file']

        # read PET raster
        pet_file = config['files']['pet_file']

        # read shape maps file path and information
        shademaps_path         = config['files']['shademaps_path']
        shade_calculation_mode = config['parameters']['shade_calculation_mode']
        time_interval          = config['shade_info_single']['time_interval']
        single_day_time_range  = config['shade_info_single']['single_day_time_range']
        morning_range          = config['shade_info_single']['morning_range']
        afternoon_range        = config['shade_info_single']['afternoon_range']
        late_afternoon_range   = config['shade_info_single']['late_afternoon_range']
        search_range           = config['shade_info_single']['search_range']
        num_days               = config['shade_info_multi']['num_days']

        # read parameters
        road_buffer_attribute  = config['parameters']['road_buffer_attribute']
        output_coolspace_type  = config['parameters']['output_coolspace_type']
        building_buffer        = config['parameters']['building_buffer']
        capacity_search_buffer = config['parameters']['capacity_search_buffer']
        progress.advance(task)

    with Progress() as progress:
        task = progress.add_task("Loading GeoPackage layers...", total=6)

        # Load each layer and update the progress
        landuse = gpd.read_file(gpkg_file, layer=landuse_layer)
        progress.advance(task)

        road = gpd.read_file(gpkg_file, layer=road_layer)
        progress.advance(task)

        building = gpd.read_file(gpkg_file, layer=building_layer)
        progress.advance(task)

        building_pop = gpd.read_file(gpkg_file, layer=building_population_layer)
        progress.advance(task)

        street_furniture = gpd.read_file(gpkg_file, layer=street_furniture_layer)
        progress.advance(task)

        heatrisk = gpd.read_file(gpkg_file, layer=heatrisk_layer)
        progress.advance(task)

    begin = time.time()
    coolspace = identification(coolspace_file=landuse,
                               road_file=road,
                               building_file=building,
                               shademaps_path=shademaps_path,
                               road_buffer_attri=road_buffer_attribute,
                               building_buffer_num=building_buffer,
                               mode=shade_calculation_mode,
                               single_day_time_range=single_day_time_range,
                               time_interval=time_interval,
                               morning_range=morning_range,
                               afternoon_range=afternoon_range,
                               late_afternoon_range=late_afternoon_range,
                               output_coolspace_type=output_coolspace_type)
    end = time.time()
    total = end - begin
    minutes, seconds = divmod(total, 60)
    print(f"Total time for identification: {int(minutes)} minutes and {seconds:.2f} seconds")

    # coolspace_output = evaluation(coolspace_file=coolspace,
    #                               building_population_file="",
    #                               bench_file="",
    #                               heatrisk_file="",
    #                               pet_file="",
    #                               search_buffer=700)

    list_to_string(coolspace)
    drop_or_wkt(coolspace, mode='to_wkt')
    coolspace.to_file(gpkg_file, layer=output_layer)

