from include.road_process import Road
from include.cool_space import CoolSpace
from include.building import Building
from include.cooleval import CoolEval
from datetime import datetime, timedelta
from typing import Union
import geopandas as gpd
import numpy as np
import rasterio
import matplotlib.pyplot as plt
import glob
import os


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
    print(start_time_dt)
    print(search_start_time_dt)
    print(search_end_time_dt)

    # Calculate the time interval as timedelta
    interval = timedelta(minutes=time_interval)

    # Compute the index range for shadows based on the input time range
    start_idx = int((search_start_time_dt - start_time_dt) / interval)
    end_idx = int((search_end_time_dt - start_time_dt) / interval)

    return [start_idx, end_idx]


def identification(landuse_file: str,
                   road_file: str,
                   building_file: str,
                   shadowmaps_file: str,
                   road_buffer_attri: str,
                   building_buffer_num: int = 4,
                   mode: str = 'single-day',
                   num_days: int = None,
                   start_time: int = 900,
                   end_time: int = 1800,
                   time_interval: int = 30,
                   search_start_time: int = 900,
                   search_end_time: int = 1800,
                   output_coolspace_type: str = 'land-use') -> gpd.geodataframe:

    # Read datasets
    road = Road(gpd.read_file(road_file))
    building = Building(gpd.read_file(building_file))
    coolSpace = CoolSpace(gpd.read_file(landuse_file))
    shadows = []
    shadow_files = glob.glob(os.path.join(shadowmaps_file, '*.TIF'))
    for shadow_file in shadow_files:
        shadow_map = rasterio.open(shadow_file, crs=coolSpace.data.crs)
        shadows.append(shadow_map)
    shadows = shadows[0:3]  # Use two for testing, this line should not be included in the final program

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
        search_range = compute_search_range(search_start_time, search_end_time, start_time, end_time, time_interval)
        morning = compute_search_range(900, 1200, start_time, end_time, time_interval)
        afternoon = compute_search_range(1200, 1600, start_time, end_time, time_interval)
        late_aftrn = compute_search_range(1600, 1800, start_time, end_time, time_interval)

        coolSpace.evaluate_shade_coverage(attri_name="Query", start=search_range[0], end=search_range[1])
        coolSpace.evaluate_shade_coverage(attri_name="Morn", start=morning[0], end=morning[1])
        coolSpace.evaluate_shade_coverage(attri_name="Aftrn", start=afternoon[0], end=afternoon[1])
        coolSpace.evaluate_shade_coverage(attri_name="LtAftrn", start=late_aftrn[0], end=late_aftrn[1])
        coolSpace.evaluate_shade_coverage(attri_name="Day", start=0, end=len(shadows) - 1)

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
    if output_coolspace_type == 'land-use':
        output_gdf = coolSpace.get_cool_spaces(geom_type='geometry')
    elif output_coolspace_type == 'public-space':
        output_gdf = coolSpace.get_cool_spaces(geom_type='clipped')
    else:
        print(f"Invalid type {output_coolspace_type}. Expected types are 'land-use' or 'public-space', "
              f"to continue the process, 'land-use' will be used for output.")
        output_gdf = coolSpace.get_cool_spaces(geom_type='geometry')

    return output_gdf


def output(gdf: gpd.geodataframe, output_file_path: str) -> None:
    for col in gdf.columns:
        if gdf[col].apply(lambda x: isinstance(x, list)).any():
            gdf[col] = gdf[col].apply(lambda x: str(x) if isinstance(x, list) else x)
    gdf.to_file(output_file_path)


def evaluation(coolspace_file: gpd.geodataframe,
               building_population_file: str,
               bench_file: str,
               heatrisk_file: str,
               pet_file: str,
               buffer_house: int = 700) -> gpd.geodataframe:

    coolspace = coolspace_file
    building_population = gpd.read_file(building_population_file)
    bench = gpd.read_file(bench_file)
    heatrisk = gpd.read_file(heatrisk_file)

    cool_eval = CoolEval(cool_places=coolspace,
                         buildings=building_population,
                         bench=bench,
                         heatrisk=heatrisk,
                         pet=pet_file,
                         buffer_house=buffer_house)

    cool_eval.evaluate_capacity()
    cool_eval.evaluate_sfurniture()
    cool_eval.evaluate_heatrisk()
    cool_eval.eval_pet()

    return cool_eval


# entry
if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"

    landuse_file = directory_win + "ams_landuse_top10NL.shp"
    road_file = directory_win + "ams_roads_top10NL.shp"
    building_file = directory_win + "ams_buildings_bagplus.shp"
    shadow_path = directory_win + "shademaps\\"

    coolspace = identification(landuse_file=landuse_file,
                               road_file=road_file,
                               building_file=building_file,
                               shadowmaps_file=shadow_path,
                               road_buffer_attri='buffer',
                               building_buffer_num=4,
                               mode='single-day',
                               output_coolspace_type='land-use')

    coolspace_output = evaluation(coolspace_file=coolspace,
                                  building_population_file="",
                                  bench_file="",
                                  heatrisk_file="",
                                  pet_file="",
                                  buffer_house=700)

    output(coolspace_output, directory_win + "test.shp")

