from rich.progress import Progress
import geopandas as gpd
import time
from include.coolspace_process import (read_config,
                                       identification,
                                       evaluation,
                                       list_to_string,
                                       drop_or_wkt,
                                       output_all_shade_geoms)


# entry
if __name__ == '__main__':

# ======== Read config =================================================================================================
    config = read_config("coolspaceConfig.json")
    with Progress() as progress:
        task = progress.add_task("Reading config parameters...", total=1)
        # read GeoPackage and layers
        folder_path                 = config['files']['folder_path']
        gpkg_file                   = config['files']['gpkg_file']
        landuse_layer               = config['files']['landuse_file']
        road_layer                  = config['files']['road_file']
        building_layer              = config['files']['building_file']
        building_population_layer   = config['files']['building_population_file']
        street_furniture_layer      = config['files']['street_furniture_file']
        heatrisk_layer              = config['files']['heatrisk_file']
        output_identification_layer = config['files']['output_identification_file']
        output_evaluation_layer     = config['files']['output_evaluation_file']

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
        hasIdentificationOutput = config['parameters']['hasIdentificationOutput']
        outputShadeGeometries   = config['parameters']['outputShadeGeometries']
        performEvaluation       = config['parameters']['performEvaluation']
        road_buffer_attribute   = config['parameters']['road_buffer_attribute']
        output_coolspace_type   = config['parameters']['output_coolspace_type']
        building_buffer         = config['parameters']['building_buffer']
        capacity_search_buffer  = config['parameters']['capacity_search_buffer']
        progress.advance(task)
# ======================================================================================================================


# ======== Identification Process ======================================================================================
    if not hasIdentificationOutput:
        with Progress() as progress:
            task = progress.add_task("Loading GeoPackage layers for identification...", total=3)

            # Load each layer and update the progress
            landuse = gpd.read_file(gpkg_file, layer=landuse_layer)
            progress.advance(task)

            road = gpd.read_file(gpkg_file, layer=road_layer)
            progress.advance(task)

            building = gpd.read_file(gpkg_file, layer=building_layer)
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

        with Progress() as progress:
            task = progress.add_task("Processing attributes for output data...", total=3)
            list_to_string(coolspace)
            progress.advance(task)
            drop_or_wkt(coolspace, mode='to_wkt')
            progress.advance(task)
            coolspace.to_file(gpkg_file, layer=output_identification_layer)
            progress.advance(task)

        end = time.time()
        total = end - begin
        minutes, seconds = divmod(total, 60)
        print(f"Total time for identification: {int(minutes)} minutes and {seconds:.2f} seconds")
# ======================================================================================================================


# ======== Output Shade Geometries if needed ===========================================================================
    if outputShadeGeometries:
        progress.add_task("Outputing shade geometries...", total=1)
        coolspace_output = gpd.read_file(gpkg_file, layer=output_identification_layer)
        output_all_shade_geoms(coolspace_output, folder_path)
# ======================================================================================================================


# ======== Evaluation Process ==========================================================================================
    if performEvaluation:
        with Progress() as progress:
            task = progress.add_task("Loading GeoPackage layers for evaluation...", total=5)

            # read PET raster
            pet_file = config['files']['pet_file']
            progress.advance(task)

            cs_output = gpd.read_file(gpkg_file, layer=output_identification_layer)
            progress.advance(task)

            building_pop = gpd.read_file(gpkg_file, layer=building_population_layer)
            progress.advance(task)

            street_furniture = gpd.read_file(gpkg_file, layer=street_furniture_layer)
            progress.advance(task)

            heatrisk = gpd.read_file(gpkg_file, layer=heatrisk_layer)
            progress.advance(task)

        begin2 = time.time()
        coolspace_output = evaluation(coolspace=cs_output,
                                      building_population_file=building_pop,
                                      bench_file=street_furniture,
                                      heatrisk_file=heatrisk,
                                      pet_file=pet_file,
                                      gpkg_file=gpkg_file,
                                      output_layer=output_evaluation_layer,
                                      search_buffer=700,
                                      single_day_time_range=single_day_time_range,
                                      time_interval=time_interval,
                                      morning_range=morning_range,
                                      afternoon_range=afternoon_range,
                                      late_afternoon_range=late_afternoon_range,
                                      search_range=search_range)

        end2 = time.time()
        total2 = end2 - begin2
        minutes2, seconds2 = divmod(total2, 60)
        print(f"Total time for evaluation: {int(minutes2)} minutes and {seconds2:.2f} seconds")
# ======================================================================================================================
