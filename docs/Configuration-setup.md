## Configuration Files
To run the code, three different configuration files have to be given as arguments for the three steps of shade map creation, cool spaces 
identification and scoring and routing. These configuration files are written in JSON.

---
## Content

**[1. Shade Calculation Configuration file](#heading--1)**

**[2. Cool Spaces Configuration file](#heading--2)**

**[3. Network Configuration file](#heading--3)**

---
## Shade Calculation Configuration file <a name="heading--1"/>
If no `config_file_shade` argument is given in *main.py*, the code shall default to .\configuration_files\shade_config.json.
One can choose to either edit the file in the repository directly, or to download the file, edit it and provide it as argument. 

The configuration file is structured as following:

```yaml
{
    "Downloading_Las_tiles": {
        "download_las": true,              # download the LAZ GeoTiles subtiles provided in the .txt file (bool)     
        "las_output_folder": "path",       # path to directory where the subtile directories will be placed containing 
                                           # the LAZ files(str)
                                           # if download_las is false, but you want to perform create_chm, 
                                           # place the path to the directory containing the directories with the LAZ tiles here.
        "subtile_list_file": "path",       # path to the .txt file containing the GeoTiles subtile names to be downloaded (str)
        "notes": {
            "output": "las_output_folder: Provide path to folder containing subtile LAZ folders if LAS tile download is not 
                       needed, but Creation of First CHM file is."
        }
    },
    "Creation_of_AHN_geotiles_DSM_and_DTM_dataset": {
        "download_dsm_dtm": true,          # download the DSM and DTM GeoTiles tiles provided in the .txt file (bool)          
        "download_buildings" : true,       # download the building geometries (bool)
        "buildings_name" : "name",         # name for the GeoPackage containing the building geometries (str)
        "buildings_output_folder": "path", # path to the directory to save the GeoPackage in (str)
        "ahn_output_folder": "path",       # path to the directory to save the merged DSM and DTM in (str)
        "tile_list_file": "path",          # path to the .txt file containing the GeoTiles main tile names to be downloaded (str)
        "merged_name": "name",             # prefix for the merged DSM and DTM (str)
        "optional": {
                                           # if download_dsm_dtm is set to false, but create_final_dsm_chm is true, provide
                                           # the paths to the .TIFF files of the merged DSM and merged DTM.
            "merged_dsm": "path",          
            "merged_dtm": "path"
        },
        "notes": {
            "optional1": "merged_dsm: Provide path to file containing the merged DSM file if creation of AHN geotiles dataset is 
                          not needed.",
            "optional2": "merged_dtm: Provide path to file containing the merged DTM file if creation of AHN geotiles dataset is 
                          not needed."
        }
    },
    "Creation_of_first_CHM_file": {
        "create_chm": true,          # perform the creation of CHM files (bool)
        "chm_output_folder": "path", # path to directory where the subtile directories will be made containing the CHM files (str)
                                     # if create_chm is false, but you want to perform create_final_dsm_chm, 
                                     # place the path to the directory containing the directories with the CHM here.
        "ndvi_threshold": 0.0,       # lowest NDVI allowed to classify points as vegetation (float)
        "resolution": 0.5,           # resolution of the output CHM in meters (float)
        "remove_las": false,         # remove the .LAZ subtile after creating the CHM (bool)
        "smooth_chm": false,         # smooth the output chm with a median filter (bool)
        "filter_size": 3,            # footprint size of the CHM filter (int)
        "pre_filter": false,         # remove all unclassified points up to 1.5 m above the lowest unclassified point, 
                                       helps speeding up processing of subtiles with a large amount of unclassified points (bool)
        "chm_max_workers": 4,        # maximum amount of parallel processes for the creation of CHM (int)
        "notes": {
            "output": "chm_output_folder: Provide path if CHM creation is not needed but want to perform Final CHM & DSM Creation"
        }
    },
    "Creation_of_final_CHM_and_DSM_with_ground_and_buildings": {    # the parameters for the final DSM & CHM creation step
        "create_final_dsm_chm": true,       # perform the creation of final DSM & CHM (bool)
        "buildings_path": "path",           # path to GeoPackage file with the building geometries, 
                                            # ONLY ADD IF NOT CREATED WITH {download_buildings} (str)
        "output_dsm_chm": "path",           # path to directory where the subtile directories will be placed containing the DSM and CHM files (str)
                                            # if create_final_dsm_chm is false, but you want to perform create_shade 
                                            # place the path to the directory containing the directories with the DSM & CHM here.
        "dsm_max_workers": 4,               # maximum amount of parallel processes for the creation of CHM & DSM (int)
        "speed_up": false,                  # speed up the processing of DTM tiles with large no data values by performing 
                                            # linear interpolation (bool) not recommended if parallel processing.
        "min_vegetation_height": 2,         # minimum allowed vegetation height in final CHM, all below will be set to 0 (int)
        "max_vegetation_height": 40,        # maximum allowed vegetation height in final CHM, all above will be set to 0 (int)
        "notes": {
            "output": "output_dsm_chm: Provide path if final CHM and DSM creation is not needed but want to perform Shade Map Creation."
        }
    },
    "Creation_of_Shade_maps": {             # the parameters for the shade map creation step
        "create_shade": true,               # perform the creation of shade maps (bool)
        "output_base_shademap": "path",     # path to directory where the subtile directories will be placed (str)
                                            # if create_shade is false, but you want to perform merging of shade maps, 
                                            # place the path to the directory containing the directories with shade maps here.
        "date": "2020-12-12",               # date to perform shade pattern calculation for, needs to be in format "YYYY-MM-DD"
        "start_time": 9,                    # start time of calculations, note: time in +2 UTC (int)
        "end_time": 20,                     # end time of calculations, note: time in +2 UTC (int)
        "interval": 30,                     # minutes between calculations
        "use_chm": true,                    # use vegetation data in the shade map calculations (bool)
        "trans": 10,                        # transmissivity of light through vegetation in percentage (int)
        "trunkheight": 25,                  # percentage of CHM height that will be seen as trunk height (int)
        "max_shade_workers": 4,             # maximum amount of parallel processes for the creation of shade maps (int)
        "notes": {
            "output": "output_base_shademap: Provide path if shade map creation is not needed but want to perform Merging of Shade Maps."
        }
    },
    "Merging_of_Shade_maps": {                      # the parameters for the shade map creation step
        "merge_shademaps": true,                    # perform the merging of shade maps (bool)
        "output_folder_merged_shademaps": "path",   # path to directory where the merged maps will be placed (str)
        "merged_name": "amsterdam",                 # prefix for the merged maps (str)
        "files_start_time": 900,                    # starting time of the subtile shade map files to be merged (int)
        "files_end_time": 2000,                     # ending time of the subtile shade map files to be merged (int)
        "delete_input_shade": false                 # delete the input seperate subtile shade maps (bool)
    }
}
```


## Cool Spaces Configuration file <a name="heading--2"/>
If no `config_file_cool_spaces` argument is given in *main.py*, the code shall default to .\configuration_files\coolspaceConfig.json.
One can choose to either edit the file in the repository directly, or to download the file, edit it and provide it as argument. 

Note that there are two parameters need to be taken carefully:
> - `shademaps_path`: the shade maps names should be as: `xxxx_time.TIF` such as `amsterdam_900.TIF`.
>   The program needs the time to sort the shade maps into correct time order, and the program only 
>   accept `.TIF` format.
> - `output_shadeGeometry_gpkg`: the file name of shade geometry Geopackage must be {name}_{date}.gpkg, 
>    the date must be YYMMDD, such as 20230621. This is because the routing part needs this format to
>   extract date.

The configuration file is structured as following:

```yaml
{
  "files": {
    "folder_path": "G:\\TUD\\Synthesis\\cool_place\\",              # the folder path for both inputs and outputs
    "gpkg_file": "G:\\TUD\\Synthesis\\cool_place\\coolspace.gpkg",  # the Geopackage contains all the vector input files
    "landuse_file": "ams_landuse_top10NL",                          # the layer name of land use data within the Geopackage
    "road_file": "ams_roads_top10NL",                               # the layer name of road data
    "building_file": "ams_buildings_bagplus",                       # the layer name of building data
    "building_population_file": "ams_bpop_2020",                    # the layer name of building data that contains population statistics
    "street_furniture_file": "ams_bench_osm",                       # the layer name of street furniture data
    "heatrisk_file": "ams_heatrisk_2023",                           # the layer name of heat risk data
    "pet_file": "G:\\TUD\\Synthesis\\cool_place\\ams_PET_average.tiff",        # the file path of PET file, raster format
    "shademaps_path": "G:\\TUD\\Synthesis\\cool_place\\shademaps_20230621\\",  # the folder path that contains all the input shade maps
    "output_identification_file": "coolspace_identification_20230621",         # the layer name of the output identification result
    "output_evaluation_file": "final_cool_space_20230621",                     # the layer name of the output evaluataion result (the final result)
    "output_shadeGeometry_gpkg": "shadeGeoms_20230621.gpkg"                    # the file name of output shade geometry Geopackage, note: this file name
                                                                               # must be {name}_{date}.gpkg, the date must be YYMMDD, such as 20230621.
  },
  "parameters": {
    "hasIdentificationOutput": false,         # if set to true, it will try to read the identification result and skip the identification process
    "outputShadeGeometries": false,           # if set to true, all shade geometries will be output as Geopackages, grouped by shade map index
    "performEvaluation": true,                # if set to false, the evaluation process will be skipped
    "useMultiProcessing": true,               # if set to false, no multi-processing will be executed
    "road_buffer_attribute": "buffer",        # the name of buffer attribute within the road data, which is used to set distance buffer for road polygons
    "building_buffer": 4,                     # the distance buffer for building polygons
    "shade_calculation_mode": "single-day",   # if set to "multi-days", the shade coverage indicator only calculates the daytime range score and a total average score for all days.
    "output_coolspace_type": "public-space",  # if set to "land-use", the output geometry will be the initial land use polygons
    "capacity_search_buffer": 700             # the radius of the circle buffer for searching population around a cool space
  },
  "shade_info_single": {
    "single_day_time_range": [900, 1800],     # specify the daytime range, by default is 900 to 1800, which means 9:00 to 18:00
    "time_interval": 30,                      # time interval of shade maps, 30 means that one shade map represents a 30-minute interval
    "morning_range": [900, 1200],             # specify the morning time range
    "afternoon_range": [1200, 1600],          # specify the afternoon time range
    "late_afternoon_range": [1600, 1800],     # specify the late afternoon time range
    "search_range": null                      # specify a user-define time range such as 9:00 to 14:00 -> [900, 1400]
  },
  "shade_info_multi": {
    "num_days": 2                             # specify the number of days for "multi-days" mode
  }
}

```


## Network Configuration file <a name="heading--3"/>
If no `config_file_network` argument is given in *main.py*, the code shall default to .\configuration_files\network_config.json. 
One can choose to either edit the file in the repository directly, or to download the file, edit it and provide it as argument. 

Note that if you want to skip certain parts, **all** values of that part should be set to empty("").

The configuration file is structured as following:

```yaml
{
    "Dataset_Preparation": {
        "graphml_file": "file_path_of_graph",                                   # File path to GraphML file for the network graph (if this is set then leave "area_name_for_OSM_network" empty)
        "area_name_for_OSM_network": "Amsterdam, Netherlands",                  # Name of the OSM area for network extraction (if this is set then leave "graphml_file" empty)
        "shade_maps_path": "directory_path_to_shade_maps",                      # Directory path to save shade maps
        "new_graphml_files_path": "directory_path_to_new_graphs",               # Directory path to save new GraphML files with shade attributes

        "cool_places_gpkg_dir": "directory_path_to_gpkg",                       # Directory path to GeoPackage files with cool places geometry data
        "cool_places_shp_dir": "directory_path_to_save_shp",                    # Directory path to save shapefiles of cool places geometry
        "cool_places_nodes_path": "directory_path_to_save_nodes"                # Directory path to save calculated cool place nodes
    },
    "Routing": {
        "graph_dir": "directory_path_to_graghs",                                # Directory path containing the main graph (with shade weight attributes) files (if "graph_dir" and "nodes_dir" are set, leave "graph_file" and "nodes_file" empty)
        "nodes_dir": "directory_path_to_nodes",                                 # Directory path containing the cool places nodes data files (if "graph_dir" and "nodes_dir" are set, leave "graph_file" and "nodes_file" empty)
        "graph_file": "file_path_of_graph",                                     # File path for the main routing graph (if "graph_file" and "nodes_file" are set, leave "graph_dir" and "nodes_dir" empty)
        "nodes_file": "file_path_of_nodes",                                     # File path for the cool places nodes data (if "graph_file" and "nodes_file" are set, leave "graph_dir" and "nodes_dir" empty)
        "route_option": "nearest_cool_place or origin_destination",             # Option for selecting the type of route calculation: "nearest_cool_place" or "origin_destination"
        "location_indication_option": "location_name or coordinates",           # Option for specifying start and end location format: "location_name" or "coordinates"
        "origin_name": "Amsterdam Central Station",                             # Name of the origin location 
                                                                                # (required if "location_indication_option" is "location_name")
        "destination_name": "Dam Square",                                       # Name of the destination location 
                                                                                # (required if "route_option" is "origin_destination" and "location_indication_option" is "location_name")
        "origin_latitude": "52.373169",                                         # Latitude of origin 
                                                                                # (required if "location_indication_option" is "coordinates")
        "origin_longitude": "4.890660",                                         # Longitude of origin
                                                                                # (required if "location_indication_option" is "coordinates")
        "destination_latitude": "52.376522",                                    # Latitude of destination 
                                                                                # (required if "route_option" is "origin_destination" and "location_indication_option" is "coordinates")
        "destination_longitude": "4.908490",                                    # Longitude of destination
                                                                                # (required if "route_option" is "origin_destination" and "location_indication_option" is "coordinates")
        "date_time": "2024-10-31 10:05:30"                                      # Date and time for route calculations, format "YYYY-MM-DD HH:MM:SS"
                                                                                # If not set, current time will be used by default
    },
    "Walking_Shed": {
        "graph_file": "file_path_of_graph",                                     # File path to GraphML with shade weight attributes
        "cool_places_path": "file_path_of_cool_places_polygons",                # File path to shapefile with cool places polygons
                                                                                # (which could be obtained from "Dataset_Preparation" -> "cool_places_shp_dir" if previous step is followed)
        "building_shapefile_path": "file_path_of_buildings",                    # File path to building shapefile
        "weight": "length or shade_weight",                                     # Attribute name for calculations: "length" or "shade_weight"
        "output_building_shapefile": "file_path_of_buildings_output",           # Output file path for buildings with distance categories
        "output_cool_place_shapefile": "file_path_of_cool_places_nodes"         # Output file path for cool place nodes shapefile
    }
}
```
