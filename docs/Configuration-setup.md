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


## Network Configuration file <a name="heading--3"/>