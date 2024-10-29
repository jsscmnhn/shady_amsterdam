## Content

**[1. Shade Calculation Functions](#heading--1)**

  * [1.1. General Functions (*functions.py*)](#heading--1-1)
     * [raster_center_coords](#heading--1-1-1)
     * [create_affine_transform](#heading--1-1-2)
     * [write_output](#heading--1-1-3)
  * [1.2. Data Acquisition (*collect_data.py*) ](#heading--1-2)
      * [download_las_tiles](#heading--1-2-1)
      * [download_and_extract](#heading--1-2-2)
      * [merge_tif_files](#heading--1-2-3)
      * [download_raster_tiles](#heading--1-2-4)
      * [download_wfs_data](#heading--1-2-5)
      * [setup_WFS_download](#heading--1-2-6)
  * [1.3. CHM Creation (*create_first_chm_parallel.py*) ](#heading--1-3)  
     * [median_filter_chm](#heading--1-3-1)
     * [extract_vegetation_points`](#heading--1-3-2)
     * [chm_creation`](#heading--1-3-3)
     * [interpolation_vegetation`](#heading--1-3-4)
     * [process_single_laz_file`](#heading--1-3-5)
     * [process_laz_files`](#heading--1-3-6)
  * [1.4. DSM and final CHM creation (*create_dsm_and_chm_parallel.py*) ](#heading--1-4)
     * [get_bbox](#heading--1-4-1) 
     * [crop_raster](#heading--1-4-2) 
     * [extract_center_cells](#heading--1-4-3) 
     * [fill_raster](#heading--1-4-4)
     * [chm_finish](#heading--1-4-5) 
     * [replace_buildings](#heading--1-4-6) 
     * [load_buildings](#heading--1-4-7) 
     * [extract_tilename](#heading--1-4-8) 
     * [process_single_file](#heading--1-4-9) 
     * [process_files](#heading--1-4-10) 
     * [process_folders](#heading--1-4-11) 
  * [1.5. Creating Shade Maps (*shade_parallel.py*) ](#heading--1-5)
     * [process_chm_dsm](#heading--1-5-1) 
     * [run_shade_calculation](#heading--1-5-2) 
     * [process_folders](#heading--1-5-3) 
  * [1.6. Merging Shade Maps (*merge_shademaps.py*) ](#heading--1-6)
     * [merge_tif_files_by_time](#heading--1-6-1)
  * [1.7. Running code (*main_shade.py*) ](#heading--1-7)
     * [flatten_dict](#heading--1-7-1)
     * [read_config](#heading--1-7-2)
	
**[2. Cool Spaces Functions](#heading--2)**

  * [2.1. Process functions](#heading--2-1)
     * [read_config](#heading--2-1-1)
     * [convert_to_datetime](#heading--2-1-2)
     * [compute_search_range](#heading--2-1-3)
     * [identification](#heading--2-1-4)
     * [evaluation](#heading--2-1-5)
     * [list_to_string](#heading--2-1-6)
     * [drop_or_wkt](#heading--2-1-7)
     * [output_all_shade_geoms](#heading--2-1-8)
  * [2.2. Cool Space Class for Identification](#heading--2-2)
     * [2.2.1. Initialization and Geometry Functions](#heading--2-2-1)
       * [__init__](#heading--2-2-1-1)
       * [clip](#heading--2-2-1-2)
     * [2.2.2. Shade Calculation Functions](#heading--2-2-2)
       * [calculate_shade](#heading--2-2-2-1)
       * [calculate_shade_for_raster](#heading--2-2-2-2)
       * [calculate_shade_multi](#heading--2-2-2-3)
     * [2.2.3. Shade Geometry Retrieval Functions](#heading--2-2-3)
       * [get_shade_geometries](#heading--2-2-3-1)
       * [get_cool_spaces](#heading--2-2-3-2)
     * [2.2.4. Evaluation Functions](#heading--2-2-4)
       * [evaluate_shade_coverage](#heading--2-2-4-1)
  * [2.3. CoolEval Class](#heading--2-3)
     * [2.3.1. Initialization Functions](#heading--2-3-1)
       * [__init__](#heading--2-3-1-1)
     * [2.3.2. Walking Shed Functions](#heading--2-3-2)
       * [calculate_walking_shed_origin](#heading--2-3-2-1)
       * [calculate_walking_shed](#heading--2-3-2-2)
       * [calculate_walking_shed_multi](#heading--2-3-2-3)
     * [2.3.3. Resident and Capacity Evaluation Functions](#heading--2-3-3)
       * [evaluate_resident](#heading--2-3-3-1)
       * [evaluate_capacity](#heading--2-3-3-2)
     * [2.3.4. Shade Furniture and Heat Risk Evaluation Functions](#heading--2-3-4)
       * [evaluate_sfurniture](#heading--2-3-4-1)
       * [evaluate_heatrisk](#heading--2-3-4-2)
     * [2.3.5. PET Evaluation Functions](#heading--2-3-5)
       * [eval_pet](#heading--2-3-5-1)
       * [eval_pet_multi](#heading--2-3-5-2)
     * [2.3.6. Aggregation and Recommendation Functions](#heading--2-3-6)
       * [aggregate_to_cool_places](#heading--2-3-6-1)
       * [final_recom](#heading--2-3-6-2)
     * [2.3.7. Export Functions](#heading--2-3-7)
       * [export_eval_gpkg](#heading--2-3-7-1)
  * [2.4. Building Class](#heading--2-4)
     * [__init__](#heading--2-4-1)
     * [create_buffer](#heading--2-4-2)
  * [2.5. Road Class](#heading--2-5)
     * [__init__](#heading--2-5-1)
     * [assign_buffer_hardSurface](#heading--2-5-2)
     * [assign_buffer_roadtype](#heading--2-5-3)
     * [create_attribute](#heading--2-5-4)
     * [create_buffer](#heading--2-5-5)

**[3. Network Functions](#heading--3)**

  * [3.1. TODO](#heading--3-1)


----


## Shade Calculation Functions <a name="heading--1"/>

### General Functions (*functions.py*) <a name="heading--1-1"/>
#### <span style="color: red;">`raster_center_coords`</span><span style="color: gray;">(min_x, max_x, min_y, max_y, resolution)</span> <a name="heading--1-1-1"/>

> Compute the center xy coordinates of all cells in a grid.
>
>**PARAMETERS**  
>- **`min_x`**, **`max_x`** — *float*: Minimum and maximum x coordinates of the grid.
>- **`min_y`**, **`max_y`** — *float*: Minimum and maximum y coordinates of the grid.
>- **`resolution`** — *float*: The length of each cell; function only works for square cells.
>
>**RETURNS**  
>- *numpy.ndarray*: `grid_center_x` - A grid where each cell contains the center x-coordinate.
>- *numpy.ndarray*: `grid_center_y` - A grid where each cell contains the center y-coordinate.

#### <span style="color: red;">`create_affine_transform`</span><span style="color: gray;">(top_left_x, top_left_y, res)</span> <a name="heading--1-1-2"/>

> Create Affine transform for the `write_output` function.
>
> **PARAMETERS**  
> - **`top_left_x`** — *float*: The x-coordinate of the top-left corner.
> - **`top_left_y`** — *float*: The y-coordinate of the top-left corner.
> - **`res`** — *float*: The resolution of the grid.
>
> **RETURNS**  
> - *Affine*: An affine transformation object.

#### <span style="color: red;">`write_output`</span><span style="color: gray;">(dataset, output, transform, name, change_nodata=False)</span> <a name="heading--1-1-3"/>

> Write a grid (numpy array) to the given path (name, string) as a .TIF file.
>
> **PARAMETERS**  
> - **`dataset`** — *object*: Can be either a rasterio dataset (for rasters) or a laspy dataset (for point clouds).
> - **`output`** — *array*: The output grid, a numpy array.
> - **`name`** — *string*: The name of the output file (will overwrite).
> - **`transform`** — *affine*: A user-defined rasterio affine object, used for transforming pixel coordinates to spatial coordinates.
> - **`change_nodata`** — *boolean*: If true, use a no data value of -9999; if false, use the dataset's no data value.
> 
> **RETURNS**  
> - (nothing): The function saves the .TIF directly at the output file path.
 
---

### Data Acquisition (*collect_data.py*) <a name="heading--1-2"/>

#### <span style="color: red;">`download_las_tiles`</span><span style="color: gray;">(tile_list_file, output_folder)</span> <a name="heading--1-2-1"/>

> Download LAZ files for each GeoTiles subtile specified in a text file. \
> *Note: the base link for downloading is hard coded as "https://geotiles.citg.tudelft.nl/AHN4_T", this function will have to be updated if GeoTiles changes*
>
> **PARAMETERS**  
> - **`tile_list_file`** — *str*: Path to the text file containing the list of subtiles to download.
> - **`output_folder`** — *str*: Path to the folder where the downloaded files should be saved.
>
> **RETURNS**  
> - (nothing): The function writes output files directly to the specified `output_folder/{subtile}`.

#### <span style="color: red;">`download_and_extract`</span><span style="color: gray;">(url, file_path, output_folder)</span> <a name="heading--1-2-2"/>

> Download a ZIP file from a URL, extract its contents, and delete the ZIP file.
>
> **PARAMETERS**  
> - **`url`** — *str*: URL of the ZIP file to download.
> - **`file_path`** — *str*: Path to save the downloaded ZIP file.
> - **`output_folder`** — *str*: Path to the folder where the extracted contents should be saved.
>
> **RETURNS**  
> - (nothing)

#### <span style="color: red;">`merge_tif_files`</span><span style="color: gray;">(input_folder, output_file, file_prefix, nodata_value=-9999)</span> <a name="heading--1-2-3"/>

> Merge TIF files in a folder with a specific prefix into a single raster file.
>
> **PARAMETERS**  
> - **`input_folder`** — *str*: Path to the folder containing the TIF files.
> - **`output_file`** — *str*: Path to the output raster file.
> - **`file_prefix`** — *str*: Prefix of files to merge ("M_" for DTM or "R_" for DSM).
> - **`nodata_value`** — *int*, *optional*: Value to replace the nodata value. Default: -9999.
>
> **RETURNS**  
> - *list*: A list containing the name of the tile and the extent of that tile, for downloading the building data.
> - The function writes the merged file directly to the specified `output_file`.

#### <span style="color: red;">`download_raster_tiles`</span><span style="color: gray;">(tile_list_file, output_folder, name)</span> <a name="heading--1-2-4"/>

> Download DSM and DTM files for each GeoTiles tile specified in a text file. \
>*Note: the base links for downloading are hard coded as "https://ns_hwh.fundaments.nl/hwh-ahn/ahn4/02a_DTM_0.5m" and "https://ns_hwh.fundaments.nl/hwh-ahn/ahn4/03a_DSM_0.5m", this function will have to be updated if GeoTiles changes*
> 
>**PARAMETERS**  
> - **`tile_list_file`** — *str*: Path to the text file containing the list of tiles to download.
> - **`output_folder`** — *str*: Directory where the downloaded and unzipped files will be saved.
> - **`name`** — *str*: Name of the output raster file.
>
> **RETURNS**  
> - *list*: A list containing the name of the tile and the extent of that tile, for downloading the building data at a later step.
> - The function writes output files directly to the specified `output_folder`.

#### <span style="color: red;">`download_wfs_data`</span><span style="color: gray;">(wfs_url, layer_name, bbox, gpkg_name, output_folder, tile_name)</span> <a name="heading--1-2-5"/>

> Download data from a WFS server in batches and save it to a GeoPackage.
>
> **PARAMETERS**  
> - **`wfs_url`** — *str*: URL of the WFS service.
> - **`layer_name`** — *str*: The layer name to download.
> - **`bbox`** — *tuple*: Bounding box as (minx, miny, maxx, maxy).
> - **`gpkg_name`** — *str*: Name for the output GeoPackage file.
> - **`tile_name`** — *str*: Layer name for saving in the GeoPackage.
>
> **RETURNS**  
> - (none): Saves a GeoPackage file to the given `{output_gpkg}` at layer `{tile_name}`.

#### <span style="color: red;">`setup_WFS_download`</span><span style="color: gray;">(gpkg_name, tile_bounds, output_folder)</span> <a name="heading--1-2-6"/>

> Collecting the needed information for downloading the WFS data.
>
> **PARAMETERS**  
> - **`gpkg_name`** — *str*: Name for the output GeoPackage file.
> - **`tile_bounds`** — *list*: List containing nested lists where at [0] the tile name is saved and at [1] the bounds of the tile.
> - **`output_folder`** — *str*: Path to the output GeoPackage file.
>
> **RETURNS**  
> - (none): Saves a GeoPackage file to the given `{gpkg_name}` at layer `{tile_name}`.
 
 ---

### CHM Creation (*create_first_chm_parallel.py*) <a name="heading--1-3"/>

#### <span style="color: red;">`median_filter_chm`</span><span style="color: gray;">(chm_array, nodata_value=-9999, size=3)</span> <a name="heading--1-3-1"/>

> Apply a median filter to a numpy array, ignoring the NoData values.
>
> **PARAMETERS**  
> - **`chm_array`** — *np.ndarray*: The array representing the height values of the CHM.
> - **`nodata_value`** — *float*: Value representing NoData in the input raster.
> - **`size`** — *int*: Size of the median filter. It defines the footprint of the filter.
>
> **RETURNS**  
> - *np.ndarray*: The smoothed CHM array.


#### <span style="color: red;">`extract_vegetation_points`</span><span style="color: gray;">(LasData, ndvi_threshold=0.1, pre_filter=False)</span> <a name="heading--1-3-2"/>

> Extract vegetation points based on AHN classification and NDVI threshold.
>
> **PARAMETERS**  
> - **`LasData`** — *laspy.LasData*: Input point cloud data in LAS format.
> - **`ndvi_threshold`** — *float*: The NDVI threshold for identifying vegetation points. NDVI values greater than this threshold are considered vegetation.
> - **`pre_filter`** — *bool*: If True, applies an additional filter to remove vegetation points below a certain height threshold (1.5 meters above the lowest vegetation point).
>
> **RETURNS**  
> - *laspy.LasData*: A new LasData object containing only the filtered vegetation points based on the specified criteria.


#### <span style="color: red;">`chm_creation`</span><span style="color: gray;">(LasData, vegetation_data, output_filename, resolution=0.5, smooth=False, nodata_value=-9999, filter_size=3)</span> <a name="heading--1-3-3"/>

> Create a CHM from LiDAR vegetation data and save it as a raster.
>
> **PARAMETERS**  
> - **`LasData`** — *laspy.LasData*: Input LiDAR point cloud data used for metadata and output CRS.
> - **`vegetation_data`** — *tuple*: A tuple containing:
>   - **`veg_raster`** — *numpy.ndarray*: The array representing the height values of vegetation.
>   - **`grid_centers`** — *tuple of numpy.ndarrays*: Contains two arrays (x, y) with the coordinates of the center points of each grid cell.
> - **`output_filename`** — *str*: The name of the output .tif file for saving the CHM.
> - **`resolution`** — *float*: The spatial resolution of the output raster in the same units as the input data (default: 0.5).
> - **`smooth`** — *bool*: If True, applies a median filter to smooth the CHM.
> - **`nodata_value`** — *float*: The value for NoData pixels (default: -9999).
> - **`filter_size`** — *int*: Size of the median filter (default: 3).
>
> **RETURNS**  
> - *(nothing)*: The function saves the CHM as a raster file (.tif) to the specified output path.

#### <span style="color: red;">`interpolation_vegetation`</span><span style="color: gray;">(LasData, veg_points, resolution, no_data_value=-9999)</span> <a name="heading--1-3-4"/>

> Create a vegetation raster from LAS vegetation points using Laplace interpolation.
>
> **PARAMETERS**  
> - **`LasData`** — *laspy.LasData*: Input LiDAR point cloud data.
> - **`veg_points`** — *laspy.LasData*: Vegetation points to be interpolated.
> - **`resolution`** — *float*: Resolution of the raster.
> - **`no_data_value`** — *int*: Value for no data (default: -9999).
>
> **RETURNS**  
> - *np.ndarray*: Generated raster for vegetation.
> - *tuple*: Grid of x, y center coordinates for each raster cell.


#### <span style="color: red;">`process_single_laz_file`</span><span style="color: gray;">(file_path, output_folder, ndvi_threshold=0.0, resolution=0.5, remove=False, smooth_chm=False, filter_size=3, pre_filter=False)</span> <a name="heading--1-3-5"/>

> Process a LAZ file to extract vegetation points and generate a CHM.
>
> **PARAMETERS**  
> - **`file_path`** — *str*: The file_path containing the folder containing the input .LAZ file.
> - **`output_folder`** — *str*: The folder where the output CHM .tif files will be saved.
> - **`ndvi_threshold`** — *float*: The NDVI threshold for classifying vegetation points.
> - **`resolution`** — *float*: The resolution of the output CHM rasters, defining the size of each pixel (default: 0.5).
> - **`remove`** — *bool*: If True, deletes the original .LAZ files after processing (default: False).
> - **`smooth_chm`** — *bool*: If True, applies smoothing to the CHM using a median filter (default: False).
> - **`filter_size`** — *int*: Size of the median filter to use if smoothing (default: 3).
> - **`pre_filter`** — *bool*: If True, applies an additional filter to remove vegetation points below a certain height threshold (1.5 meters above the lowest vegetation point).
>
> **RETURNS**  
> - *(Nothing)*: The function processes a .LAZ file, creates a corresponding CHM .tif file, and saves it to the output folder. Optionally deletes the original .LAZ files if `remove` is set to True.


#### <span style="color: red;">`process_laz_files`</span><span style="color: gray;">(input_folder, output_folder, ndvi_threshold=0.0, resolution=0.5, remove=False, smooth_chm=False, filter_size=3, pre_filter=False, max_workers=4)</span> <a name="heading--1-3-6"/>

> Process a folder of LAZ files in parallel to extract vegetation points and generate CHMs.
>
> **PARAMETERS**  
> - **`input_folder`** — *str*: The folder containing the input .LAZ files.
> - **`output_folder`** — *str*: The folder where the output CHM .tif files will be saved.
> - **`ndvi_threshold`** — *float*: The NDVI threshold for classifying vegetation points.
> - **`resolution`** — *float*: The resolution of the output CHM rasters, defining the size of each pixel (default: 0.5).
> - **`remove`** — *bool*: If True, deletes the original .LAZ files after processing (default: False).
> - **`smooth_chm`** — *bool*: If True, applies smoothing to the CHM using a median filter (default: False).
> - **`filter_size`** — *int*: Size of the median filter to use if smoothing (default: 3).
> - **`pre_filter`** — *bool*: If True, applies an additional filter to remove vegetation points below a certain height threshold (1.5 meters above the lowest vegetation point).
> - **`max_workers`** — *int*: The number of parallel processes to process the LAZ files.
>
> **RETURNS**  
> - *(Nothing)*: The function processes the .LAZ files, creates corresponding CHM .tif files, and saves them to the output folder. Optionally deletes the original .LAZ files if `remove` is set to True.

---

### DSM and Final CHM Creation (*create_dsm_and_chm_parallel.py*) <a name="heading--1-4"/>

#### <span style="color: red;">`get_bbox`</span><span style="color: gray;">(raster_paths)</span> <a name="heading--1-4-1"/>

> Compute the overlapping bounding box of the CHM, DSM, and DTM file.
>
> **Input:**  
> - `raster_paths` (list of strings): Paths to the CHM, DSM, and DTM .TIF files.
>
> **Output:**  
> - *rasterio.coords.BoundingBox*: Bounding box of the raster.


#### <span style="color: red;">`crop_raster`</span><span style="color: gray;">(raster_path, bbox, no_data=-9999)</span> <a name="heading--1-4-2"/>

> Crop the input rasters to the size of the overlapping bounding box.
>
> **Input:**  
> - `raster_path` (string): Path to the .TIF file.  
> - `bbox` (4-tuple): Overlapping bounding box of the CHM, DSM, and DTM file.  
> - `no_data` (int, optional): NoData value to replace source NoData value with.
>
> **Output:**  
> - *2D numpy array*: Cropped raster data.  
> - *src.window_transform(window)*: Affine transform matrix for the given window.  
> - *src.crs*: A PROJ4 dict representation of the CRS of the input raster.


#### <span style="color: red;">`extract_center_cells`</span><span style="color: gray;">(cropped_data, no_data=-9999)</span> <a name="heading--1-4-3"/>

> Extract the values of each cell in the input data and save these with the x and y (row and col) indices.
>
> **Input:**  
> - `cropped_data` (2D numpy array): Cropped raster data.  
> - `no_data` (int, optional): NoData value to replace source NoData value with.
>
> **Output:**  
> - *List*: Containing x, y, and z coordinates of the cells.


#### <span style="color: red;">`fill_raster`</span><span style="color: gray;">(cropped_data, nodata_value, transform, speed_up=False)</span> <a name="heading--1-4-4"/>

> Fill the NoData values of a given raster using Laplace interpolation.
>
> **Input:**  
> - `cropped_data` (2D numpy array): Cropped raster data.  
> - `nodata_value` (int): NoData value to replace NaN after interpolation with.  
> - `transform` (rasterio transform): Affine transform matrix.  
> - `speed_up` (boolean): If True, checks for a large NoData area and uses linear interpolation if so. Default is False.
>
> **Output:**  
> - *2D numpy array*: Filled raster data with first and last rows and columns removed to ensure there are no NoData values from Laplace interpolation.  
> - *rasterio transform*: Affine transform matrix reflecting the one-column one-row removal shift.


#### <span style="color: red;">`chm_finish`</span><span style="color: gray;">(chm_array, dtm_array, transform, min_height=2, max_height=40)</span> <a name="heading--1-4-5"/>

> Finish the CHM file by first removing the ground height, then removing vegetation height below and above a certain range.
>
> **Input:**  
> - `chm_array` (2D numpy array): Cropped raster array of the CHM.  
> - `dtm_array` (2D numpy array): Cropped raster array of the filled DSM.  
> - `transform` (rasterio transform): Affine transform matrix.  
> - `min_height` (float, optional): Minimal height for vegetation to be included.  
> - `max_height` (float, optional): Maximum height for vegetation to be included.
>
> **Output:**  
> - *2D numpy array*: Array of the CHM with normalized height and min and max heights removed.  
> - *rasterio transform*: Affine transform matrix reflecting the one-column one-row removal shift.


#### <span style="color: red;">`replace_buildings`</span><span style="color: gray;">(filled_dtm, dsm_buildings, buildings_geometries, transform)</span> <a name="heading--1-4-6"/>

> Replace the values of the filled DTM with the values of the filled DSM, if there is a building at the value location.
>
> **Input:**  
> - `filled_dtm` (2D np array): Filled array of the cropped AHN DTM.  
> - `dsm_buildings` (2D np array): Filled array of the cropped AHN DSM.  
> - `building_geometries` (list): A list of the building geometries.  
> - `transform` (rasterio transform): Affine transform matrix.
>
> **Output:**  
> - *2D numpy array*: A numpy array representing the final DSM, containing only ground and building heights.


#### <span style="color: red;">`load_buildings`</span><span style="color: gray;">(buildings_path, layer)</span> <a name="heading--1-4-7"/>

> Load the building shapes from a geopackage file.
>
> **Input:**  
> - `buildings_path` (string): Path to the geopackage file.  
> - `layer` (string): (Tile) name of the layer of buildings to be used.
>
> **Output:**  
> - *List of dictionaries*: A list of geometries in GeoJSON-like dictionary format. Each dictionary represents a building geometry with its spatial coordinates.


#### <span style="color: red;">`extract_tilename`</span><span style="color: gray;">(filename)</span> <a name="heading--1-4-8"/>

> Extract the name of the AHN tile from the file name.
>
> **Input:**  
> - `filename` (string): The name of the input chm.tif.
>
> **Output:**  
> - *string*: The name of the AHN tile.


#### <span style="color: red;">`process_single_file`</span><span style="color: gray;">(chm_path, dtm_path, dsm_path, building_geometries, output_base_folder, nodata_value=-9999, speed_up=False, min_height=2, max_height=40)</span> <a name="heading--1-4-9"/>

> Process from one CHM file to a final DSM and CHM.
>
> **Input:**  
> - `chm_path` (string): Path to the input CHM .tif file.  
> - `dtm_path` (string): Path to the (merged DTM) .tif file.  
> - `dsm_path` (string): Path to the (merged DSM) .tif file.  
> - `building_geometries` (list): Path to the geopackage file containing building geometries.  
> - `output_base_folder` (string): Base path for saving the output DSM and CHM files.  
> - `nodata_value` (int): NoData value for raster processing (default: -9999).  
> - `speed_up` (boolean): If True, checks for a large NoData area in the DTM file and uses a different interpolation method if so. Default is False.  
> - `min_height` (float, optional): Minimal height for vegetation to be included in final CHM.  
> - `max_height` (float, optional): Maximum height for vegetation to be included in final CHM.
>
> **Output:**  
> - None: The function writes output files directly to the specified `output_base_folder`.  
> - Output files: For each input CHM file, the function generates:  
>   - A DSM file with building data incorporated: `DSM_<tile>_<subtile_number>.tif`  
>   - A processed CHM file: `CHM_<tile>_<subtile_number>.tif`  
>   - These files are saved in folders named after the tile in `output_base_folder`.


#### <span style="color: red;">`process_files`</span><span style="color: gray;">(chm_folder, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value=-9999, max_workers=4, speed_up=False, min_height=2, max_height=40)</span> <a name="heading--1-4-10"/>

> Run the whole process of creating the final DSM and CHM.
>
> **Input:**  
> - `chm_folder` (string): Path to the folder containing input CHM .tif files.  
> - `dtm_path` (string): Path to the DTM .tif file.  
> - `dsm_path` (string): Path to the DSM .tif file.  
> - `buildings_path` (string): Path to the geopackage file containing building geometries.  
> - `output_base_folder` (string): Base path for saving the output DSM and CHM files.  
> - `nodata_value` (int): NoData value for raster processing (default: -9999).  
> - `max_workers` (int): Number of parallel processes.  
> - `speed_up` (boolean): If True, checks for a large NoData area in the DTM file and uses a different interpolation method if so. Default is False.  
> - `min_height` (float, optional): Minimal height for vegetation to be included in final CHM.  
> - `max_height` (float, optional): Maximum height for vegetation to be included in final CHM.
>
> **Output:**  
> - (Nothing): The function writes output files directly to the specified `output_base_folder`.  
> - Output files: For each input CHM file, the function generates:  
>   - A DSM file with building data incorporated: `DSM_<tile>_<subtile_number>.tif`  
>   - A processed CHM file: `CHM_<tile>_<subtile_number>.tif`  
>   - These files are saved in folders named after the tile in `output_base_folder`.


#### <span style="color: red;">`process_folders`</span><span style="color: gray;">(base_chm_folder, dtm_path, dsm_path, buildings_path, output_base_folder, max_workers=4, speed_up=False, min_height=2, max_height=40)</span> <a name="heading--1-4-11"/>

> Process each folder containing CHM files concurrently.
>
> **Input:**  
> - `base_chm_folder` (string): Path to the base folder containing subfolders of CHM files.  
> - `dtm_path` (string): Path to the DTM .tif file.  
> - `dsm_path` (string): Path to the DSM .tif file.  
> - `buildings_path` (string): Path to the geopackage file containing building geometries.  
> - `output_base_folder` (string): Base path for saving the output DSM and CHM files.  
> - `nodata_value` (int): NoData value for raster processing (default: -9999).  
> - `max_workers` (int): Number of parallel processes.
>
> **Output:**  
> - (Nothing): The function writes output files directly to the specified `output_base_folder`.

---

### Creating Shade Maps (*shade_parallel.py*) <a name="heading--1-5"/>

#### <span style="color: red;">`process_chm_dsm`</span><span style="color: gray;">(chm_filename, dsm_filename, folder_path, output_dir, date, start_time=10, end_time=21, interval=30, use_chm=True, trans=10, trunkheight=25)</span> <a name="heading--1-5-1"/>

> Function to process a single DSM and CHM file pair and create shade maps from them for the given time frame and interval.
>
> **Input:**  
> - `chm_filename` (str): Name of the CHM file to be processed.  
> - `dsm_filename` (str): Name of the DSM file to be processed.  
> - `folder_path` (str): Path to the folder containing the CHM and DSM files.  
> - `output_dir` (str): Directory where the output results will be saved.  
> - `date` (str): Date parameter used in the shade calculation.  
> - `start_time` (int, optional): Starting hour for shade calculation (default: 10).  
> - `end_time` (int, optional): Ending hour for shade calculation (default: 21).  
> - `interval` (int, optional): Time interval for the calculations in minutes (default: 30).  
> - `use_chm` (bool, optional): Whether to use the CHM file in the shade calculation (default: True).  
> - `trans` (int, optional): Transmissivity value for calculations, representing the percentage for tree shade (default: 10).  
> - `trunkheight` (int, optional): Trunk height for vegetation calculations, representing the percentage of CHM height for trunk height (default: 25).
>
> **Output:**  
> - (Nothing): This function prints a message upon completion of shade calculation.


#### <span style="color: red;">`run_shade_calculation`</span><span style="color: gray;">(file_pairs, output_base_folder, date, start_time, end_time, interval, use_chm=True, trans=10, trunkheight=25, max_workers=4)</span> <a name="heading--1-5-2"/>

> Process CHM and DSM file pairs in parallel to create shade maps.
>
> **Input:**  
> - `file_pairs` (list): List of tuples containing `(CHM file, DSM file, folder_path, tile)`.  
> - `output_base_folder` (str): Base path for saving the output files.  
> - `date` (str): Date parameter used in the shade calculation.  
> - `max_workers` (int, optional): Maximum number of worker threads for concurrent processing (default: 4).  
> - `start_time` (int): Starting hour for shade calculation.  
> - `end_time` (int): Ending hour for shade calculation.  
> - `interval` (int): Time interval for calculations in minutes.  
> - `trans` (int, optional): Transmissivity value for calculations, representing the percentage for tree shade (default: 10).  
> - `trunkheight` (int, optional): Trunk height for vegetation calculations, representing the percentage of CHM height for trunk height (default: 25).
>
> **Output:**  
> - (Nothing): This function prints messages regarding processing status and creates output directories.


#### <span style="color: red;">`process_folders`</span><span style="color: gray;">(base_folder, output_base_folder, date, start_time=9, end_time=20, interval=30, use_chm=True, trans=10, trunkheight=25, max_workers=4)</span> <a name="heading--1-5-3"/>

> Function to process all subfolders in a base folder for shade calculations.
>
> **Input:**  
> - `base_folder` (str): Path to the base folder containing subfolders with CHM and DSM files.  
> - `output_base_folder` (str): Base path for saving the output files.  
> - `date` (str): Date parameter used in the shade calculation.  
> - `max_workers` (int, optional): Maximum number of worker threads for concurrent processing (default: 4).  
> - `start_time` (int, optional): Starting hour for shade calculation (default: 9).  
> - `end_time` (int, optional): Ending hour for shade calculation (default: 20).  
> - `interval` (int, optional): Time interval for calculations in minutes (default: 30).  
> - `trans` (int, optional): Transmissivity value for calculations, representing the percentage for tree shade (default: 10).  
> - `trunkheight` (int, optional): Trunk height for vegetation calculations, representing the percentage of CHM height for trunk height (default: 25).
>
> **Output:**  
> - (Nothing): This function prints messages regarding the processing status for each folder.

---

### Merging Shade Maps (*merge_shademaps.py*) <a name="heading--1-6"/>

#### <span style="color: red;">`merge_tif_files_by_time`</span><span style="color: gray;">(main_folder, output_folder, merged_name, start_time=900, end_time=2000, delete_input_files=False, nodata_value=-9999)</span> <a name="heading--1-6-1"/>

> Merge multiple TIF files from specified subfolders into a single file for each specified time step.
>
> **Input:**  
> - `main_folder` (str): Path to the main folder containing subfolders with TIF files to be merged.  
> - `output_folder` (str): Path to the folder where the merged TIF files will be saved.  
> - `merged_name` (str): Base name for the merged output files. Each file will be suffixed with the date and time.  
> - `nodata_value` (int, optional): Value to represent 'no data' in the output TIF files (default: -9999).  
> - `start_time` (int, optional): Starting time for filtering TIF files based on their filenames (default: 900).  
> - `end_time` (int, optional): Ending time for filtering TIF files based on their filenames (default: 2000).  
> - `delete_input_files` (bool, optional): If True, deletes the original TIFF files after merging (default: False).
>
> **Output:**  
> - (Nothing): The function saves merged TIFF files directly to the specified `output_folder` with the naming format `{merged_name}_{date_str}_{time}.TIF`.

---
### Running Code (*main_shade.py*) <a name="heading--1-7"/>

#### <span style="color: red;">`flatten_dict`</span><span style="color: gray;">(nested_dict, parent_key='', sep='_')</span> <a name="heading--1-7-1"/>

> Flatten a nested dictionary, ignoring the first layer.
>
> **Input:**  
> - `nested_dict` (dict): The dictionary to flatten.  
> - `parent_key` (str): The base key string for the flattened keys.  
> - `sep` (str): Separator to use between concatenated keys.
>
> **Output:**  
> - *dictionary*: A flat dictionary.


#### <span style="color: red;">`read_config`</span><span style="color: gray;">(file_path)</span> <a name="heading--1-7-2"/>

> Read the contents of a JSON configuration file and save them to a dictionary.
>
> **Input:**  
> - `file_path` (str): The path to the configuration file.
>
> **Output:**  
> - *dictionary*: The configuration dictionary containing all required parameters.
 
 ---

## Cool Spaces Functions <a name="heading--2"/>
### Process Functions (*coolspace_process.py*) <a name="heading--2-1"/>

#### <span style="color: red;">`read_config`</span><span style="color: gray;">(filename)</span> <a name="heading--2-1-1"/>

> Load configuration settings from a JSON file.
>
>**PARAMETERS**  
>- **`filename`** — *str*: The path to the JSON configuration file.
>
>**RETURNS**  
>- *dict*: `config` - A dictionary containing the parsed configuration settings.

#### <span style="color: red;">`convert_to_datetime`</span><span style="color: gray;">(time)</span> <a name="heading--2-1-2"/>

> Convert an integer in HHMM format to a `datetime` object.
>
>**PARAMETERS**  
>- **`time`** — *int*: The time in HHMM format, where the first two digits represent hours and the last two represent minutes.
>
>**RETURNS**  
>- *datetime*: A `datetime` object set to October 16, 2024, with the specified hour and minute. The date components (year, month, day) are placeholders and do not impact the function's purpose.

#### <span style="color: red;">`compute_search_range`</span><span style="color: gray;">(search_start_time, search_end_time, start_time, end_time, time_interval)</span> <a name="heading--2-1-3"/>

> Calculate the index range for shade maps based on a specified time interval within a given range.
>
>**PARAMETERS**  
>- **`search_start_time`**, **`search_end_time`** — *int*: The starting and ending times in HHMM format for the desired search range.
>- **`start_time`**, **`end_time`** — *int*: The beginning and ending times in HHMM format for the available shade maps.
>- **`time_interval`** — *int*: The interval in minutes between each shade map.
>
>**RETURNS**  
>- *list*: A list containing `start_idx` and `end_idx`, which are the computed index range for the shade maps based on the specified search range. If the search range exceeds the time bounds of the shade maps, returns `[None, None]` and prints an error message.

#### <span style="color: red;">`identification`</span><span style="color: gray;">(coolspace_file, road_file, building_file, shademaps_path, road_buffer_attri, building_buffer_num=4, mode='single-day', num_days=None, single_day_time_range=None, time_interval=None, search_range=None, morning_range=None, afternoon_range=None, late_afternoon_range=None, output_coolspace_type='land-use', useMultiProcessing=False)</span> <a name="heading--2-1-4"/>

> Identify and evaluate cool spaces based on shade coverage across different time ranges, area ratios and output configurations.
>
>**PARAMETERS**  
>- **`coolspace_file`** — *GeoDataFrame*: GeoDataFrame containing cool space data.
>- **`road_file`** — *GeoDataFrame*: GeoDataFrame containing road data.
>- **`building_file`** — *GeoDataFrame*: GeoDataFrame containing building data.
>- **`shademaps_path`** — *str*: Path to the directory containing shade map files.
>- **`road_buffer_attri`** — *str*: Attribute name in `road_file` used to create a buffer.
>- **`building_buffer_num`** — *int, optional*: Buffer size for buildings; defaults to 4.
>- **`mode`** — *str, optional*: Calculation mode, either `'single-day'` or `'multi-days'`; defaults to `'single-day'`.
>- **`num_days`** — *int, optional*: Number of days, required if `mode` is `'multi-days'`.
>- **`single_day_time_range`** — *list, optional*: Start and end times in HHMM format for a single-day range.
>- **`time_interval`** — *int, optional*: Time interval in minutes for shade map calculations.
>- **`search_range`**, **`morning_range`**, **`afternoon_range`**, **`late_afternoon_range`** — *list, optional*: Specific time ranges in HHMM format for different periods within the single-day mode.
>- **`output_coolspace_type`** — *str, optional*: Type of output geometry, either `'land-use'` or `'public-space'`; defaults to `'land-use'`.
>- **`useMultiProcessing`** — *bool, optional*: Enables multiprocessing for shade calculation if set to `True`; defaults to `False`.
>
>**RETURNS**  
>- *GeoDataFrame*: `output_gdf` - A GeoDataFrame containing cool space polygons with shade evaluations. The output geometry type depends on `output_coolspace_type`, with other geometries converted to WKT format and stored in columns.

#### <span style="color: red;">`evaluation`</span><span style="color: gray;">(coolspace, building_population_file, bench_file, heatrisk_file, pet_file, gpkg_file, output_layer, search_buffer=700, single_day_time_range=None, time_interval=None, search_range=None, morning_range=None, afternoon_range=None, late_afternoon_range=None, useMultiProcessing=False)</span> <a name="heading--2-1-5"/>

> Evaluate cool spaces by analyzing their accessibility, capacity, street furniture, heat risk, and PET scores based on spatial and temporal parameters.
>
>**PARAMETERS**  
>- **`coolspace`** — *GeoDataFrame*: GeoDataFrame containing cool space data to evaluate.
>- **`building_population_file`** — *GeoDataFrame*: GeoDataFrame with building population data.
>- **`bench_file`** — *GeoDataFrame*: GeoDataFrame with data on available benches.
>- **`heatrisk_file`** — *GeoDataFrame*: GeoDataFrame containing heat risk data.
>- **`pet_file`** — *str*: Path to the PET (Physiological Equivalent Temperature) data file.
>- **`gpkg_file`** — *str*: File path for saving the evaluation results in a GeoPackage format.
>- **`output_layer`** — *str*: Name of the output layer in the GeoPackage.
>- **`search_buffer`** — *int, optional*: Buffer distance in meters for the walking shed calculation; defaults to 700.
>- **`single_day_time_range`** — *list, optional*: Start and end times in HHMM format for a single-day evaluation range.
>- **`time_interval`** — *int, optional*: Time interval in minutes for shade map evaluations.
>- **`search_range`**, **`morning_range`**, **`afternoon_range`**, **`late_afternoon_range`** — *list, optional*: Specific time ranges in HHMM format for different periods within the single-day mode.
>- **`useMultiProcessing`** — *bool, optional*: Enables multiprocessing for shade calculations if set to `True`; defaults to `False`.
>
>**RETURNS**  
>- *GeoDataFrame*: A GeoDataFrame with cool space evaluations, including assessed shade accessibility, resident capacity, street furniture presence, heat risk, and PET scores.

#### <span style="color: red;">`list_to_string`</span><span style="color: gray;">(gdf)</span> <a name="heading--2-1-6"/>

> Convert list values in columns of a GeoDataFrame to strings.
>
>**PARAMETERS**  
>- **`gdf`** — *GeoDataFrame*: The GeoDataFrame in which columns containing lists will have their lists converted to string format.
>
>**RETURNS**  
>- *None*: This function modifies the GeoDataFrame `gdf` in place, converting any list-type elements in its columns to strings.

#### <span style="color: red;">`drop_or_wkt`</span><span style="color: gray;">(gdf, mode='to_wkt')</span> <a name="heading--2-1-7"/>

> Convert or drop geometry columns in a GeoDataFrame.
>
>**PARAMETERS**  
>- **`gdf`** — *GeoDataFrame*: The GeoDataFrame where geometry columns other than the active geometry may be converted or removed.
>- **`mode`** — *str, optional*: Determines action for geometry columns; `'to_wkt'` converts columns to WKT format, while other values drop them. Defaults to `'to_wkt'`.
>
>**RETURNS**  
>- *None*: This function modifies the GeoDataFrame `gdf` in place, converting specified geometry columns to WKT format or dropping them based on the `mode` setting.

#### <span style="color: red;">`output_all_shade_geoms`</span><span style="color: gray;">(gdf, folder_path, output_gpkg)</span> <a name="heading--2-1-8"/>

> Export shade geometries from a GeoDataFrame to a GeoPackage, filtering out non-public areas.
>
>**PARAMETERS**  
>- **`gdf`** — *GeoDataFrame*: Input GeoDataFrame containing shade geometries and area attributes.
>- **`folder_path`** — *str*: Path to the directory where the GeoPackage will be saved.
>- **`output_gpkg`** — *str*: Name of the GeoPackage file (e.g., "shadeGeoms.gpkg") for storing output layers.
>
>**RETURNS**  
>- *None*: This function filters out rows with `typelandge` values of `'overig'` and `'bebouwd gebied'` and exports all shade geometry columns starting with `'sdGeom'` to individual layers in a GeoPackage. The function runs iteratively, using a progress bar to track the output process for each geometry column.

---

### Cool Space Class for Identification (*identification.py*) <a name="heading--2-2"/>

#### **Initialization and Geometry Functions** <a name="heading--2-2-1"/>

##### <span style="color: red;">`__init__`</span><span style="color: gray;">  <a name="heading--2-2-1-1"/>
> Initializes a `CoolSpace` object.
> - **`data`**: GeoDataFrame with spatial data; initializes with a `clipped` column for storing modified geometries.

##### <span style="color: red;">`clip`</span><span style="color: gray;">  <a name="heading--2-2-1-2"/>
> Clips the geometries in `data` using a specified clipping method.
> - **`clipper`**: A GeoDataFrame used to clip `data`.
> - **`how`**: Specifies the clipping method, e.g., `'difference'`, `'intersection'`, etc.
> - **`use_clip`**: If `True`, uses the `clipped` geometry in `data`.
> - **`filter_thin`**: If `True`, applies additional filtering to exclude thin geometries.
>
>**Process Summary**:
> - Validates `how` method and coordinate systems.
> - Computes the clipping operation and optionally explodes multi-polygons.
> - Filters thin geometries by buffering and computes area-to-perimeter ratios to filter small or thin geometries.
> - Updates `self.data` with the filtered geometries in the `clipped` column.

#### **Shade Calculation Functions** <a name="heading--2-2-2"/>

##### <span style="color: red;">`calculate_shade`</span><span style="color: gray;"> <a name="heading--2-2-2-1"/> 
> Calculates shade metrics (average shade, shade area, and geometry) for each raster.
> - **`rasters`**: List of rasters representing shade data.
> - **`area_thres`**: Minimum area threshold for shade regions.
> - **`shade_thres`**: Maximum shade value for valid shaded regions.
> - **`ratio_thres`**: Minimum area-to-perimeter ratio for shade polygons.
> - **`use_clip`**: If `True`, uses `clipped` geometry.
>
>**Process Summary**:
> - For each raster, checks if geometries intersect the raster bounds.
> - Clips geometries to the raster, extracts shaded regions, and filters based on area and shape.
> - Computes average shade values, areas, and geometries, storing them as new columns in `data`.
> - Updates the `intervals` attribute with the number of rasters processed.

##### <span style="color: red;">`calculate_shade_for_raster`</span><span style="color: gray;">  <a name="heading--2-2-2-2"/>
> Processes a single raster and calculates shade metrics.
> - Supports multi-processing by handling a single raster in isolation.
> - **`raster_idx`**, **`raster`**, **`area_thres`**, **`shade_thres`**, **`ratio_thres`**: Similar to `calculate_shade`.
>
>**Process Summary**:
> - Intersects each geometry with raster bounds, extracts valid shaded regions.
> - Aggregates shade metrics and returns them for integration into the main dataset.

##### <span style="color: red;">`calculate_shade_multi`</span><span style="color: gray;">  <a name="heading--2-2-2-3"/>
> Multi-processing version of `calculate_shade`.
> - **`rasters`**, **`area_thres`**, **`shade_thres`**, **`ratio_thres`**, **`use_clip`**: Similar to `calculate_shade`.
>
>**Process Summary**:
> - Initializes parallel processes to handle rasters individually.
> - Each process uses `calculate_shade_for_raster`, collects the results, and integrates them into `data`.
> - Sets `intervals` to the number of processed rasters.

#### **Shade Geometry Retrieval Functions** <a name="heading--2-2-3"/>

##### <span style="color: red;">`get_shade_geometries`</span><span style="color: gray;">  <a name="heading--2-2-3-1"/>
> Retrieves specific shade geometries for a given raster.
> - **`raster_idx`**: Index of the raster to retrieve shade geometries from.
>
>**Process Summary**:
> - Retrieves shade geometry, area, and average shade values for each geometry in the specified raster, returning them as a GeoDataFrame.

##### <span style="color: red;">`get_cool_spaces`</span><span style="color: gray;">  <a name="heading--2-2-3-2"/>
> Retrieves cool spaces that contain shade geometries within a specified range.
> - **`start`**: Starting raster index.
> - **`end`**: Ending raster index.
> - **`geom_type`**: Geometry type to output (`'geometry'` or `'clipped'`).
>
>**Process Summary**:
> - Counts the number of shade geometries for each geometry within the specified range.
> - Filters based on the `count` of valid shade geometries, setting `clipped` or `geometry` as output based on `geom_type`.

#### **Evaluation Functions** <a name="heading--2-2-4"/>

##### <span style="color: red;">`evaluate_shade_coverage`</span><span style="color: gray;">  <a name="heading--2-2-4-1"/>
> Evaluates shade coverage over a specific time range.
> - **`attri_name`**: Label for the output attributes.
> - **`start`** and **`end`**: Raster indices defining the evaluation time range.
> - **`geom_type`**: Geometry type to use (`'geometry'` or `'clipped'`).
>
>**Process Summary**:
> - Aggregates shade scores across rasters in the specified range.
> - Classifies shade coverage into four categories: **0** (<50%), **1** (50%-70%), **2** (70%-90%), **3** (90%-100%).
> - Adds new columns to `data` to store classified shade coverage based on time and area.

---

### Cool Space Class for Evaluation (*evaluation.py*) <a name="heading--2-3"/>

#### **Initialization Functions** <a name="heading--2-3-1"/>

##### <span style="color: red;">`__init__`</span><span style="color: gray;"> <a name="heading--2-3-1-1"/>
> Initializes a `CoolEval` object with data about cool places, nearby buildings, benches, heat risk, and PET (Physiological Equivalent Temperature) values.
> - **Parameters**:
>   - `cool_places`: GeoDataFrame containing polygons of cool places.
>   - `buildings`: GeoDataFrame with building data, including attributes like resident population.
>   - `bench`: GeoDataFrame with bench locations.
>   - `heatrisk`: GeoDataFrame containing heat risk information.
>   - `pet`: Path to the raster file for PET data.
>   - `search_buffer`: Distance buffer used to determine building proximity to cool places.
> - **Attributes**:
>   - `eval_shades`: List for storing evaluation results from each shade geometry.

#### **Walking Shed Functions** <a name="heading--2-3-2"/>

##### <span style="color: red;">`calculate_walking_shed_origin`</span><span style="color: gray;"> <a name="heading--2-3-2-1"/>
> Calculates the walking shed by assigning each building to the nearest cool place within a specified buffer distance.
> - **Process**:
>   - Projects geometries to a suitable coordinate system.
>   - Buffers each cool place, finds intersecting buildings, and calculates the nearest cool place for each building.
>   - Stores the nearest cool place ID and distance in the `buildings` GeoDataFrame.

##### <span style="color: red;">`calculate_walking_shed`</span><span style="color: gray;"> <a name="heading--2-3-2-2"/>
> Similar to `calculate_walking_shed_origin`, but processes buildings and cool places in batches for efficiency.
> - **Process**:
>   - Divides `cool_places` into batches, applies the same distance calculation, and stores results for each building.
>   - Returns the `buildings` GeoDataFrame with assigned cool place IDs.

##### <span style="color: red;">`calculate_walking_shed_multi`</span><span style="color: gray;"> <a name="heading--2-3-2-3"/>
> Multi-processing version of `calculate_walking_shed`, designed for parallel execution.
> - **Process**:
>   - Uses multiple processes to compute the nearest cool places for each batch of buildings.
>   - Aggregates results to create the final `buildings` dataset with walking shed assignments.

#### **3. Resident and Capacity Evaluation Functions** <a name="heading--2-3-3"/>

##### <span style="color: red;">`evaluate_resident`</span><span style="color: gray;"> <a name="heading--2-3-3-1"/>
> Aggregates the number of residents, elderly residents, and children within a specified distance to each cool place.
> - **Process**:
>   - Groups building attributes (residents, elderly, kids) by cool place ID and joins them to the `cool_places` GeoDataFrame.
>   - Returns the updated `cool_places` with aggregated resident information.

##### <span style="color: red;">`evaluate_capacity`</span><span style="color: gray;"> <a name="heading--2-3-3-2"/>
> Evaluates the capacity of shaded areas based on area and nearby residents.
> - **Process**:
>   - Computes the area of each shaded polygon, calculates the capacity (e.g., people per square meter), and assigns capacity status based on nearby residents.
>   - Adds columns for capacity attributes and returns the updated shade GeoDataFrame.

#### **4. Shade Furniture and Heat Risk Evaluation Functions** <a name="heading--2-3-4"/>

##### <span style="color: red;">`evaluate_sfurniture`</span><span style="color: gray;"> <a name="heading--2-3-4-1"/>
> Checks for benches within each shaded area.
> - **Process**:
>   - Uses a spatial join to count benches in each shaded polygon, assigning counts to each shade area.
>   - Adds bench availability information to the shade GeoDataFrame.

##### <span style="color: red;">`evaluate_heatrisk`</span><span style="color: gray;"> <a name="heading--2-3-4-2"/>
> Calculates and classifies the heat risk within shaded areas.
> - **Process**:
>   - Joins `heatrisk` data to shaded areas based on intersections, calculates average heat risk, and classifies risk levels.
>   - Adds heat risk scores and classifications to the shade GeoDataFrame.

#### **PET Evaluation Functions** <a name="heading--2-3-5"/>

##### <span style="color: red;">`eval_pet`</span><span style="color: gray;"> <a name="heading--2-3-5-1"/>
> Computes the average PET (Physiological Equivalent Temperature) values for shaded areas and assigns recommendations.
> - **Process**:
>   - Divides shaded areas into chunks, uses zonal statistics to calculate mean PET for each polygon, and applies classifications.
>   - Adds PET values and recommendations to the shade GeoDataFrame.

##### <span style="color: red;">`eval_pet_multi`</span><span style="color: gray;"> <a name="heading--2-3-5-2"/>
> Multi-processing version of `eval_pet`, leveraging parallel computation for large datasets.
> - **Process**:
>   - Splits shade data across available processors, computes PET values, and aggregates results.
>   - Adds PET values and recommendations in parallel, enhancing performance for large datasets.

#### **Aggregation and Recommendation Functions** <a name="heading--2-3-6"/>

##### **<span style="color: red;">`aggregate_to_cool_places`</span><span style="color: gray;">** <a name="heading--2-3-6-1"/>
> Aggregates evaluation results from shaded areas back to cool places.
> - **Process**:
>   - Aggregates key attributes like capacity, benches, and PET across all shade layers and calculates average values.
>   - Joins these averages to the `cool_places` GeoDataFrame for a comprehensive summary.

##### <span style="color: red;">`final_recom`</span><span style="color: gray;"> <a name="heading--2-3-6-2"/>
> Calculates a recommendation score for each cool place based on capacity, benches, heat risk, PET, and shade metrics.
> - **Process**:
>   - Normalizes features, assigns weights, and classifies final scores into "Not recommended," "Recommended," or "Highly Recommended."
>   - Adds final recommendation classifications to the `cool_places` dataset.

#### **7. Export Functions** <a name="heading--2-3-7"/>

##### <span style="color: red;">`export_eval_gpkg`</span><span style="color: gray;"> <a name="heading--2-3-7-1"/>
> Exports the final evaluation of cool places to a GeoPackage.
> - **Process**:
>   - Converts geometries to WKT format if necessary, and saves the evaluated `cool_places` data as a new layer in a specified GeoPackage file.

---

### Building Class (*building.py*) <a name="heading--2-4">

#### <span style="color: red;">`__init__`</span><span style="color: gray;"> <a name="heading--2-4-1"/>
> Initializes a `Building` object with building geometry data.
> - **Parameters**: `data` - A GeoDataFrame containing building geometries.
> - **Attributes**: Adds a `"buffered"` column to store buffered geometries.

#### <span style="color: red;">`create_buffer`</span><span style="color: gray;"> <a name="heading--2-4-2"/>
> Creates buffer geometries around each building.
> - **Parameters**: `buffer_size` - The buffer distance in meters.
> - **Process**: Applies the buffer around each building's geometry and stores it in the `"buffered"` column, confirming creation status with a message.

---

### Road Class (*road_process.py*) <a name="heading--2-5">

**Note**: The <span style="color: red;">`assign_buffer_hardSurface`</span><span style="color: gray;">, <span style="color: red;">`assign_buffer_roadtype`</span><span style="color: gray;">, and <span style="color: red;">`create_attribute`</span><span style="color: gray;"> methods are designed for a specific dataset and are used primarily during development. In practical applications, the input road dataset **must** include an attribute specifying the buffer distance for different road types, allowing direct use of the <span style="color: red;">`create_buffer`</span><span style="color: gray;"> method.

#### <span style="color: red;">`__init__`</span><span style="color: gray;"> <a name="heading--2-5-1"/>
> Initializes a `Road` object with road geometry data.
> - **Parameters**: `data` - A GeoDataFrame containing road geometries.
> - **Attributes**: Adds a `"buffered"` column to store buffered geometries.

#### <span style="color: red;">`assign_buffer_hardSurface`</span><span style="color: gray;"> <a name="heading--2-5-2"/>
> Assigns buffer size based on the width of hard surfaces.
> - **Parameters**: `roadtype` - A string specifying road width range (e.g., "> 7 meter").
> - **Returns**: Buffer size in meters.

#### <span style="color: red;">`assign_buffer_roadtype`</span><span style="color: gray;"> <a name="heading--2-5-3"/>
> Assigns buffer size based on road type.
> - **Parameters**: `roadtype` - A string specifying the type of road (e.g., "autosnelweg").
> - **Returns**: Buffer size in meters.

#### <span style="color: red;">`create_attribute`</span><span style="color: gray;"> <a name="heading--2-5-4"/>
> Creates a new attribute for buffer size based on road surface type or road type.
> - **Parameters**:
>   - `attri_in` - The existing attribute name (e.g., "verharding").
>   - `new_attri` - The name for the new attribute to store buffer sizes.

#### <span style="color: red;">`create_buffer`</span><span style="color: gray;"> <a name="heading--2-5-5"/>
> Creates buffer geometries for roads based on a specified buffer attribute.
> - **Parameters**: `buffer_attri` - The column name holding buffer sizes.
> - **Process**: Buffers each road geometry according to the buffer size in `buffer_attri` and stores the result in `"buffered"`. 

