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
     * [median_filter_chm`](#heading--1-3-1)
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

  * [2.1. TODO](#heading--2-1)


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
> - (nothing)
 
### Data Acquisition (*collect_data.py*) <a name="heading--1-2"/>

#### <span style="color: red;">`download_las_tiles`</span><span style="color: gray;">(tile_list_file, output_folder)</span> <a name="heading--1-2-1"/>

> Download LAZ files for each subtile specified in a text file.
>
> **PARAMETERS**  
> - **`tile_list_file`** — *str*: Path to the text file containing the list of subtiles to download.
> - **`output_folder`** — *str*: Path to the folder where the downloaded files should be saved.
>
> **RETURNS**  
> - (none): The function writes output files directly to the specified `output_folder/{subtile}`.

#### <span style="color: red;">`download_and_extract`</span><span style="color: gray;">(url, file_path, output_folder)</span> <a name="heading--1-2-2"/>

> Download a ZIP file from a URL, extract its contents, and delete the ZIP file.
>
> **PARAMETERS**  
> - **`url`** — *str*: URL of the ZIP file to download.
> - **`file_path`** — *str*: Path to save the downloaded ZIP file.
> - **`output_folder`** — *str*: Path to the folder where the extracted contents should be saved.
>
> **RETURNS**  
> - (none)

#### <span style="color: red;">`merge_tif_files`</span><span style="color: gray;">(input_folder, output_file, file_prefix, nodata_value=-9999)</span> <a name="heading--1-2-3"/>

> Merge TIF files in a folder with a specific prefix into a single raster file.
>
> **PARAMETERS**  
> - **`input_folder`** — *str*: Path to the folder containing the TIF files.
> - **`output_file`** — *str*: Path to the output raster file.
> - **`file_prefix`** — *str*: Prefix of files to merge ("M_" for DTM or "R_" for DSM).
> - **`nodata_value`** — *int*, *optional*: Value to replace the nodata value. Defaults to -9999.
>
> **RETURNS**  
> - *list*: A list containing the name of the tile and the extent of that tile, for downloading the building data.
> - (none): The function writes the merged file directly to the specified `output_file`.

#### <span style="color: red;">`download_raster_tiles`</span><span style="color: gray;">(tile_list_file, output_folder, name)</span> <a name="heading--1-2-4"/>

> Download DSM and DTM files for each tile specified in a text file.
>
> **PARAMETERS**  
> - **`tile_list_file`** — *str*: Path to the text file containing the list of tiles to download.
> - **`output_folder`** — *str*: Directory where the downloaded and unzipped files will be saved.
> - **`name`** — *str*: Name of the output raster file.
>
> **RETURNS**  
> - *list*: A list containing the name of the tile and the extent of that tile, for downloading the building data.
> - (none): The function writes output files directly to the specified `output_folder`.

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

### CHM Creation (*create_first_chm_parallel.py*) <a name="heading--1-3"/>

#### <span style="color: red;">`median_filter_chm`</span><span style="color: gray;">(chm_array, nodata_value=-9999, size=3)</span> <a name="heading--1-3-1"/>

> Apply a median filter to a CHM, handling NoData values.
>
> **PARAMETERS**  
> - **`chm_array`** — *np.ndarray*: The array representing the height values of the CHM.
> - **`nodata_value`** — *float*: Value representing NoData in the input raster.
> - **`size`** — *int*: Size of the median filter. It defines the footprint of the filter.
>
> **RETURNS**  
> - *np.ndarray*: The smoothed CHM array.

---

#### <span style="color: red;">`extract_vegetation_points`</span><span style="color: gray;">(LasData, ndvi_threshold=0.1, pre_filter=False)</span> <a name="heading--1-3-2"/>

> Extract vegetation points based on classification and NDVI threshold.
>
> **PARAMETERS**  
> - **`LasData`** — *laspy.LasData*: Input point cloud data in LAS format.
> - **`ndvi_threshold`** — *float*: The NDVI threshold for identifying vegetation points. NDVI values greater than this threshold are considered vegetation.
> - **`pre_filter`** — *bool*: If True, applies an additional filter to remove vegetation points below a certain height threshold (1.5 meters above the lowest vegetation point).
>
> **RETURNS**  
> - *laspy.LasData*: A new LasData object containing only the filtered vegetation points based on the specified criteria.

---

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
> - *(none)*: The function saves the CHM as a raster file (.tif) to the specified output path.

---

#### <span style="color: red;">`interpolation_vegetation`</span><span style="color: gray;">(LasData, veg_points, resolution, no_data_value=-9999)</span> <a name="heading--1-3-4"/>

> Create a vegetation raster using Laplace interpolation.
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

---

#### <span style="color: red;">`process_single_laz_file`</span><span style="color: gray;">(file_path, output_folder, ndvi_threshold=0.0, resolution=0.5, remove=False, smooth_chm=False, filter_size=3, pre_filter=False)</span> <a name="heading--1-3-5"/>

> Process a LAZ file to extract vegetation points and generate a CHM.
>
> **PARAMETERS**  
> - **`file_path`** — *str*: The file_path containing the folders containing the input .LAZ file.
> - **`output_folder`** — *str*: The folder where the output CHM .tif files will be saved.
> - **`ndvi_threshold`** — *float*: The NDVI threshold for classifying vegetation points.
> - **`resolution`** — *float*: The resolution of the output CHM rasters, defining the size of each pixel (default: 0.5).
> - **`remove`** — *bool*: If True, deletes the original .LAZ files after processing (default: False).
> - **`smooth_chm`** — *bool*: If True, applies smoothing to the CHM using a median filter (default: False).
> - **`filter_size`** — *int*: Size of the median filter to use if smoothing (default: 3).
> - **`pre_filter`** — *bool*: If True, applies an additional filter to remove vegetation points below a certain height threshold (1.5 meters above the lowest vegetation point).
>
> **RETURNS**  
> - *(none)*: The function processes a .LAZ file, creates a corresponding CHM .tif file, and saves it to the output folder. Optionally deletes the original .LAZ files if `remove` is set to True.

---

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
> - *(none)*: The function processes the .LAZ files, creates corresponding CHM .tif files, and saves them to the output folder. Optionally deletes the original .LAZ files if `remove` is set to True.

### DSM and Final CHM Creation (*create_dsm_and_chm_parallel.py*) <a name="heading--1-4"/>

#### <span style="color: red;">`get_bbox`</span><span style="color: gray;">(raster_paths)</span> <a name="heading--1-4-1"/>

> Compute the overlapping bounding box of the CHM, DSM, and DTM file.
>
> **Input:**  
> - `raster_paths` (list of strings): Paths to the CHM, DSM, and DTM .TIF files.
>
> **Output:**  
> - *rasterio.coords.BoundingBox*: Bounding box of the raster.

---

#### <span style="color: red;">`crop_raster`</span><span style="color: gray;">(raster_path, bbox, no_data=-9999)</span> <a name="heading--1-4-2"/>

> Crop the input rasters to the size of the overlapping bounding box.
>
> **Input:**  
> - `raster_path` (string): Path to the .TIF file.  
> - `bbox` (4-tuple): Overlapping bounding box of the CHM, DSM, and DTM file.  
> - `no_data` (int, optional): NoData value to replace source NoData value with.
>
> **Output:**  
> - `cropped_data` (2D numpy array): Cropped raster data.  
> - `src.window_transform(window)` (affine transform matrix): For the given window.  
> - `src.crs` (rasterio src): A PROJ4 dict representation of the CRS of the input raster.

---

#### <span style="color: red;">`extract_center_cells`</span><span style="color: gray;">(cropped_data, no_data=-9999)</span> <a name="heading--1-4-3"/>

> Extract the values of each cell in the input data and save these with the x and y (row and col) indices.
>
> **Input:**  
> - `cropped_data` (2D numpy array): Cropped raster data.  
> - `no_data` (int, optional): NoData value to replace source NoData value with.
>
> **Output:**  
> - `xyz_filled` (list): List containing x, y, and z coordinates of the cells.

---

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
> - `new_data[0, 1:-1, 1:-1]` (2D numpy array): Filled raster data with first and last rows and columns removed to ensure there are no NoData values from Laplace interpolation.  
> - `new_transform` (rasterio transform): Affine transform matrix reflecting the one-column one-row removal shift.

---

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
> - `result_array` (2D numpy array): Array of the CHM with normalized height and min and max heights removed.  
> - `new_transform` (rasterio transform): Affine transform matrix reflecting the one-column one-row removal shift.

---

#### <span style="color: red;">`replace_buildings`</span><span style="color: gray;">(filled_dtm, dsm_buildings, buildings_geometries, transform)</span> <a name="heading--1-4-6"/>

> Replace the values of the filled DTM with the values of the filled DSM, if there is a building.
>
> **Input:**  
> - `filled_dtm` (2D np array): Filled array of the cropped AHN DTM.  
> - `dsm_buildings` (2D np array): Filled array of the cropped AHN DSM.  
> - `building_geometries` (list): A list of the building geometries.  
> - `transform` (rasterio transform): Affine transform matrix.
>
> **Output:**  
> - `final_dsm` (2D numpy array): A numpy array representing the final DSM, containing only ground and building heights.

---

#### <span style="color: red;">`load_buildings`</span><span style="color: gray;">(buildings_path, layer)</span> <a name="heading--1-4-7"/>

> Load the building shapes from a geopackage file.
>
> **Input:**  
> - `buildings_path` (string): Path to the geopackage file.  
> - `layer` (string): (Tile) name of the layer of buildings to be used.
>
> **Output:**  
> - List of dictionaries: A list of geometries in GeoJSON-like dictionary format. Each dictionary represents a building geometry with its spatial coordinates.

---

#### <span style="color: red;">`extract_tilename`</span><span style="color: gray;">(filename)</span> <a name="heading--1-4-8"/>

> Extract the name of the AHN tile from the file name.
>
> **Input:**  
> - `filename` (string): The name of the input chm.tif.
>
> **Output:**  
> - `match.group(1)`: The name of the AHN tile.

---

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

---

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
> - None: The function writes output files directly to the specified `output_base_folder`.  
> - Output files: For each input CHM file, the function generates:  
>   - A DSM file with building data incorporated: `DSM_<tile>_<subtile_number>.tif`  
>   - A processed CHM file: `CHM_<tile>_<subtile_number>.tif`  
>   - These files are saved in folders named after the tile in `output_base_folder`.

---

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
> - None: The function writes output files directly to the specified `output_base_folder`.

### Creating Shade Maps (*shade_parallel.py*) <a name="heading--1-5"/>

#### <span style="color: red;">`process_chm_dsm`</span><span style="color: gray;">(chm_filename, dsm_filename, folder_path, output_dir, date, start_time=10, end_time=21, interval=30, use_chm=True, trans=10, trunkheight=25)</span> <a name="heading--1-5-1"/>

> Function to process a single DSM and CHM file pair.
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
> - None: This function prints a message upon completion of shade calculation.

---

#### <span style="color: red;">`run_shade_calculation`</span><span style="color: gray;">(file_pairs, output_base_folder, date, start_time, end_time, interval, use_chm=True, trans=10, trunkheight=25, max_workers=4)</span> <a name="heading--1-5-2"/>

> Process CHM and DSM file pairs in parallel.
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
> - None: This function prints messages regarding processing status and creates output directories.

---

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
> - None: This function prints messages regarding the processing status for each folder.

### Merging Shade Maps (*merge_shademaps.py*) <a name="heading--1-6"/>

#### <span style="color: red;">`merge_tif_files_by_time`</span><span style="color: gray;">(main_folder, output_folder, merged_name, start_time=900, end_time=2000, delete_input_files=False, nodata_value=-9999)</span> <a name="heading--1-6-1"/>

> Merge multiple TIF files from specified subfolders into a single mosaic file for each specified time step.
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
> - None: The function saves merged TIFF files directly to the specified `output_folder` with the naming format `{merged_name}_{date_str}_{time}.TIF`.

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
> - `items`: A flat dictionary.

---

#### <span style="color: red;">`read_config`</span><span style="color: gray;">(file_path)</span> <a name="heading--1-7-2"/>

> Read the contents of a JSON configuration file and save them to a dictionary.
>
> **Input:**  
> - `file_path` (str): The path to the configuration file.
>
> **Output:**  
> - `params` (dict): The configuration dictionary containing all required parameters.
## Cool Spaces Functions <a name="heading--2"/>