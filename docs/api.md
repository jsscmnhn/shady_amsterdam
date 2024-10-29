## Content

**[1. Shade Calculation Functions](#heading--1)**

  * [1.1. General Functions (*functions.py*)](#heading--1-1)
  * [1.2. Data Acquisition (*collect_data.py*) ](#heading--1-2)
  * [1.3. CHM Creation (*create_first_chm_parallel.py*) ](#heading--1-3)  
  * [1.4. DSM and final CHM creation (*create_dsm_and_chm_parallel.py*) ](#heading--1-4)
  * [1.5. Creating Shade Maps (*shade_parallel.py*) ](#heading--1-5)
  * [1.6. Merging Shade Maps (*merge_shademaps.py*) ](#heading--1-6)
  * [1.7. Running code (*main_shade.py*) ](#heading--1-7)
	
**[2. Cool Spaces Functions](#heading--2)**

  * [2.1. TODO](#heading--2-1)


**[3. Network Functions](#heading--3)**

  * [3.1. TODO](#heading--3-1)


----


## Shade Calculation Functions <a name="heading--1"/>

### General Functions (*functions.py*) <a name="heading--1-1"/>
#### <span style="color: red;">`raster_center_coords`</span><span style="color: gray;">(min_x, max_x, min_y, max_y, resolution)</span>

<div style="margin-left: 20px;">
Compute the center xy coordinates of all cells in a grid.

**PARAMETERS**  
- **`min_x`**, **`max_x`** — *float*: Minimum and maximum x coordinates of the grid.
- **`min_y`**, **`max_y`** — *float*: Minimum and maximum y coordinates of the grid.
- **`resolution`** — *float*: The length of each cell; function only works for square cells.

**RETURNS**  
- *numpy.ndarray*: `grid_center_x` - A grid where each cell contains the center x-coordinate.
- *numpy.ndarray*: `grid_center_y` - A grid where each cell contains the center y-coordinate.
- 
</div>

#### <span style="color: red;">`raster_center_coords`</span><span style="color: lightgray;">(min_x, max_x, min_y, max_y, resolution)</span>


Compute the center xy coordinates of a grid.

**PARAMETERS**  
- **`min_x`**, **`max_x`** — *float*: Minimum and maximum x coordinates of the grid.
- **`min_y`**, **`max_y`** — *float*: Minimum and maximum y coordinates of the grid.
- **`resolution`** — *float*: The length of each cell; function only works for square cells.

**RETURNS**  
- *numpy.ndarray*: `grid_center_x` - A grid where each cell contains the center x-coordinate.
- *numpy.ndarray*: `grid_center_y` - A grid where each cell contains the center y-coordinate.

### Data Acquisition (*collect_data.py*) <a name="heading--1-2"/>

* Add relevant content here.

### CHM Creation (*create_first_chm_parallel.py*) <a name="heading--1-3"/>

* Add relevant content here.

### DSM and Final CHM Creation (*create_dsm_and_chm_parallel.py*) <a name="heading--1-4"/>

* Add relevant content here.

### Creating Shade Maps (*shade_parallel.py*) <a name="heading--1-5"/>

* Add relevant content here.

### Merging Shade Maps (*merge_shademaps.py*) <a name="heading--1-6"/>

* Add relevant content here.

### Running Code (*main_shade.py*) <a name="heading--1-7"/>

* Add relevant content here.
* 

## Cool Spaces Functions <a name="heading--2"/>