## Content

**[1. Shade Calculation Functions](#heading--1)**

  * [1.1. functions.py](#heading--1-1)
  * [1.2. Markdown formatting details](#heading--1-2)

**[2. Cool Spaces Functions](#heading--2)**

  * [2.1. Basic text formatting](#heading--2-1)


**[3. Network Functions](#heading--3)**

  * [2.1. Basic text formatting](#heading--2-1)


----


## Shade Calculation Functions <a name="heading--1"/>

### functions.py
#### <span style="color: red;">`raster_center_coords`</span><span style="color: lightgray;">(min_x, max_x, min_y, max_y, resolution)</span>

<ul>
Compute the center xy coordinates of a grid.

**PARAMETERS**  
- **`min_x`**, **`max_x`** — *float*: Minimum and maximum x coordinates of the grid.
- **`min_y`**, **`max_y`** — *float*: Minimum and maximum y coordinates of the grid.
- **`resolution`** — *float*: The length of each cell; function only works for square cells.

**RETURNS**  
- *numpy.ndarray*: `grid_center_x` - A grid where each cell contains the center x-coordinate.
- *numpy.ndarray*: `grid_center_y` - A grid where each cell contains the center y-coordinate.
</ul>

#### <span style="color: red;">`raster_center_coords`</span><span style="color: lightgray;">(min_x, max_x, min_y, max_y, resolution)</span>

<ul>
Compute the center xy coordinates of a grid.

**PARAMETERS**  
- **`min_x`**, **`max_x`** — *float*: Minimum and maximum x coordinates of the grid.
- **`min_y`**, **`max_y`** — *float*: Minimum and maximum y coordinates of the grid.
- **`resolution`** — *float*: The length of each cell; function only works for square cells.

**RETURNS**  
- *numpy.ndarray*: `grid_center_x` - A grid where each cell contains the center x-coordinate.
- *numpy.ndarray*: `grid_center_y` - A grid where each cell contains the center y-coordinate.
</ul>



## Cool Spaces Functions <a name="heading--2"/>