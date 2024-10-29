# Shade Calculation Functions

### functions.py
### <span style="color: red;">`raster_center_coords(min_x, max_x, min_y, max_y, resolution)`</span>

Compute the center xy coordinates of a grid.
**PARAMETERS**  
- **`min_x`**, **`max_x`** — *float*: Minimum and maximum x coordinates of the grid.
- **`min_y`**, **`max_y`** — *float*: Minimum and maximum y coordinates of the grid.
- **`resolution`** — *float*: The length of each cell; function only works for square cells.

**RETURNS**  
- *numpy.ndarray*: `grid_center_x` - A grid where each cell contains the center x-coordinate.
- *numpy.ndarray*: `grid_center_y` - A grid where each cell contains the center y-coordinate.

**EXAMPLE**  
```python
>>> grid_center_x, grid_center_y = raster_center_coords(0, 100, 0, 100, 10)
