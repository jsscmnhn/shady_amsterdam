import geopandas as gpd
import rasterio
from rasterio.mask import mask
from scipy.ndimage import label
import numpy as np


class CoolSpace:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.clipped = None
        self.shade_avg = None
        self.shadeAreas = []

    def clip(self, to_clip: gpd.geodataframe, clipper: gpd.geodataframe, how='difference') -> None:
        if how not in ['difference', 'intersection', 'union', 'symmetric_difference']:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Choose from "
                f"'difference', 'intersection', 'union', 'symmetric_difference'.")

        self.clipped = gpd.overlay(to_clip, clipper, how=how)

    def calculate_shade(self, raster: rasterio.io.DatasetReader, area_thres=200, use_clip=False) -> None:

        geometries = self.clipped.geometry if use_clip and self.clipped is not None else self.data.geometry
        self.shadeAreas = []  # Reset shade areas for each clipped geometry
        self.shade_avg = None

        for geom in geometries:
            geom_geojson = [geom.__geo_interface__]  # Convert single geometry to geojson format
            out_image, out_transform = mask(raster, geom_geojson, crop=True)
            out_image = out_image[0]  # Assume shade value is stored in the first band

            # Calculate average shading value for this clipped area
            valid_data = out_image[out_image >= 0]
            avg = valid_data.mean() if valid_data.size > 0 else 0

            # Calculate continuous shaded area which are larger than the area threshold, by default is 200m2
            labeled_array, num_features = label(out_image >= 0)
            pixel_size = out_transform[0] * (-out_transform[4])  # x_length * y_length for each pixel
            pixel_areas = []
            for region_label in range(1, num_features + 1):
                region_area = np.sum(labeled_array == region_label) * pixel_size
                if region_area >= area_thres:
                    pixel_areas.append(region_area)

            # Store the shade areas and average for each geometry
            self.shadeAreas.append({
                'geom': geom,
                'shade_areas': pixel_areas,
                'avg_shade': avg
            })

    def is_all_qualified(self, ratio: float) -> bool:
        if not self.shadeAreas:
            return False

        for entry in self.shadeAreas:
            total_area = entry['geom'].area  # Get total area of the current clipped geometry
            if not entry['shade_areas']:
                return False  # If there are no shade areas, return False

            max_shade_area = max(entry['shade_areas'])  # Get the largest shade area
            if max_shade_area / total_area < ratio:
                return False  # If the largest shade area is less than 50% of total area, return False

        return True  # If all areas pass the 50% test, return True

    def is_qualified(self, index: int, ratio: float) -> bool:
        # check whether shadeAreas has been calculated or the index is legit
        if len(self.shadeAreas) == 0:
            raise IndexError(f"Haven't calculate shade areas yet.")
        elif index >= len(self.shadeAreas) or index < 0:
            raise IndexError(f"Invalid index: {index}. Must be between 0 and {len(self.shadeAreas) - 1}.")

        # get current geom's area and its shaded area
        entry = self.shadeAreas[index]
        total_area = entry['geom'].area
        if not entry['shade_areas']:
            return False

        max_shade_area = max(entry['shade_areas'])
        return max_shade_area / total_area >= ratio


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"
    filename = "ams_landuse_top10NL.shp"
    filepath_mac = directory_mac + filename
    filepath_win = directory_win + filename

    coolSpace = CoolSpace(gpd.read_file(filepath_win))
    print(coolSpace.data.head(5))

