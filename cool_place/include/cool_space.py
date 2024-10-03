import geopandas as gpd
import rasterio
from rasterio.mask import mask
from scipy.ndimage import label
import numpy as np


class CoolSpace:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.data['clipped'] = None
        self.data['shade_avg'] = None
        self.data['shade_areas'] = None

    def clip(self, to_clip: gpd.geodataframe, clipper: gpd.geodataframe, how='difference') -> None:
        if how not in ['difference', 'intersection', 'union', 'symmetric_difference']:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Choose from "
                f"'difference', 'intersection', 'union', 'symmetric_difference'."
            )

        clipped = gpd.overlay(to_clip, clipper, how=how)
        self.data['clipped'] = clipped.geometry

    def calculate_shade(self, raster: rasterio.io.DatasetReader, area_thres=200, use_clip=False) -> None:
        self.data['shade_avg'] = None
        self.data['shade_areas'] = None

        geometries = self.data['clipped'] if use_clip and 'clipped' in self.data else self.data.geometry

        for idx, geom in geometries.items():
            geom_geojson = [geom.__geo_interface__]  # Convert single geometry to geojson format
            out_image, out_transform = mask(raster, geom_geojson, crop=True)
            out_image = out_image[0]  # Assume shade value is stored in the first band

            # Calculate average shading value for this clipped area
            valid_data = out_image[out_image >= 0]
            avg = valid_data.mean() if valid_data.size > 0 else 0

            # Calculate continuous shaded area which are larger than the area threshold
            labeled_array, num_features = label(out_image >= 0)
            pixel_size = out_transform[0] * (-out_transform[4])  # x_length * y_length for each pixel
            pixel_areas = []
            for region_label in range(1, num_features + 1):
                region_area = np.sum(labeled_array == region_label) * pixel_size
                if region_area >= area_thres:
                    pixel_areas.append(region_area)

            # Store the results in the GeoDataFrame columns
            self.data.at[idx, 'shade_avg'] = avg
            self.data.at[idx, 'shade_areas'] = pixel_areas

    def is_all_qualified(self, ratio: float) -> bool:
        if self.data['shade_areas'].isnull().all():
            return False

        for idx, row in self.data.iterrows():
            total_area = row['geometry'].area  # Get total area of the current clipped geometry
            if not row['shade_areas']:
                return False  # If there are no shade areas, return False

            max_shade_area = max(row['shade_areas'])  # Get the largest shade area
            if max_shade_area / total_area < ratio:
                return False  # If the largest shade area is less than the ratio, return False

        return True  # If all areas pass the ratio test, return True

    def is_qualified(self, index: int, ratio: float) -> bool:
        if index >= len(self.data) or index < 0:
            raise IndexError(f"Invalid index: {index}. Must be between 0 and {len(self.data) - 1}.")

        row = self.data.iloc[index]
        total_area = row['geometry'].area
        if not row['shade_areas']:
            return False

        max_shade_area = max(row['shade_areas'])
        return max_shade_area / total_area >= ratio

    def get_qualified_geometries(self, min_shade_area=200) -> gpd.GeoDataFrame:
        """
        return geometries that have shaded areas larger than 200m2
        """
        qualified_geometries = self.data[self.data['shade_areas'].apply(
            lambda x: any(area >= min_shade_area for area in x) if x else False)]
        return qualified_geometries



if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"
    filename = "ams_landuse_top10NL.shp"
    filepath_mac = directory_mac + filename
    filepath_win = directory_win + filename

    coolSpace = CoolSpace(gpd.read_file(filepath_win))
    print(coolSpace.data.head(5))

