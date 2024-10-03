import geopandas as gpd
import rasterio
from rasterio.mask import mask
from scipy.ndimage import label
import numpy as np
import matplotlib.pyplot as plt


class CoolSpace:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.data["clipped"] = None
        self.data["shadeAvg"] = None
        self.data["shadeAreas"] = None
        self.data["shadeGeoms"] = None
        self.data["shade_qualified"] = None
        self.data["shape_qualified"] = None
        self.data["qualiGeoms"] = None

    def clip(self, to_clip: gpd.geodataframe, clipper: gpd.geodataframe, how='difference') -> None:
        if how not in ['difference', 'intersection', 'union', 'symmetric_difference']:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Choose from "
                f"'difference', 'intersection', 'union', 'symmetric_difference'."
            )

        clipped = gpd.overlay(to_clip, clipper, how=how)
        self.data['clipped'] = clipped.geometry

    def calculate_shade(self, raster: rasterio.io.DatasetReader, area_thres=200, use_clip=False) -> None:
        self.data["shadeAvg"] = None
        self.data["shadeAreas"] = None

        if use_clip and 'clipped' in self.data:
            geometries = self.data['clipped']
        else:
            print("No clipped geometry, default geometry will be used.")
            geometries = self.data.geometry

        for idx, geom in geometries.items():
            geom_geojson = [geom.__geo_interface__]  # Convert single geometry to geojson format
            out_image, out_transform = mask(raster, geom_geojson, crop=True)
            out_image = out_image[0]  # Assume shade value is stored in the first band

            # Calculate average shading value for this clipped area (just a test)
            valid_data = out_image[out_image >= 0]
            avg = valid_data.mean() if valid_data.size > 0 else 0

            # Calculate continuous shaded area which are larger than the area threshold
            # This step should put another threshold for each pixel, right now is >= 0
            labeled_array, num_features = label(out_image >= 0)
            pixel_size = out_transform[0] * (-out_transform[4])  # x_length * y_length for each pixel
            pixel_areas = []
            for region_label in range(1, num_features + 1):
                region_area = np.sum(labeled_array == region_label) * pixel_size
                if region_area >= area_thres:
                    pixel_areas.append(region_area)

            # Store the results in the GeoDataFrame columns
            self.data.at[idx, "shadeAvg"]: float = avg
            self.data.at[idx, "shadeAreas"]: list = pixel_areas

    def get_qualified_shaded_geometries(self, min_shade_area=200) -> None:
        """
        get geometries that have continuous shaded areas larger than 200m2
        """
        self.data["shade_qualified"]: bool = self.data["shadeAreas"].apply(
            lambda x: any(area >= min_shade_area for area in x) if x else False
        )

        if self.data["shade_qualified"].any():
            print(f"There are {self.data['shade_qualified'].sum()}/{self.data.shape[0]} shade-qualified geometries.")
        else:
            print("No shade-qualified geometries found.")

    def get_qualified_shape_geometries(self, ratio_thres=0.35, use_clip=False) -> None:
        """
        get geometries that have the ratio of area to perimeter > 0.35
        """
        geometries = None
        if use_clip and 'clipped' in self.data:
            geometries = self.data['clipped']
        else:
            print("No clipped geometry, default geometry will be used.")
            geometries = self.data.geometry

        self.data["clipped_perimeter"] = geometries.length
        self.data["clipped_area"] = geometries.area
        self.data["clipped_perimeter_area_ratio"] = self.data["clipped_area"] / self.data["clipped_perimeter"]
        self.data["shape_qualified"]: bool = self.data["clipped_perimeter_area_ratio"] > ratio_thres

        if self.data["shape_qualified"].any():
            print(f"There are {self.data['shape_qualified'].sum()}/{self.data.shape[0]} shape-qualified geometries.")
        else:
            print("No shape-qualified geometries found.")

    def get_qualified_geometries(self) -> None:
        """
        get geometries that fulfill all the qualification checks.
        """
        if not self.data["shade_qualified"].any() or not self.data["shape_qualified"].any():
            raise ValueError("Please check shade and shape qualification first.")

        self.data["qualiGeoms"] = None
        self.data["qualiGeoms"] = self.data[self.data["shade_qualified"] and self.data["shape_qualified"]]


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"
    filename = "ams_landuse_top10NL.shp"
    filepath_mac = directory_mac + filename
    filepath_win = directory_win + filename

    coolSpace = CoolSpace(gpd.read_file(filepath_win))
    coolSpace.get_qualified_shaded_geometries()

