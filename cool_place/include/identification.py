import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from scipy.ndimage import label
from shapely.geometry import shape, box
from shapely.ops import unary_union
import numpy as np
import matplotlib.pyplot as plt

pd.set_option('future.no_silent_downcasting', True)


class CoolSpace:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.data["clipped"] = None
        self.intervals = 0

    def clip(self, clipper: gpd.geodataframe, how='difference', use_clip=False) -> None:
        if how not in ['difference', 'intersection', 'union', 'symmetric_difference']:
            raise ValueError(
                f"Invalid 'how' parameter: {how}. Choose from "
                f"'difference', 'intersection', 'union', 'symmetric_difference'."
            )
        if clipper.crs is None:
            clipper = clipper.set_crs(self.data.crs)
        if self.data.crs != clipper.crs:
            clipper = clipper.to_crs(self.data.crs)
        if use_clip and 'clipped' in self.data.columns:
            self.data.set_geometry('clipped', inplace=True)
        else:
            self.data.set_geometry('geometry', inplace=True)

        self.data['orig_id'] = self.data.index

        clipped = gpd.overlay(self.data, clipper, how=how, keep_geom_type=True)

        # Explode multipolygons into individual polygons
        exploded = clipped.explode(index_parts=False).reset_index(drop=True)

        # Calculate area to perimeter ratio
        exploded['area'] = exploded.geometry.area
        exploded['perimeter'] = exploded.geometry.length
        exploded['area_to_perimeter'] = exploded['area'] / (exploded['perimeter'] + 0.0001)

        # Filter polygons based on area to perimeter ratio and area
        filtered = exploded[(exploded['area'] >= 200) & (exploded['area_to_perimeter'] >= 0.35)]

        # Group by original geometry and merge polygons back to multipolygon if needed
        filtered = filtered.dissolve(by='orig_id', as_index=False)

        # Set back to original index
        filtered.set_index('orig_id', inplace=True)

        # Update self.data with the filtered multipolygon geometries
        self.data['clipped'] = filtered.geometry.reindex(self.data.index)
        self.data.drop(columns='orig_id', inplace=True)

    def calculate_shade(self,
                        rasters: list[rasterio.io.DatasetReader],
                        area_thres=200,
                        shade_thres=0.5,
                        ratio_thres=0.35,
                        use_clip=False) -> None:
        """
        input a series of rasters, use previous geometries to clip the rasters, and calculate each
        raster's average shade, areas of shade, and shade geometries.

        - The average shade value is calculated as the average of all pixels with 0 <=shade value <= 0.5
          (0 means maximum shade, 1 means sun).
        - The shade area is calculated as the area of all continuous areas larger than 200m2.
        - The shade geometry is calculated as the union of all shade polygons with area-to-perimeter ratio >= 0.35.

        - The average shade value will be added as a new column to "self.data" as "sdAvg{raster_idx}"
        - The shade area will be added as a new column to "self.data" as "sdArea{raster_idx}"
        - The shade geometry will be added as a new column to "self.data" as "sdGeom{raster_idx}"

        - An extra column "intervals" will be added to self (NOT "self.data"!) to indicate the number of rasters
          used for shade calculation.

        :param rasters: input list of rasters
        :param area_thres: minimum threshold of shade area, by default is 200m2
        :param shade_thres: maximum threshold of shade value, by default is 0.5
        :param ratio_thres: minimum threshold of area-to-perimeter ratio, by default is 0.35
        :param use_clip: use clipped geometry or not, if there is no clipped geometry, original will be used
        """

        if use_clip and 'clipped' in self.data.columns and self.data['clipped'].any():
            self.data.set_geometry('clipped', inplace=True)
        else:
            print("No clipped geometry, default geometry will be used.")
            self.data.set_geometry('geometry', inplace=True)

        # same calculation for each raster
        for raster_idx, raster in enumerate(rasters):
            print(f"Processing raster {raster_idx + 1}/{len(rasters)}")

            # initialize attributes for current raster
            self.data[f"sdAvg{raster_idx}"] = None
            self.data[f"sdArea{raster_idx}"] = None
            self.data[f"sdGeom{raster_idx}"] = None

            minx, miny, maxx, maxy = raster.bounds
            raster_bounds = box(minx, miny, maxx, maxy)
            print(raster_bounds)

            # create GeoSeries to store all the sdGeom (polygon transformed from pixels)
            all_shade_geoms = gpd.GeoSeries(index=self.data.index)

            for idx, row in self.data.iterrows():
                geom = row[self.data.geometry.name]
                if geom is None:
                    # print(f"Geometry {idx} is None, skipping.")
                    all_shade_geoms.at[idx] = None
                    continue

                # check if geometry intersects the raster or not
                if geom.intersects(raster_bounds):
                    clipped_geom = geom.intersection(raster_bounds)

                    # if intersection fail (empty) or the clipped result is too small, skip it
                    if clipped_geom.is_empty or clipped_geom.area < 1e-6:
                        print(f"Geometry {idx} is too small ({clipped_geom.area}) "
                              f"or empty after intersection, skipping.")
                        all_shade_geoms.at[idx] = None
                        continue
                else:
                    continue  # if geometry not intersects with raster, skip it

                geom_geojson = [clipped_geom.__geo_interface__]  # transform geom to GeoJSON
                out_image, out_transform = mask(raster, geom_geojson, crop=True)
                out_image = out_image[0]  # assume shade value is stored in band 1

                # for all the pixels have shade value <= 0.5 (0 means maximum shade, 1 means sun),
                # calculate the continuous area
                labeled_array, num_features = label((out_image >= 0) & (out_image <= shade_thres))
                pixel_size = out_transform[0] * (-out_transform[4])  # area of a pixel
                pixel_areas = []
                shade_polygons = []

                # only select those continuous area larger than 200m2
                for region_label in range(1, num_features + 1):
                    region_area = np.sum(labeled_array == region_label) * pixel_size
                    if region_area >= area_thres:
                        pixel_areas.append(np.round(region_area, 4))

                        # create region mask (bool)
                        region_mask = (labeled_array == region_label).astype(np.uint8)

                        # use rasterio.features.shapes to transform the mask to geometry
                        for geom, value in shapes(region_mask, transform=out_transform):
                            if value == 1:  # reserve the "true" part
                                # check the ratio of area to perimeter
                                if shape(geom).area / shape(geom).length >= ratio_thres:
                                    shade_polygons.append(shape(geom))

                if shade_polygons:
                    merged_polygon = unary_union(shade_polygons)
                    all_shade_geoms.at[idx] = merged_polygon
                else:
                    all_shade_geoms.at[idx] = None

                self.data.at[idx, f"sdArea{raster_idx}"] = pixel_areas

                # calculate average shade value
                valid_data = out_image[(out_image >= 0) & (out_image <= shade_thres)]
                if len(pixel_areas) == 0:
                    avg = 1
                else:
                    avg = valid_data.mean() if valid_data.size > 0 else 1

                self.data.at[idx, f"sdAvg{raster_idx}"] = np.float64(avg)

            self.data[f"sdGeom{raster_idx}"] = all_shade_geoms
            self.intervals = len(rasters)

    def get_shade_geometries(self, raster_idx: int) -> gpd.geodataframe:
        """
        get shade geometries for a specific raster

        :param raster_idx: index of raster
        :return: shade geometries for a specific raster
        """
        if self.intervals == 0:
            print("No shade calculation has been done, please run calculate_shade() first.")
            return None
        
        shade_geom_col = f"sdGeom{raster_idx}"
        shade_avg_col = f"sdAvg{raster_idx}"
        shade_area_col = f"sdArea{raster_idx}"
        if self.data[shade_geom_col].isnull().all():
            print(f"No shade geometry for raster {raster_idx} in data.")
            return None

        output = self.data[self.data[shade_geom_col].notnull()][['id',
                                                                 shade_avg_col,
                                                                 shade_area_col,
                                                                 shade_geom_col]].copy()
        output = output.set_geometry(shade_geom_col, crs=self.data.crs)
        return output

    def get_cool_spaces(self, start: int = None, end: int = None, geom_type: str = 'geometry') -> gpd.geodataframe:
        """
        get cool spaces geometries (the geometries that have shade geometries from the rasters within the search range).

        - if a cool space has at least one shade geometry from all searching rasters, it will be return, otherwise it
          will be discarded.

        :param geom_type: by default, the output cool space polygon will be the filtered input land-use polygons.
                    If user wants to output the clipped geometry (public space) as output, set 'type' to 'clipped'
        :param start: the index of the first raster of the search range
        :param end: the index of the last raster of the search range
        :return: cool spaces geometries for a specific raster
        """
        if self.intervals == 0:
            print("No shade calculation has been done, please run calculate_shade() first.")
            return None
        
        raster_nums = self.intervals

        if geom_type not in ['geometry', 'clipped']:
            raise ValueError("Invalid type. Expected 'geometry' or 'clipped'.")
        cool_geom_col = f"{geom_type}"

        self.data["count"] = 0

        if start is not None and end is not None:
            if start < 0 or end > raster_nums - 1:
                print(f"The search range is: {start} - {end}, which exceeds the range of"
                      f"shade maps: 0 - {raster_nums - 1}, the shade maps range will be used.")
                search_range = range(raster_nums)
            else:
                search_range = range(start, end + 1)
        else:
            search_range = range(raster_nums)

        for raster_idx in search_range:
            shade_geom_col = f"sdGeom{raster_idx}"
            self.data["count"] += self.data[shade_geom_col].notnull().astype(int)

        if cool_geom_col == "clipped" and self.data[cool_geom_col].isnull().all():
            print("No clipped geometry in data, the original geometry will be returned.")
            cool_geom_col = "geometry"

        output = self.data[self.data["count"] > 0].copy()
        output.set_geometry(cool_geom_col, inplace=True)
        output.drop(columns=["count"], inplace=True)
        # if cool_geom_col == "geometry" and "clipped" in output.columns:
        #     output["clipped"] = output["clipped"].to_wkt()
        # elif cool_geom_col == "clipped":
        #     output["geometry"] = output["geometry"].to_wkt()
        #
        # # transform all shade geometries into WKT and store as new columns
        # for i in search_range:
        #     output[f"sdGeom{i}"] = output[f"sdGeom{i}"].to_wkt()

        self.data.drop(columns=["count"], inplace=True)
        return output

    def evaluate_shade_coverage(self, attri_name: str = "Query", start: int = None, end: int = None) -> None:
        """
        calculate the shade coverage based on the average shade value of all rasters.

        - The shade coverage is classified into 4 categories: 0 (<50%), 1 (50% - 70%), 2 (70% - 90%), 3 (90% - 100%).
        - The shade coverage will be added as a new column to "self.data" as "sc{attri_name}".
        """
        raster_nums = self.intervals
        if raster_nums == 0:
            print("No shade calculation has been done, please run calculate_shade() first.")
            return None

        if start is not None and end is not None:
            try:
                if start < 0 or end > raster_nums - 1:
                    raise ValueError(f"{attri_name}: The search range is:  "
                                     f"[{start} - {end}], {end + 1 - start} in total, which exceeds the range of"
                                     f"shade maps: [0 - {raster_nums - 1}], {raster_nums} in total")
            except ValueError as e:
                print(f"Error: {e}, this evaluation will be ignored.")
                return None

            search_range = range(start, end + 1)
        else:
            search_range = range(raster_nums)

        self.data["tol_shade_avg"] = 0
        for i in search_range:
            shade_avg_col = f"sdAvg{i}"
            self.data["tol_shade_avg"] += self.data[shade_avg_col].fillna(1)

        self.data["tol_shade_avg"] /= raster_nums
        self.data["tol_shade_avg"] = self.data["tol_shade_avg"].round(4)

        def classify_shade_coverage(avg) -> int:
            if 0 <= avg <= 0.1:
                return 3
            elif avg <= 0.3:
                return 2
            elif avg <= 0.5:
                return 1
            else:
                return 0
        attribute_name = f"sc{attri_name}"
        self.data[attribute_name] = self.data["tol_shade_avg"].apply(classify_shade_coverage).astype(int)
        self.data.drop(columns=["tol_shade_avg"], inplace=True)