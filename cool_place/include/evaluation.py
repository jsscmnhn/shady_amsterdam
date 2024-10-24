import geopandas as gpd
import jenkspy
from rasterstats import zonal_stats
import pandas as pd
import rasterio
from shapely.geometry import shape, box
import numpy as np
import fiona
import os
from shapely.geometry import base
from shapely import wkt
import re
import osmnx as ox
from shapely.geometry import Point, shape
import time
from rich.progress import Progress


class CoolEval:
    def __init__(self, cool_places: gpd.GeoDataFrame, buildings: gpd.GeoDataFrame, bench: gpd.GeoDataFrame,
                 heatrisk: gpd.GeoDataFrame, pet, search_buffer: float) -> None:
        self.cool_places = cool_places  # This will hold the final output
        self.buildings = buildings
        self.bench = bench
        self.heatrisk = heatrisk
        self.pet = pet
        self.search_buffer = search_buffer
        self.eval_shades = []  # Collect evaluation results for each shade geometry column (sdgeom1, sdgeom2, etc.)

    def calculate_walking_shed(self):
        with Progress() as progress:
            task = progress.add_task("Calculate walking shed...", total=3)
            start_time = time.time()
            # Ensure buildings are in a projected CRS
            if self.buildings.crs.is_geographic:
                self.buildings = self.buildings.to_crs(epsg=28992)
            if self.cool_places.crs.is_geographic:
                self.cool_places = self.cool_places.to_crs(epsg=28992)

            buildings['c_id'] = None
            building_sindex = self.buildings.sindex

            for cool_place in self.cool_places.itertuples():
                current_buffer = cool_place.geometry.buffer(self.search_buffer)
                possible_matches_index = list(building_sindex.intersection(current_buffer.bounds))
                possible_matches = self.buildings.iloc[possible_matches_index]
                buildings_within = possible_matches[possible_matches.geometry.intersects(current_buffer)]

                for idx in buildings_within.index:
                    actual_distance = self.buildings.loc[idx].geometry.distance(cool_place.geometry)

                    if self.buildings.loc[idx, 'c_id'] is None or actual_distance < self.buildings.loc[idx, 'dist']:
                        self.buildings.loc[idx, 'dist'] = actual_distance
                        self.buildings.loc[
                            idx, 'c_id'] = cool_place.id

            missing_cool_place = buildings['c_id'].isna().sum()
            if missing_cool_place > 0:
                print(f"Warning: {missing_cool_place} buildings were not assigned to any cool place.")

            progress.advance(task)
            end_time = time.time()
            duration = end_time - start_time
            print(f"Calculating walking shed took {duration:.2f} seconds")
            self.buildings.to_file("build_ws.shp")

            return self.buildings

    def evaluate_resident(self) -> gpd.GeoDataFrame:
        with Progress() as progress:
            task = progress.add_task("Evaluate resident...", total=4)
            start_time = time.time()

            if self.cool_places.crs != self.buildings.crs:
                self.buildings = self.buildings.to_crs(self.cool_places.crs)

            self.buildings['c_id'] = self.buildings['c_id'].fillna(0).astype(float).astype(int).astype(str).str.strip()
            self.cool_places['id'] = self.cool_places['id'].astype(int).astype(str).str.strip()

            building_grouped = self.buildings.groupby('c_id').agg({
                'resident': 'sum',
                'elder_resi': 'sum',
                'kid': 'sum'
            }).reset_index()
            progress.advance(task)

            self.cool_places = self.cool_places.merge(building_grouped, left_on='id', right_on='c_id', how='left')
            progress.advance(task)

            self.cool_places.to_file("cp_ws.shp")
            progress.advance(task)

            end_time = time.time()
            duration = end_time - start_time
            print(f"Calculating walking shed took {duration:.2f} seconds")

            return self.cool_places

    def evaluate_capacity(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:

        shade['Area'] = shade['geometry'].area.round(2)
        shade['cap_area'] = shade.apply(
            lambda row: int(row['Area'] / 10) if row['Area'] > 0 else 0,
            axis=1
        )
        shade = gpd.sjoin(shade, self.cool_places[['geometry', 'resident']], how='left', op='intersects')
        shade['resident'] = shade['resident'].fillna(0)
        shade['cap_status'] = shade.apply(
            lambda row: int(row['cap_area'] - row['resident']) if pd.notnull(row['cap_area']) and pd.notnull(
                row['resident']) else 0,
            axis=1
        )
        # shade['Cap_status'] = shade.apply(
        #     lambda row: (
        #         (lambda x: f"More capacity for {x:.0f} person(s) available" if x > 0
        #         else f"Require more space for {abs(x):.0f} person(s)")(
        #             row['Cap_area'] - row['resident'])
        #     ) if row['resident'] and row['Area'] > 0 else "No data",
        #     axis=1
        # )
        return shade

    def evaluate_sfurniture(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        start_time = time.time()

        bench_join = gpd.sjoin(self.bench, shade, how="inner", predicate='intersects')
        count_benches = bench_join.groupby('id').size().reset_index(name='benches_count')
        shade = shade.merge(count_benches, on='id', how='left')
        shade['Benches'] = shade['benches_count'].fillna(0).apply(
            lambda x: f"{int(x)} bench" if x > 0 else "No available bench"
        )
        end_time = time.time()
        duration = end_time - start_time
        print(f"Calculating benches took {duration:.2f} seconds")
        return shade

    def evaluate_heatrisk(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        start_time = time.time()

        join = gpd.sjoin(shade, self.heatrisk[['geometry', 'HI_TOTAAL_S']], how="inner", predicate='intersects')
        avg_hr = join.groupby('id').agg({'HI_TOTAAL_S': 'mean'}).reset_index()
        shade = shade.merge(avg_hr, on='id', how='left')
        shade['heat_rs'] = shade['HI_TOTAAL_S'].round(2).fillna(0)

        # Classify heat risk score
        classify_hit = jenkspy.JenksNaturalBreaks(n_classes=5)
        classify_hit.fit(shade['heat_rs'].dropna())
        shade['heat_rlv'] = classify_hit.labels_

        classify_level = {
            0: 'Lowest',
            1: 'Below Average',
            2: 'Average',
            3: 'Above Average',
            4: 'Highest'
        }
        shade['heat_rlv'] = shade['heat_rlv'].map(classify_level)
        end_time = time.time()
        duration = end_time - start_time
        print(f"Evaluate heatrisk took {duration:.2f} seconds")
        return shade

    def eval_pet(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        start_time = time.time()
        with rasterio.open(self.pet) as src:
            stats = zonal_stats(shade, self.pet, stats=['mean'], affine=src.transform)
        avg_pet = [stat['mean'] for stat in stats]
        shade['PET'] = avg_pet
        shade['PET'] = shade['PET'].round()
        shade['PET Recom'] = shade['PET'].apply(lambda x: "not recommended" if x > 35 else "recommended")

        end_time = time.time()
        duration = end_time - start_time
        print(f"Eval pet took {duration:.2f} seconds")
        return shade

    def aggregate_to_cool_places(self) -> None:
        """
        After evaluating the shade geometries, aggregate the results to cool places
        based on the common 'id' field between cool places and shade geometries.
        """
        start_time = time.time()
        for i, eval_shade in enumerate(self.eval_shades):
            eval_shade_renamed = eval_shade.add_suffix(f'_{i}')
            columns_to_add = []
            for col in eval_shade_renamed.columns:
                if col in self.cool_places.columns:
                    new_col_name = f"{col}_l{i}"
                    eval_shade_renamed = eval_shade_renamed.rename(columns={col: new_col_name})
                    columns_to_add.append(new_col_name)
                else:
                    columns_to_add.append(col)
            self.cool_places = self.cool_places.join(eval_shade_renamed[columns_to_add], how='left')
            print(f"Processed layer {i + 1} and added columns: {columns_to_add}")

        columns_to_average = ['cap_area', 'cap_status,' 'benches_count', 'heat_rs', 'PET']

        for feature in columns_to_average:
            feature_columns = [col for col in self.cool_places.columns if re.match(f'{feature}_\d+$', col)]
            if feature_columns:
                self.cool_places[feature_columns] = self.cool_places[feature_columns].apply(pd.to_numeric,
                                                                                            errors='coerce')
                self.cool_places[f'{feature}_avg'] = self.cool_places[feature_columns].mean(axis=1)
        end_time = time.time()
        duration = end_time - start_time
        print(f"Integrating attribute to cool space took {duration:.2f} seconds")

    def final_recom(self):
        w_capacity = 0.2
        w_bench = 0.1
        w_hr = 0.2
        w_pet = 0.2
        w_shade = 0.3

        columns_to_fill = ['cap_status_avg', 'benches_count_avg', 'heat_rs_avg', 'PET_avg', 'scDay']

        for col in columns_to_fill:
            if col in self.cool_places.columns:
                self.cool_places[col].fillna(0, inplace=True)

        def classify_pet(pet):
            if pet > 35:
                return 0  # Lowest class
            elif 30 < pet <= 35:
                return 0.8
            else:
                return 1

        def classify_shade(shade):
            if shade == 1:
                return 0.6
            elif shade == 2:
                return 0.8
            elif shade == 3:
                return 1
            else:
                return 0

        def min_max_normalize(column):
            if column.max() == column.min():
                return np.zeros_like(column)
            return (column - column.min()) / (column.max() - column.min())

        self.cool_places['capst_norm'] = min_max_normalize(self.cool_places['cap_status_avg']).round(2)
        self.cool_places['hrs_norm'] = min_max_normalize(self.cool_places['heat_rs_avg']).round(2)
        self.cool_places['bc_norm'] = min_max_normalize(self.cool_places['benches_count_avg']).round(2)

        self.cool_places['PET_c'] = self.cool_places['PET_avg'].apply(classify_pet).round(2)
        self.cool_places['scDay_c'] = self.cool_places['scDay'].apply(classify_shade).round(2)

        self.cool_places['score'] = ((self.cool_places['capst_norm'] * w_capacity) +
                                     (self.cool_places['bc_norm'] * w_bench) +
                                     (self.cool_places['hrs_norm'] * w_hr) +
                                     (self.cool_places['PET_c'] * w_pet) +
                                     (self.cool_places['scDay_c'] * w_shade)).round(2)

        def classify(row):
            if row['scDay_c'] == 0 or row['PET_c'] == 0:  # If shade is 0, automatically "Not Recommended"
                return 'Not recommended'
            elif row['score'] <= 0.4:
                return 'Not recommended'
            elif 0.4 < row['score'] <= 0.6:
                return 'Recommended'
            else:
                return 'Highly Recommended'

        self.cool_places['recom'] = self.cool_places.apply(classify, axis=1)

    def export_eval_gpkg(self, output_gpkg_path: str, layer_name: str) -> None:
        """
        Export the evaluation results to a GeoPackage as a specific layer.
        """

        print("Columns before dropping:", self.cool_places.columns.tolist())

        geometry_columns = [col for col in self.cool_places.columns if col != self.cool_places.geometry.name]
        for col in geometry_columns:
            if isinstance(self.cool_places[col].iloc[0], base.BaseGeometry):
                self.cool_places[col] = gpd.GeoSeries(self.cool_places[col]).to_wkt()

        if 'geometry' in self.cool_places.columns:
            self.cool_places.set_geometry('geometry', inplace=True)
        else:
            raise ValueError("The 'geometry' column is missing in the GeoDataFrame.")

        print("GeoDataFrame structure before export:")
        print(self.cool_places.head())
        print("Remaining columns after dropping:", self.cool_places.columns.tolist())

        self.cool_places.to_file(output_gpkg_path, layer=layer_name, driver='GPKG')
        print(f"Layer '{layer_name}' saved to GeoPackage: {output_gpkg_path}")


