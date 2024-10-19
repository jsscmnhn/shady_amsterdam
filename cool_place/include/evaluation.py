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


    def evaluate_capacity(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        # Ensure CRS are the same
        if shade.crs != self.buildings.crs:
            self.buildings = self.buildings.to_crs(shade.crs)

        # Buffer the shade geometries
        shade['buffer'] = shade.geometry.buffer(self.search_buffer)

        # Spatial join to find buildings within the buffers
        buffers = shade.set_geometry('buffer')
        buildings_within_buffers = gpd.sjoin(self.buildings, buffers, how='inner', predicate='intersects')

        # Initialize columns to hold residents and area
        shade['Area'] = shade['geometry'].area.round(2)

        # Sum residents from the buildings within each buffer
        resident_sum = buildings_within_buffers.groupby('id').agg(
            {'resident': 'sum', 'elder_resi': 'sum', 'kid': 'sum'}).reset_index()

        # Merge resident sum back to the shade geometries
        shade = shade.merge(resident_sum, on='id', how='left', suffixes=('', f'_{layer_name}'))

        # Step: Calculate capacity per area
        shade['Capacity per area'] = shade.apply(
            lambda row: int(row['Area'] / 10) if row['Area'] > 0 else 0,
            axis=1
        )

        # Add capacity status
        shade['Capacity status'] = shade.apply(
            lambda row: (
                (lambda x: f"More capacity for {x:.0f} person(s) available" if x > 0
                else f"Require more space for {abs(x):.0f} person(s)")(
                    row['Capacity per area'] - row['resident'])
            ) if row['resident'] and row['Area'] > 0 else "No data",
            axis=1
        )

        # Add descriptions for residents, elderly, and children
        shade['Residents'] = shade['resident'].fillna(0).apply(
            lambda x: f"{int(x)} residents" if x > 0 else "No nearest building"
        )
        shade['Elders'] = shade['elder_resi'].fillna(0).apply(
            lambda x: f"{int(x)} elderly residents" if x > 0 else "No elderly resident"
        )
        shade['Children'] = shade['kid'].fillna(0).apply(
            lambda x: f"{int(x)} children" if x > 0 else "No children"
        )

        return shade

    def evaluate_sfurniture(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        # Make sure the CRS matches
        self.bench = self.bench.to_crs(shade.crs)

        # Perform spatial join between benches and the shade layer
        bench_join = gpd.sjoin(self.bench, shade, how="inner", predicate='intersects')

        # Count the number of benches per shade area
        count_benches = bench_join.groupby('id').size().reset_index(name='benches_count')

        # Merge the bench counts back into the shade geometries
        shade = shade.merge(count_benches, on='id', how='left')
        shade['Benches'] = shade['benches_count'].fillna(0).apply(
            lambda x: f"{int(x)} bench" if x > 0 else "No available bench"
        )
        return shade

    def evaluate_heatrisk(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        # Spatial join between the heat risk and shade layers
        join = gpd.sjoin(shade, self.heatrisk[['geometry', 'HI_TOTAAL_S']], how="inner", predicate='intersects')

        # Calculate the mean heat risk score per shade geometry
        avg_hr = join.groupby('id').agg({'HI_TOTAAL_S': 'mean'}).reset_index()

        # Merge heat risk score back into the shade geometries
        shade = shade.merge(avg_hr, on='id', how='left')
        shade['Heat Risk Score'] = shade['HI_TOTAAL_S'].round(2).fillna(0)

        # Classify heat risk score
        classify_hit = jenkspy.JenksNaturalBreaks(n_classes=5)
        classify_hit.fit(shade['Heat Risk Score'].dropna())
        shade['Heat Risk Level'] = classify_hit.labels_

        classify_level = {
            0: 'Lowest',
            1: 'Below Average',
            2: 'Average',
            3: 'Above Average',
            4: 'Highest'
        }
        shade['Heat Risk Level'] = shade['Heat Risk Level'].map(classify_level)

        return shade

    def eval_pet(self, shade: gpd.GeoDataFrame, layer_name: str) -> gpd.GeoDataFrame:
        # Zonal statistics for PET (mean PET per shade geometry)
        with rasterio.open(self.pet) as src:
            stats = zonal_stats(shade, self.pet, stats=['mean'], affine=src.transform)
        avg_pet = [stat['mean'] for stat in stats]
        shade['PET'] = avg_pet
        shade['PET'] = shade['PET'].round()
        shade['PET Recom'] = shade['PET'].apply(lambda x: "not recommended" if x > 35 else "recommended")

        return shade

    def aggregate_to_cool_places(self) -> None:
        """
        After evaluating the shade geometries, aggregate the results to cool places
        based on the common 'id' field between cool places and shade geometries.
        """
        for i, eval_shade in enumerate(self.eval_shades):
            # Rename columns with suffix to avoid duplication
            eval_shade_renamed = eval_shade.add_suffix(f'_{i}')

            # Ensure no duplicate column names
            columns_to_add = []
            for col in eval_shade_renamed.columns:
                if col in self.cool_places.columns:
                    # Append a unique suffix if column already exists
                    new_col_name = f"{col}_l{i}"
                    eval_shade_renamed = eval_shade_renamed.rename(columns={col: new_col_name})
                    columns_to_add.append(new_col_name)
                else:
                    columns_to_add.append(col)


            # Merge only the new/renamed columns to avoid duplicates
            self.cool_places = self.cool_places.join(eval_shade_renamed[columns_to_add], how='left')



            print(f"Processed layer {i + 1} and added columns: {columns_to_add}")

        # Step 2: Identify the columns for each feature group (e.g., resident_0, resident_1, etc.)
        columns_to_average = ['resident', 'elder_resi', 'kid', 'Capacity per area', 'benches_count',
                                  'Heat Risk Score', 'PET']

        # Iterate over each feature group and calculate the average
        for feature in columns_to_average:
        # Use regular expression to match columns exactly like resident_0, resident_1, etc.
            feature_columns = [col for col in self.cool_places.columns if re.match(f'{feature}_\d+$', col)]

            if feature_columns:  # Check if there are matching columns
                # Step 3: Convert columns to numeric, replacing non-numeric values with NaN
                self.cool_places[feature_columns] = self.cool_places[feature_columns].apply(pd.to_numeric,
                                                                                                errors='coerce')

                # Step 4: Calculate the average across these columns, ignoring NaNs
                self.cool_places[f'{feature}_avg'] = self.cool_places[feature_columns].mean(axis=1)

    def export_eval_gpkg(self, output_gpkg_path: str, layer_name: str) -> None:
        """
        Export the evaluation results to a GeoPackage as a specific layer.
        """
        # Print the columns before any manipulation
        print("Columns before dropping:", self.cool_places.columns.tolist())

        geometry_columns = [col for col in self.cool_places.columns if col != self.cool_places.geometry.name]
        for col in geometry_columns:
            if isinstance(self.cool_places[col].iloc[0], base.BaseGeometry):
                self.cool_places[col] = gpd.GeoSeries(self.cool_places[col]).to_wkt()



        # Ensure the 'geometry' column is set as the active geometry
        if 'geometry' in self.cool_places.columns:
            self.cool_places.set_geometry('geometry', inplace=True)
        else:
            raise ValueError("The 'geometry' column is missing in the GeoDataFrame.")

        # Print the structure of the GeoDataFrame before exporting
        print("GeoDataFrame structure before export:")
        print(self.cool_places.head())
        print("Remaining columns after dropping:", self.cool_places.columns.tolist())

        # Export to the GeoPackage
        self.cool_places.to_file(output_gpkg_path, layer=layer_name, driver='GPKG')
        print(f"Layer '{layer_name}' saved to GeoPackage: {output_gpkg_path}")


