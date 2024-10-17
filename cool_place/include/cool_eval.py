import geopandas as gpd
import jenkspy
from rasterstats import zonal_stats
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.features import shapes
from scipy.ndimage import label
from shapely.geometry import shape, box
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import fiona
import os

class CoolEval:
    def __init__(self, cool_places: gpd.GeoDataFrame, buildings: gpd.GeoDataFrame, bench: gpd.GeoDataFrame, heatrisk: gpd.GeoDataFrame, pet, buffer_house: float) -> None:
        self.cool_places = cool_places
        self.buildings = buildings
        self.bench = bench
        self.heatrisk = heatrisk
        self.pet = pet
        self.buffer_house = buffer_house
        self.eval = None

    def evaluate_capacity(self) -> None:
        # Ensure CRS are the same
        if self.cool_places.crs != self.buildings.crs:
            self.buildings = self.buildings.to_crs(self.cool_places.crs)

        # Buffer the cool places
        self.cool_places['buffer'] = None

        self.cool_places['buffer'] = self.cool_places.geometry.buffer(self.buffer_house)
        if 'buffer' in self.cool_places.columns:
            print("Buffered geometry column created.")
        else:
            print("Buffered geometry column failed to be created.")
        buffers = self.cool_places.set_geometry('buffer')
        buffer_shp_path = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\buffers.shp"
        # Export the buffer if a path is provided
        if buffer_shp_path:
            buffer_gdf = self.cool_places[['id', 'buffer']]  # Select relevant columns, including buffer
            buffer_gdf = buffer_gdf.set_geometry('buffer')  # Set 'buffer' as the active geometry column

            # Export buffer to a shapefile
            buffer_gdf.to_file(buffer_shp_path, driver="ESRI Shapefile")
            print(f"Buffer saved to shapefile: {buffer_shp_path}")


        # Spatial join to find buildings within the buffers
        buildings_within_buffers = gpd.sjoin(self.buildings, buffers,  how='inner', predicate='intersects')

        # Initialize a capacity status column
        self.cool_places['capacity_status'] = "No nearest building"
        self.cool_places['Area'] = self.cool_places['geometry'].area.round(2)

        # buildings_shp_path = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\buildings_within_buffers.shp"
        # buildings_within_buffers.to_file(buildings_shp_path, driver="ESRI Shapefile")

        # Sum residents from the buildings within each buffer
        # resident_sum = buildings_within_buffers.groupby(buildings_within_buffers.index).agg(
        #     {'resident': 'sum', 'elder_resi': 'sum'}).reset_index()
        resident_sum = buildings_within_buffers.groupby('index_right')['resident'].sum().reset_index()

        # Rename the resident column to resident_total
        resident_sum = resident_sum.rename(columns={'resident': 'resident_total'})

        # Reset index only if 'index' column doesn't already exist
        if 'index' in self.cool_places.columns:
            self.cool_places.drop(columns='index', inplace=True)
        self.cool_places.reset_index(inplace=True)
        if 'index' in resident_sum.columns:
            resident_sum.drop(columns='index', inplace=True)
        resident_sum.reset_index(inplace=True)
        print(resident_sum.columns)

        # Join the sum back to the cool places
        self.eval = self.cool_places.merge(resident_sum,  left_index=True, right_on='index', how='left')

        # Step: Calculate the number of residents per square meter (density)
        self.eval['Capacity per area'] = self.eval.apply(
            lambda row: int(row['Area'] / 10) if row['Area'] > 0 else 0,
            axis=1
        )

        # Check for NaN values and update capacity status accordingly
        self.eval['Capacity status'] = self.eval.apply(
            lambda row: (
                (lambda x: f"More capacity for {x:.0f} person(s) available" if x > 0
                else f"Require more space for {abs(x):.0f} person(s)")(
                    row['Capacity per area'] - row['resident_total'])
            ) if row['resident_total'] and row['Area'] > 0 else "No data",
            axis=1
        )

        self.eval['Resident'] = self.eval['resident_total'].fillna(0).apply(
            lambda x: f"{int(x)} residents" if x > 0 else "No nearest building"
        )
        # self.eval['Elderly Resident'] = self.eval['elder_resi'].fillna(0).apply(
        #     lambda x: f"{int(x)} elderly residents" if x > 0 else "No elderly resident"
        # )

    def evaluate_sfurniture(self) -> None:

        self.bench = self.bench.to_crs(self.eval.crs)
        bench_join = gpd.sjoin(self.bench, self.eval, how="inner", predicate='intersects')
        # Debugging: Check the columns of the resulting join
        print("bench_join columns:", bench_join.columns)
        count_benches = bench_join.groupby('index_right0').size()
        self.eval['Benches'] = self.eval.index.map(count_benches).fillna(0).astype(int)
        self.eval['Benches'] = self.eval['Benches'].fillna(0).apply(
                lambda x: f"{int(x)} bench" if x > 0 else "No available bench"
            )

    def evaluate_heatrisk(self) -> None:
        # Check and drop any existing 'index_right' column to prevent conflicts
        if 'index_right' in self.eval.columns:
            self.eval = self.eval.drop(columns='index_right')
        join = gpd.sjoin(self.eval, self.heatrisk[['geometry', 'HI_TOTAAL_S']], how="inner", predicate='intersects')
        print("join columns:", join.columns)
        avg_hr = join.groupby('index').agg({'HI_TOTAAL_S': 'mean'}).reset_index()
        self.eval = self.eval.merge(avg_hr, how='left', left_index=True, right_on='index', suffixes=('_left', '_right'))
        self.eval['Heat Risk Score'] = self.eval['HI_TOTAAL_S'].round(2).fillna(0)

        classify_hit = jenkspy.JenksNaturalBreaks(n_classes=5)
        classify_hit.fit(self.eval['Heat Risk Score'].dropna())
        self.eval['Heat Risk Level'] = classify_hit.labels_
        classify_level = {
            0: 'Lowest',
            1: 'Below Average',
            2: 'Average',
            3: 'Above Average',
            4: 'Highest'
        }
        self.eval['Heat Risk Level'] = self.eval['Heat Risk Level'].map(classify_level)

    def eval_pet(self) -> None:
        with rasterio.open(self.pet) as src:
            stats = zonal_stats(self.eval, self.pet, stats=['mean'], affine=src.transform)
        avg_pet = [stat['mean'] for stat in stats]
        self.eval['PET'] = avg_pet
        self.eval['PET'] = self.eval['PET'].round()
        # Add recommendation based on PET values
        self.eval['PET Recom'] = self.eval['PET'].apply(lambda x: "not recommended" if x > 35 else "recommended")
    def export_eval_gpkg(self, output_gpkg_path: str, layer_name: str) -> None:
        """
        Export the evaluation results to a GeoPackage as a specific layer.
        """
        # Dynamically find the 'shadeArea' column (could be shadeArea0, shadeArea1, etc.)
        shade_area_col = [col for col in self.eval.columns if col.startswith('shadeArea')][0]

        # Drop buffer if it exists
        self.eval = self.eval.drop(columns='buffer', errors='ignore')

        # Dynamically include the identified shadeArea column
        self.eval = self.eval[['geometry', 'id', shade_area_col, 'Capacity per area',
                               'Capacity status', 'Resident','resident_total' ,
                               'Area', 'Benches', 'Heat Risk Score', 'Heat Risk Level',
                               'PET', 'PET Recom']]
        self.eval.to_file(output_gpkg_path, layer=layer_name, driver='GPKG')
        print(f"Layer '{layer_name}' saved to GeoPackage: {output_gpkg_path}")

if __name__ == '__main__':


    pop = "bpop20.shp"
    bench = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\benchamsosm.shp"
    coolplace_gpkg = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\datasets\\shadeGeoms.gpkg"
    heatrisk = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\Hitte kaarten gemeente Amsterdam\\Klimaatrisicokaarten QGIS\\Risicokaarten definitief.gdb"
    pet = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\Hitte kaarten gemeente Amsterdam\\Kaarten door TAUW ontwikkeld (2022)\\PET_average - Gemeente Amsterdam.tiff"

    # Read the data

    buildings = gpd.read_file(pop)
    bench = gpd.read_file(bench)
    heatrisk = gpd.read_file(heatrisk, layer="Risico_per_buurt_20231009_enkel_thema_scores")

    # Output GeoPackage
    output_gpkg = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\eval_cool_placesleft.gpkg"

    # Get the list of layers in the input GeoPackage
    layers = fiona.listlayers(coolplace_gpkg)

    # Ask the user to provide a range of layers to process
    start_layer = int(input("Enter the starting layer index (e.g., 0 for shadeGeom0): "))
    end_layer = int(input("Enter the ending layer index (e.g., 5 for shadeGeom5): "))

    # Filter the layers to process based on user input
    selected_layers = [layer for layer in layers if layer.startswith('shadeGeom')
                       and int(layer.replace('shadeGeom', '')) in range(start_layer, end_layer + 1)]

    # Process each selected layer
    for layer in selected_layers:
        cool_places = gpd.read_file(coolplace_gpkg, layer=layer)

        # Initialize the CoolEval class for the current cool space layer
        cool_eval = CoolEval(cool_places, buildings, bench, heatrisk, pet, 800)

        # Perform evaluations
        cool_eval.evaluate_capacity()
        cool_eval.evaluate_sfurniture()
        cool_eval.evaluate_heatrisk()
        cool_eval.eval_pet()

        # Export results as a new layer in the output GeoPackage
        cool_eval.export_eval_gpkg(output_gpkg, layer_name=f"eval_{layer}")

    print("Selected layers processed and saved.")







