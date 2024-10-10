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

        # Spatial join to find buildings within the buffers
        buildings_within_buffers = gpd.sjoin(self.buildings, buffers, how="inner", predicate='intersects')

        # Debugging output
        print("Buildings within buffers:")
        print(buildings_within_buffers[['resident']].head())

        # Initialize a capacity status column
        self.cool_places['capacity_status'] = "No nearest building"
        self.cool_places['Area'] = self.cool_places['geometry'].area.round(2)


        # Check if there are buildings found
        if not buildings_within_buffers.empty:
            # Sum residents from the buildings within each buffer
            resident_sum = buildings_within_buffers.groupby(buildings_within_buffers.index).agg(
                {'resident': 'sum', 'elder_resi': 'sum'}).reset_index()

            # Rename the resident column to resident_total
            resident_sum = resident_sum.rename(columns={'resident': 'resident_total'})

            # Debugging output to check resident_sum
            print("Resident sum after grouping:")
            print(resident_sum)

            # Reset index only if 'index' column doesn't already exist
            if 'index' in self.cool_places.columns:
                self.cool_places.drop(columns='index', inplace=True)
            self.cool_places.reset_index(inplace=True)

            if 'index' in resident_sum.columns:
                resident_sum.drop(columns='index', inplace=True)
            resident_sum.reset_index(inplace=True)
            # Join the sum back to the cool places
            self.eval = self.cool_places.join(resident_sum.set_index('index'), on='index', how='left')

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
                         row['Capacity per area']- row['resident_total'])
                ) if row['resident_total'] and row['Area'] > 0 else "No data",
                axis=1
            )

            self.eval['Resident'] = self.eval['resident_total'].fillna(0).apply(
                lambda x: f"{int(x)} residents" if x > 0 else "No nearest building"
            )
            self.eval['Elderly Resident'] = self.eval['elder_resi'].fillna(0).apply(
                lambda x: f"{int(x)} elderly residents" if x > 0 else "No elderly resident"
            )
        else:
            self.eval = self.cool_places.copy()  # Just keep the original cool places if no buildings

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
        join = gpd.sjoin(self.eval, self.heatrisk[['geometry', 'HI_TOTAAL_S']], how="inner", predicate='intersects')
        print("join columns:", join.columns)
        avg_hr = join.groupby('index').agg({'HI_TOTAAL_S': 'mean'}).reset_index()
        self.eval = self.eval.merge(avg_hr, how='left', left_index=True, right_on='index')
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
    def export_eval(self, output_path: str) -> None:
        # Specify the output path for the shapefile
        self.eval = self.eval.drop(columns='buffer')
        self.eval = self.eval[['geometry','Capacity per area', 'Capacity status', 'Resident', 'Elderly Resident', 'Area', 'Benches','Heat Risk Score', 'Heat Risk Level', 'PET']]
        # Export the GeoDataFrame to a shapefile
        self.eval.to_file(output_path, driver='ESRI Shapefile')

        print(f"Output shapefile saved at: {output_path}")

if __name__ == '__main__':

    # lu = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\datasets\\ams_landuse_top10NL.shp"
    pop = "bpop20.shp"
    bench = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\benchamsosm.shp"
    # coolplace = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\qualified_geometriess.shp"
    coolplace = "lu_space.shp"
    heatrisk = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\Hitte kaarten gemeente Amsterdam\\Klimaatrisicokaarten QGIS\\Risicokaarten definitief.gdb"
    pet = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\Hitte kaarten gemeente Amsterdam\\Kaarten door TAUW ontwikkeld (2022)\\PET_average - Gemeente Amsterdam.tiff"

    # Read the data
    cool_places = gpd.read_file(coolplace)
    buildings = gpd.read_file(pop)
    bench = gpd.read_file(bench)
    heatrisk = gpd.read_file(heatrisk, layer="Risico_per_buurt_20231009_enkel_thema_scores")

    output_path = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\eval_cool_places20.shp"

    # Initialize the CoolEval class
    cool_eval = CoolEval(cool_places, buildings, bench, heatrisk, pet, 700)
    cool_eval.evaluate_capacity()
    cool_eval.evaluate_sfurniture()
    cool_eval.evaluate_heatrisk()
    cool_eval.eval_pet()
    cool_eval.export_eval(output_path)








