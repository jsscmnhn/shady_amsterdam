from include.road_process import Road
from include.cool_space import CoolSpace
from include.building import Building
import geopandas as gpd
import numpy as np
import rasterio
import matplotlib.pyplot as plt
import glob
import os

def print_hi(name):
    print(f'Hi, {name}')


# entry
if __name__ == '__main__':
    print_hi('PyCharm')
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

