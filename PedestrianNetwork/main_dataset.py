import shade_weight_calculation
import cool_places_nodes_calculation
import os
import json
import osmnx as ox

print('Provide the configuration json file with paths for your data sources and outputs.')

def load_config(config_path="C:/Github_synthesis/configuration/config_pedestrian_network_analysis.json"):
    if os.path.exists(config_path):
        print("Loading config...")
        with open(config_path, 'r') as file:
            config = json.load(file)
    else:
        raise FileNotFoundError(f"Configuration file not found at {config_path}.")
    return config

if __name__ == "__main__":
    # Load the configuration
    config = load_config()

    # Access paths from the structured JSON file
    graphml_file = config["shade_weights"].get("graphml_file")
    area_name = config["shade_weights"].get("area_name")
    shade_maps_path = config["shade_weights"].get("shade_maps_path")
    new_graphml_files_path = config["shade_weights"].get("new_graphml_files_path")

    cool_places_path = config["cool_places"].get("cool_places_path")
    cool_places_nodes_path = config["cool_places"].get("cool_places_nodes_path")

    # Determine which network configuration to use
    if graphml_file == "None":
        network_WGS84 = ox.graph_from_place(area_name, network_type="walk")
        network = ox.project_graph(network_WGS84, to_crs='EPSG:28992')
        print(f"Using OSM network for area: {area_name}")
    elif area_name == "None":
        network_WGS84 = ox.load_graphml(graphml_file)
        network = ox.project_graph(network_WGS84, to_crs='EPSG:28992')
        print(f"Using provided GraphML file: {graphml_file}")
    else:
        raise ValueError(
            "Both 'graphml_file' and 'area_name' are missing or set to 'None'. Please provide one.")

    # Pass parameters to `process_multiple_shade_maps`
    shade_weight_calculation.process_multiple_shade_maps(
        graph_file=network,
        raster_dir=shade_maps_path,
        output_dir=new_graphml_files_path
    )

    cool_places_nodes_calculation.process_all_shapefiles(
        polygon_directory=cool_places_path,
        graph_directory='C:/Github_synthesis/AMS/graphs_with_shade_output',
        output_directory=cool_places_nodes_path
    )
