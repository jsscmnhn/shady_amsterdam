import shade_weight_calculation
import cool_places_nodes_calculation
import os
import json
import osmnx as ox

print('Provide the configuration json file. In this file the following data should be included:')
print('1) The path of the graphml file of the network in your study area or the name of the area from which you wish to acquire the OSM network.')
print('2) The path of the file in which you store the shade maps.')
print('3) The path of the file in which you wish to store the new graphml files after incorporating the shade weights into it.')
print('4) The path of the shapefile of the cool place polygons in your study area.')
print('5) The path of the file in which you wish to store the shapefile of the cool places nodes.')


def load_config(config_path="C:/Github_synthesis/configuration/config_dataset.json"):
    if os.path.exists(config_path):
        print("Loading config...")
        with open(config_path, 'r') as file:
            config = json.load(file)
    else:
        raise FileNotFoundError(f"Configuration file not found at {config_path}.")
    return config


if __name__ == "__main__":
    # Configuration
    config = load_config()

    # Load configuration values with default handling
    graphml_file = config.get('graphml_file')
    area_name = config.get('area name for OSM network')
    shade_maps_path = config.get('shade maps path')
    new_graphml_files_path = config.get('new graphml files path')

    cool_places_path = config.get('cool places path')
    cool_places_nodes_path = config.get('cool places nodes path')

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
            "Both 'graphml_file' and 'area name for OSM network' are missing or set to 'None'. Please provide one.")

    # Pass parameters directly to `process_multiple_shade_maps`
    shade_weight_calculation.process_multiple_shade_maps(graph_file=network,raster_dir=shade_maps_path,output_dir=new_graphml_files_path)
    cool_places_nodes_calculation.process_all_shapefiles(polygon_directory=cool_places_path, graph_directory='C:/Github_synthesis/AMS/graphs_with_shade_output',output_directory=cool_places_nodes_path)


    


