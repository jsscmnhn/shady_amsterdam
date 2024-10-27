import pickle
import osmnx as ox
import networkx as nx
import geopandas as gpd
import json
import os

# Default paths
default_config = {
    "polygon_path": "C:/Androniki/pythonProject1/ams_public_space.shp",
    "graph_file_path": "C:/Androniki/pythonProject1/AMS/Amsterdam_pedestrian_network.graphml",
    "output_cool_place_nodes": "C:/Androniki/pythonProject1/outputs/cool_place_nodes.pkl"
}


# Load configuration from JSON file or use defaults
def load_config(config_path="C:/Androniki/pythonProject1/config100.json"):
    if os.path.exists(config_path):
        print("Loading config...")
        with open(config_path, 'r') as file:
            config = json.load(file)
            # Fill missing values with defaults
            for key, value in default_config.items():
                config.setdefault(key, value)
    else:
        config = default_config
        print(f"Configuration file not found. Using default paths.")
    return config


# Configuration
config = load_config()


def load_graph_from_file(graph_file_path):
    graph = ox.load_graphml(graph_file_path)

    # Ensure shade_weight is numeric (float) for shortest path calculations
    for u, v, key, data in graph.edges(keys=True, data=True):
        if 'shade_weight' in data:
            data['shade_weight'] = float(data['shade_weight'])
        else:
            data['shade_weight'] = 1000.0  # Default high value if shade_weight is missing
    return graph


def load_cool_place_polygons(polygon_path):
    polygons = gpd.read_file(polygon_path)
    return polygons


def calculate_and_save_cool_place_nodes(graph, cool_place_polygons, output_path):
    nodes = ox.graph_to_gdfs(graph, nodes=True, edges=False)

    if nodes.crs != cool_place_polygons.crs:
        cool_place_polygons = cool_place_polygons.to_crs(nodes.crs)

    spatial_index = cool_place_polygons.sindex

    def is_within_polygon(node_geometry):
        possible_matches_index = list(spatial_index.intersection(node_geometry.bounds))
        possible_matches = cool_place_polygons.iloc[possible_matches_index]
        return possible_matches.contains(node_geometry).any()

    cool_place_nodes = nodes[nodes['geometry'].apply(is_within_polygon)].index.tolist()

    with open(output_path, 'wb') as f:
        pickle.dump(cool_place_nodes, f)

    print(f"Cool place nodes calculated and saved to {output_path}.")

# Load data based on config
polygon_path = config["polygon_path"]
graph_file_path = config["graph_file_path"]
pre_calculated_nodes_path = config["output_cool_place_nodes"]

# Load the cool place polygons and calculate cool place nodes
G = ox.load_graphml(graph_file_path)
graph = ox.project_graph(G, to_crs='EPSG:28992')
cool_place_polygons = load_cool_place_polygons(polygon_path)
calculate_and_save_cool_place_nodes(graph, cool_place_polygons, pre_calculated_nodes_path)
