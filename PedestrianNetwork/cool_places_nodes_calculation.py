import pickle
import osmnx as ox
import networkx as nx
import geopandas as gpd

def load_graph_from_file(graph_file_path):
    # Load the graph from a GraphML file
    graph = ox.load_graphml(graph_file_path)

    # Ensure shade_weight is numeric (float) for shortest path calculations
    for u, v, key, data in graph.edges(keys=True, data=True):
        if 'shade_weight' in data:
            data['shade_weight'] = float(data['shade_weight'])  # Convert to float if necessary
        else:
            data['shade_weight'] = 1000.0  # Assign a high default value if shade_weight is missing

    return graph

def load_cool_place_polygons(polygon_path):
    # Load polygons representing cool places from a shapefile or GeoJSON
    polygons = gpd.read_file(polygon_path)
    return polygons


def calculate_and_save_cool_place_nodes(graph, cool_place_polygons, output_path):
    # Get nodes as GeoDataFrame
    nodes = ox.graph_to_gdfs(graph, nodes=True, edges=False)

    # Ensure that both geometries are in the same CRS
    if nodes.crs != cool_place_polygons.crs:
        cool_place_polygons = cool_place_polygons.to_crs(nodes.crs)

    # Use spatial indexing to find nodes within polygons
    spatial_index = cool_place_polygons.sindex

    def is_within_polygon(node_geometry):
        possible_matches_index = list(spatial_index.intersection(node_geometry.bounds))
        possible_matches = cool_place_polygons.iloc[possible_matches_index]
        return possible_matches.contains(node_geometry).any()

    # Apply the function to each node geometry
    cool_place_nodes = nodes[nodes['geometry'].apply(is_within_polygon)].index.tolist()

    # Save the cool place nodes to a file (using pickle for simplicity)
    with open(output_path, 'wb') as f:
        pickle.dump(cool_place_nodes, f)

    print(f"Cool place nodes calculated and saved to {output_path}.")


polygon_path = 'C:/pedestrian_demo_data/public_spaces/ams_public_space.shp'
graph_file_path = 'C:/pedestrian_demo_data/ams_graph_with_shade_900_cropped.graphml'
pre_calculated_nodes_path = 'C:/pedestrian_demo_data/cool_place_nodes.pkl'

graph = load_graph_from_file(graph_file_path)
graph = ox.project_graph(graph, to_crs='EPSG:28992')
cool_place_polygons = load_cool_place_polygons(polygon_path)
calculate_and_save_cool_place_nodes(graph, cool_place_polygons, pre_calculated_nodes_path)