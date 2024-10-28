import geopandas as gpd
import networkx as nx
import osmnx as ox
from shapely.geometry import Point
import time
from concurrent.futures import ThreadPoolExecutor, as_completed



def load_building_polygons(building_shapefile_path):
    # Load building polygons from a shapefile on your local machine
    buildings = gpd.read_file(building_shapefile_path)
    return buildings


def load_graph_from_osm(place):
    # Load the graph from OpenStreetMap using OSMnx for a given place
    graph = ox.graph_from_place(place, network_type='walk')
    return graph


def load_graph_from_file(graph_file_path):
    # Load the graph from a GraphML file
    graph = ox.load_graphml(graph_file_path)
    return graph


def process_graph(graph):
    # Ensure the graph is in the correct CRS (EPSG:28992)
    if graph.graph['crs'] != 'epsg:28992':
        graph = ox.project_graph(graph, to_crs='EPSG:28992')
    return graph


def load_cool_place_polygons(polygon_path):
    # Load polygons representing cool places from a shapefile or GeoJSON
    polygons = gpd.read_file(polygon_path)
    return polygons


def find_cool_place_nodes(graph, cool_place_polygons):
    start_time = time.time()

    # Get nodes as GeoDataFrame
    nodes = ox.graph_to_gdfs(graph, nodes=True, edges=False)

    # Ensure that both geometries are in the same CRS (Coordinate Reference System)
    if nodes.crs != cool_place_polygons.crs:
        cool_place_polygons = cool_place_polygons.to_crs(nodes.crs)

    # Use spatial indexing to find nodes within polygons
    spatial_index = cool_place_polygons.sindex

    def is_within_polygon(node_geometry):
        possible_matches_index = list(spatial_index.intersection(node_geometry.bounds))
        possible_matches = cool_place_polygons.iloc[possible_matches_index]
        return possible_matches.contains(node_geometry).any()

    cool_place_nodes = nodes[nodes['geometry'].apply(is_within_polygon)].index

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding cool place nodes took {duration:.2f} seconds")

    return cool_place_nodes.tolist()


def find_nodes_within_distances(graph, cool_place_nodes, distances=[200, 300, 400, 500]):
    start_time = time.time()
    print("Starting to find nodes within specified distance thresholds...")

    # Dictionary to store nodes within each distance threshold for each cool place
    distance_nodes = {cool_place_node: {distance: set() for distance in distances} for cool_place_node in
                      cool_place_nodes}

    # For each cool place node, perform Dijkstra to find nodes within the distance thresholds
    for i, cool_place_node in enumerate(cool_place_nodes):
        if i % 10 == 0:  # Print every 10 cool place nodes to track progress
            print(f"Processing cool place node {i + 1}/{len(cool_place_nodes)}")

        for distance in distances:
            lengths = nx.single_source_dijkstra_path_length(graph, cool_place_node, cutoff=distance, weight='length')
            distance_nodes[cool_place_node][distance].update(lengths.keys())
            print(f"  Found {len(lengths)} nodes within {distance} meters of cool place node {cool_place_node}")

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding nodes within distance took {duration:.2f} seconds")

    return distance_nodes


def assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes):
    start_time = time.time()
    print("Starting to assign distance categories to buildings...")

    # Calculate nearest nodes only once
    building_centroids = buildings.geometry.centroid
    nearest_nodes = ox.distance.nearest_nodes(graph, X=building_centroids.x, Y=building_centroids.y)
    buildings['nearest_node'] = nearest_nodes

    # Define the modified distance category assignment function
    def get_distance_category(building_index, node_id):
        centroid_node = nearest_nodes[building_index]  # Use precomputed nearest node

        for cool_place_node, distances_dict in distance_nodes.items():
            for distance in sorted(distances_dict.keys()):
                if node_id in distances_dict[distance]:
                    try:
                        path_length = nx.shortest_path_length(graph, source=centroid_node, target=cool_place_node,
                                                              weight='length')
                        if path_length <= distance:
                            return distance
                    except nx.NetworkXNoPath:
                        continue
        return '>500m'

    # Use ThreadPoolExecutor to parallelize building assignment
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(get_distance_category, idx, node_id) for idx, node_id in
                   enumerate(buildings['nearest_node'])]
        buildings['distance_category'] = [future.result() for future in as_completed(futures)]

    end_time = time.time()
    duration = end_time - start_time
    print(f"Assigning buildings to distance category took {duration:.2f} seconds")

    return buildings


def assign_building_colors(buildings):
    # Define a color map for different distance categories
    color_map = {200: 'green', 300: 'blue', 400: 'yellow', 500: 'orange'}

    # Assign colors based on the distance category
    buildings['color'] = buildings['distance_category'].map(color_map)

    # Fill NaN values with a default color
    buildings['color'] = buildings['color'].fillna('red')

    return buildings


def walking_shed_calculation(place=None, graph_file_path=None, polygon_path=None, building_shapefile_path=None):
    if graph_file_path:
        graph = load_graph_from_file(graph_file_path)
    elif place:
        graph = load_graph_from_osm(place)
    else:
        raise ValueError("You must provide either a place name or a graph file path.")

    graph = process_graph(graph)

    if polygon_path:
        cool_place_polygons = load_cool_place_polygons(polygon_path)
    else:
        raise ValueError("You must provide the path to a cool place polygon file.")

    cool_place_nodes = find_cool_place_nodes(graph, cool_place_polygons)
    print(f"Found {len(cool_place_nodes)} cool place nodes on the graph.")

    buildings = load_building_polygons(building_shapefile_path)

    print("Calculating walking shed with distance categories...")
    # Calculate nodes within distance categories for each cool place
    distance_nodes = find_nodes_within_distances(graph, cool_place_nodes)

    # Assign distance category to buildings based on proximity to cool place nodes
    buildings = assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes)

    print("Assigning colors to the buildings...")
    # Assign colors to the buildings based on their proximity to cool places
    buildings = assign_building_colors(buildings)

    # Plot the walking shed with colored buildings
    plot_colored_walking_shed(buildings)


if __name__ == "__main__":
    walking_shed_calculation(
        place="Amsterdam, Netherlands",
        polygon_path="C:/Androniki/pythonProject1/ams_public_space.shp",
        building_shapefile_path="C:/Androniki/pythonProject1/merged_buildings.shp"
        # polygon_path="C:/pedestrian_demo_data/public_spaces/ams_public_space.shp",
        # building_shapefile_path="C:/pedestrian_demo_data/ams_buildings/ams_buildings.shp"
    )