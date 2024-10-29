import geopandas as gpd
import networkx as nx
import osmnx as ox
from shapely.geometry import Point
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import time
from scipy.spatial import KDTree
import numpy as np
from multiprocessing import Pool, cpu_count
from osmnx.distance import nearest_nodes
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Manager
from scipy.spatial import KDTree
import osmnx as ox



from shapely.geometry import MultiPolygon



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


def find_cool_place_nodes(graph, cool_place_polygons, max_distance=100):
    start_time = time.time()

    # Convert graph nodes to a GeoDataFrame
    nodes_gdf = ox.graph_to_gdfs(graph, nodes=True, edges=False)

    # Ensure CRS matches between nodes and polygons
    if nodes_gdf.crs != cool_place_polygons.crs:
        cool_place_polygons = cool_place_polygons.to_crs(nodes_gdf.crs)

    # Filter out rows with None geometries
    cool_place_polygons = cool_place_polygons[cool_place_polygons.geometry.notnull()]

    # Extract node coordinates and create KDTree for quick nearest neighbor search
    node_coords = [(geom.x, geom.y) for geom in nodes_gdf.geometry]
    kd_tree = KDTree(node_coords)

    # Map from KDTree result indices to graph node IDs
    node_indices = list(nodes_gdf.index)

    # Find the nearest network node for each cool place polygon or MultiPolygon centroid
    cool_place_nodes = []
    for _, polygon in cool_place_polygons.iterrows():
        if polygon.geometry is None:
            print("Skipping None geometry in cool_place_polygons.")
            continue

        # Check if the geometry is a MultiPolygon
        if isinstance(polygon.geometry, MultiPolygon):
            # Process each polygon in the MultiPolygon
            for sub_polygon in polygon.geometry.geoms:
                centroid = sub_polygon.centroid
                distance, nearest_idx = kd_tree.query((centroid.x, centroid.y))
                if distance <= max_distance:
                    nearest_node = node_indices[nearest_idx]
                    cool_place_nodes.append(nearest_node)
                else:
                    print(f"Skipping centroid at {centroid} with nearest node distance {distance:.2f}m (greater than {max_distance}m)")
        else:
            # Process single Polygon
            centroid = polygon.geometry.centroid
            distance, nearest_idx = kd_tree.query((centroid.x, centroid.y))
            if distance <= max_distance:
                nearest_node = node_indices[nearest_idx]
                cool_place_nodes.append(nearest_node)
            else:
                print(f"Skipping centroid at {centroid} with nearest node distance {distance:.2f}m (greater than {max_distance}m)")

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding cool place nodes took {duration:.2f} seconds")

    return cool_place_nodes



def get_nearest_node(graph, point):
    # Find the nearest graph node to a given point
    return ox.distance.nearest_nodes(graph, X=point.x, Y=point.y)


def find_nodes_within_distances(graph, cool_place_nodes, distances=[200, 300, 400, 500]):
    start_time = time.time()
    print("Starting to find nodes within specified distance thresholds...")

    # Dictionary to store nodes within each distance threshold
    distance_nodes = {distance: set() for distance in distances}

    # For each cool place node, perform Dijkstra to find nodes within the distance thresholds
    for i, cool_place_node in enumerate(cool_place_nodes):
        if i % 10 == 0:  # Print every 10 cool place nodes to track progress
            print(f"Processing cool place node {i+1}/{len(cool_place_nodes)}")

        for distance in distances:
            lengths = nx.single_source_dijkstra_path_length(graph, cool_place_node, cutoff=distance, weight='length')
            distance_nodes[distance].update(lengths.keys())
            print(f"  Found {len(lengths)} nodes within {distance} meters of cool place node {cool_place_node}")

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding nodes within distance took {duration:.2f} seconds")

    return distance_nodes


def assign_buildings_to_distance_category(buildings, graph, distance_nodes):
    start_time = time.time()
    print("Starting to assign distance categories to buildings based on proximity...")

    # Initialize distance category for each building
    buildings['distance_category'] = None

    # Assign each building to a distance category based on its nearest graph node
    for idx, building in buildings.iterrows():
        if idx % 10 == 0:  # Print every 10 buildings to track progress
            print(f"Assigning distance category for building {idx+1}/{len(buildings)}")

        building_centroid = building.geometry.centroid
        nearest_building_node = get_nearest_node(graph, building_centroid)

        # Assign the building to the closest distance category where its node falls
        assigned_category = None
        for distance in sorted(distance_nodes.keys()):
            if nearest_building_node in distance_nodes[distance]:
                buildings.at[idx, 'distance_category'] = distance
                assigned_category = distance
                break

        print(f"  Building {idx+1} assigned to distance category: {assigned_category or '> 500m'}")

    end_time = time.time()
    duration = end_time - start_time
    print(f"Assigning buildings to distance category took {duration:.2f} seconds")

    return buildings


def assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes):
    start_time = time.time()
    print("Starting to assign distance categories to buildings...")

    # Find nearest nodes for all building centroids
    building_centroids = buildings.geometry.centroid
    nearest_nodes = ox.distance.nearest_nodes(graph, X=building_centroids.x, Y=building_centroids.y)
    buildings['nearest_node'] = nearest_nodes

    # Assign distance category based on precomputed distance nodes
    def get_distance_category(node_id):
        for distance in sorted(distance_nodes.keys()):
            if node_id in distance_nodes[distance]:
                return distance
        return '>500m'  # If the node isn't within the maximum distance, assign ">500m"

    buildings['distance_category'] = buildings['nearest_node'].apply(get_distance_category)

    end_time = time.time()
    duration = end_time - start_time
    print(f"Assigning buildings to distance category took {duration:.2f} seconds")

    return buildings


def assign_building_distance_category(building, graph, distance_nodes):
    # Calculate the nearest node for the building and assign the distance category
    building_centroid = building.geometry.centroid
    nearest_building_node = get_nearest_node(graph, building_centroid)

    # Assign the building to the closest distance category where its node falls
    assigned_category = '> 500m'
    for distance in sorted(distance_nodes.keys()):
        if nearest_building_node in distance_nodes[distance]:
            assigned_category = distance
            break

    building['distance_category'] = assigned_category

    # Print the result for this building for monitoring purposes
    print(f"Processed building with ID {building.name}: Assigned distance category {assigned_category}")
    return building

def parallel_assign_buildings_to_distance_category(buildings, graph, distance_nodes):
    start_time = time.time()
    print("Starting parallel assignment of distance categories to buildings with ThreadPoolExecutor...")

    # List to store the processed buildings
    processed_buildings = []

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor() as executor:
        # Submit tasks to the executor for each building
        futures = [executor.submit(assign_building_distance_category, building, graph, distance_nodes) for _, building in buildings.iterrows()]

        # Track the progress as each future completes
        for i, future in enumerate(as_completed(futures), 1):
            processed_buildings.append(future.result())

            # Print progress for every 10 buildings processed
            if i % 10 == 0:
                print(f"Processed {i}/{len(buildings)} buildings...")

    # Combine processed buildings back into a single GeoDataFrame
    buildings = gpd.GeoDataFrame(processed_buildings)

    end_time = time.time()
    duration = end_time - start_time
    print(f"ThreadPoolExecutor assignment of buildings to distance category took {duration:.2f} seconds")
    print("All buildings have been assigned distance categories.")

    return buildings


def calculate_walking_shed_with_distance_categories(buildings, cool_place_nodes, graph, distances=[200, 300, 400, 500]):
    # Find nodes within each distance threshold for cool place nodes
    distance_nodes = find_nodes_within_distances(graph, cool_place_nodes, distances)

    # Assign distance category to buildings based on proximity to cool place nodes
    # buildings = assign_buildings_to_distance_category(buildings, graph, distance_nodes)
    # buildings  = parallel_assign_buildings_to_distance_category(buildings, graph, distance_nodes)
    buildings = assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes)

    return buildings


def assign_building_colors(buildings):
    # Define a color map for different distance categories
    color_map = {200: 'green', 300: 'blue', 400: 'yellow', 500: 'orange'}

    # Assign colors based on the distance category
    buildings['color'] = buildings['distance_category'].map(color_map)

    # Fill NaN values with a default color
    buildings['color'] = buildings['color'].fillna('red')

    return buildings


def plot_colored_walking_shed(buildings):
    fig, ax = plt.subplots(figsize=(10, 10))

    buildings_polygons = buildings[buildings.geometry.type == 'Polygon']
    buildings_polygons.plot(ax=ax, facecolor=buildings_polygons['color'], edgecolor='none', alpha=0.7)

    legend_elements = [
        Patch(facecolor='green', edgecolor='none', label='< 200m'),
        Patch(facecolor='blue', edgecolor='none', label='200m - 300m'),
        Patch(facecolor='yellow', edgecolor='none', label='300m - 400m'),
        Patch(facecolor='orange', edgecolor='none', label='400m - 500m'),
        Patch(facecolor='red', edgecolor='none', label='> 500m')
    ]

    ax.legend(handles=legend_elements, loc='upper right', title='Distance to Cool Place')

    plt.title("Walking Shed with Distance-based Building Colors")
    plt.show()


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

    print('number of cool spaces:',len(cool_place_polygons))

    cool_place_nodes = find_cool_place_nodes(graph, cool_place_polygons)
    print(f"Found {len(cool_place_nodes)} cool place nodes on the graph.")

    buildings = load_building_polygons(building_shapefile_path)

    print("Calculating walking shed with distance categories...")
    # Calculate walking shed using network distances
    buildings = calculate_walking_shed_with_distance_categories(buildings, cool_place_nodes, graph)

    print("Assigning colors to the buildings...")
    # Assign colors to the buildings based on their proximity to cool places
    buildings = assign_building_colors(buildings)


    # Plot the walking shed with colored buildings
    plot_colored_walking_shed(buildings)


if __name__ == "__main__":
    walking_shed_calculation(
        place="Amsterdam, Netherlands",
        polygon_path="C:/Github_synthesis/AMS/cool_places_polygons/cool_places_polygons_20230816_900.shp",
        building_shapefile_path="C:/Androniki/pythonProject1/merged_buildings.shp")
        # polygon_path="C:/pedestrian_demo_data/public_spaces/ams_public_space.shp",
        # building_shapefile_path="C:/pedestrian_demo_data/ams_buildings/ams_buildings.shp"
