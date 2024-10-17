import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from rasterstats import zonal_stats
import osmnx as ox
import networkx as nx
from shapely.geometry import Point, shape
import matplotlib.pyplot as plt
from rasterio.plot import show
import numpy as np
from shapely.geometry import box
import time
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor


def load_graph_from_osm(place):
    # Load graph from OpenStreetMap using osmnx
    graph = ox.graph_from_place(place, network_type='walk')
    return graph

def load_graph_from_file(graph_file_path):
    # Load the graph from a GraphML file
    graph = ox.load_graphml(graph_file_path)
    return graph


def load_combined_graph(places):
    # Load graph for each place and merge them
    combined_graph = None

    for place in places:
        # Load graph for each place
        graph = ox.graph_from_place(place, network_type='walk')

        # Combine the graph with the previous one
        if combined_graph is None:
            combined_graph = graph
        else:
            combined_graph = nx.compose(combined_graph, graph)

    return combined_graph

def load_raster(raster_path):
    # Load the raster file
    raster = rasterio.open(raster_path)
    return raster


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
    # First, create the spatial index for polygons
    spatial_index = cool_place_polygons.sindex

    # Define a lambda function to check if a point is within any cool place polygon
    def is_within_polygon(node_geometry):
        possible_matches_index = list(spatial_index.intersection(node_geometry.bounds))
        possible_matches = cool_place_polygons.iloc[possible_matches_index]
        return possible_matches.contains(node_geometry).any()

    # Apply the function to each node geometry in a vectorized way
    cool_place_nodes = nodes[nodes['geometry'].apply(is_within_polygon)].index

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding cool place nodes took {duration:.2f} seconds")

    return cool_place_nodes.tolist()


def compute_shade_weights(edges, raster, affine, nodata_value=None):
    start_time = time.time()

    # Print the raster CRS and graph CRS for confirmation
    print(f"Graph CRS: {edges.crs}")
    print(f"Raster CRS: {raster.crs}")

    # Compute zonal statistics to get the count of shaded pixels (0 values) for each edge
    shade_stats = zonal_stats(
        vectors=edges['geometry'],
        raster=raster.read(1),
        affine=affine,
        stats=['sum', 'count'],  # Get sum and total pixel count (for normalization)
        nodata=nodata_value,
        all_touched=True  # Include all pixels that touch the geometry
    )

    raster_data = raster.read(1)

    print(f"Raster min value: {np.min(raster_data)}")
    print(f"Raster max value: {np.max(raster_data)}")
    print(f"Raster unique values: {np.unique(raster_data)}")

    # Calculate the proportion of shaded pixels (0 values)
    # Sum of pixel values will give the total "unshaded" pixels (since 1 represents no shade).
    edges['total_unshaded'] = [stat['sum'] if stat['sum'] is not None else 0 for stat in shade_stats]
    edges['total_pixels'] = [stat['count'] if stat['count'] is not None else 1 for stat in
                             shade_stats]  # Avoid division by zero

    # Calculate the proportion of shaded pixels (where 0 indicates shade)
    edges['shade_proportion'] = [(stat['count'] - stat['sum']) / stat['count'] if stat['count'] > 0 else 0
                                 for stat in shade_stats]

    print("Shade Proportion (First 10 edges):", edges['shade_proportion'].head(10))
    print("Shade Proportion Description:", edges['shade_proportion'].describe())

    # Inverse of the shade proportion to calculate weight (more shaded = lower weight)
    edges['shade_weight'] = [1 / (shade + 1e-5) if shade > 0 else 1000 for shade in edges['shade_proportion']]

    print("Shade Weight Statistics (First 10 edges):", edges['shade_weight'].head(10))
    print("Shade Weight Description:", edges['shade_weight'].describe())

    end_time = time.time()
    duration = end_time - start_time
    print(f"Shade weight calculation took {duration:.2f} seconds")

    if 'shade_weight' in edges.columns:
        print("The 'shade_weight' column exists in the edges DataFrame.")
    else:
        print("The 'shade_weight' column is missing after the function.")

    return edges

###########################################################
def process_zonal_stats_chunk(edges_chunk, raster_data, affine):
    # Perform zonal stats for a chunk of edges
    shade_stats_chunk = zonal_stats(
        vectors=edges_chunk['geometry'],
        raster=raster_data,
        affine=affine,
        stats=['sum', 'count'],
        all_touched=True
    )
    return shade_stats_chunk


def compute_shade_weights_parallel(edges, raster, affine, chunk_size=1000, nodata_value=None):
     # Read raster data
    raster_data = raster.read(1)

    # Split edges into chunks for parallel processing
    edges_chunks = [edges.iloc[i:i + chunk_size] for i in range(0, len(edges), chunk_size)]

    # Use ThreadPoolExecutor for parallel computation
    with ThreadPoolExecutor() as executor:
        shade_stats_results = list(executor.map(
            process_zonal_stats_chunk,
            edges_chunks,
            [raster_data] * len(edges_chunks),
            [affine] * len(edges_chunks)
        ))

    # Flatten the results back into a single list
    shade_stats = [item for sublist in shade_stats_results for item in sublist]

    # Add the results to the edges DataFrame
    edges['total_unshaded'] = [stat['sum'] if stat['sum'] is not None else 0 for stat in shade_stats]
    edges['total_pixels'] = [stat['count'] if stat['count'] is not None else 1 for stat in
                             shade_stats]  # Avoid division by zero

    # Calculate the proportion of shaded pixels (0 indicates shade)
    edges['shade_proportion'] = [(stat['count'] - stat['sum']) / stat['count'] if stat['count'] > 0 else 0 for stat in
                                 shade_stats]

    # Inverse of the shade proportion to calculate weight (more shaded = lower weight)
    edges['shade_weight'] = [1 / (shade + 1e-5) if shade > 0 else 1000 for shade in edges['shade_proportion']]

    return edges


def compute_shade_for_edges(edges, raster_path, chunk_size=1000):
    start_time = time.time()
    # Load raster data
    with rasterio.open(raster_path) as raster:
        affine = raster.transform

        # Compute shade weights in parallel
        updated_edges = compute_shade_weights_parallel(edges, raster, affine, chunk_size)

    end_time = time.time()
    duration = end_time - start_time
    print(f"Shade weight calculation took {duration:.2f} seconds")

    return updated_edges

########################################################################################

def add_edge_identifiers_to_gdf(graph, edges_gdf):
    start_time = time.time()
    # Extract edge tuples from the graph (u, v, key)
    edge_tuples = list(graph.edges(keys=True))  # Extract all (u, v, key) tuples

    # Extract u, v, and key from the edge tuples
    u_values = [u for u, v, key in edge_tuples]
    v_values = [v for u, v, key in edge_tuples]
    key_values = [key for u, v, key in edge_tuples]

    # Add these as new columns to the GeoDataFrame
    edges_gdf['u'] = u_values
    edges_gdf['v'] = v_values
    edges_gdf['key'] = key_values

    end_time = time.time()
    duration = end_time - start_time
    print(f"Adding edge identifiers took {duration:.2f} seconds")

    return edges_gdf


def update_graph_with_shade_weight(graph, edges_gdf):
    start_time = time.time()

    # Create a multi-index on 'u', 'v', and 'key' for faster lookups
    edges_gdf.set_index(['u', 'v', 'key'], inplace=True)

    # Create a dictionary for fast lookup of 'shade_weight'
    shade_weight_dict = edges_gdf['shade_weight'].to_dict()

    # Iterate through all edges in the graph and update 'shade_weight'
    for u, v, key, data in graph.edges(keys=True, data=True):
        # Create a tuple key for this edge
        edge_key = (u, v, key)

        # Lookup shade_weight in the dictionary; assign default if not found
        data['shade_weight'] = shade_weight_dict.get(edge_key, 1000)  # Default to 1000 if not found

    # Reset index on edges_gdf after we're done
    edges_gdf.reset_index(inplace=True)

    end_time = time.time()
    duration = end_time - start_time
    print(f"Applying shade weight to the graph edges took {duration:.2f} seconds")

    return graph


# def compute_combined_weights(edges, user_shade_preference):
#     start_time = time.time()
#
#     # Check if 'shade_weight' exists
#     if 'shade_weight' not in edges.columns:
#         print("'shade_weight' column is missing in edges!")
#         return
#
#     # Normalize the length of each edge to be between 0 and 1
#     edges['length_weight'] = edges['length'] / edges['length'].max()
#
#     print("Length Weight Statistics (First 10 edges):", edges['length_weight'].head(10))
#     print("Length Weight Description:", edges['length_weight'].describe())
#
#     # User preference for shade (e.g., 0.7 means 70% weight on shade, 30% on distance)
#     shade_factor = user_shade_preference
#     distance_factor = 1 - shade_factor
#
#     # Calculate the combined weight using the user-defined preference
#     edges['combined_weight'] = (edges['shade_weight'] * shade_factor) + (edges['length_weight'] * distance_factor)
#
#     print("Combined Weight Statistics (First 10 edges):", edges['combined_weight'].head(10))
#     print("Combined Weight Description:", edges['combined_weight'].describe())
#
#     end_time = time.time()
#     duration = end_time - start_time
#     print(f"Combined weight calculation took {duration:.2f} seconds")
#
#     return edges


def find_nearest_cool_place(graph, node, cool_place_nodes, max_distance):
    start_time = time.time()
    # Find cool place nodes within a certain distance
    nearby_nodes = []

    for target in cool_place_nodes:
        try:
            # Only append target if there is a path and it is within max_distance
            path_length = nx.shortest_path_length(graph, node, target)
            if path_length < max_distance:
                nearby_nodes.append(target)
        except nx.NetworkXNoPath:
            # If no path exists between node and target, continue
            continue

    if nearby_nodes:
        # If there are nearby cool place nodes, find the closest one
        nearest_cool_place = min(nearby_nodes, key=lambda target: nx.shortest_path_length(graph, node, target))

        end_time = time.time()
        duration = end_time - start_time
        print(f"Finding the nearest cool place took {duration:.2f} seconds")

        return nearest_cool_place

    return None


def find_nearest_cool_place_dijkstra(graph, start_node, cool_place_nodes, max_distance):
    start_time = time.time()

    try:
        # Perform Dijkstra search from all cool place nodes to the target (start_node)
        distance, path = nx.multi_source_dijkstra(graph, cool_place_nodes, target=start_node, weight='length',
                                                  cutoff=max_distance)

        # Find the node in `cool_place_nodes` that was the source of the shortest path to the start_node
        for cool_place_node in cool_place_nodes:
            if cool_place_node in path:
                # Return the nearest cool place node and its path

                end_time = time.time()
                duration = end_time - start_time
                print(f"Finding the nearest cool place took {duration:.2f} seconds")

                return cool_place_node

    except nx.NetworkXNoPath:
        # If no path exists between the start node and any cool place node, return None
        return None, None



def calculate_routes_to_cool_places(graph, start_node, cool_place_nodes, max_distance=1000, weight='length'):
    # Try finding routes to nearby cool places within the max distance
    nearest_cool_place = find_nearest_cool_place(graph, start_node, cool_place_nodes, max_distance)

    if nearest_cool_place is None:
        print(f"No cool places found within {max_distance} meters.")
        return None, None

    # Calculate shortest path
    shortest_path = nx.shortest_path(graph, start_node, nearest_cool_place, weight='length')

    # Calculate shadiest path
    shadiest_path = nx.shortest_path(graph, start_node, nearest_cool_place, weight='shade_weight')

    return shortest_path, shadiest_path


def calculate_balanced_route(graph, start_node, destination_node):
    start_time = time.time()

    # Calculate the total sum of all lengths and shade weights in the graph
    total_length = sum(data['length'] for u, v, key, data in graph.edges(keys=True, data=True))
    total_shade_weight = sum(data['shade_weight'] for u, v, key, data in graph.edges(keys=True, data=True))

    # Step 1: Normalize the length and assign to "normalized_length"
    for u, v, key, data in graph.edges(keys=True, data=True):
        data['normalized_length'] = data['length'] / total_length

    # Step 2: Normalize the shade weight and assign to "normalized_shade_weight"
    for u, v, key, data in graph.edges(keys=True, data=True):
        data['normalized_shade_weight'] = data['shade_weight'] / total_shade_weight

    # Step 3: Calculate the weighted sum and assign to "weighted_sum_weight"
    for u, v, key, data in graph.edges(keys=True, data=True):
        data['weighted_sum_weight'] = (data['normalized_length'] / 2) + (data['normalized_shade_weight'] / 2)

    # Step 4: Calculate the shortest path using "weighted_sum_weight" as the weight
    balanced_route = nx.shortest_path(graph, start_node, destination_node, weight='weighted_sum_weight')

    end_time = time.time()
    duration = end_time - start_time
    print(f"Balanced route calculation took {duration:.2f} seconds")

    return balanced_route


# def calculate_balanced_route(graph, start_node, destination_node, user_shade_preference, max_distance=1000):
#     for u, v, key, data in graph.edges(keys=True, data=True):
#         if 'shade_weight' not in data:
#             print("shade_weight is missing in the edges attribute!")
#             return
#
#     start_time = time.time()
#
#     # Get the maximum edge length for normalization
#     max_length = max(data['length'] for u, v, key, data in graph.edges(keys=True, data=True))
#
#     # Combine weights for both distance and shade according to user preference
#     for u, v, key, data in graph.edges(keys=True, data=True):
#         length_weight = data['length'] / max_length  # Normalize length
#         shade_weight = data['shade_weight']
#
#         # Calculate combined weight based on user preference
#         combined_weight = (user_shade_preference * shade_weight) + ((1 - user_shade_preference) * length_weight)
#         data['combined_weight'] = combined_weight
#
#     # Calculate the balanced route using the combined weight
#     # balanced_route = nx.shortest_path(graph, start_node, nearest_cool_place, weight='combined_weight')
#     balanced_route = nx.shortest_path(graph, start_node, destination_node, weight='combined_weight')
#
#     end_time = time.time()
#     duration = end_time - start_time
#     print(f"Combined route calculation took {duration:.2f} seconds")
#
#     return balanced_route


# Walking shed calculation
def load_building_polygons(place):
    # Load building polygons from OSM for the specified place
    buildings = ox.geometries_from_place(place, tags={'building': True})
    return buildings


def calculate_walking_shed(buildings, cool_place_nodes, graph, distances=[200, 400, 600]):
    # Ensure buildings and nodes are in a projected CRS
    if buildings.crs.is_geographic:
        buildings = buildings.to_crs(epsg=28992)
    if graph.graph['crs'] != 'epsg:28992':
        graph = ox.project_graph(graph, to_crs='EPSG:28992')

    # Convert nodes to GeoDataFrame
    cool_place_points = gpd.GeoDataFrame(
        {'geometry': [Point(graph.nodes[node]['x'], graph.nodes[node]['y']) for node in cool_place_nodes]},
        crs='EPSG:28992'
    )

    # Create buffers for each distance and assign a 'distance_category' to each building
    buildings['distance_category'] = None  # Initialize column

    for distance in distances:
        # Buffer around cool place nodes by the current distance
        cool_place_buffer = cool_place_points.buffer(distance)

        # Find buildings within this buffer but not already assigned a closer category
        within_distance = buildings[buildings['distance_category'].isnull()]
        buildings_within = within_distance[within_distance.geometry.centroid.within(cool_place_buffer.unary_union)]

        # Assign the distance category (e.g., 200m, 400m) to the buildings within this buffer
        buildings.loc[buildings_within.index, 'distance_category'] = distance

    return buildings


def assign_building_colors(buildings):
    color_map = {200: 'red', 400: 'orange', 600: 'yellow'}

    # Assign colors based on the distance category
    buildings['color'] = buildings['distance_category'].map(color_map)

    # Fill NaN values with a default color, e.g., 'gray'
    buildings['color'] = buildings['color'].fillna('gray')

    return buildings


def plot_colored_walking_shed(buildings):
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot all buildings with their assigned color
    buildings.plot(ax=ax, facecolor=buildings['color'], edgecolor='none', alpha=0.7)

    plt.title("Walking Shed with Distance-based Building Colors")
    plt.show()


def demo_shade_route_calculation(places, raster_path, polygon_path):
    # Load graph, raster, and polygons
    # graph = load_graph_from_osm(place)
    # graph = load_graph_from_file(graph_path)
    graph = load_combined_graph(places)
    raster = load_raster(raster_path)
    cool_place_polygons = load_cool_place_polygons(polygon_path)

    graph = ox.project_graph(graph, to_crs='EPSG:28992')

    # Identify cool place nodes
    cool_place_nodes = find_cool_place_nodes(graph, cool_place_polygons)
    print(f"Found {len(cool_place_nodes)} cool place nodes on the graph.")

    # Calculate shade weight for each edge
    nodes, edges = ox.graph_to_gdfs(graph)

    edges = add_edge_identifiers_to_gdf(graph, edges)

    affine = raster.transform
    # edges = compute_shade_weights(edges, raster, affine)
    # Test batch processing
    edges = compute_shade_for_edges(edges, raster_path, chunk_size=1000)

    graph = update_graph_with_shade_weight(graph, edges)

    # Find the shortest, shadiest, and combined routes for a sample node
    sample_node = list(graph.nodes())[756]  # Choose a sample node as the start point
    # nearest_cool_place = find_nearest_cool_place(graph, sample_node, cool_place_nodes, 1000)
    #  Use Dijkstra's algorithm to find the nearest cool place in terms of shortest path length
    nearest_cool_place = find_nearest_cool_place_dijkstra(graph, sample_node, cool_place_nodes, 1000)
    sample_node_destination = list(graph.nodes())[2389]

    start_time = time.time()
    # Load building polygons
    buildings = load_building_polygons('Amsterdam, Netherlands')

    # Calculate walking shed with different distances
    buildings = calculate_walking_shed(buildings, cool_place_nodes, graph)

    # Assign colors based on distance
    buildings = assign_building_colors(buildings)

    # Plot the buildings with walking shed colors
    plot_colored_walking_shed(buildings)

    end_time = time.time()
    duration = end_time - start_time
    print(f"Walking shed calculation took {duration:.2f} seconds")

    # if nearest_cool_place:
    #     # Shortest path based on distance
    #     shortest_route = nx.shortest_path(graph, sample_node, sample_node_destination, weight='length')
    #
    #     # Shadiest path based on shade weight
    #     shadiest_route = nx.shortest_path(graph, sample_node, sample_node_destination, weight='shade_weight')
    #
    #     # Combined route based on user preference for shade vs. distance
    #     balanced_route = calculate_balanced_route(graph, sample_node, sample_node_destination)
    #
    #
    #     # Plot all routes on the same graph
    #     fig, ax = plt.subplots(figsize=(10, 10))
    #
    #     # Plot the graph
    #     ox.plot_graph(graph, ax=ax, show=False, close=False)
    #
    #     # Plot the shortest route
    #     ox.plot_graph_route(graph, shortest_route, ax=ax, route_color='blue', route_linewidth=2,
    #                         orig_dest_node_color='red', alpha=0.7, show=False, close=False)
    #
    #     # Plot the shadiest route
    #     ox.plot_graph_route(graph, shadiest_route, ax=ax, route_color='green', route_linewidth=3,
    #                         orig_dest_node_color='yellow', alpha=0.7, show=False, close=False)
    #
    #     # Plot the combined route
    #     ox.plot_graph_route(graph, balanced_route, ax=ax, route_color='purple', route_linewidth=6,
    #                         orig_dest_node_color='orange', alpha=0.7, show=True, close=False)
    #
    #     plt.title(f"Shortest (blue), Shadiest (green), and Combined (purple) Routes")
    # else:
    #     print("Could not find a nearest cool place.")


# place = 'Amsterdam, Netherlands'
places = ['Amsterdam, Netherlands', 'Diemen, Netherlands', 'Ouder-Amstel, Netherlands']
graph_path = 'C:/pedestrian_demo_data/ams.graphml'
# raster_path = 'C:/Androniki/pythonProject1/amsterdam_20150701_1630.TIF'
# polygon_path = 'C:/Androniki/pythonProject1/ams_public_space.shp'
raster_path = 'C:/pedestrian_demo_data/amsterdam_time_900.tif'
polygon_path = 'C:/pedestrian_demo_data/public_spaces/ams_public_space.shp'

# User-defined preference
# user_shade_preference = 10

# demo_shade_route_calculation(place, raster_path, polygon_path)
# demo_shade_route_calculation(graph_path, raster_path, polygon_path, user_shade_preference)
demo_shade_route_calculation(places, raster_path, polygon_path)