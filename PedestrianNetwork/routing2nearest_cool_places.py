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

def load_graph_from_osm(place):
    # Load graph from OpenStreetMap using osmnx
    graph = ox.graph_from_place(place, network_type='walk')
    return graph


def load_raster(raster_path):
    # Load the raster file
    raster = rasterio.open(raster_path)
    return raster


def load_cool_place_polygons(polygon_path):
    # Load polygons representing cool places from a shapefile or GeoJSON
    polygons = gpd.read_file(polygon_path)
    return polygons


def find_cool_place_nodes(graph, cool_place_polygons):
    # Get nodes as GeoDataFrame
    nodes, edges = ox.graph_to_gdfs(graph)

    print(f"Graph CRS: {edges.crs}")
    print(f"Polygon CRS: {cool_place_polygons.crs}")

    fig, ax = plt.subplots(figsize=(10, 10))
    nodes.plot(ax=ax, color='yellow', markersize=10, alpha = 0.7, label='Graph Nodes')
    cool_place_polygons.plot(ax=ax, color='red', alpha=0.5, label='Cool Place Polygons')

    # Create a spatial index for the nodes
    cool_place_nodes = []

    for idx, node in nodes.iterrows():
        node_point = Point(node['x'], node['y'])
        # Check if node is within any cool place polygon
        if cool_place_polygons.contains(node_point).any():
            cool_place_nodes.append(idx)

    return cool_place_nodes


def compute_shade_weights(edges, raster, affine, nodata_value=None):
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

    return edges


def compute_combined_weights(edges, user_shade_preference):
    # Check if 'shade_weight' exists
    if 'shade_weight' not in edges.columns:
        print("'shade_weight' column is missing in edges!")
        return

    # Normalize the length of each edge to be between 0 and 1
    edges['length_weight'] = edges['length'] / edges['length'].max()

    print("Length Weight Statistics (First 10 edges):", edges['length_weight'].head(10))
    print("Length Weight Description:", edges['length_weight'].describe())

    # User preference for shade (e.g., 0.7 means 70% weight on shade, 30% on distance)
    shade_factor = user_shade_preference
    distance_factor = 1 - shade_factor

    # Calculate the combined weight using the user-defined preference
    edges['combined_weight'] = (edges['shade_weight'] * shade_factor) + (edges['length_weight'] * distance_factor)

    print("Combined Weight Statistics (First 10 edges):", edges['combined_weight'].head(10))
    print("Combined Weight Description:", edges['combined_weight'].describe())

    return edges


def find_nearest_cool_place(graph, node, cool_place_nodes, max_distance):
    # Find cool place nodes within a certain distance
    nearby_nodes = [target for target in cool_place_nodes if
                    nx.shortest_path_length(graph, node, target) < max_distance]

    if nearby_nodes:
        # If there are nearby cool place nodes, find the closest one
        nearest_cool_place = min(nearby_nodes, key=lambda target: nx.shortest_path_length(graph, node, target))
        return nearest_cool_place
    return None


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


def calculate_balanced_route(graph, start_node, cool_place_nodes, user_shade_preference, max_distance=1000):
    # Try finding routes to nearby cool places within the max distance
    nearest_cool_place = find_nearest_cool_place(graph, start_node, cool_place_nodes, max_distance)

    if nearest_cool_place is None:
        print(f"No cool places found within {max_distance} meters.")
        return None

    # Get the maximum edge length for normalization
    max_length = max(data['length'] for u, v, key, data in graph.edges(keys=True, data=True))

    # Combine weights for both distance and shade according to user preference
    for u, v, key, data in graph.edges(keys=True, data=True):
        length_weight = data['length'] / max_length  # Normalize length
        shade_weight = data['shade_weight']

        # Calculate combined weight based on user preference
        combined_weight = (user_shade_preference * shade_weight) + ((1 - user_shade_preference) * length_weight)
        data['combined_weight'] = combined_weight

    # Calculate the balanced route using the combined weight
    balanced_route = nx.shortest_path(graph, start_node, nearest_cool_place, weight='combined_weight')

    return balanced_route


def demo_shade_route_calculation(place, raster_path, polygon_path, user_shade_preference):
    # Load graph, raster, and polygons
    graph = load_graph_from_osm(place)
    raster = load_raster(raster_path)
    cool_place_polygons = load_cool_place_polygons(polygon_path)

    graph = ox.project_graph(graph, to_crs='EPSG:28992')

    # Identify cool place nodes
    cool_place_nodes = find_cool_place_nodes(graph, cool_place_polygons)
    print(f"Found {len(cool_place_nodes)} cool place nodes on the graph.")

    # Calculate shade weight for each edge
    nodes, edges = ox.graph_to_gdfs(graph)

    edge_tuples = list(graph.edges(keys=True))  # Extract all (u, v, key) tuples from the graph
    u_values = [u for u, v, key in edge_tuples]
    v_values = [v for u, v, key in edge_tuples]
    key_values = [key for u, v, key in edge_tuples]

    edges['u'] = u_values
    edges['v'] = v_values
    edges['key'] = key_values

    affine = raster.transform
    edges = compute_shade_weights(edges, raster, affine)

    # Apply the shade weight to the graph edges
    for u, v, key, data in graph.edges(keys=True, data=True):
        # Match the graph edges with the edges GeoDataFrame using u, v, and key
        edge_idx = edges[(edges['u'] == u) & (edges['v'] == v) & (edges['key'] == key)].index
        if len(edge_idx) > 0:
            data['shade_weight'] = edges.loc[edge_idx[0], 'shade_weight']
        else:
            data['shade_weight'] = 1000  # Assign a large weight if no shade weight is found for safety

    # Find the shortest, shadiest, and combined routes for a sample node
    sample_node = list(graph.nodes())[555]  # Choose a sample node as the start point
    nearest_cool_place = find_nearest_cool_place(graph, sample_node, cool_place_nodes, 1000)

    if nearest_cool_place:
        # Shortest path based on distance
        shortest_route = nx.shortest_path(graph, sample_node, nearest_cool_place, weight='length')

        # Shadiest path based on shade weight
        shadiest_route = nx.shortest_path(graph, sample_node, nearest_cool_place, weight='shade_weight')

        # Combined route based on user preference for shade vs. distance
        balanced_route = calculate_balanced_route(graph, sample_node, cool_place_nodes, user_shade_preference)

        # Plot all routes on the same graph
        fig, ax = plt.subplots(figsize=(10, 10))

        # Plot the graph
        ox.plot_graph(graph, ax=ax, show=False, close=False)

        # Plot the shortest route
        ox.plot_graph_route(graph, shortest_route, ax=ax, route_color='blue', route_linewidth=2,
                            orig_dest_node_color='red', alpha=0.7, show=False, close=False)

        # Plot the shadiest route
        ox.plot_graph_route(graph, shadiest_route, ax=ax, route_color='green', route_linewidth=3,
                            orig_dest_node_color='yellow', alpha=0.7, show=False, close=False)

        # Plot the combined route
        ox.plot_graph_route(graph, balanced_route, ax=ax, route_color='purple', route_linewidth=6,
                            orig_dest_node_color='orange', alpha=0.7, show=True, close=False)

        plt.title(f"Shortest (blue), Shadiest (green), and Combined (purple) Routes")
    else:
        print("Could not find a nearest cool place.")


place = 'Amsterdam, Netherlands'
raster_path = 'C:/pedestrian_demo_data/amsterdam_time_900.tif'
polygon_path = 'C:/pedestrian_demo_data/public_spaces/ams_public_space.shp'

# User-defined preference
user_shade_preference = 0.0

# demo_shade_route_calculation(place, raster_path, polygon_path)
demo_shade_route_calculation(place, raster_path, polygon_path, user_shade_preference)