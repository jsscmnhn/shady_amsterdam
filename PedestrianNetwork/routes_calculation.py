import osmnx as ox
import networkx as nx
import geopandas as gpd
import matplotlib.pyplot as plt
import time
from geopy.geocoders import Nominatim  # For location name geocoding
from pyproj import Transformer  # For coordinate transformation
import pickle
from datetime import datetime, timedelta
import os
import walking_shed_network


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


def find_nearest_node(graph, lat, lon):
    # Convert WGS84 (lat, lon) to EPSG:28992 for the graph
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:28992", always_xy=True)  # WGS84 to EPSG:28992
    x, y = transformer.transform(lon, lat)  # Note: pyproj takes (lon, lat), hence the reverse order
    print(f"lon: {x}, lat:{y}")
    nearest_node = ox.nearest_nodes(graph, x, y)  # OSMnx expects projected coordinates in the graph's CRS
    return nearest_node


def geocode_location(location_name):
    geolocator = Nominatim(user_agent="route_planner")
    location = geolocator.geocode(location_name)
    if location:
        return location.latitude, location.longitude
    else:
        raise ValueError(f"Location '{location_name}' not found.")


def find_cool_place_nodes(graph, cool_place_polygons):
    start_time = time.time()

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
    cool_place_nodes = nodes[nodes['geometry'].apply(is_within_polygon)].index

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding cool place nodes took {duration:.2f} seconds")

    return cool_place_nodes.tolist()


def load_cool_place_nodes(pre_calculated_path):
    # Load pre-calculated cool place nodes from a file (using pickle)
    with open(pre_calculated_path, 'rb') as f:
        cool_place_nodes = pickle.load(f)

    print(f"Loaded {len(cool_place_nodes)} cool place nodes from {pre_calculated_path}.")
    return cool_place_nodes


def find_nearest_cool_place_dijkstra(graph, start_node, cool_place_nodes, max_distance):
    start_time = time.time()

    try:
        # Perform Dijkstra search from all cool place nodes to the target (start_node)
        distance, path = nx.multi_source_dijkstra(graph, cool_place_nodes, target=start_node, weight='length',
                                                  cutoff=max_distance)
        # print("Performed Dijkstra search from all cool place nodes to the target.")

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
        print(f"No path found from start node {start_node} to any cool place node.")
        return None


def calculate_routes(graph, start_node, end_node=None, cool_place_nodes=None, max_distance=5000):
    # If end_node is None, calculate routes to the nearest cool place
    if end_node is None and cool_place_nodes is not None:
        nearest_cool_place = find_nearest_cool_place_dijkstra(graph, start_node, cool_place_nodes, max_distance)
        if nearest_cool_place is None:
            return None, None, None
        end_node = nearest_cool_place

    # Shortest path based on distance
    shortest_path = nx.shortest_path(graph, start_node, end_node, weight='length')

    # Shadiest path based on shade weight
    shadiest_path = nx.shortest_path(graph, start_node, end_node, weight='shade_weight')

    # Combined route 1: 70/30 shade/length
    balanced_route_1 = calculate_balanced_route(graph, start_node, end_node, shade_weight_ratio=70,
                                                length_weight_ratio=30)

    # Combined route 2: 30/70 shade/length
    balanced_route_2 = calculate_balanced_route(graph, start_node, end_node, shade_weight_ratio=30,
                                                length_weight_ratio=70)

    return shortest_path, shadiest_path, balanced_route_1, balanced_route_2


def calculate_balanced_route(graph, start_node, destination_node, shade_weight_ratio=70, length_weight_ratio=30):
    start_time = time.time()

    # Calculate the total sum of all lengths and shade weights in the graph
    total_length = sum(data['length'] for u, v, key, data in graph.edges(keys=True, data=True))
    total_shade_weight = sum(data['shade_weight'] for u, v, key, data in graph.edges(keys=True, data=True))

    # Normalize the length and assign to "normalized_length"
    for u, v, key, data in graph.edges(keys=True, data=True):
        data['normalized_length'] = data['length'] / total_length

    # Normalize the shade weight and assign to "normalized_shade_weight"
    for u, v, key, data in graph.edges(keys=True, data=True):
        data['normalized_shade_weight'] = data['shade_weight'] / total_shade_weight

    # Calculate the weighted sum and assign to "weighted_sum_weight"
    for u, v, key, data in graph.edges(keys=True, data=True):
        # data['weighted_sum_weight'] = (data['normalized_length'] * (30 / 100)) + (data['normalized_shade_weight'] * (70 / 100))
        data['weighted_sum_weight'] = (data['normalized_length'] * (length_weight_ratio / 100)) + \
                                      (data['normalized_shade_weight'] * (shade_weight_ratio / 100))

    # Calculate the shortest path using "weighted_sum_weight" as the weight
    balanced_route = nx.shortest_path(graph, start_node, destination_node, weight='weighted_sum_weight')

    end_time = time.time()
    duration = end_time - start_time
    print(f"Balanced route calculation with {shade_weight_ratio}% shade and {length_weight_ratio}% length took {duration:.2f} seconds")

    return balanced_route


def plot_routes(graph, shortest_route, shadiest_route, balanced_route_1, balanced_route_2):
    # Extract the coordinates of all nodes involved in the routes
    route_nodes = set(shortest_route + shadiest_route + balanced_route_1 + balanced_route_2)

    # Get the geometry of each node in the route
    node_coordinates = [(graph.nodes[node]['x'], graph.nodes[node]['y']) for node in route_nodes]
    xs, ys = zip(*node_coordinates)  # Split into x and y coordinates

    # Define the bounding box with a margin
    margin = 50  # Set a margin (in meters) around the routes for better visualization
    x_min, x_max = min(xs) - margin, max(xs) + margin
    y_min, y_max = min(ys) - margin, max(ys) + margin

    # Plot all routes on the same graph
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot the full graph in the background
    ox.plot_graph(graph, ax=ax, show=False, close=False, node_color='none', edge_color='gray')

    # Plot the shortest route
    ox.plot_graph_route(graph, shortest_route, ax=ax, route_color='blue', route_linewidth=2,
                        orig_dest_node_color='red', alpha=0.7, show=False, close=False)

    # Plot the shadiest route
    ox.plot_graph_route(graph, shadiest_route, ax=ax, route_color='green', route_linewidth=3,
                        orig_dest_node_color='yellow', alpha=0.7, show=False, close=False)

    # Plot the combined (balanced) route 1: 70/30 shade/length
    ox.plot_graph_route(graph, balanced_route_1, ax=ax, route_color='purple', route_linewidth=6,
                        orig_dest_node_color='orange', alpha=0.7, show=False, close=False)

    # Plot the combined (balanced) route 2: 30/70 shade/length
    ox.plot_graph_route(graph, balanced_route_2, ax=ax, route_color='pink', route_linewidth=4,
                        orig_dest_node_color='cyan', alpha=0.7, show=False, close=False)

    # Set the plot limits to zoom in on the routes
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    # Create proxy artists for the legend
    from matplotlib.lines import Line2D
    legend_lines = [Line2D([0], [0], color='blue', lw=2),
                    Line2D([0], [0], color='green', lw=3),
                    Line2D([0], [0], color='purple', lw=6),
                    Line2D([0], [0], color='pink', lw=4)]

    # Add a legend
    legend_labels = ['Shortest Route (blue)', 'Shadiest Route (green)',
                     'Combined1 Route (70% Shade, 30% Length) (purple)',
                     'Combined2 Route (30% Shade, 70% Length) (pink)']

    ax.legend(legend_lines, legend_labels, loc='upper left')

    plt.title(f"Shortest (blue), Shadiest (green), Combined1 (purple), Combined2 (pink) Routes")
    plt.show()


def demo_shade_route_calculation(graph_file_path, pre_calculated_nodes_path, user_input, input_type, mode="nearest_cool_place"):
    # Load precomputed graph with shade weights
    print("Loading precomputed graph with shade weights...")
    graph = load_graph_from_file(graph_file_path)

    graph = ox.project_graph(graph, to_crs='EPSG:28992')

    # Load pre-calculated cool place nodes
    cool_place_nodes = load_cool_place_nodes(pre_calculated_nodes_path)
    if not cool_place_nodes:
        print("No cool place nodes were found. Exiting...")
        return

    print(f"Found {len(cool_place_nodes)} cool place nodes on the graph.")

    # Handle user input (coordinates or location name)
    if input_type == "coordinates":
        if mode == "origin_destination":
            # Unpack the two coordinates for origin and destination
            origin_lat, origin_lon, dest_lat, dest_lon = user_input
            start_node = find_nearest_node(graph, origin_lat, origin_lon)
            end_node = find_nearest_node(graph, dest_lat, dest_lon)
        else:
            lat, lon = user_input
            start_node = find_nearest_node(graph, lat, lon)
            end_node = None
    elif input_type == "location_name":
        if mode == "origin_destination":
            origin_name, destination_name = user_input
            origin_lat, origin_lon = geocode_location(origin_name)
            dest_lat, dest_lon = geocode_location(destination_name)
            start_node = find_nearest_node(graph, origin_lat, origin_lon)
            end_node = find_nearest_node(graph, dest_lat, dest_lon)
        else:
            lat, lon = geocode_location(user_input)
            start_node = find_nearest_node(graph, lat, lon)
            end_node = None
    else:
        raise ValueError("Invalid input type. Use 'coordinates' or 'location_name'.")

    # If mode is 'nearest_cool_place', calculate routes to nearest cool place
    if mode == "nearest_cool_place":
        shortest_route, shadiest_route, balanced_route_1, balanced_route_2 = calculate_routes(graph, start_node,
                                                                          cool_place_nodes=cool_place_nodes)
    # If mode is 'origin_destination', expect user_input to contain both origin and destination
    elif mode == "origin_destination":
        shortest_route, shadiest_route, balanced_route_1, balanced_route_2 = calculate_routes(graph, start_node, end_node)
    else:
        raise ValueError("Invalid mode. Use 'nearest_cool_place' or 'origin_destination'.")

    # Plot the routes if found
    if shortest_route and shadiest_route and balanced_route_1 and balanced_route_2:
        print(f"Found routes.")
        plot_routes(graph, shortest_route, shadiest_route, balanced_route_1, balanced_route_2)
    else:
        print("Could not find any route.")


def find_nearest_timestamp_files(date_time, graph_dir, nodes_dir):
    # Parse available files in the directories
    graph_files = [f for f in os.listdir(graph_dir) if f.startswith("ams_graph_with_shade")]
    nodes_files = [f for f in os.listdir(nodes_dir) if f.startswith("cool_places_nodes")]

    # Extract dates and times, ensuring they are zero-padded and parsed as datetime objects
    graph_timestamps = []
    for f in graph_files:
        parts = f.split('_')
        date_str = parts[4]
        time_str = parts[5].split('.')[0].zfill(4)  # Zero-pad time to ensure HHMM format
        timestamp = datetime.strptime(date_str + time_str, "%Y%m%d%H%M")
        graph_timestamps.append((timestamp, f))  # Store both the datetime and filename

    # Find closest date to user input
    closest_date = min(graph_timestamps, key=lambda x: abs(x[0].date() - date_time.date()))[0].date()

    # Filter files to keep only those with the closest date
    closest_date_files = [t for t in graph_timestamps if t[0].date() == closest_date]

    # From files with closest date, find the closest time
    nearest_timestamp, nearest_graph_file = min(
        closest_date_files,
        key=lambda x: abs(
            x[0].time().hour * 60 + x[0].time().minute - (date_time.time().hour * 60 + date_time.time().minute))
    )
    print(f"Nearest Timestamp: {nearest_timestamp} for Graph File: {nearest_graph_file}")

    # Generate the expected nodes filename for the nearest timestamp
    nearest_date = nearest_timestamp.strftime("%Y%m%d")
    nearest_time = nearest_timestamp.strftime("%H%M")
    nodes_file = f"cool_places_nodes_{nearest_date}_{nearest_time}.pkl"

    # Verify if the nodes file exists
    if nodes_file in nodes_files:
        print(f"Selected Graph File: {nearest_graph_file}")
        print(f"Selected Cool Places Nodes File: {nodes_file}")
        return os.path.join(graph_dir, nearest_graph_file), os.path.join(nodes_dir, nodes_file)
    else:
        print("Error: No matching nodes file for the nearest timestamp found.")
        return None, None


def demo_shade_route_calculation_with_time(graph_dir, nodes_dir, user_input, input_type, mode="nearest_cool_place",
                                           date_time=None):
    # Determine the datetime if not provided
    if date_time is None:
        date_time = datetime.now()

    # Find the nearest files based on date and time
    graph_file_path, pre_calculated_nodes_path = find_nearest_timestamp_files(date_time, graph_dir, nodes_dir)

    # Proceed with the original function using the selected files
    demo_shade_route_calculation(graph_file_path, pre_calculated_nodes_path, user_input, input_type, mode)


# # Example Usage
# graph_dir = 'C:/pedestrian_demo_data/graphs_with_shade/'
# nodes_dir = 'C:/pedestrian_demo_data/cool_place_nodes/'

# demo_shade_route_calculation_with_time(graph_dir, nodes_dir, user_input=("Amsterdam Central Station", "Dam Square"),
#                                        input_type="location_name", mode="origin_destination")

# # Option 1: Routes to nearest cool place with coordinates input
# demo_shade_route_calculation(graph_file_path, polygon_path, user_input=(52.373169, 4.890660), input_type="coordinates",
#                              mode="nearest_cool_place")

# # Option 2: Routes between two locations with location name input
# demo_shade_route_calculation(graph_file_path, pre_calculated_nodes_path, user_input=("Amsterdam Central Station", "Dam Square"),
#                              input_type="location_name", mode="origin_destination")

# # Option 3: Routes to nearest cool place with location name input
# demo_shade_route_calculation(graph_file_path, polygon_path, user_input=("Amsterdam Central Station"),
#                              input_type="location_name", mode="nearest_cool_place")

# # Option 4: Routes between two locations with coordinates input
# demo_shade_route_calculation(graph_file_path, pre_calculated_nodes_path, user_input=((52.373169, 4.890660, 52.376522, 4.908490)), input_type="coordinates",
#                              mode="origin_destination")