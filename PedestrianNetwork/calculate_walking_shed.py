import geopandas as gpd
import matplotlib.pyplot as plt
import osmnx as ox
import networkx as nx
from shapely.geometry import Point
import time
from matplotlib.patches import Patch


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


def calculate_walking_shed(buildings, cool_place_nodes, graph, distances=[200, 300, 400, 500]):
    # Ensure buildings and nodes are in a projected CRS
    if buildings.crs.is_geographic:
        buildings = buildings.to_crs(epsg=28992)  # Convert to EPSG:28992 for consistency

    if graph.graph['crs'] != 'epsg:28992':
        graph = ox.project_graph(graph, to_crs='EPSG:28992')

    # Convert cool place nodes to a GeoDataFrame
    cool_place_points = gpd.GeoDataFrame(
        {'geometry': [Point(graph.nodes[node]['x'], graph.nodes[node]['y']) for node in cool_place_nodes]},
        crs='EPSG:28992'
    )

    # Create buffers for each distance and assign a 'distance_category' to each building
    buildings['distance_category'] = None  # Initialize the distance category column

    for distance in distances:
        # Buffer around cool place nodes by the current distance
        cool_place_buffer = cool_place_points.buffer(distance)

        # Find buildings within this buffer but not already assigned a closer distance category
        within_distance = buildings[buildings['distance_category'].isnull()]
        buildings_within = within_distance[within_distance.geometry.centroid.within(cool_place_buffer.unary_union)]

        # Assign the distance category to the buildings within this buffer
        buildings.loc[buildings_within.index, 'distance_category'] = distance

    return buildings


def assign_building_colors(buildings):
    # Define a color map for different distance categories
    color_map = {200: 'green', 300: 'blue', 400: 'yellow', 500:'orange'}

    # Assign colors based on the distance category
    buildings['color'] = buildings['distance_category'].map(color_map)

    # Fill NaN values with a default color, e.g., 'gray' for buildings not within any distance category
    buildings['color'] = buildings['color'].fillna('red')

    return buildings


def plot_colored_walking_shed(buildings):
    fig, ax = plt.subplots(figsize=(10, 10))

    buildings_polygons = buildings[buildings.geometry.type == 'Polygon']

    buildings_polygons.plot(ax=ax, facecolor=buildings_polygons['color'], edgecolor='none', alpha=0.7)

    # # Plot all buildings with their assigned color
    # buildings.plot(ax=ax, facecolor=buildings['color'], edgecolor='none', alpha=0.7)

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
        # Load from a file
        graph = load_graph_from_file(graph_file_path)
    elif place:
        # Load from OpenStreetMap (OSM) using the specified place
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

    # Load building polygons from a shapefile
    buildings = load_building_polygons(building_shapefile_path)

    # Calculate walking shed with different distances
    buildings = calculate_walking_shed(buildings, cool_place_nodes, graph)

    # Assign colors to the buildings based on their proximity to cool places
    buildings = assign_building_colors(buildings)

    # Plot the walking shed with colored buildings
    plot_colored_walking_shed(buildings)


# Example usage
walking_shed_calculation(
    place="Amsterdam, Netherlands",
    polygon_path="C:/Androniki/pythonProject1/ams_public_space.shp",
    building_shapefile_path="C:/Androniki/pythonProject1/new_graph/ams_buildings.shp"
)

# walking_shed_calculation(place="Amsterdam, Netherlands", polygon_path="C:/pedestrian_demo_data/public_spaces/ams_public_space.shp")
