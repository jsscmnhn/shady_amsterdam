import geopandas as gpd
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import time
from osmnx.distance import nearest_nodes
from shapely.geometry import MultiPolygon, Point

def load_building_polygons(building_shapefile_path):
    """
    Load building polygons from a shapefile.

    Parameters:
    - building_shapefile_path: Path to the shapefile containing building polygons.

    Returns:
    - GeoDataFrame containing building geometries.
    """

    # Load building polygons from a shapefile on your local machine
    buildings = gpd.read_file(building_shapefile_path)
    return buildings


def load_graph_from_osm(place):
    """
    Load a graph from OpenStreetMap for a specified place using OSMnx.

    Parameters:
    - place: The name of the place to load the graph for (e.g., "Amsterdam, Netherlands").

    Returns:
    - Graph representing the walking network for the specified place.
    """

    # Load the graph from OpenStreetMap using OSMnx for a given place
    graph = ox.graph_from_place(place, network_type='walk')
    return graph


def load_graph_from_file(graph_file_path):
    """
    Load a pre-saved graph from a GraphML file.

    Parameters:
    - graph_file_path: Path to the GraphML file containing the graph.

    Returns:
    - NetworkX graph loaded from the file.
    """

    # Load the graph from a GraphML file
    graph = ox.load_graphml(graph_file_path)
    return graph


def process_graph(graph):
    """
    Ensure the graph is projected to the correct CRS (EPSG:28992).

    Parameters:
    - graph: NetworkX graph to process.

    Returns:
    - Projected graph in EPSG:28992 coordinate reference system.
    """

    # Ensure the graph is in the correct CRS (EPSG:28992)
    if graph.graph['crs'] != 'epsg:28992':
        graph = ox.project_graph(graph, to_crs='EPSG:28992')
    return graph


def load_cool_place_polygons(polygon_path):
    """
    Load polygons representing cool places from a shapefile or GeoJSON.

    Parameters:
    - polygon_path: Path to the file containing cool place polygons.

    Returns:
    - GeoDataFrame of cool place polygons.
    """

    # Load polygons representing cool places from a shapefile or GeoJSON
    polygons = gpd.read_file(polygon_path)
    return polygons

def find_cool_place_nodes(graph, cool_place_polygons, max_distance=500, known_crs="EPSG:28992"):
    """
    Find the nearest graph nodes for each cool place polygon centroid if within max_distance.

    Parameters:
    - graph: NetworkX graph with nodes representing the walkable area.
    - cool_place_polygons: GeoDataFrame of cool place polygons.
    - max_distance: Maximum distance in meters for a node to be considered as a cool place node.
    - known_crs: The known CRS for the data to avoid CRS estimation issues.

    Returns:
    - List of unique nearest nodes to cool place centroids within max_distance.
    """
    start_time = time.time()

    # Convert graph nodes to GeoDataFrame and set known CRS
    nodes_gdf = ox.graph_to_gdfs(graph, nodes=True, edges=False)
    nodes_gdf = nodes_gdf.set_crs(known_crs, allow_override=True)

    # Ensure cool_place_polygons also has the known CRS
    cool_place_polygons = cool_place_polygons.to_crs(known_crs)

    cool_place_nodes = []

    for _, polygon in cool_place_polygons.iterrows():
        if polygon.geometry is None:
            continue

        # Handle MultiPolygon by iterating over each part
        if isinstance(polygon.geometry, MultiPolygon):
            polygons = polygon.geometry.geoms  # Extract individual polygons from MultiPolygon
        else:
            polygons = [polygon.geometry]

        for poly in polygons:
            centroid = poly.centroid

            # Find the nearest node to this centroid
            nearest_node = ox.nearest_nodes(graph, X=centroid.x, Y=centroid.y)

            # Check the distance manually since we can't use automatic estimation
            nearest_node_data = nodes_gdf.loc[nearest_node]
            nearest_node_point = Point(nearest_node_data['geometry'].x, nearest_node_data['geometry'].y)
            distance = centroid.distance(nearest_node_point)

            if distance <= max_distance:
                cool_place_nodes.append(nearest_node)
            else:
                print(f"Skipping centroid at {centroid} with nearest node distance {distance:.2f}m (greater than {max_distance}m)")

    end_time = time.time()
    duration = end_time - start_time
    print(f"Finding cool place nodes based on centroids took {duration:.2f} seconds")

    # Return the list of unique nearest nodes
    return list(set(cool_place_nodes))  # Remove duplicates, if any


def get_nearest_node(graph, point):
    """
    Find the nearest graph node to a given point.

    Parameters:
    - graph: NetworkX graph.
    - point: Point geometry to find the nearest node for.

    Returns:
    - Node ID of the nearest node to the point.
    """

    # Find the nearest graph node to a given point
    return ox.distance.nearest_nodes(graph, X=point.x, Y=point.y)


def find_nodes_within_distances(graph, cool_place_nodes, distances=[200, 300, 400, 500]):
    """
    Find nodes within specified distance thresholds from each cool place node.

    Parameters:
    - graph: NetworkX graph of the area.
    - cool_place_nodes: List of nodes identified as cool places.
    - distances: List of distance thresholds to search within (in meters).

    Returns:
    - Dictionary where keys are distances and values are sets of nodes within each distance threshold.
    """

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


def assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes):
    """
    Assign each building to a distance category based on the proximity of its nearest node to cool place nodes.

    Parameters:
    - buildings: GeoDataFrame containing building polygons.
    - graph: NetworkX graph of the area.
    - distance_nodes: Dictionary where each key is a distance threshold (e.g., 200, 300) and each value is
                         a set of nodes within that distance from a cool place node.

    Returns:
    - Updated GeoDataFrame of buildings with a 'distance_category' column indicating the closest distance category.
    """

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


def calculate_walking_shed_with_distance_categories(buildings, cool_place_nodes, graph, distances=[200, 300, 400, 500]):
    """
    Calculate walking shed with distance categories for each building based on proximity to cool places.

    Parameters:
    - buildings: GeoDataFrame of building polygons.
    - cool_place_nodes: List of nodes representing cool places in the graph.
    - graph: NetworkX graph of the area.
    - distances: List of distance thresholds to categorize nodes/buildings by proximity (default: [200, 300, 400, 500] meters).

    Returns:
    - Updated GeoDataFrame of buildings with assigned distance categories based on proximity to cool places.
    """

    # Find nodes within each distance threshold for cool place nodes
    distance_nodes = find_nodes_within_distances(graph, cool_place_nodes, distances)

    # Assign distance category to buildings based on proximity to cool place nodes
    # buildings = assign_buildings_to_distance_category(buildings, graph, distance_nodes)
    # buildings  = parallel_assign_buildings_to_distance_category(buildings, graph, distance_nodes)
    buildings = assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes)

    return buildings


def assign_building_colors(buildings):
    """
    Assign colors to buildings based on distance category.

    Parameters:
    - buildings: GeoDataFrame of buildings with distance categories.

    Returns:
    - Updated buildings GeoDataFrame with color assigned for each category.
    """

    # Define a color map for different distance categories
    color_map = {200: 'green', 300: 'blue', 400: 'yellow', 500: 'orange'}

    # Assign colors based on the distance category
    buildings['color'] = buildings['distance_category'].map(color_map)

    # Fill NaN values with a default color
    buildings['color'] = buildings['color'].fillna('red')

    return buildings


def plot_colored_walking_shed(buildings):
    """
    Plot the walking shed with color-coded buildings.

    Parameters:
    - buildings: GeoDataFrame of buildings with color codes.
    """

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


def walking_shed_calculation(graph=None, polygon_path=None, building_shapefile_path=None):
    """
    Calculate and plot a walking shed with distance categories for a given place.

    Parameters:
    - place: Name of the place to calculate the walking shed for (optional).
    - graph_file_path: Path to a pre-saved graph (optional).
    - polygon_path: Path to the cool place polygon file.
    - building_shapefile_path: Path to the building polygons file.
    """
    graph = load_graph_from_file(graph)
    graph = process_graph(graph)
    cool_place_polygons = load_cool_place_polygons(polygon_path)
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
        graph="C:/Github_synthesis/AMS/graphs_with_shade/ams_graph_with_shade_20150701_1800_cropped.graphml",
        polygon_path = "C:/Github_synthesis/AMS/cool_places_polygons/cool_places_polygons_20230816_900.shp",
        building_shapefile_path = "C:/Androniki/pythonProject1/merged_buildings.shp",
        # polygon_path="C:/pedestrian_demo_data/cool_places_polygons/cool_places_polygons_20230816_1300.shp",
        # building_shapefile_path="C:/pedestrian_demo_data/merged_buildings/merged_buildings.shp"
    )