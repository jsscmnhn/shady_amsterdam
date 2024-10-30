import geopandas as gpd
import networkx as nx
import osmnx as ox
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import time
from osmnx.distance import nearest_nodes
from shapely.geometry import MultiPolygon, Polygon, Point
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.spatial import KDTree


def load_building_polygons(building_shapefile_path):
    buildings = gpd.read_file(building_shapefile_path)
    return buildings


def load_graph_from_file(graph_file_path):
    graph = ox.load_graphml(graph_file_path)
    return graph


def process_graph(graph):
    if graph.graph['crs'] != 'epsg:28992':
        graph = ox.project_graph(graph, to_crs='EPSG:28992')
    return graph


def load_cool_place_polygons(polygon_path):
    polygons = gpd.read_file(polygon_path)
    return polygons


def find_cool_place_nodes(graph, cool_place_polygons, max_distance=50, known_crs="EPSG:28992"):
    start_time = time.time()

    nodes_gdf = ox.graph_to_gdfs(graph, nodes=True, edges=False)
    nodes_gdf = nodes_gdf.set_crs(known_crs, allow_override=True)
    cool_place_polygons = cool_place_polygons.to_crs(known_crs).dropna(subset=['geometry'])

    # Create KDTree for quick distance lookup
    node_coords = [(geom.x, geom.y) for geom in nodes_gdf.geometry]
    kd_tree = KDTree(node_coords)

    def get_extent_points(poly):
        bounds = poly.bounds
        return [Point(bounds[2], bounds[3]), Point(bounds[0], bounds[3]), Point(bounds[2], bounds[1]), Point(bounds[0], bounds[1])]

    def process_polygon_extents(feature):
        geometry = feature.geometry
        extent_points = []

        if isinstance(geometry, MultiPolygon):
            for poly in geometry.geoms:
                extent_points.extend(get_extent_points(poly))
        elif isinstance(geometry, Polygon):
            extent_points = get_extent_points(geometry)

        # Perform the nearest node lookup
        query_points = [(point.x, point.y) for point in extent_points]
        distances, indices = kd_tree.query(query_points)

        # Collect nodes if within the max distance
        nearby_nodes = [
            nodes_gdf.index[nearest_idx] for distance, nearest_idx in zip(distances, indices) if distance <= max_distance
        ]
        return nearby_nodes

    cool_place_nodes = []
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_polygon_extents, feature) for _, feature in cool_place_polygons.iterrows()]
        for future in as_completed(futures):
            cool_place_nodes.extend(future.result())

    print(f"Finding cool place nodes took {time.time() - start_time:.2f} seconds")
    return list(set(cool_place_nodes))


def find_nodes_within_distances(graph, cool_place_nodes, distances=[200, 300, 400, 500], weight="length"):
    start_time = time.time()
    print("Starting to find nodes within specified distance thresholds...")

    distance_nodes = {distance: set() for distance in distances}

    for i, cool_place_node in enumerate(cool_place_nodes):
        if i % 10 == 0:
            print(f"Processing cool place node {i+1}/{len(cool_place_nodes)}")

        for distance in distances:
            lengths = nx.single_source_dijkstra_path_length(graph, cool_place_node, cutoff=distance, weight=weight)
            distance_nodes[distance].update(lengths.keys())
            print(f"  Found {len(lengths)} nodes within {distance} meters of cool place node {cool_place_node}")

    print(f"Finding nodes within distance took {time.time() - start_time:.2f} seconds")
    return distance_nodes


def assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes):
    start_time = time.time()
    print("Starting to assign distance categories to buildings...")

    building_centroids = buildings.geometry.centroid
    nearest_nodes = ox.distance.nearest_nodes(graph, X=building_centroids.x, Y=building_centroids.y)
    buildings['nearest_node'] = nearest_nodes

    def get_distance_category(node_id):
        for distance in sorted(distance_nodes.keys()):
            if node_id in distance_nodes[distance]:
                return distance
        return '>500m'

    buildings['distance_category'] = buildings['nearest_node'].apply(get_distance_category)

    print(f"Assigning buildings to distance category took {time.time() - start_time:.2f} seconds")
    return buildings


def calculate_walking_shed_with_distance_categories(buildings, cool_place_nodes, graph, distances=[200, 300, 400, 500], weight="length"):
    distance_nodes = find_nodes_within_distances(graph, cool_place_nodes, distances, weight)
    buildings = assign_buildings_to_distance_category_with_precomputed(buildings, graph, distance_nodes)
    return buildings


def assign_building_colors(buildings):
    color_map = {200: 'green', 300: 'blue', 400: 'yellow', 500: 'orange'}
    buildings['color'] = buildings['distance_category'].map(color_map)
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


def preprocess_graph_weights(graph, weight):
    for u, v, data in graph.edges(data=True):
        if weight in data:
            try:
                data[weight] = float(data[weight])
            except (ValueError, TypeError):
                print(f"Warning: Edge ({u}, {v}) has an invalid '{weight}' value: {data[weight]}, setting to infinity.")
                data[weight] = float('inf')
        else:
            data[weight] = float('inf')


def walking_shed_calculation(graph=None, polygon_path=None, building_shapefile_path=None, weight="length"):
    graph = load_graph_from_file(graph)
    graph = process_graph(graph)
    preprocess_graph_weights(graph, weight)

    cool_place_polygons = load_cool_place_polygons(polygon_path)
    cool_place_nodes = find_cool_place_nodes(graph, cool_place_polygons)
    print(f"Found {len(cool_place_nodes)} cool place nodes on the graph.")

    buildings = load_building_polygons(building_shapefile_path)
    print("Calculating walking shed with distance categories...")
    buildings = calculate_walking_shed_with_distance_categories(buildings, cool_place_nodes, graph, weight=weight)

    print("Assigning colors to the buildings...")
    buildings = assign_building_colors(buildings)

    plot_colored_walking_shed(buildings)


if __name__ == "__main__":
    walking_shed_calculation(
        graph="C:/Github_synthesis/AMS/graphs_with_shade/ams_graph_with_shade_20150701_1800_cropped.graphml",
        polygon_path="C:/Github_synthesis/AMS/cool_places_polygons/cool_places_polygons_20230816_900.shp",
        building_shapefile_path="C:/Androniki/pythonProject1/merged_buildings.shp",
        weight="shade_weight"
    )
