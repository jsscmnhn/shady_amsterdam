import pickle
import osmnx as ox
import networkx as nx
import geopandas as gpd
import os
import re
import fiona
from shapely.geometry import MultiPolygon, Point, Polygon
from scipy.spatial import KDTree
from concurrent.futures import ThreadPoolExecutor, as_completed
import pickle
import time


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

    # Identify and fix invalid geometries
    if not polygons.is_valid.all():
        print("Some geometries are invalid. Attempting to fix them...")
        polygons['geometry'] = polygons['geometry'].buffer(0)

        # Check validity again after repair
        invalid_geometries = polygons[~polygons.is_valid]
        if not invalid_geometries.empty:
            print(f"Unable to fix {len(invalid_geometries)} geometries. They will be skipped.")

    # Filter out any remaining invalid geometries if necessary
    polygons = polygons[polygons.is_valid]

    return polygons


def calculate_and_save_cool_place_nodes(graph, cool_place_polygons, output_path, max_distance=50, known_crs="EPSG:28992"):
    start_time = time.time()

    # Convert graph nodes to GeoDataFrame
    nodes_gdf = ox.graph_to_gdfs(graph, nodes=True, edges=False)
    nodes_gdf = nodes_gdf.set_crs(known_crs, allow_override=True)

    # Ensure cool place polygons are in the same CRS
    cool_place_polygons = cool_place_polygons.to_crs(known_crs).dropna(subset=['geometry'])

    # Create KDTree for efficient distance lookup
    node_coords = [(geom.x, geom.y) for geom in nodes_gdf.geometry]
    kd_tree = KDTree(node_coords)

    def get_extent_points(poly):
        # Get the corner points of the bounding box of the polygon
        bounds = poly.bounds
        return [
            Point(bounds[2], bounds[3]),  # Top-right
            Point(bounds[0], bounds[3]),  # Top-left
            Point(bounds[2], bounds[1]),  # Bottom-right
            Point(bounds[0], bounds[1])   # Bottom-left
        ]

    def process_polygon_extents(feature):
        geometry = feature.geometry
        extent_points = []

        if isinstance(geometry, MultiPolygon):
            for poly in geometry.geoms:
                extent_points.extend(get_extent_points(poly))
        elif isinstance(geometry, Polygon):
            extent_points = get_extent_points(geometry)

        # Perform nearest node lookup for the extent points
        query_points = [(point.x, point.y) for point in extent_points]
        distances, indices = kd_tree.query(query_points)

        # Collect nodes that are within the maximum distance
        nearby_nodes = [
            nodes_gdf.index[nearest_idx] for distance, nearest_idx in zip(distances, indices) if distance <= max_distance
        ]
        return nearby_nodes

    cool_place_nodes = []
    # Use ThreadPoolExecutor for parallel processing of polygons
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_polygon_extents, feature) for _, feature in cool_place_polygons.iterrows()]
        for future in as_completed(futures):
            cool_place_nodes.extend(future.result())

    # Remove duplicates from cool place nodes
    cool_place_nodes = list(set(cool_place_nodes))

    # Save the cool place nodes to a file
    with open(output_path, 'wb') as f:
        pickle.dump(cool_place_nodes, f)

    print(f"Cool place nodes calculated and saved to {output_path}.")
    print(f"Finding cool place nodes took {time.time() - start_time:.2f} seconds")

    return cool_place_nodes


def export_layers_to_shapefiles(geopackage_path, output_directory, date_str):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    layer_names = fiona.listlayers(geopackage_path)

    for layer_name in layer_names:
        match = re.search(r'sdGeom(\d+)', layer_name)
        if match:
            time_index = int(match.group(1))
            time_str = f"{9 + time_index // 2}{30 if time_index % 2 else '00'}"
            output_path = os.path.join(output_directory, f'cool_places_polygons_{date_str}_{time_str}.shp')
            gdf = gpd.read_file(geopackage_path, layer=layer_name)[['geometry']].to_crs('EPSG:28992')
            gdf.to_file(output_path, driver='ESRI Shapefile')
            print(f"Exported {layer_name} to {output_path}")


import os
import re

def process_all_geopackages_in_directory(gpkg_directory, graph_directory, shapefile_output_directory, output_directory):
    # List all .gpkg files in the specified directory
    gpkg_files = [f for f in os.listdir(gpkg_directory) if f.endswith('.gpkg')]

    if not gpkg_files:
        print("No GeoPackage files found in the specified directory.")
        return

    # Ensure the shapefile output directory exists
    if not os.path.exists(shapefile_output_directory):
        os.makedirs(shapefile_output_directory)

    for gpkg_file in gpkg_files:
        geopackage_path = os.path.join(gpkg_directory, gpkg_file)

        # Extract date from GeoPackage filename
        match = re.search(r'shadeGeoms_(\d{8})', gpkg_file)
        if not match:
            print(f"Date not found in filename {gpkg_file}. Skipping this file.")
            continue

        date_str = match.group(1)

        # Check if all potential .pkl output files for this GeoPackage date already exist
        all_layers_exported = True
        for layer_index in range(19):  # Assuming up to 19 layers (0-18)
            time_suffix = f"{9 + layer_index // 2}{30 if layer_index % 2 else '00'}"
            output_path = os.path.join(output_directory, f'cool_places_nodes_{date_str}_{time_suffix}.pkl')
            if not os.path.exists(output_path):
                all_layers_exported = False
                break

        # Skip this GeoPackage if all .pkl files for its layers already exist
        if all_layers_exported:
            print(f"All cool places nodes files for date {date_str} already exist. Skipping {gpkg_file}.")
            continue

        # Export layers to shapefiles in the shared shapefile output directory
        export_layers_to_shapefiles(geopackage_path, shapefile_output_directory, date_str)

        # Process each shapefile in the shapefile output directory
        for filename in os.listdir(shapefile_output_directory):
            # Check if filename matches the required shapefile naming pattern
            if filename.endswith('.shp') and re.match(r'cool_places_polygons_\d{8}_\d{3,4}\.shp', filename):
                print(f"Processing shapefile: {filename}")
                date_time_match = re.search(r'cool_places_polygons_(\d{8})_(\d{3,4})', filename)
                if not date_time_match:
                    print(f"Filename {filename} does not match expected date and time pattern. Skipping...")
                    continue

                date, time = date_time_match.groups()
                polygon_path = os.path.join(shapefile_output_directory, filename)
                graph_filename = f'ams_graph_with_shade_{date}_{time}_cropped.graphml'
                graph_file_path = os.path.join(graph_directory, graph_filename)
                output_path = os.path.join(output_directory, f'cool_places_nodes_{date}_{time}.pkl')

                # Skip processing if the .pkl file already exists
                if os.path.exists(output_path):
                    print(f"Output file {output_path} already exists. Skipping processing for this file.")
                    continue

                if not os.path.exists(graph_file_path):
                    print(f"Graph file {graph_filename} not found for shapefile {filename}. Skipping...")
                    continue

                try:
                    graph = load_graph_from_file(graph_file_path)
                    print(f"Loaded graph from {graph_file_path}")
                    graph = ox.project_graph(graph, to_crs='EPSG:28992')
                    cool_place_polygons = load_cool_place_polygons(polygon_path)
                    print(f"Loaded polygons from {polygon_path}")
                    calculate_and_save_cool_place_nodes(graph, cool_place_polygons, output_path)
                    print(f"Cool place nodes saved to {output_path}")

                except Exception as e:
                    print(f"An error occurred while processing {filename}: {e}")