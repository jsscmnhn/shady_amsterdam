import pickle
import osmnx as ox
import networkx as nx
import geopandas as gpd
import os
import re
import fiona


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


def export_layers_to_shapefiles(geopackage_path, output_directory, date_str):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Get all layer names from the GeoPackage
    layer_names = fiona.listlayers(geopackage_path)

    for layer_name in layer_names:
        # Match the layer name pattern to extract the time index
        match = re.search(r'sdGeom(\d+)', layer_name)
        if match:
            # Convert the index to time format (e.g., sdGeom0 -> 900, sdGeom1 -> 930)
            time_index = int(match.group(1))
            time_str = f"{9 + time_index // 2}{30 if time_index % 2 else '00'}"

            # Construct the output shapefile path
            output_path = os.path.join(output_directory, f'cool_places_polygons_{date_str}_{time_str}.shp')

            # Read the layer and save as a shapefile
            gdf = gpd.read_file(geopackage_path, layer=layer_name)
            gdf = gdf[['geometry']].to_crs('EPSG:28992')  # Keep only geometry and reproject
            gdf.to_file(output_path, driver='ESRI Shapefile')
            print(f"Exported {layer_name} to {output_path}")


def process_all_shapefiles(polygon_directory, graph_directory, output_directory):
    # Loop through each shapefile in the directory with the specified pattern
    for filename in os.listdir(polygon_directory):
        # Check if the filename matches the required pattern
        if filename.endswith('.shp') and re.match(r'cool_places_polygons_\d{8}_\d{3,4}\.shp', filename):
            print(f"Processing shapefile: {filename}")  # Debug statement to confirm each file being processed

            # Extract date and time from filename
            date_time_match = re.search(r'cool_places_polygons_(\d{8})_(\d{3,4})', filename)
            if not date_time_match:
                print(f"Filename {filename} does not match expected date and time pattern. Skipping...")
                continue

            date, time = date_time_match.groups()
            print(f"Extracted date: {date}, time: {time}")  # Debug statement for date and time extraction

            # Construct file paths for the polygon and corresponding graph
            polygon_path = os.path.join(polygon_directory, filename)
            graph_filename = f'ams_graph_with_shade_{date}_{time}_cropped.graphml'
            graph_file_path = os.path.join(graph_directory, graph_filename)
            output_path = os.path.join(output_directory, f'cool_places_nodes_{date}_{time}.pkl')

            # Check if the corresponding graph file exists
            if not os.path.exists(graph_file_path):
                print(f"Graph file {graph_filename} not found for shapefile {filename}. Skipping...")
                continue

            try:
                # Load the graph
                graph = load_graph_from_file(graph_file_path)
                print(f"Loaded graph from {graph_file_path}")  # Confirm graph loaded

                # Project graph to correct CRS
                graph = ox.project_graph(graph, to_crs='EPSG:28992')

                # Load polygons
                cool_place_polygons = load_cool_place_polygons(polygon_path)
                print(f"Loaded polygons from {polygon_path}")  # Confirm polygons loaded

                # Calculate and save cool place nodes
                calculate_and_save_cool_place_nodes(graph, cool_place_polygons, output_path)
                print(f"Cool place nodes saved to {output_path}")  # Confirm nodes saved

            except Exception as e:
                print(f"An error occurred while processing {filename}: {e}")


# polygon_path = 'C:/pedestrian_demo_data/public_spaces/ams_public_space.shp'
# graph_file_path = 'C:/pedestrian_demo_data/ams_graph_with_shade_900_cropped.graphml'

# graph = load_graph_from_file(graph_file_path)
# graph = ox.project_graph(graph, to_crs='EPSG:28992')
# cool_place_polygons = load_cool_place_polygons(polygon_path)
# calculate_and_save_cool_place_nodes(graph, cool_place_polygons, pre_calculated_nodes_path)

# Example use
polygon_directory = 'C:/pedestrian_demo_data/cool_places_polygons/'
graph_directory = 'C:/pedestrian_demo_data/graphs_with_shade/'
output_directory = 'C:/pedestrian_demo_data/cool_place_nodes/'
process_all_shapefiles(polygon_directory,
                       graph_directory,
                       output_directory)

# # From gpkg to shapefiles
# geopackage_path = 'C:/pedestrian_demo_data/shadeGeoms.gpkg'
# output_directory = 'C:/pedestrian_demo_data/cool_places_polygons/'
# date_str = '20230816'  # Specify the date string for filenames
# export_layers_to_shapefiles(geopackage_path, output_directory, date_str)