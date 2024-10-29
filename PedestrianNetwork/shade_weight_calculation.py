import osmnx as ox
import networkx as nx
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import time
import os
import re




def add_edge_identifiers_to_gdf(graph, edges_gdf):
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

    return edges_gdf


def compute_shade_weights(edges, raster, affine):
    start_time = time.time()

    # Compute zonal statistics for shade weights on each edge
    shade_stats = zonal_stats(
        vectors=edges['geometry'],
        raster=raster.read(1),
        affine=affine,
        stats=['sum', 'count'],  # Get sum and total pixel count
        all_touched=True
    )

    # Calculate the proportion of shaded pixels (0 values in the raster)
    edges['shade_proportion'] = [(stat['count'] - stat['sum']) / stat['count'] if stat['count'] > 0 else 0
                                 for stat in shade_stats]

    # Inverse of shade proportion to calculate weight (more shade = lower weight)
    edges['shade_weight'] = [1 / (shade + 1e-5) if shade > 0 else 1000 for shade in edges['shade_proportion']]

    print("Shade weights calculated.")
    end_time = time.time()
    print(f"Shade weight calculation took {end_time - start_time:.2f} seconds")

    return edges


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


def store_graph_with_shade_weights(graph, edges_gdf, output_path):
    # Save the updated graph to a file
    ox.save_graphml(graph, output_path)
    print(f"Graph with shade weights saved to {output_path}")


def precalculate_and_store_shade_weights(graph, raster_path, output_graph_path):

    # Convert the graph edges to a GeoDataFrame
    _, edges = ox.graph_to_gdfs(graph)

    # Add edge identifiers (u, v, key)
    edges = add_edge_identifiers_to_gdf(graph, edges)

    # Load the raster for shade data
    print("Loading raster data...")
    with rasterio.open(raster_path) as raster:
        affine = raster.transform

        # Calculate shade weights for the edges
        print("Calculating shade weights...")
        edges_with_weights = compute_shade_weights(edges, raster, affine)

    # Update the graph with shade weights
    print("Updating graph with shade weights...")
    graph = update_graph_with_shade_weight(graph, edges_with_weights)

    # Store the graph with pre-calculated shade weights
    print("Storing graph with shade weights...")
    store_graph_with_shade_weights(graph, edges_with_weights, output_graph_path)


def generate_graph_name_from_raster(raster_filename):
    # Extract date and time from raster filename using regular expression
    match = re.search(r'amsterdam_(\d{8})_(\d{3,4})', raster_filename)
    if match:
        date_str = match.group(1)
        time_str = match.group(2)
        return f'ams_graph_with_shade_{date_str}_{time_str}_cropped.graphml'
    else:
        raise ValueError(f"Filename {raster_filename} does not match expected pattern.")

def process_multiple_shade_maps(graph_file, raster_dir, output_dir):
    # Get a list of all TIF files in the raster directory
    graph = graph_file
    raster_files = [f for f in os.listdir(raster_dir) if f.endswith('.TIF')]

    for raster_file in raster_files:
        # Construct full path to the raster file
        raster_path = os.path.join(raster_dir, raster_file)

        # Generate the output graph name based on the raster file name
        output_graph_name = generate_graph_name_from_raster(raster_file)
        output_graph_path = os.path.join(output_dir, output_graph_name)

        # Process the shade map and generate the graph with shade weights
        print(f"Processing {raster_file}...")
        precalculate_and_store_shade_weights(graph=graph, raster_path=raster_path,output_graph_path=output_graph_path)

# place = 'Amsterdam, Netherlands'
# graph_file_path = 'C:/Github_synthesis/AMS/Amsterdam_pedestrian_network.graphml'
# raster_path = 'C:/pedestrian_demo_data/amsterdam_time_900.tif'
# output_graph_path = 'C:/pedestrian_demo_data/ams_graph_with_shade_900_cropped.graphml'

# raster_dir = 'C:/Github_synthesis/AMS/shade_maps/'
# output_dir = 'C:/Github_synthesis/AMS/graphs_with_shade/'

# precalculate_and_store_shade_weights(place, raster_path, polygon_path, output_graph_path)
# process_multiple_shade_maps(main_dataset.graph_file_path, main_dataset.raster_dir, main_dataset.output_dir)