import shade_weight_calculation
import cool_places_nodes_calculation
import routes_calculation
import walking_shed_network
import os
import json
import osmnx as ox
from datetime import datetime


def load_config(config_path="network_config.json"):
    if os.path.exists(config_path):
        print("Loading combined configuration...")
        with open(config_path, 'r') as file:
            config = json.load(file)
        print("Configuration loaded successfully.")
    else:
        raise FileNotFoundError(f"Configuration file not found at {config_path}.")
    return config


def check_required_inputs(config_section, required_keys):
    empty_keys = [key for key in required_keys if not config_section.get(key)]

    # If all keys are empty or None, skip the section
    if len(empty_keys) == len(required_keys):
        print("Skipping section: All required inputs are either None or empty.")
        return False
    return True


def dataset_preparation(config):
    dataset_config = config["Dataset_Preparation"]
    required_inputs = ["graphml_file", "area_name_for_OSM_network", "cool_places_gpkg_dir", "cool_places_shp_dir", "cool_places_nodes_path"]

    # Check if dataset preparation is needed based on required inputs
    if not check_required_inputs(dataset_config, required_inputs):
        return

    # Extract values
    graphml_file = dataset_config.get('graphml_file', None)
    area_name = dataset_config.get('area_name_for_OSM_network', None)
    shade_maps_path = dataset_config.get('shade_maps_path', None)
    new_graphml_files_path = dataset_config.get('new_graphml_files_path', None)
    cool_places_gpkg_dir = dataset_config.get('cool_places_gpkg_dir')
    cool_places_shp_dir = dataset_config.get('cool_places_shp_dir')
    cool_places_nodes_path = dataset_config.get('cool_places_nodes_path')

    # Determine network configuration
    if graphml_file:
        network_WGS84 = ox.load_graphml(graphml_file)
        network = ox.project_graph(network_WGS84, to_crs='EPSG:28992')
        print(f"Using provided GraphML file: {graphml_file}")
    elif area_name:
        network_WGS84 = ox.graph_from_place(area_name, network_type="walk")
        network = ox.project_graph(network_WGS84, to_crs='EPSG:28992')
        print(f"Using OSM network for area: {area_name}")
    else:
        raise ValueError("Please provide either 'graphml_file' or 'area_name_for_OSM_network'.")

    # Calculate shade weights if paths are provided
    if shade_maps_path and new_graphml_files_path:
        print("Calculating shade weights...")
        shade_weight_calculation.process_multiple_shade_maps(
            graph_file=network, raster_dir=shade_maps_path, output_dir=new_graphml_files_path)
    else:
        print("Skipping shade weight calculation as paths are not provided.")

    # Calculate cool places nodes
    print("Calculating cool places nodes...")
    cool_places_nodes_calculation.process_all_geopackages_in_directory(
        gpkg_directory=cool_places_gpkg_dir,
        graph_directory=new_graphml_files_path if new_graphml_files_path else '',
        shapefile_output_directory=cool_places_shp_dir,
        output_directory=cool_places_nodes_path)


def routing(config):
    routing_config = config["Routing"]
    required_inputs = ["graph_dir", "nodes_dir", "graph_file", "nodes_file", "route_option", "location_indication_option", "origin_name",
                       "destination_name", "origin_latitude", "origin_longitude", "destination_latitude", "destination_longitude", "date_time"]

    # Check if routing is needed based on required inputs
    if not check_required_inputs(routing_config, required_inputs):
        return

    # Extract configuration values
    graph_dir = routing_config.get('graph_dir')
    nodes_dir = routing_config.get('nodes_dir')
    graph_file = routing_config.get('graph_file')  # New optional field
    nodes_file = routing_config.get('nodes_file')  # New optional field
    route_option = routing_config.get('route_option')
    location_indication_option = routing_config.get('location_indication_option')
    origin_name = routing_config.get('origin_name')
    destination_name = routing_config.get('destination_name')
    origin_latitude = routing_config.get('origin_latitude')
    origin_longitude = routing_config.get('origin_longitude')
    destination_latitude = routing_config.get('destination_latitude')
    destination_longitude = routing_config.get('destination_longitude')
    date_time_str = routing_config.get('date_time')

    # Convert date_time to datetime if specified
    date_time = None
    if date_time_str and date_time_str != "None":
        date_time = datetime.strptime(date_time_str, "%Y-%m-%d %H:%M:%S")

    # Determine user_input based on location indication
    if location_indication_option == 'location_name':
        user_input = (origin_name, destination_name) if route_option == 'origin_destination' else (origin_name,)
    else:
        if route_option == 'origin_destination':
            user_input = (origin_latitude, origin_longitude, destination_latitude, destination_longitude)
        else:
            user_input = (origin_latitude, origin_longitude)

    # Decide which function to call based on whether specific files are provided
    if graph_file and nodes_file:
        # Use specific files
        print("Using specified graph and nodes files for route calculation.")
        routes_calculation.demo_shade_route_calculation(
            graph_file_path=graph_file,
            pre_calculated_nodes_path=nodes_file,
            user_input=user_input,
            input_type=location_indication_option,
            mode=route_option
        )
    else:
        # Use directories to find nearest timestamped files
        print("Using directory search to find nearest timestamped graph and nodes files.")
        routes_calculation.demo_shade_route_calculation_with_time(
            graph_dir=graph_dir,
            nodes_dir=nodes_dir,
            user_input=user_input,
            input_type=location_indication_option,
            mode=route_option,
            date_time=date_time
        )


def walking_shed(config):
    walking_shed_config = config["Walking_Shed"]
    required_inputs = ["graph_file", "cool_places_path", "building_shapefile_path", "weight", "output_building_shapefile", "output_cool_place_shapefile"]

    # Check if walking shed is needed based on required inputs
    if not check_required_inputs(walking_shed_config, required_inputs):
        return

    # Extract configuration values
    graph_file = walking_shed_config.get('graph_file')
    cool_places_path = walking_shed_config.get('cool_places_path')
    building_shapefile_path = walking_shed_config.get('building_shapefile_path')
    weight = walking_shed_config.get('weight')
    output_building_shapefile = walking_shed_config.get('output_building_shapefile')
    output_cool_place_shapefile = walking_shed_config.get('output_cool_place_shapefile')

    # Perform walking shed calculations
    print("Calculating walking shed network...")
    walking_shed_network.walking_shed_calculation(
        graph=graph_file,
        polygon_path=cool_places_path,
        building_shapefile_path=building_shapefile_path,
        weight=weight,
        output_building_shapefile=output_building_shapefile,
        output_cool_place_shapefile=output_cool_place_shapefile
    )


if __name__ == "__main__":
    config_path = "C:/Users/17731/Documents/GitHub/shady_amsterdam/PedestrianNetwork/network_config.json"

    # Load the combined configuration
    config = load_config(config_path)

    # Run dataset preparation if necessary
    print("Checking dataset preparation requirements...")
    dataset_preparation(config)

    # Run routing if necessary
    print("Checking routing requirements...")
    routing(config)

    # Run walking shed calculation if necessary
    print("Checking walking shed requirements...")
    walking_shed(config)