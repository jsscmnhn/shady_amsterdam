import routes_calculation
import walking_shed_network
import os
import json
from datetime import datetime

# Instructions for users on what to include in the configuration JSON file
print("Please provide a configuration JSON file with the following data:")
print("1) 'graph_dir': Path to the directory containing GraphML files for the network in your study area, with shade weights.")
print("2) 'nodes_dir': Path to the directory where the cool place nodes are stored in .pkl files.")
print("3) Optional - 'graph_file': Specify an exact GraphML file (including shade weights) if you want to bypass timestamp-based selection.")
print("4) Optional - 'nodes_file': Specify an exact .pkl file with pre-calculated cool place nodes if you want to bypass timestamp-based selection.")
print("5) 'route option': Route type for calculation; use either 'origin_destination' for a specified destination or 'nearest_cool_place' for the closest cool place.")
print("6) 'location indication option': Specify 'coordinates' or 'location_name' to indicate if you are using exact coordinates or place names.")
print("7) If using 'location_name':")
print("    - 'origin name': Name of the origin location (e.g., 'Amsterdam Central Station').")
print("    - 'destination name': Name of the destination location (if using 'origin_destination').")
print("8) If using 'coordinates':")
print("    - 'origin latitude' and 'origin longitude': Coordinates for the origin location.")
print("    - 'destination latitude' and 'destination longitude': Coordinates for the destination (if using 'origin_destination').")
print("9) 'date time': (Optional) A specific date and time for routing, formatted as 'YYYY-MM-DD HH:MM:SS'. If omitted or set to 'None', the current date and time will be used.")


def load_config(config_path="C:/Users/17731/Documents/GitHub/shady_amsterdam/PedestrianNetwork/config_map_routing.json"):
    if os.path.exists(config_path):
        print("Loading config...")
        with open(config_path, 'r') as file:
            config = json.load(file)
    else:
        raise FileNotFoundError(f"Configuration file not found at {config_path}.")
    return config

if __name__ == "__main__":
    # Configuration
    config = load_config()

    # Extract configuration values
    graph_dir = config.get('graph_dir')
    nodes_dir = config.get('nodes_dir')
    graph_file = config.get('graph_file')  # New optional field
    nodes_file = config.get('nodes_file')  # New optional field
    route_option = config.get('route option')
    location_indication_option = config.get('location indication option')
    origin_name = config.get('origin name')
    destination_name = config.get('destination name')
    origin_latitude = config.get('origin latitude')
    origin_longitude = config.get('origin longitude')
    destination_latitude = config.get('destination latitude')
    destination_longitude = config.get('destination longitude')
    date_time_str = config.get('date time')

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