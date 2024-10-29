import routes_calculation
import walking_shed_network
import os
import json

print('Provide the configuration json file. In this file the following data should be included:')
print('1)The path of the graphml file of the network in your study area. It should contain shade weights.')
print('2)The path of the file in which you store the cool place polygons in a shapefile.')
print('3)The path of the shapefile of the buildings in your study area')

print('4)The path of the pkl file of the cool places nodes in your study area')
print('5)Your option regarding the type of the route you wish to examine between a location indicated as the origin and another indicated as the destination or between the origin and the nearest coolest place. Please in the json configuration file type one of these two options:origin_destination or nearest_cool_place')
print('6)The coordinates of the origin and the destination of the route or their names to implement geocoding.')
print('7)Your option regarding indicating these two locations. Please in the json configuration file type one of these two options:: coordinates or location_name')



def load_config(config_path="C:/Github_synthesis/configuration/config_map_routing.json"):
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

    # Load configuration values with default handling
    graphml_file = config.get('graphml_file')
    cool_places_polygons = config.get('cool places path')
    buildings_path = config.get('buildings')

    cool_places_nodes = config.get('cool places nodes path')
    route_option = config.get('route option')
    location_indication_option = config.get('location indication option')
    origin_name = config.get('origin name')
    destination_name = config.get('destination name')
    origin_latitude = config.get('origin latitude')
    origin_longitude = config.get('origin longitude')
    destination_latitude = config.get('destination latitude')
    destination_longitude = config.get('destination longitude')
    date_time = config.get('date time')




    walking_shed_network.walking_shed_calculation(graph=graphml_file, polygon_path=cool_places_polygons, building_shapefile_path=buildings_path)
    # check polygon_path again
    if location_indication_option == 'location_name':
        if route_option == 'origin_destination':
            routes_calculation.demo_shade_route_calculation_with_timedemo_shade_route_calculation_with_time(graph_dir=graphml_file, nodes_dir=cool_places_nodes, user_input=(origin_name,destination_name), input_type=location_indication_option, mode=route_option)
        else:
            routes_calculation.demo_shade_route_calculation_with_timedemo_shade_route_calculation_with_time(graph_dir=graphml_file, nodes_dir=cool_places_nodes, user_input=(origin_name),input_type=location_indication_option, mode=route_option)
    else:
        if route_option == 'origin_destination':
            routes_calculation.demo_shade_route_calculation_with_timedemo_shade_route_calculation_with_time(graph_dir=graphml_file, nodes_dir=cool_places_nodes, user_input=(origin_latitude, origin_longitude, destination_latitude, destination_longitude), input_type=location_indication_option, mode=route_option)
        else:
            routes_calculation.demo_shade_route_calculation_with_timedemo_shade_route_calculation_with_time(graph_dir=graphml_file, nodes_dir=cool_places_nodes, user_input=(origin_latitude,origin_longitude), input_type=location_indication_option, mode=route_option)



