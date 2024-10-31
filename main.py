import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'shade_calculation')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'cool_place')))

from shade_calculation.main_shade import shade_main
# from cool_place.main import coolspace_main



if __name__ == '__main__':
    # TODO: FILL IN DESCRIPTION
    parser = argparse.ArgumentParser(description='TO DO FILL IN THIS ')
    parser.add_argument(
        'config_file_shade',
        nargs='?',
        # default="configuration_files/shade_config.json",
        default="example_run\config_files\shade_config.json",
        type=str,
        help='Path to the configuration file of the shade maps (default: configuration_files/shade_config.json)'
    )
    parser.add_argument(
        'config_file_cool_spaces',
        nargs='?',
        # default="configuration_files/coolspaceConfig.json",
        default="example_run\config_files\coolspaceConfig.json" ,
        type=str,
        help='Path to the configuration file for cool spaces (default: configuration_files/coolspaceConfig.json)'
    )

    #TODO: FILL IN FOR PEDESTRIAN NETWORK
    parser.add_argument(
        'config_file_network',
        nargs='?',
        default="configuration_files/coolspaceConfig.json",
        type=str,
        help='Path to the configuration file for the network (default: )'
    )
    args = parser.parse_args()

    config_file_shade = args.config_file_shade
    config_file_cool_spaces = args.config_file_cool_spaces
    config_file_network = args.config_file_network

    # Check if the config file exists
    if not os.path.exists(config_file_shade):
        raise FileNotFoundError(f"Config file not found: {config_file_shade}")

    if not os.path.exists(config_file_cool_spaces):
        raise FileNotFoundError(f"Config file not found: {config_file_cool_spaces}")

    if not os.path.exists(config_file_network):
        raise FileNotFoundError(f"Config file not found: {config_file_network}")

########################################### Functions for creating the shade maps ######################################
    shade_main(str(config_file_shade))

########################################### Functions for creating the cool spaces #####################################
   # coolspace_main(config_file_cool_spaces)

########################################### Functions for creating the pedestrian network ##############################