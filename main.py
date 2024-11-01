import argparse
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'shade_calculation')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'cool_place')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'PedestrianNetwork')))

from shade_calculation.main_shade import shade_main
from cool_place.main import coolspace_main
from PedestrianNetwork.main_network import main_network


if __name__ == '__main__':
    # TODO: FILL IN DESCRIPTION
    parser = argparse.ArgumentParser(description='TO DO FILL IN THIS ')
    parser.add_argument(
        '--config_file_shade',
        # default="configuration_files/shade_config.json",
        default="example_run\config_files\shade_config.json",
        type=str,
        help='Path to the configuration file of the shade maps (default: configuration_files/shade_config.json)'
    )
    parser.add_argument(
        '--config_file_cool_spaces',
        # default="configuration_files/coolspaceConfig.json",
        default="example_run\config_files\coolspaceConfig.json" ,
        type=str,
        help='Path to the configuration file for cool spaces (default: configuration_files/coolspaceConfig.json)'
    )

    parser.add_argument(
        '--config_file_network',
        # default="configuration_files/network_config.json",
        default="example_run\config_files/network_config.json",
        type=str,
        help='Path to the configuration file for the network (default: configuration_files/network_config.json)'
    )


    parser.add_argument(
        '--run_shade',
        type=bool,
        default=False,
        help='Set to True to run shade calculation, False to skip (default: True)'
    )
    parser.add_argument(
        '--run_cool_spaces',
        type=bool,
        default=False,
        help='Set to True to run cool spaces calculation, False to skip (default: True)'
    )
    parser.add_argument(
        '--run_network',
        type=bool,
        default=True,
        help='Set to True to run pedestrian network creation, False to skip (default: True)'
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

    ########################################### Run selected functions ######################################
    if args.run_shade:
        print("Running shade calculation...")
        shade_main(config_file_shade)

    if args.run_cool_spaces:
        print("Running cool spaces calculation...")
        coolspace_main(config_file_cool_spaces)

    if args.run_network:
        print("Running pedestrian network creation...")
        main_network(config_file_network)