from include.road_process import Road
from include.cool_space import CoolSpace
import geopandas as gpd
import matplotlib.pyplot as plt


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"

    landuse_file = directory_win + "ams_landuse_top10NL.shp"
    road_file = directory_win + "ams_roads_top10NL.shp"

    coolSpace = CoolSpace(gpd.read_file(landuse_file))
    road = Road(gpd.read_file(road_file))

    road.create_attribute('verharding', 'buffer')
    road.create_buffer('buffer')
    coolSpace.clip(road.buffered)

    coolSpace.is_qualified(10, 0.5)

    # coolSpace.clipped.to_file(directory_win + "test_coolSpace.shp")
    # coolSpace.clipped.plot()
    # plt.show()
