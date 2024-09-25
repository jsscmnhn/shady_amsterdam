from include.road_process import Road
from include.cool_space import CoolSpace
from include.building import Building
import geopandas as gpd
import matplotlib.pyplot as plt


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"

    landuse_file = directory_mac + "ams_landuse_top10NL.shp"
    road_file = directory_mac + "ams_roads_top10NL.shp"
    building_file = directory_mac + "ams_buildings_bagplus.shp"

    coolSpace = CoolSpace(gpd.read_file(landuse_file))
    road = Road(gpd.read_file(road_file))
    building = Building(gpd.read_file(building_file))

    road.create_attribute('typeweg', 'buffer')
    road.create_buffer('buffer')

    building.create_buffer(4)

    coolSpace.clip(coolSpace.data, road.buffered)
    coolSpace.clip(coolSpace.clipped, building.buffered)

    coolSpace.clipped.to_file(directory_mac + "ams_public_space.shp")
    # coolSpace.clipped.plot()
    # plt.show()
