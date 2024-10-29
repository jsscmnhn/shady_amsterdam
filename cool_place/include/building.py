import geopandas as gpd
import matplotlib.pyplot as plt


class Building:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.data["buffered"] = None

    def create_buffer(self, buffer_size: int) -> None:
        self.data["buffered"] = None
        self.data["buffered"] = self.data['geometry'].buffer(buffer_size)
        if 'buffered' in self.data.columns:
            print("Building buffer geometries created.")
        else:
            print("Fail to create building buffer geometries.")


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"

    building_file = directory_win + "ams_buildings_bagplus.shp"

    building = Building(gpd.read_file(building_file))
    building.create_buffer(4)
    building.data['buffered'].plot()
    plt.show()

