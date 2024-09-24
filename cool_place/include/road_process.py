import geopandas as gpd


class Road:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data

    def assign_buffer_hardSurface(self, roadtype: str) -> int:
        if roadtype is None:
            return 0
        elif "> 7 meter" in roadtype:
            return 10
        elif "4 - 7 meter" in roadtype:
            return 7
        elif "2 - 4 meter" in roadtype:
            return 4

    def assign_buffer_roadtype(self, roadtype: str) -> int:
        if roadtype is None:
            return 0
        elif "autosnelweg" in roadtype:
            return 10
        elif "hoofdweg" in roadtype:
            return 8
        elif "lokale weg" in roadtype:
            return 6
        elif "overig" in roadtype:
            return 4
        elif "parkeerplaats" in roadtype:
            return 2
        elif "regionale weg" in roadtype:
            return 2
        elif "straat" in roadtype:
            return 0

    def create_attribute(self, attri_in: str, new_attri: str) -> None:
        if attri_in == 'verharding':
            self.data[new_attri] = self.data[attri_in].apply(self.assign_buffer_hardSurface)
        if attri_in == 'typeweg':
            self.data[new_attri] = self.data[attri_in].apply(self.assign_buffer_roadtype)


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"
    filename = "ams_roads_top10NL.shp"
    filepath_mac = directory_mac + filename
    filepath_win = directory_win + filename

    road = Road(gpd.read_file(filepath_mac))
    road.create_attribute('verharding', 'buffer')
    # road.data.to_file(directory_mac + "ams_roads_top10NL_buffered.shp")
    print(road.data.head(5))
