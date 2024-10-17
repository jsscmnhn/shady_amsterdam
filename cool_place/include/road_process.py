import geopandas as gpd


class Road:
    """
    Note: the assign methods and the create_attribute method only works for ONE specific dataset,
    they are only used in development. In practice, the input road dataset MUST have an attribute
    which specify the buffer distance for different types of road, so that the create_buffer method
    can be called.
    """
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.data["buffered"] = None

    def assign_buffer_hardSurface(self, roadtype: str) -> int:
        if roadtype is None:
            return 0
        elif "> 7 meter" in roadtype:
            return 20
        elif "4 - 7 meter" in roadtype:
            return 15
        elif "2 - 4 meter" in roadtype:
            return 10

    def assign_buffer_roadtype(self, roadtype: str) -> int:
        if roadtype is None:
            return 0
        elif "autosnelweg" in roadtype:
            return 20
        elif "hoofdweg" in roadtype:
            return 15
        elif "lokale weg" in roadtype:
            return 10
        elif "overig" in roadtype:
            return 0
        elif "parkeerplaats" in roadtype:
            return 5
        elif "regionale weg" in roadtype:
            return 10
        elif "straat" in roadtype:
            return 5

    def create_attribute(self, attri_in: str, new_attri: str) -> None:
        if attri_in == 'verharding':
            self.data[new_attri] = self.data[attri_in].apply(self.assign_buffer_hardSurface)
        if attri_in == 'typeweg':
            self.data[new_attri] = self.data[attri_in].apply(self.assign_buffer_roadtype)

    def create_buffer(self, buffer_attri: str):

        if buffer_attri not in self.data.columns:
            raise ValueError(f"Column {buffer_attri} does not exist in the data.")
        self.data["buffered"] = None
        self.data["buffered"] = self.data.apply(lambda row: row.geometry.buffer(row[buffer_attri]), axis=1)
        if 'buffered' in self.data.columns:
            print("Road buffer geometries created.")
        else:
            print("Fail to create road buffer geometries.")


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"
    filename = "ams_roads_top10NL.shp"
    filepath_mac = directory_mac + filename
    filepath_win = directory_win + filename

    road = Road(gpd.read_file(filepath_win))
    road.create_attribute('typeweg', 'buffer')
    road.create_buffer('buffer')
    road.data['buffered'].to_file(directory_win + "test.shp")
