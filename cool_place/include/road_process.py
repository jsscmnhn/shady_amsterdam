import geopandas as gpd


class Road:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.buffered = None

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
            return 5
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
        self.buffered = None
        self.data['buffered_geom'] = self.data.apply(lambda row: row.geometry.buffer(row[buffer_attri]), axis=1)
        self.buffered = gpd.GeoDataFrame(self.data, geometry='buffered_geom', crs=self.data.crs)
        self.data.drop(columns=['buffered_geom'], inplace=True)
        self.buffered.drop(columns=['geometry'], inplace=True)


if __name__ == '__main__':
    directory_mac = "/Volumes/T7 Shield/TUD/Synthesis/cool_place/"
    directory_win = "G:\\TUD\\Synthesis\\cool_place\\"
    filename = "ams_roads_top10NL.shp"
    filepath_mac = directory_mac + filename
    filepath_win = directory_win + filename

    road = Road(gpd.read_file(filepath_mac))
    road.create_attribute('typeweg', 'buffer')
    road.create_buffer('buffer')
    road.buffered.to_file(directory_mac + "ams_roads_top10NL_buffered.shp")
