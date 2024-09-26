import geopandas as gpd


class Building:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data
        self.buffered = None

    def create_buffer(self, buffer_size: int) -> None:
        self.buffered = None
        self.data['buffered_geom'] = self.data['geometry'].buffer(buffer_size)
        self.buffered = gpd.GeoDataFrame(self.data, geometry='buffered_geom', crs=self.data.crs)
        self.data.drop(columns=['buffered_geom'], inplace=True)
        self.buffered.drop(columns=['geometry'], inplace=True)

