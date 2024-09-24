import geopandas as gpd


class Building:
    def __init__(self, data: gpd.geodataframe) -> None:
        self.data = data

