import geopandas as gpd
import matplotlib.pyplot as plt

# Step 1: Load the GeoPackage and specify the layer you are interested in
gpkg_file = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\top10nl_Compleet.gpkg"
road_layer = 'top10nl_wegdeel_vlak'  # The specific layer name in the GeoPackage
adms_layer = 'top10nl_registratief_gebied_multivlak'
lu_layer = 'top10nl_terrein_vlak'
publicshp = "D:\\OneDrive - Delft University of Technology\\Synthesis project\\cool place\\datasets\\ams_publicspace_bagplus.shp"

# Load the layer into a GeoDataFrame
road = gpd.read_file(gpkg_file, layer=road_layer)
adms = gpd.read_file(gpkg_file, layer=adms_layer)
lu = gpd.read_file(gpkg_file, layer=lu_layer)
public = gpd.read_file(publicshp)

# Step 2: Filter the data by the attribute

ams = adms[adms['naamofficieel'] == 'Amsterdam']

lu_ams = gpd.clip(lu, ams)
lu_green = lu_ams[lu_ams['typelandgebruik'].isin(['bos: gemengd bos', 'bos:griend', 'bos:loofbos'])]
publicspace = public[public['type_omsch']== 'Landschappelijk gebied' ]

geodata = [lu_green, publicspace]
intersection = geodata[0]
for gdata in geodata[1:]:
    intersection = gpd.overlay(intersection, gdata, how='intersection')
# Print the GeoDataFrame to see its attributes
print(intersection)


road_ams = gpd.clip(road, ams)
road_filter = road_ams[road_ams['hoofdverkeersgebruik'] == 'gemengd verkeer']
road_width = road_ams[road_ams['verhardingsbreedteklasse'] == 'gemengd verkeer']
# Define a function to determine buffer size based on an attribute value
def get_buffer_size(row):
    if row['verhardingsbreedteklasse'] == '> 7 meter':
        return 20  # Buffer size for value1
    elif row['verhardingsbreedteklasse'] == '2 - 4 meter':
        return 5  # Buffer size for value2
    elif row['verhardingsbreedteklasse'] == '4 - 7 meter':
        return 10  # Buffer size for value2
    else:
        return 0   # Default buffer size

# Apply the function to create a new column for buffer sizes
road_ams['buffer_size'] = road_ams.apply(get_buffer_size, axis=1)
# Create buffers based on the calculated buffer sizes
road_ams['geometry'] = road_ams.geometry.buffer(road_ams['buffer_size'])
# Filter only POLYGON geometries
road_ams = road_ams[road_ams.geometry.type == 'Polygon']
# Now try saving it to a shapefile again
road_ams.to_file('road_ams.shp')
# # Check buffer sizes
# print(road_ams['buffer_size'])
# # Check the geometry after buffering
# print(road_ams.geometry)
# Plot the original and buffered geometries for comparison
fig, ax = plt.subplots(figsize=(10, 10))
road_ams.plot(ax=ax, color='red', alpha=0.5, label='Buffered Roads')
#road_ams.geometry.boundary.plot(ax=ax, color='black', label='Original Roads')
plt.title('Road Buffers by Attribute Value')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()



# intersection = gpd.overlay(lu_green, publicspace, how='intersection')
intersection.plot()
# Customize the plot (optional)
plt.title("Intersection Map")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
# Show the map
plt.show()



# Plot both layers on the map: clipped data and boundary data
fig, ax = plt.subplots(figsize=(10, 10))
# Plot the boundary (as a background)
ams.plot(ax=ax, color='lightgray', edgecolor='black')
# Plot the clipped data (as the foreground)
road_filter.plot(ax=ax, color='blue', edgecolor='black', alpha=0.7)
# Add titles and labels (optional)
ax.set_title("Clipped Landschappelijk gebied with Boundary", fontsize=16)
ax.set_xlabel("Longitude", fontsize=12)
ax.set_ylabel("Latitude", fontsize=12)
# Show the map
plt.show()

# Plot the filtered data
publicspace.plot()
# Customize the plot (optional)
plt.title("Public space Map")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
# Show the map
plt.show()

# Plot the filtered data
lu_green.plot()
# Customize the plot (optional)
plt.title("lu_green Map")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
# Show the map
plt.show()