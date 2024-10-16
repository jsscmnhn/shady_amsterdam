import osmnx as ox
import networkx as nx
import geopandas as gpd
from shapely.geometry import Polygon
import random
import pandas as pd
from pyproj import Transformer

# Load the walkable street network for Amsterdam
G = ox.graph_from_place("Amsterdam, Netherlands", network_type="walk")

# # Project the graph for distance-based calculations
# G_proj = ox.project_graph(G)

# Get the edges as a GeoDataFrame
edges = ox.convert.graph_to_gdfs(G, nodes=False, edges=True)  # Get the edges GeoDataFrame (streets)

edges['u'] = edges.index.get_level_values('u')
edges['v'] = edges.index.get_level_values('v')
edges['key'] = edges.index.get_level_values('key')

# Generate Random Polygons (these represent shaded areas)
# Example: Generating 5 random shaded polygons
def generate_random_polygons(num_polygons, bounds):
    polygons = []
    for _ in range(num_polygons):
        # Generate random coordinates within the bounds of the graph
        minx, miny, maxx, maxy = bounds
        x1, x2 = sorted([random.uniform(minx, maxx) for _ in range(2)])
        y1, y2 = sorted([random.uniform(miny, maxy) for _ in range(2)])

        # Create a polygon from the random coordinates
        polygon = Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
        polygons.append(polygon)
    return polygons


# Get the bounds of the graph (minx, miny, maxx, maxy)
bounds = edges.total_bounds  # Get the spatial bounds of the edges GeoDataFrame
random_polygons = generate_random_polygons(5, bounds)  # Generate 5 random polygons

# Convert polygons to a GeoDataFrame
gdf_polygons = gpd.GeoDataFrame(geometry=random_polygons, crs=edges.crs)

# Calculate the shade index for each edge
edges['shade_index'] = 0.0  # Initialize the shade index column with 0

# Merge all polygons to avoid overlapping areas being counted multiple times
merged_polygons = gdf_polygons.unary_union

# For each edge, calculate how much of it falls inside the shaded polygons
for idx, edge in edges.iterrows():
    # Find the intersection of the edge with the shaded polygons
    # intersection = gdf_polygons.geometry.intersection(edge.geometry)
    intersection = merged_polygons.intersection(edge.geometry)

    # Calculate the length of the part of the edge that is in the shade
    # shaded_length = intersection.length.sum()
    shaded_length = intersection.length

    # Calculate the total length of the edge
    total_length = edge.geometry.length

    # Calculate the shade index as the proportion of the edge that is shaded
    if total_length > 0:
        shade_index = shaded_length / total_length
        edges.at[idx, 'shade_index'] = shade_index  # Update the shade index

# Create a dictionary to map (u, v, key) to shade_index
edge_attr_dict = {}
epsilon = 0.001  # Small value to avoid division by zero

for idx, edge in edges.iterrows():
    u = edge['u']
    v = edge['v']
    key = edge['key'] if 'key' in edge else 0  # In case of multi-edges
    shade_index = edge['shade_index'] if pd.notnull(edge['shade_index']) else 0.0

    # if shade_index == 0:  # Replace 0 with a small positive value to ensure it's traversable
    #     shade_index = 0.001

    # Invert the shade index so that more shaded areas have lower weights
    inverted_weight = 1 / (shade_index + epsilon)

    # Add to the dictionary where the key is (u, v, key)
    # edge_attr_dict[(u, v, key)] = {'shade_weight': shade_index}
    edge_attr_dict[(edge['u'], edge['v'], edge['key'])] = {'shade_weight': inverted_weight}

# Use NetworkX to set the edge attribute in the graph
nx.set_edge_attributes(G, edge_attr_dict)

# Now the graph has 'shade_weight' as an edge attribute
# Access the new 'shade_weight' attribute
for u, v, key, data in G.edges(keys=True, data=True):
    if 'shade_weight' in data:
        print(f"Edge ({u}, {v}) has shade_weight: {data['shade_weight']}")
    else:
        print(f"Edge ({u}, {v}) is missing 'shade_weight' attribute.")

# Define origin and destination nodes for routing (EPSG:4326)
orig_lat, orig_lon = 52.389865, 4.875868  # Example coordinates
# orig_lat, orig_lon = 52.377425, 4.880914
dest_lat, dest_lon = 52.324293, 4.887134  # Example coordinates
# dest_lat, dest_lon = 52.365972, 4.865461

# # Define a transformer to convert from WGS84 (EPSG:4326) to UTM (EPSG:32631)
# transformer = Transformer.from_crs("epsg:4326", "epsg:32631")

# # Transform the origin and destination coordinates to UTM (EPSG:32631)
# orig_x, orig_y = transformer.transform(orig_lat, orig_lon)
# dest_x, dest_y = transformer.transform(dest_lat, dest_lon)

# Get the nearest nodes to the given lat/lon
orig_node = ox.distance.nearest_nodes(G, X=orig_lon, Y=orig_lat)
dest_node = ox.distance.nearest_nodes(G, X=dest_lon, Y=dest_lat)

# Print the nearest nodes to the transformed coordinates
print(f"Nearest node to the origin: {orig_node}")
print(f"Nearest node to the destination: {dest_node}")

# Get the largest strongly connected component
G_largest_component = max(nx.strongly_connected_components(G), key=len)

# Check if the origin and destination nodes are in the same connected component
if orig_node in G_largest_component and dest_node in G_largest_component:
    print("Origin and destination nodes are in the same connected component.")
else:
    print("Origin and destination nodes are not in the same connected component.")

# Find the shadiest path (minimizing the shade_weight)
# shadiest_route = ox.shortest_path(G_proj, orig_node, dest_node, weight='shade_weight')
# Try to compute the shadiest path using Dijkstra's algorithm
try:
    shadiest_route = ox.routing.shortest_path(G, orig_node, dest_node, weight='shade_weight')
    print(f"Shadiest route found with {len(shadiest_route)} nodes.")
except nx.NetworkXNoPath:
    print("No path could be found between the origin and destination.")

# Check the CRS of the graph (edges)
print(f"Graph CRS: {edges.crs}")

# Check the CRS of the polygons
print(f"Polygon CRS: {gdf_polygons.crs}")

# Plot the shadiest route
fig, ax = ox.plot_graph_route(G, shadiest_route, route_linewidth=4, route_color='r',
                              figsize=(12, 12))

# # Plot the random polygons on the graph (shaded areas)
# gdf_polygons.plot(ax=ax, facecolor='none', edgecolor='green', linewidth=2)
