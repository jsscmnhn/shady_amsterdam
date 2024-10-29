# Load the graph from the GraphML file
start = datetime.now()
print(f"Starting time for walkable network: {start}")
graph_path = "C:/Androniki/pythonProject1/AMS/metropolitan_region_amsterdam_walk.graphml"
G_walk = ox.load_graphml(graph_path)
end = datetime.now()
print(f"Ending time for walkable network: {end}")

fig, ax = ox.plot_graph(G_walk, figsize=(10, 10))

# Define the coordinates (latitude, longitude) for the origin and destination
orig_lat, orig_lon = 52.389865, 4.875868
dest_lat, dest_lon = 52.324293, 4.887134

# Get the nearest network nodes to these lat/lng points
orig_node = ox.distance.nearest_nodes(G_walk , X=orig_lon, Y=orig_lat)
dest_node = ox.distance.nearest_nodes(G_walk , X=dest_lon, Y=dest_lat)

# Get the coordinates of the origin and destination nodes
orig_x = G_walk.nodes[orig_node]['x']
orig_y = G_walk.nodes[orig_node]['y']
dest_x = G_walk.nodes[dest_node]['x']
dest_y = G_walk.nodes[dest_node]['y']

# find the shortest path between nodes, minimizing travel time, then plot it
route = ox.shortest_path(G_walk, orig_node, dest_node, weight='length')
ox.plot_graph_route(G_walk ,route,figsize=(10,10))

# Load the graph from the GraphML file
print(f"Starting time for bike network: {start}")
graph_path = "C:/Androniki/pythonProject1/AMS/metropolitan_region_amsterdam_bike.graphml"
G_bike = ox.load_graphml(graph_path)
end = datetime.now()
print(f"Ending time for bike network: {end}")
fig, ax = ox.plot_graph(G_bike, figsize=(10, 10))

# Get the nearest network nodes to these lat/lng points
orig_node = ox.distance.nearest_nodes(G_bike , X=orig_lon, Y=orig_lat)
dest_node = ox.distance.nearest_nodes(G_bike , X=dest_lon, Y=dest_lat)

# Get the coordinates of the origin and destination nodes
orig_x = G_bike.nodes[orig_node]['x']
orig_y = G_bike.nodes[orig_node]['y']
dest_x = G_bike.nodes[dest_node]['x']
dest_y = G_bike.nodes[dest_node]['y']

# find the shortest path between nodes, minimizing travel time, then plot it
route = ox.shortest_path(G_bike, orig_node, dest_node, weight='length')
ox.plot_graph_route(G_bike ,route,figsize=(10,10))