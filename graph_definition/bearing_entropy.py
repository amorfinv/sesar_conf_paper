import osmnx as ox
import geopandas as gpd
from os import path
import numpy as np
from osmnx.bearing import orientation_entropy

gis_data_path = 'gis'

# load cleaned graph
G = ox.load_graphml(filepath=path.join(gis_data_path, 'streets', 'directed_groups.graphml'))

# get undirected graph
G_un = ox.get_undirected(G)

## get correct edge bearings
ox.add_edge_bearings(G_un)

# get bearing entropy
print(f'Vienna entropy {ox.orientation_entropy(G_un)}')

### MARTA's paper
boundaries = [37.7648, 37.738, -122.474, -122.51]
G = ox.graph_from_bbox(boundaries[0], boundaries[1], boundaries[2], boundaries[3], network_type='drive')

# get undirected graph
G_un = ox.get_undirected(G)

## get correct edge bearings
ox.add_edge_bearings(G_un)

# get bearing entropy
print(f'Marta entropy {ox.orientation_entropy(G_un)}')
