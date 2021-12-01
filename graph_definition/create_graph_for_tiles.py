import osmnx as ox
import geopandas as gpd
from os import path

# code to create the map tile geopackage

# load final graph
G = ox.load_graphml(filepath=path.join('gis', 'streets', 'directed_groups.graphml'))

ox.save_graph_geopackage(G, filepath=path.join('gis', 'streets', 'sesar_map_tile.gpkg'), directed=True)



