import osmnx as ox
import geopandas as gpd
from os import path
import graph_funcs

def main():
    # gis data pah
    gis_data_path = 'gis'

    # load cleaned graph
    G = ox.load_graphml(filepath=path.join(gis_data_path, 'streets', 'cleaned_streets.graphml'))

    # convert to gdf
    nodes, edges = ox.graph_to_gdfs(G)

    # split groups at degree-2 nodes with 90 degree
    nodes, edges= graph_funcs.new_groups_90(nodes, edges)
    
    # convert back to graph and save
    G = ox.graph_from_gdfs(nodes, edges)
    ox.save_graphml(G, filepath=path.join(gis_data_path, 'streets', 'split_groups.graphml'))

    ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'split_group.gpkg'))



if __name__ == '__main__':
    main()

