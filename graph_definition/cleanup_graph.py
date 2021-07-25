import os
from platform import node
from networkx.classes.function import degree
import osmnx as ox
import geopandas as gpd
import networkx as nx
from osmnx.projection import project_gdf
from osmnx.utils_graph import graph_from_gdfs
import graph_funcs
from os import path
import momepy

# use osmnx environment here

def main():

    # working path
    gis_data_path = 'gis'

    # get airspace polygon from geopackage
    edges = gpd.read_file(path.join(gis_data_path, 'streets', 'initial_clean_edges.gpkg'), layer='initial_clean_edges')
    nodes = gpd.read_file(path.join(gis_data_path, 'streets', 'initial_clean_nodes.gpkg'))

    # recreate edges in format for osmnx
    edges = graph_funcs.edge_gdf_format_from_gpkg(edges)
    nodes = graph_funcs.node_gdf_format_from_gpkg(nodes)
    G_dummy = ox.graph_from_gdfs(nodes, edges)
    
    # run coins to reorder edges
    coins_obj = momepy.COINS(edges)
    edges['stroke_group'] = coins_obj.stroke_attribute()
    group_gdf = coins_obj.stroke_gdf()

    #set intial group directions to get everything back in order
    init_edge_directions = graph_funcs.get_first_group_edges(G_dummy, group_gdf, edges)
    edges = graph_funcs.set_direction2(edges, init_edge_directions)

    # add edge interior angles due to edge changes of COINS
    edges = graph_funcs.add_edge_interior_angles(edges)

    # reomove degree 2 edges with int angle greater than 120, requires fresh add_edge_interior angles
    nodes, edges = graph_funcs.simplify_graph(nodes, edges, angle_cut_off=120)

    # # manually adapt some edges
    nodes, edges = graph_funcs.manual_edits(nodes, edges)

    # rerun coins to group new edges and interior angles
    coins_obj = momepy.COINS(edges)
    edges['stroke_group'] = coins_obj.stroke_attribute()
    edges = graph_funcs.add_edge_interior_angles(edges)

    # get edge geomtery (linestring info) from geodataframe and calculate integral bearing differenc
    edges_geometry = edges['geometry'].to_numpy()
    integral_bearings_diff = graph_funcs.calculate_integral_bearing_difference(edges_geometry)

    # add bearing difference column to geodataframe
    edges['int_bearing_diff'] = integral_bearings_diff

    # find degree node 2 require fresh add interior angles TODO: migrate from group_split
    # _, _ = graph_funcs.new_groups_90(nodes, edges)

    # create graph and save edited
    G = ox.graph_from_gdfs(nodes, edges)
    ox.add_edge_bearings(G)
    ox.distance.add_edge_lengths(G)

    # save as osmnx graph
    ox.save_graphml(G, filepath=path.join(gis_data_path, 'streets', 'cleaned_streets.graphml'))
    
    # # Save geopackage for import to QGIS and momepy
    ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'cleaned_streets.gpkg'))
    edges.to_csv(path.join(gis_data_path, 'streets', 'edges_clean.csv'))
    nodes.to_csv(path.join(gis_data_path, 'streets', 'nodes_clean.csv'))

if __name__ == '__main__':
    main()