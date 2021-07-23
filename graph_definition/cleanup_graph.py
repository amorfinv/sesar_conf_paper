from platform import node
from networkx.classes.function import degree
import osmnx as ox
import geopandas as gpd
import networkx as nx
from osmnx.utils_graph import graph_from_gdfs
import graph_funcs
from os import path
import momepy

# use osmnx environment here

def main():

    # working path
    gis_data_path = 'gis'

    # get airspace polygon from geopackage
    edges = gpd.read_file(path.join(gis_data_path, 'streets', 'pre_cleaned_streets.gpkg'), layer='pre_cleaned_streets')
    nodes = gpd.read_file(path.join(gis_data_path, 'streets', 'pre_cleaned_intersections.gpkg'))

    # recreate edges in format for osmnx
    edges = graph_funcs.edge_gdf_format_from_gpkg(edges)
    nodes = graph_funcs.node_gdf_format_from_gpkg(nodes)

    # add edge interior anfles
    edges = graph_funcs.add_edge_interior_angles(edges)

    # simplify graph by removing intersitial nodes
    edges = graph_funcs.simplify_graph(nodes, edges)

    # create osmnx graph
    G = ox.graph_from_gdfs(nodes, edges)

    # find single degree nodes
    node_ids = list(nodes.index.values)

    for osmid, node_deg in G.degree(node_ids):
        if node_deg == 1:
            print(osmid)

    ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'cleaned_streets.gpkg'))

    # # create MultiDigraph from polygon
    # G = ox.graph_from_polygon(airspace_poly, network_type='drive', simplify=True)
    
    # # save as osmnx graph
    # ox.save_graphml(G, filepath=path.join(gis_data_path, 'streets', 'raw_streets.graphml'))
    # ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'raw_streets.gpkg'))
    
    # # remove unconnected streets and add edge bearing attrbute 
    # G = graph_funcs.remove_unconnected_streets(G)
    # ox.add_edge_bearings(G)
    
    # #ox.plot.plot_graph(G)
    
    # # manually remove some nodes and edges to clean up graph
    # nodes_to_remove = [423816204, 423816203, 25267592, 25267593, 33079226, 33079227, 25267602, 33183703, 33183702, 
    #                     25267606, 3704365809, 33182085, 1370935902, 1381660583, 1381660593, 30696020, 34166967,
    #                     30696021, 33242270]
    # # edges_to_remove = [(291088171, 3155094143), (60957703, 287914700), (2451285012, 287914700),
    # #                    (25280685, 30696019), (30696019, 25280685), (392251, 25280685), 
    # #                    (25280685, 392251), (33301346, 1119870220),  
    # #                    (33345331, 33345333), (378699, 378696), (378696, 33143911), 
    # #                    (33143911, 33144821), (264061926, 264055537), (33144706, 33144712),
    # #                    (33144712, 33174086), (33174086, 33144719), (33144719, 92739749),
    # #                    (33345319, 29048469), (287914700, 60957703), (213287623, 251207325),
    # #                    (251207325, 213287623)]
    # G.remove_nodes_from(nodes_to_remove)
    # # G.remove_edges_from(edges_to_remove)
    
    # #ox.plot.plot_graph(G)
    
    # # convert graph to geodataframe
    # g = ox.graph_to_gdfs(G)
    
    # # # get node and edge geodataframe
    # nodes = g[0]
    # edges = g[1]
    
    # # remove double two way edges
    # edges = graph_funcs.remove_two_way_edges(edges)
    
    # # remove non parallel opposite edges (or long way)
    # edges = graph_funcs.remove_long_way_edges(edges)
    
    # # allocated edge height based on cardinal method (TODO: do with groups)
    # layer_allocation, _ = graph_funcs.allocate_edge_height(edges, 0)
    # edges['layer_height'] = layer_allocation
    
    # # Perform COINS algorithm to add stroke groups
    # coins_obj = momepy.COINS(edges)
    # edges['stroke_group'] = coins_obj.stroke_attribute()
    # group_gdf = coins_obj.stroke_gdf()
    
    # #init_edge_directions = graph_funcs.get_first_group_edges(G, group_gdf, edges)
    
    # # Apply direction algorithm
    # # edge_directions = calcDirectionality(group_gdf, nodes, init_edge_directions)
    
    # # reoroder edge geodatframe
    # # edges = graph_funcs.set_direction(edges, init_edge_directions)

    # # # add interior angles at all intersections
    # edges = graph_funcs.add_edge_interior_angles(edges)
    
    # # create graph and save edited
    # G = ox.graph_from_gdfs(nodes, edges)
    
    # # save as osmnx graph
    # ox.save_graphml(G, filepath=path.join(gis_data_path, 'streets', 'processed_streets.graphml'))
    
    # # Save geopackage for import to QGIS and momepy
    # ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'processed_streets.gpkg'))
    # nx.all_pairs_shortest_path
    # # save csv for reference
    # edges.to_csv(path.join(gis_data_path, 'streets', 'edges.csv'))
    # nodes.to_csv(path.join(gis_data_path, 'streets', 'nodes.csv'))

if __name__ == '__main__':
    main()