import osmnx as ox
import geopandas as gpd
import networkx as nx
import graph_funcs
from os import path
import momepy

# use osmnx environment here

def main():

    # working path
    gis_data_path = 'gis'

    # get airspace polygon from geopackage
    airspace_gdf = gpd.read_file(path.join(gis_data_path, 'airspace', 'initial_airspace_unprojected.gpkg', ))
    airspace_poly = airspace_gdf.geometry.iloc[0]

    # create MultiDigraph from polygon
    G = ox.graph_from_polygon(airspace_poly, network_type='drive', simplify=True)
    
    # save as osmnx graph
    ox.save_graphml(G, filepath=path.join(gis_data_path, 'streets', 'raw_streets.graphml'))
    ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'raw_streets.gpkg'))
    
    # remove unconnected streets and add edge bearing attrbute 
    G = graph_funcs.remove_unconnected_streets(G)
    ox.add_edge_bearings(G)
    
    #ox.plot.plot_graph(G)
    
    # manually remove some nodes and edges to clean up graph
    nodes_to_remove = [423816204, 423816203, 25267592, 25267593, 33079226, 33079227, 25267602, 33183703, 33183702, 
                        25267606, 3704365809, 33182085, 1370935902, 1381660583, 1381660593, 30696020, 34166967,
                        30696021, 33242270]

    G.remove_nodes_from(nodes_to_remove)
    # G.remove_edges_from(edges_to_remove)
    
    #ox.plot.plot_graph(G)
    
    # convert graph to geodataframe
    g = ox.graph_to_gdfs(G)
    
    # # get node and edge geodataframe
    nodes = g[0]
    edges = g[1]
    # remove double two way edges
    edges = graph_funcs.remove_two_way_edges(edges)
    
    # remove non parallel opposite edges (or long way)
    edges = graph_funcs.remove_long_way_edges(edges)

    # # add interior angles at all intersections
    edges = graph_funcs.add_edge_interior_angles(edges)
    
    # allocated edge height based on cardinal method (TODO: do with groups)
    layer_allocation, _ = graph_funcs.allocate_edge_height(edges, 0)
    edges['layer_height'] = layer_allocation
    
    # Perform COINS algorithm to add stroke groups
    coins_obj = momepy.COINS(edges)
    edges['stroke_group'] = coins_obj.stroke_attribute()
    group_gdf = coins_obj.stroke_gdf()
    
    # set initial edge_directions
    init_edge_directions = graph_funcs.get_first_group_edges(G, group_gdf, edges)
    
    # reoroder edge geodatframe
    edges = graph_funcs.set_direction(edges, init_edge_directions)
    
    # create graph and save edited
    G = ox.graph_from_gdfs(nodes, edges)
    
    # save as osmnx graph
    ox.save_graphml(G, filepath=path.join(gis_data_path, 'streets', 'processed_streets.graphml'))
    
    # Save geopackage for import to QGIS and momepy
    ox.save_graph_geopackage(G, filepath=path.join(gis_data_path, 'streets', 'processed_streets.gpkg'))
    # nx.all_pairs_shortest_path
    # save csv for reference
    edges.to_csv(path.join(gis_data_path, 'streets', 'edges.csv'))
    nodes.to_csv(path.join(gis_data_path, 'streets', 'nodes.csv'))

if __name__ == '__main__':
    main()