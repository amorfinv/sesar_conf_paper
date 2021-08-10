import osmnx as ox
import geopandas as gpd
from os import path
import graph_funcs
import copy
import numpy as np

# def main():
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

# Set directionality
best_solution = [1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 
                0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 
                0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 
                1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 
                0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 
                0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 
                1, 1, 0, 1, 0, 1, 1, 0, 1, 0]

best_solution2 = [1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1,
                  0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 
                  0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 
                  1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 
                  0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0,
                  1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 
                  0, 1, 0, 1, 0, 0, 1, 0, 1, 0]

init_edge_directions = graph_funcs.get_first_group_edges(G, edges)
edge_directions = copy.copy(init_edge_directions)

for i in range(len(best_solution)):
        if best_solution[i] == 1 or best_solution[i] == True:
            direction = copy.copy(edge_directions[i])
            edge_directions[i] = (direction[1], direction[0], direction[2])

# edit genetic results for some small groups 80, and 114 (TODO: FIX METHOD THIS IN FUTURE)
# IDEA IS THAT DEGREE 90 split is done after so these things can be avoided
edge_directions[edge_directions.index((3704365814, 1114680094, 0))] = (1114680094, 3704365814, 0) # 80
edge_directions[edge_directions.index((249607807, 33182067, 0))] = (33182067, 249607807, 0) # 114

new_edges = graph_funcs.set_direction2(edges, edge_directions)

new_edges = graph_funcs.add_edge_interior_angles(new_edges)

### final edits missed on past edits
nodes, new_edges = graph_funcs.manual_edits_after_genetic(nodes, new_edges)

### allocate group heights
new_edges['layer_height'] = graph_funcs.allocate_group_height(nodes, new_edges, rotation_val=0)
# set group 80 to height 2
new_edges.at[(1114680094, 3704365814, 0), 'layer_height'] = 'height 2'

G_final = ox.graph_from_gdfs(nodes, new_edges)
ox.distance.add_edge_lengths(G_final)
ox.add_edge_bearings(G_final)

ox.save_graphml(G_final, filepath=path.join(gis_data_path, 'streets', 'directed_groups.graphml'))

ox.save_graph_geopackage(G_final, filepath=path.join(gis_data_path, 'streets', 'directed_group.gpkg'))



# if __name__ == '__main__':
#     main()

