'''
FAST PATH PLANNER
'''
import osmnx as ox
from shapely.geometry import LineString
import ast

class PathPlanner():
    def __init__(self, G, nodes, edges, angle_cutoff=30):
        self.G = G

        # get edge geodataframe
        self.node_gdf = nodes
        self.edge_gdf = edges

        # get edge indices
        self.edge_idx = list(self.edge_gdf.index.values)

        # get angle cutoff to label turns as turnbool
        self.angle_cutoff = angle_cutoff

    def route(self, origin_node, dest_node):

        # get route as a list of osmids
        osmid_route = ox.shortest_path(self.G, origin_node, dest_node)

        # get_correct_order of edges inside graph and reverese linestring geometry if necessary
        edge_geom_list = []

        if osmid_route:
            print('Path found!')
            for idx in range(len(osmid_route) - 1):

                edge = (osmid_route[idx], osmid_route[idx + 1], 0)
                
                # get geometry of line
                line_geom = list(self.edge_gdf.loc[edge, 'geometry'].coords)

                # append edge and geometry for later use
                edge_geom_list.append((edge, line_geom))

            # calculate succesive interior angles and see which nodes are turn nodes
            int_angle_list = []
            turn_node_list = []
            for idx in range(len(edge_geom_list) - 1):
                current_edge = edge_geom_list[idx][0]
                next_edge = edge_geom_list[idx + 1][0]

                int_angle_dict = ast.literal_eval(self.edge_gdf.loc[current_edge, 'edge_interior_angle'])
                
                # get interior angle. search in current_edge
                interior_angle = int_angle_dict[next_edge]
                
                # get osmids of turn nodes
                if interior_angle < 180 - self.angle_cutoff:
                    node_to_append = current_edge[1]
                    turn_node_list.append(node_to_append)

                int_angle_list.append(interior_angle)     

            # create list of lat lon for path finding
            lat_list = []
            lon_list = []
            lon_lat_list = []   # this is used for searching for turn nodes
            for edge_geo in edge_geom_list:
                edge = edge_geo[0]
                geom = edge_geo[1]

                # add all geometry info. adds the first node and second to last for lon/lat
                for idx in range(len(geom) - 1):
                    lon = geom[idx][0]
                    lat = geom[idx][1]

                    lon_list.append(lon)
                    lat_list.append(lat)
                    lon_lat_list.append(f'{lon}-{lat}')

            # add destination node to lists because for loop above does not
            lon_dest = self.node_gdf.loc[dest_node, 'x']
            lat_dest = self.node_gdf.loc[dest_node, 'y']
            lon_list.append(lon_dest)
            lat_list.append(lat_dest)
            lon_lat_list.append(f'{lon_dest}-{lat_dest}')

            # find indices of turn_nodes
            turn_indices = []
            for turn_node in turn_node_list:
                
                # Find lat lon of current turn node
                lat_node = self.node_gdf.loc[turn_node, 'y']
                lon_node = self.node_gdf.loc[turn_node, 'x']

                index_turn = lon_lat_list.index(f'{lon_node}-{lat_node}')
                turn_indices.append(index_turn)

            # create turnbool. true if waypoint is a turn waypoint, else false
            turnbool = []
            for idx in range(len(lat_list)):
                if idx in turn_indices:
                    turn_flag = 1
                else:
                    turn_flag = 0

                turnbool.append(turn_flag)

            # add height to route (50 for now) TODO: make nicely dynamic
            height_list = [50 for idx in range(len(lat_list))] 
            
            route = list(zip(lon_list, lat_list, height_list))
        
        else:
            route = []
            turnbool = []
            print('No path Found!')

        return route, turnbool
