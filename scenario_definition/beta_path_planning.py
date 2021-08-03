'''
FAST PATH PLANNER
'''
from os import lseek
import osmnx as ox
import math
from shapely.geometry import MultiLineString
from shapely import ops
class PathPlanner():
    def __init__(self, G, nodes, edges, angle_cutoff=20):
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
        try:
            length_route = len(osmid_route)
        except TypeError:
            length_route = 0

        if length_route > 1:
            print('Path found!')
            # get all edges and their geometry
            edge_list = []
            geom_list = []

            for idx in range(len(osmid_route) - 1):
                
                # get edge
                edge = (osmid_route[idx], osmid_route[idx + 1], 0)
                
                # get geometry of edge
                line_geom = self.edge_gdf.loc[edge, 'geometry']

                # get all edges
                edge_list.append(edge)

                # get geometry list
                geom_list.append(line_geom)
            
            # combine geometry them into one linestring list
            merged_line = ops.linemerge(MultiLineString(geom_list))
            merged_line_list = list(merged_line.coords)

            # create route and add height to route (50 for now) TODO: make nicely dynamic
            height = 50
            route = [(merged_line_list[idx][0], merged_line_list[idx][1], height) for idx in range(len(merged_line_list))]

            # create turnbool based on interior angle between edges, origin is zero
            turnbool = [0]
        
            for idx in range(len(merged_line_list)-2):
                line_string_1 = [merged_line_list[idx], merged_line_list[idx + 1]]
                line_string_2 = [merged_line_list[idx+1], merged_line_list[idx + 2]]
                
                angle = 180 - angleBetweenTwoLines(line_string_1,line_string_2)

                if angle > self.angle_cutoff:
                    turn_flag = 1
                else:
                    turn_flag = 0
                
                turnbool.append(turn_flag)

            # add a zero for turnbool of destination
            turnbool.append(0)
        else:
            route = []
            turnbool = []
            print('No path Found!')

        return route, turnbool

"""
The below function calculates the joining angle between
two line segments. FROM COINS
"""
def angleBetweenTwoLines(line1, line2):
    l1p1, l1p2 = line1
    l2p1, l2p2 = line2
    l1orien = computeOrientation(line1)
    l2orien = computeOrientation(line2)
    """
    If both lines have same orientation, return 180
    If one of the lines is zero, exception for that
    If both the lines are on same side of the horizontal plane, calculate 180-(sumOfOrientation)
    If both the lines are on same side of the vertical plane, calculate pointSetAngle
    """
    if (l1orien==l2orien): 
        angle = 180
    elif (l1orien==0) or (l2orien==0): 
        angle = pointsSetAngle(line1, line2)
        
    elif l1p1 == l2p1:
        if ((l1p1[1] > l1p2[1]) and (l1p1[1] > l2p2[1])) or ((l1p1[1] < l1p2[1]) and (l1p1[1] < l2p2[1])):
            angle = 180 - (abs(l1orien) + abs(l2orien))
        else:
            angle = pointsSetAngle([l1p1, l1p2], [l2p1,l2p2])
    elif l1p1 == l2p2:
        if ((l1p1[1] > l2p1[1]) and (l1p1[1] > l1p2[1])) or ((l1p1[1] < l2p1[1]) and (l1p1[1] < l1p2[1])):
            angle = 180 - (abs(l1orien) + abs(l2orien))
        else:
            angle = pointsSetAngle([l1p1, l1p2], [l2p2,l2p1])
    elif l1p2 == l2p1:
        if ((l1p2[1] > l1p1[1]) and (l1p2[1] > l2p2[1])) or ((l1p2[1] < l1p1[1]) and (l1p2[1] < l2p2[1])):
            angle = 180 - (abs(l1orien) + abs(l2orien))
        else:
            angle = pointsSetAngle([l1p2, l1p1], [l2p1,l2p2])
    elif l1p2 == l2p2:
        if ((l1p2[1] > l1p1[1]) and (l1p2[1] > l2p1[1])) or ((l1p2[1] < l1p1[1]) and (l1p2[1] < l2p1[1])):
            angle = 180 - (abs(l1orien) + abs(l2orien))
        else:
            angle = pointsSetAngle([l1p2, l1p1], [l2p2,l2p1])
    return(angle)

"""
This below function calculates the acute joining angle between
two given set of points. FROM COINS
"""
def pointsSetAngle(line1, line2):
    l1orien = computeOrientation(line1)
    l2orien = computeOrientation(line2)
    if ((l1orien>0) and (l2orien<0)) or ((l1orien<0) and (l2orien>0)):
        return(abs(l1orien)+abs(l2orien))
    elif ((l1orien>0) and (l2orien>0)) or ((l1orien<0) and (l2orien<0)):
        theta1 = abs(l1orien) + 180 - abs(l2orien)
        theta2 = abs(l2orien) + 180 - abs(l1orien)
        if theta1 < theta2:
            return(theta1)
        else:
            return(theta2)
    elif (l1orien==0) or (l2orien==0):
        if l1orien<0:
            return(180-abs(l1orien))
        elif l2orien<0:
            return(180-abs(l2orien))
        else:
            return(180 - (abs(computeOrientation(line1)) + abs(computeOrientation(line2))))
    elif (l1orien==l2orien):
        return(180)

# FROM COINS

def computeOrientation(line):
    point1 = line[1]
    point2 = line[0]
    """
    If the latutide of a point is less and the longitude is more, or
    If the latitude of a point is more and the longitude is less, then
    the point is oriented leftward and wil have negative orientation.
    """
    if ((point2[0] > point1[0]) and (point2[1] < point1[1])) or ((point2[0] < point1[0]) and (point2[1] > point1[1])):
        return(-computeAngle(point1, point2))
    #If the latitudes are same, the line is horizontal
    elif point2[1] == point1[1]:
        return(0)
    #If the longitudes are same, the line is vertical
    elif point2[0] == point1[0]:
        return(90)
    else:
        return(computeAngle(point1, point2))

"""
The function below calculates the angle between two points in space. FROM COINS
"""

def computeAngle(point1, point2):
    height = abs(point2[1] - point1[1])
    base = abs(point2[0] - point1[0])
    angle = round(math.degrees(math.atan(height/base)), 3)
    return(angle)
