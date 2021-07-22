#################### airspace folder #####################
Initial airspace defined 21-07-2021: 

-initial_airspace.gpkg---> projected crs EPSG:32633 (meters)
-initial_airspace_unprojected---> unprojected crs EPSG:4326

#################### streets folder #####################
Data created via graph_definition/create_graph.py

-raw_streets.gpkg---> unprojected crs
-raw_streets.graphml---> graph object same as raw_streets.gpkg

some operations are performed on raw graph
1) remove unconnected streets
2) add edge bearings
3) manually remove some nodes/edges
4) convert to directed from multidirected
5) allocate edge heights based on cardinal method
6) perform COINS
7) give directionality to groups
8) edge interior angle calculation

-processed_streets.gpkg---> unprojected crs
-processed_streets.graphml---> graph object same as processed_streets.gpkg
