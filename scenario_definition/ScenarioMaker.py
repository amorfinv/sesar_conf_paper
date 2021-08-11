# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:47:30 2021

@author: andub
"""
import osmnx as ox
import numpy as np
import BlueskySCNTools
from beta_path_planning import PathPlanner
import os
import pickle

# Initialize stuff
bst = BlueskySCNTools.BlueskySCNTools()

# Step 1: Import the graph we will be using
dir_path = os.path.dirname(os.path.realpath(__file__))
graph_path = dir_path.replace('scenario_definition',
    'graph_definition/gis/streets/directed_groups.graphml')
G = ox.io.load_graphml(graph_path)
nodes, edges = ox.graph_to_gdfs(G)
print('Graph loaded!')

# Step 2: Initalize the Path Planner class
path_planner = PathPlanner(G, nodes, edges)
scennames = []

for idx, concurrent_ac in enumerate([ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                                      8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                      11,11,11,11,11,11,11,11,11,11]):
# ----------------------------------------------------------------------------
# If using existing pickles, comment out everything between the long lines (below)

    #Step 2: Generate traffic from it
    # aircraft_vel = 15 # [m/s]
    # max_time = 3600 # [s]
    # dt = 5
    # min_dist = 1000 # [m]
    # turn_factor = 0
    
    # # Origins are a number of points around the perimeter of the area
    # orig_coords = np.array([[16.34563374, 48.2064161],
    #                 [16.35480063, 48.22504831],
    #                 [16.35011347, 48.22652922],
    #                 [16.35212045, 48.22601715],
    #                 [16.34217249, 48.22816281],
    #                 [16.32695281, 48.20949966],
    #                 [16.32179963, 48.21180067],
    #                 [16.33112073, 48.22751132],
    #                 [16.33199617, 48.22883781],
    #                 [16.35388269, 48.2203232],
    #                 [16.33143868, 48.2279931],
    #                 [16.34145194, 48.2283065],
    #                 [16.35527402, 48.22376019],
    #                 [16.32805092, 48.22285926],
    #                 [16.35321755, 48.22562058],
    #                 [16.32265359, 48.21345807],
    #                 [16.35233842, 48.2165078],
    #                 [16.32726236, 48.22166411],
    #                 [16.32557943, 48.21911323],
    #                 [16.34794355, 48.22698312],
    #                 [16.34657655, 48.22726904],
    #                 [16.32458531, 48.21720681],
    #                 [16.33131364, 48.20878015],
    #                 [16.33022432, 48.22615298],
    #                 [16.35434871, 48.22147447],
    #                 [16.32155546, 48.21132673],
    #                 [16.34901953, 48.20830585],
    #                 [16.34021758, 48.20731047],
    #                 [16.33848969, 48.22889717],
    #                 [16.35403911, 48.22532359],
    #                 [16.34444181, 48.22771024],
    #                 [16.33057391, 48.20890221],
    #                 [16.33427001, 48.22973841],
    #                 [16.32502973, 48.21806918],
    #                 [16.34531602, 48.22753267],
    #                 [16.32225884, 48.21269193],
    #                 [16.35163361, 48.21476624],
    #                 [16.34982089, 48.21028649],
    #                 [16.3512195, 48.22629785],
    #                 [16.32330003, 48.21471264],
    #                 [16.32665069, 48.22073702],
    #                 [16.32358728, 48.21527009],
    #                 [16.34832287, 48.20658385],
    #                 [16.35498457, 48.22304522],
    #                 [16.34051602, 48.22849313],
    #                 [16.3362755, 48.22933861],
    #                 [16.34124711, 48.20714048],
    #                 [16.32919031, 48.22458603],
    #                 [16.32850517, 48.2235477],
    #                 [16.32761704, 48.22220167],
    #                 [16.33337821, 48.20843944],
    #                 [16.3485993, 48.20726714],
    #                 [16.33591201, 48.20802124],
    #                 [16.32921229, 48.20912689],
    #                 [16.34282683, 48.22803232],
    #                 [16.34353675, 48.22789074],
    #                 [16.33767197, 48.20773072],
    #                 [16.34889288, 48.22678455],
    #                 [16.35115844, 48.21359206],
    #                 [16.32132718, 48.21042759],
    #                 [16.32615409, 48.20963143],
    #                 [16.32487971, 48.20984165],
    #                 [16.3229974, 48.21412532],
    #                 [16.33330818, 48.22993013],
    #                 [16.33955482, 48.22868479],
    #                 [16.35072783, 48.21252792],
    #                 [16.32423891, 48.21653463],
    #                 [16.35338749, 48.21909977],
    #                 [16.3281711, 48.20929867],
    #                 [16.35042154, 48.21177098],
    #                 [16.325275, 48.2185451],
    #                 [16.33559567, 48.22947414],
    #                 [16.32320989, 48.21011708],
    #                 [16.35545007, 48.22481353],
    #                 [16.35305335, 48.21827425]])
    # generated_traffic, routes, turnslist, edge_ids = bst.Fast2Scn(G, concurrent_ac, aircraft_vel, max_time, 
    #                                    dt, min_dist, turn_factor, path_planner, orig_coords)
    # print('Traffic generated!')

    # # Step 3.1: Loop through traffic, find path, add to dictionary
    # scenario_dict = dict()
    # for i, flight in enumerate(generated_traffic):
    #     # First get the route and turns
    #     origin = flight[2]
    #     destination = flight[3]
    #     #plan = PathPlanning(G,edges, origin[1], origin[0], destination[1], destination[0])
    #     route,turns= routes[i], turnslist[i]
    #     route = np.array(route)
    #     # Create dictionary
    #     scenario_dict[flight[0]] = dict()
    #     # Add start time
    #     scenario_dict[flight[0]]['start_time'] = flight[1]
    #     #Add lats
    #     scenario_dict[flight[0]]['lats'] = route[:,1]
    #     #Add lons
    #     scenario_dict[flight[0]]['lons'] = route[:,0]
    #     #Add turnbool
    #     scenario_dict[flight[0]]['turnbool'] = turns
    #     #Add alts
    #     scenario_dict[flight[0]]['alts'] = None
    #     # Add edge_ids
    #     scenario_dict[flight[0]]['edge_ids'] = edge_ids[i]

    # print('All paths created!')
    
    # # Step 3.2: Pickle the traffic dictionary and save it in case we need it
    # # later on
    # with open(f'Pickles/Test_Scenario_{concurrent_ac}_{idx+1}.pickle', 'wb') as f:
    #     pickle.dump(scenario_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

# If using the existing pickles, comment out everything between the long lines (above)
# ----------------------------------------------------------------------------
    # Step 3.3: Load pickle dump to use for each file. If this is uncommented, then
    # comment out all previous lines inside the for loop, as we don't need to generate
    # traffic again.
    with open(f'Pickles/Test_Scenario_{concurrent_ac}_{idx+1}.pickle', 'rb') as f:
        scenario_dict = pickle.load(f)
    
    for resometh in ['NONE', 'MVP', 'VO', 'ORCA']:
        # Step 4a: Create scenario file from dictionary
        bst.Dict2Scn(f'Scenarios/Test_Scenario_{concurrent_ac}_{idx+1}_{resometh}.scn', 
                     scenario_dict, resometh=resometh)
        scennames.append(f'Test_Scenario_{concurrent_ac}_{idx+1}_{resometh}.scn')
    
    # airspace structure scenarios
    # Step 4b: Create scenario file from dictionary for airspace experiment
    bst.Dict2Scn(f'Scenarios/Test_Scenario_{concurrent_ac}_{idx+1}_NONE_air.scn', 
                    scenario_dict, resometh='NONE', airspace=True)
    scennames.append(f'Test_Scenario_{concurrent_ac}_{idx+1}_NONE_air.scn')
    print('Scenario file created!')
    
# Step 5: Create batch simulation file
with open('Scenarios/batch_SESAR.scn', 'w') as f:
    for scenname in scennames:
        simplename = scenname.replace('.scn','')
        f.write(f'00:00:00>SCEN {simplename}\n')
        f.write(f'00:00:00>PCALL SESAR2021/{scenname}\n\n')
        
    
    