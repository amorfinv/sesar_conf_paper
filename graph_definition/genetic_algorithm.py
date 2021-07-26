import osmnx as ox
from deap import base
from deap import creator
from deap import tools
import numpy as np
import graph_funcs
from evaluate import dijkstra_search_multiple
import random
import os
import copy

class Node:
    def __init__(self,key,x,y):
        self.key=key
        self.x=x
        self.y=y
        self.children={}# each element of neigh is like [ox_key,edge_cost] 

class GeneticAlgorithm:
    def __init__(self):
        nodes_lonlat = np.array([[16.33316675,48.22759827],
                                [16.33720514,48.22755137],
                                [16.34124351,48.22750432],
                                [16.34528188,48.22745714],
                                [16.33309666,48.22489996],
                                [16.33713484,48.22485306],
                                [16.341173,48.22480602],
                                [16.34521116,48.22475883],
                                [16.34924931,48.22471151],
                                [16.35328744,48.22466404],
                                [16.32898861,48.22224839],
                                [16.33302658,48.22220164],
                                [16.33706455,48.22215474],
                                [16.3411025,48.2221077],
                                [16.34514045,48.22206053],
                                [16.34917838,48.22201321],
                                [16.3532163,48.22196574],
                                [16.32891875,48.21955007],
                                [16.33295652,48.21950332],
                                [16.33699427,48.21945643],
                                [16.34103201,48.21940939],
                                [16.34506974,48.21936222],
                                [16.34910747,48.2193149],
                                [16.35314518,48.21926744],
                                [16.32481135,48.21689834],
                                [16.32884891,48.21685174],
                                [16.33288646,48.21680499],
                                [16.336924,48.21675811],
                                [16.34096153,48.21671108],
                                [16.34499905,48.21666391],
                                [16.34903656,48.2166166],
                                [16.32474172,48.21420001],
                                [16.32877907,48.21415341],
                                [16.33281641,48.21410667],
                                [16.33685374,48.21405979],
                                [16.34089106,48.21401277],
                                [16.34492837,48.2139656],
                                [16.34896567,48.21391829],
                                [16.32467211,48.21150168],
                                [16.32870925,48.21145508],
                                [16.33274638,48.21140835],
                                [16.33678349,48.21136147],
                                [16.3408206,48.21131445],
                                [16.3448577,48.21126729],
                                [16.34889478,48.21121999],
                                [16.33267635,48.20871002],
                                [16.33671326,48.20866315],
                                [16.34075015,48.20861613],
                                [16.34478704,48.20856898],
                                [16.34882391,48.20852168]])

        ################################## LOAD #####################################
        dir_path = os.path.dirname(os.path.realpath(__file__))
        graph_path = dir_path.replace('graph_definition',
            'graph_definition/gis/streets/split_groups.graphml')
        self.G_init = ox.io.load_graphml(graph_path)
        
        node_numbers = ox.nearest_nodes(self.G_init,nodes_lonlat[:,0] ,nodes_lonlat[:,1] )
        
        self.evaluate_nodes = list(np.unique(node_numbers))
            
        g = ox.graph_to_gdfs(self.G_init)
        self.edges = g[1]
        self.nodes = g[0]
    
        # set directionality of groups with one edge
        self.edge_directions = graph_funcs.get_first_group_edges(self.G_init, self.edges)


    def GeneticAlgorithm(self):
        toolbox = base.Toolbox()
        
        # Start genetic algorithm
        creator.create("FitnessMin", base.Fitness, weights = (-1.0,))
        creator.create("Individual", list, fitness = creator.FitnessMin)
        
        toolbox.register("boollist", self.boollist)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.boollist, n=1)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        
        toolbox.register("evaluate", self.CostFunction)
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutFlipBit, indpb = 0.05)
        
        # Do algorithm
        pop = toolbox.population(n = 200)
        fitnesses = list(map(toolbox.evaluate, pop))
        print(fitnesses)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit
            
        CXPB, MUTPB = 0.5, 0.2
        
        fits = [ind.fitness.values[0] for ind in pop]
        
        generation = 0
        prevmin = float('inf')
        globalmin = [float('inf'), None]
        combo = 0
        max_combo = 50
        
        while combo <= max_combo:
            print(f'-------------- Generation {generation} --------------')
            print(f'######## Current combo is {combo}/{max_combo} ########')
            generation = generation + 1
            offspring = toolbox.select(pop, len(pop))
            offspring = list(map(toolbox.clone, offspring))
            
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < CXPB:
                    toolbox.mate(child1[0], child2[0])
                    del child1.fitness.values
                    del child2.fitness.values
                    
            for mutant in offspring:
                if random.random() < MUTPB:
                    toolbox.mutate(mutant[0])
                    del mutant.fitness.values
                    
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
                
            pop[:] = offspring
            fits = [ind.fitness.values[0] for ind in pop]
            
            # Get minimum of this generation
            thismin_val = min(fits)
            thismin_genome = pop[fits.index(thismin_val)]
            if thismin_val <= prevmin:
                combo += 1
                if thismin_val < globalmin[0]:
                    globalmin = [thismin_val, thismin_genome]
            else:
                combo = 0
            
            prevmin = thismin_val
        
        best = pop[np.argmin([toolbox.evaluate(x) for x in pop])]
        return best, globalmin
            
        
        
    def boollist(self):
        return [random.randint(0,1) for i in range(len(self.edge_directions))]
        
    def CostFunction(self, bool_list):
        # Make a copy of edge directions
        directions = copy.copy(self.edge_directions)
        dest_nodes = copy.copy(self.evaluate_nodes)
        # Change direction in function of bool_list
        for i in range(len(bool_list[0])):
            if bool_list[0][i] == 1 or bool_list[0][i] == True:
                direction = copy.copy(directions[i])
                directions[i] = (direction[1], direction[0], direction[2])
                
        # reoroder edge geodatframe
        edges_new = graph_funcs.set_direction2(self.edges, directions)
        print(bool_list)
        G = ox.graph_from_gdfs(self.nodes, edges_new)
        
        # Process nodes to put em in the right form
        omsnx_keys_list=list(G._node.keys())
        
        ###Initialise the graph for the search
        graph={}
        for i in range(len(omsnx_keys_list)):
            key=omsnx_keys_list[i]
            x=G._node[key]['x']
            y=G._node[key]['y']
            node=Node(key,x,y)
            children=list(G._succ[key].keys())
            for ch in children:
                cost=G[key][ch][0]['length']
                node.children[ch]=cost
            
            graph[key]=node
        orig_nodes = []
        for i, node in enumerate(self.evaluate_nodes):
            orig_nodes.append(graph[node])
        
        # Get cost
        cost = dijkstra_search_multiple(graph, orig_nodes, dest_nodes)
        return cost

def main():
    # Let's do genetics
    GA = GeneticAlgorithm()
    #print(GA.CostFunction([[0]*45]))
    
    print('Best solution:', GA.GeneticAlgorithm())
    return

if __name__ == "__main__":
    main()