# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 02:16:33 2016

@author: melni
"""
import sys
import os
from os import listdir
from os.path import isfile, join
import network_creator
import itertools
from graphviz import Digraph
import networkx as nx
from networkx import graph
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import to_agraph

data = sys.argv[1]
outbreaks = [dI for dI in os.listdir(data) if os.path.isdir(os.path.join(data,dI))]

for outbreak in outbreaks:
    current_path = os.path.join(data, outbreak)
    current_samples = [f for f in listdir(current_path) if isfile(join(current_path, f))]
    current_pairs = list(itertools.combinations(current_samples, 2))
    current_outbreak_pairs_min_distances = {}
    for pair in current_pairs:
        first_sequence_set = network_creator.parse_fasta(os.path.join(data, outbreak, pair[0]))
        second_sequence_set = network_creator.parse_fasta(os.path.join(data, outbreak,pair[1]))
        current_pair_distances = []
        for sequence_1 in first_sequence_set:
            for sequence_2 in second_sequence_set:
                d = network_creator.hamming_distance(sequence_1, sequence_2)
                current_pair_distances.append(d)
        current_outbreak_pairs_min_distances[pair] = min(current_pair_distances)
    # Graphviz graph (one per outbreak)
    current_outbreak_graph = Digraph(comment=outbreak)
    for pair in current_pairs:
        current_outbreak_graph.node(pair[0])
        current_outbreak_graph.node(pair[1])
        current_outbreak_graph.edge(pair[0], pair[1])
    current_outbreak_graph.render(os.path.join('outbreak_graphs', outbreak), view=False)
    # Create the same graph in networkx to find the MST
    G = nx.Graph()
    for pair in current_pairs:
        G.add_node(pair[0])
        G.add_node(pair[1])
        G.add_edge(pair[0], pair[1], weight=current_outbreak_pairs_min_distances[pair])
        
    dict_labels = {}
    for i in range(0, len(current_samples)):
          #dict_labels[str(i)] = current_samples[i]
          dict_labels[G.nodes()[i]] = current_samples[i]
    
    MST = nx.minimum_spanning_tree(G, weight='weight')
    MST_agraph = nx.drawing.nx_agraph.to_agraph(MST)
    MST_agraph.render(os.path.join('outbreak_graphs', outbreak+'_mst.png'), view=False)
    #nx.draw(MST)
    #nx.draw_graphviz(MST, cmap=plt.get_cmap('Reds'), labels = dict_labels, with_labels = True)
    
    #plt.savefig(os.path.join('outbreak_graphs', outbreak+'_mst.png'))
    #plt.clf()
