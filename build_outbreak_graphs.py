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
    current_outbreak_graph = Digraph(comment=outbreak)
    for pair in current_pairs:
        current_outbreak_graph.node(pair[0])
        current_outbreak_graph.node(pair[1])
        current_outbreak_graph.edge(pair[0], pair[1])
    current_outbreak_graph.render(os.path.join('outbreak_graphs', outbreak), view=False)