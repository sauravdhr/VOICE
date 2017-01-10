#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/8/17
"""

import networkx as nx
import graph_utils
import hamming_dist_graph


class SimulationDistGraph(object):
    BASE_BORDER_VICINITY = 1

    def __init__(self, seqs1, seqs2):
        self.seqs = [seqs1[:], seqs2[:]]
        closest_seqs = hamming_dist_graph.find_closest_sequences(seqs1, seqs2, self.BASE_BORDER_VICINITY)
        self.borders_seqs = [[x[0] for x in closest_seqs], [x[1] for x in closest_seqs]]
        if False:
            self.seqs = self.borders_seqs
        self.borders_consensuses = [hamming_dist_graph.find_consensus(s) for s in self.borders_seqs]
        for i in [0, 1]:
            if self.borders_consensuses[i] not in self.seqs[i]:
                self.seqs[i].append(self.borders_consensuses[i])
        self.hamming_dist_graphs = [hamming_dist_graph.HammingDistGraph(s) for s in self.seqs]
        consensuses_dist = hamming_dist_graph.hamming_distance(self.borders_consensuses[0], self.borders_consensuses[1])
        for (i, j) in [(0, 1), (1, 0)]:
            if self.borders_consensuses[j] not in self.seqs[i]:
                self.hamming_dist_graphs[i].graph.add_edge(self.borders_consensuses[i], self.borders_consensuses[j])
                self.hamming_dist_graphs[i].graph[self.borders_consensuses[i]][self.borders_consensuses[j]]['weight'] =\
                    consensuses_dist
            self.hamming_dist_graphs[i].graph.node[self.borders_consensuses[j]]['source'] = True
        self.simulation_graphs = [graph_utils.ShortestPathTree.find_shortest_path_tree(
            self.hamming_dist_graphs[i].graph, self.borders_consensuses[j], nx.Graph)
                                  for (i, j) in [(0, 1), (1, 0)]]
