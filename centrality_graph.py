#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/7/16
"""

import os
import Bio.SeqIO
import glob
import graph_utils
import hamming_dist_graph
import networkx as nx


OUTBREAK_DIR = 'data/all_clipped'
OUT_EDGE_LIST_FILE_NAME = 'results/centrality_graph.txt'


class CentralityGraph(graph_utils.Graph):
#    BACK_EDGE_PENALTY = 5

    def __init__(self, fastas):
        super(CentralityGraph, self).__init__(nx.DiGraph)
        self.graph.add_nodes_from([x.split('.')[0] for x in fastas.keys()])
        self.hamming_dist_graphs = dict()
        self.build_min_dist_graphs(fastas)
        self.build_edges()

    def build_min_dist_graphs(self, fastas):
        for host in fastas.keys():
            self.hamming_dist_graphs[host.split('.')[0]] = hamming_dist_graph.HammingDistGraph([str(x.seq) for x in fastas[host]])

    def build_edges(self):
        nodes = self.graph.nodes()
        shortest_paths_trees_bunch = [self.hamming_dist_graphs[x].find_shortest_path_trees() for x in nodes]
        min_shortests_path_trees_weights = [{k: graph_utils.ShortestPathTree.find_shortest_path_tree_weight(v)
                                             for k, v in x.items()}
                                            for x in shortest_paths_trees_bunch]
        for i in range(len(nodes) - 1):
            for j in range(i, len(nodes)):
                hosts_indices = [i, j]
                hosts = [nodes[i], nodes[j]]
                neighbour_sequnces = hamming_dist_graph.find_closest_sequences(
                    self.hamming_dist_graphs[hosts[0]].graph.nodes(),
                    self.hamming_dist_graphs[hosts[1]].graph.nodes())
                neighbours = [set(x[i] for x in neighbour_sequnces) for i in [0, 1]]
                spt_border_weights = [max([min_shortests_path_trees_weights[hosts_indices[i]][v]
                                           for v in neighbours[i]])
                                      for i in [0, 1]]
                spt_central_weights = [min([min_shortests_path_trees_weights[hosts_indices[i]][v]
                                           for v in neighbours[i]])
                                       for i in [0, 1]]
                centrality_weights = [spt_border_weights[i]/spt_central_weights[i] if spt_central_weights[i] else 100
                                      for i in [0, 1]]

                self.graph.add_edges_from([(hosts[0], hosts[1]), (hosts[1], hosts[0])])
                self.graph[hosts[0]][hosts[1]]['weight'] = centrality_weights[0]
                self.graph[hosts[1]][hosts[0]]['weight'] = centrality_weights[1]


#                if centrality_weights[0] < centrality_weights[1]:
#                    self.graph[hosts[1]][hosts[0]]['weight'] = centrality_weights[0]
#                    self.graph[hosts[0]][hosts[1]]['weight'] = centrality_weights[0] * self.BACK_EDGE_PENALTY
#                elif centrality_weights[0] == centrality_weights[1]:
#                    self.graph[hosts[1]][hosts[0]]['weight'] = centrality_weights[0]
#                    self.graph[hosts[0]][hosts[1]]['weight'] = centrality_weights[1]
#                else:
#                    self.graph[hosts[1]][hosts[0]]['weight'] = centrality_weights[1] * self.BACK_EDGE_PENALTY
#                    self.graph[hosts[0]][hosts[1]]['weight'] = centrality_weights[1]


def main():
    fastas = {os.path.basename(f): list(Bio.SeqIO.parse(f, 'fasta'))
              for f in glob.glob("%s/*.fas" % OUTBREAK_DIR)}
    centrality_graph = CentralityGraph(fastas)
    nx.write_edgelist(centrality_graph.graph, OUT_EDGE_LIST_FILE_NAME)


if __name__ == '__main__':
    main()
