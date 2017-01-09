#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/9/17
"""

import networkx as nx
import graph_utils


MUTATION_PROBABILITY = 0.01


class ProbabilityGraph(graph_utils.Graph):
    e = MUTATION_PROBABILITY
    s = e/(1-3*e)

    @staticmethod
    def loop_probability(L):
        return (1 - 3 * ProbabilityGraph.e) ** L

    def edge_probability(self, m):
        return self.c*(self.s ** m)

    def __init__(self, simulation_graph, L=20):
        super(ProbabilityGraph, self).__init__(nx.Graph)
        self.L = L
        self.c = self.loop_probability(self.L)
        self.unit_edge_probability = self.c * self.s
        self.source = None
        self.artificial_nodes_counter = 0
        self.build_edges(simulation_graph)
        self.mark_nodes(simulation_graph)


    def build_edges(self, simulation_dist_graph):
#            self.graph.add_edge(node, node)
#            self.graph[node][node]['weight'] = self.c
        for edge in simulation_dist_graph.edges():
            self.add_artificial_edges(edge[0], edge[1], simulation_dist_graph[edge[0]][edge[1]]['weight'])

    def add_artificial_edges(self, u, v, weight):
        if weight == 1:
            self.graph.add_edge(u, v)
            self.graph[u][v]['weight'] = self.unit_edge_probability
            return
        new_artificial_node_counter = self.artificial_nodes_counter + weight - 1
        for i in range(self.artificial_nodes_counter, new_artificial_node_counter + 1):
            if i == self.artificial_nodes_counter:
                self.graph.add_edge(u, i)
                self.graph[u][i]['weight'] = self.unit_edge_probability
            elif i == new_artificial_node_counter:
                self.graph.add_edge(i-1, v)
                self.graph[i-1][v]['weight'] = self.unit_edge_probability
            else:
                self.graph.add_edge(i-1, i)
                self.graph[i-1][i]['weight'] = self.unit_edge_probability
        self.artificial_nodes_counter = new_artificial_node_counter

    def mark_nodes(self, simulation_dist_graph):
        for node in simulation_dist_graph.nodes():
            self.graph.node[node]['original'] = True
            if 'source' in simulation_dist_graph.node[node]:
                self.graph.node[node]['source'] = True
                self.source = node

