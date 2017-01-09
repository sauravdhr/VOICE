#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12/01/2016
"""

import json
import networkx as nx
import networkx.readwrite.json_graph
import enum


def import_graph_from_json(file_name):
    with open(file_name) as f:
        data = json.load(f)
        return nx.readwrite.json_graph.adjacency_graph(data)


class Graph(object):
    def __init__(self, graph_constructor):
        self.graph_constructor = graph_constructor
        self.graph = graph_constructor()

    def find_shortest_path_trees(self):
        nodes = self.graph.nodes()
        return dict(zip(nodes, [ShortestPathTree.find_shortest_path_tree(self.graph, v, self.graph_constructor)
                                for v in nodes]))


class ShortestPathTree(object):
    @staticmethod
    def find_shortest_path_tree(graph, vertex, graph_constructor):
        edges_dict = dict()
        dijkstra_paths = nx.single_source_dijkstra_path(graph, vertex)
        for dest in dijkstra_paths:
            v = vertex
            for w in dijkstra_paths[dest]:
                if v != w and (v, w) not in edges_dict:
                    edges_dict[(v, w)] = graph[v][w]['weight']
                v = w
        shortest_path_tree = graph_constructor()
        shortest_path_tree.add_weighted_edges_from([(x[0][0], x[0][1], x[1]) for x in edges_dict.items()])
        return shortest_path_tree

    @staticmethod
    def find_shortest_path_tree_weight(tree):
        return sum([tree[e[0]][e[1]]['weight'] for e in tree.edges()])
