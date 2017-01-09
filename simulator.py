#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/9/17
"""
import random


class Simulator(object):
    def __init__(self, probability_graph):
        self.probability_graph = probability_graph
        self.counter = 0
        self.visited = {k: False for k in self.probability_graph.nodes()}
        self.active_nodes = dict()
        self.source = self.find_source()
        random.seed()

    def find_source(self):
        for v in self.probability_graph.nodes():
            if 'source' in self.probability_graph.node[v]:
                return v

    def simulate(self):
        self.visited[self.source] = True
        self.active_nodes[self.source] = None
        while True:
            self.counter += 1
            visited_nodes = []
            for i, active_node in enumerate(self.active_nodes):
                d = random.random()
                p = 0
                for (u, v, w) in self.probability_graph.edges(active_node, data='weight'):
                    p += w
                    if p > d:
                        visited_nodes.append(v)
                        break

            new_active_nodes = []
            for v in visited_nodes:
                if not self.visited[v]:
                    self.visited[v] = True
                    new_active_nodes.append(v)

            for v in new_active_nodes:
                self.active_nodes[v] = None
            self.update_active_nodes()

            if self.are_all_nodes_visited():
                return self.counter

    def update_active_nodes(self):
        nodes_for_inactivation = []
        for i, node in enumerate(self.active_nodes):
            if self.are_all_neighbours_visited(node):
                nodes_for_inactivation.append(node)
        for node in nodes_for_inactivation:
            self.active_nodes.pop(node, None)

    def are_all_neighbours_visited(self, node):
        for (_, v) in self.probability_graph.edges(node):
            if not self.visited[v]:
                return False
        return True

    def are_all_nodes_visited(self):
        return all(self.visited.values())
