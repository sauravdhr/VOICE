#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created:
"""
#import argparse
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


SOURCES_FILE_NAME = 'data/sources.txt'
MIN_DIST_EDGE_LIST = 'data/all_clipped_min_dist.txt'
SIMULATION_EDGE_LIST = 'data/all_clipped_simulation.txt'


class GraphAnalyzer(object):
    UNRELATED = 'XX'

    def __init__(self, edges_list_file_name, graph):
        self.edges_list_file_name = edges_list_file_name
        self.graph = graph
        self.outbreaks_nodes = self.infer_outbreaks_nodes(self.graph)

    @staticmethod
    def infer_outbreaks_nodes(graph):
        outbreaks_nodes = dict()
        for node in graph.nodes():
            outbreak_name = node[:2]
            if outbreak_name not in outbreaks_nodes:
                outbreaks_nodes[outbreak_name] = [node]
            else:
                outbreaks_nodes[outbreak_name].append(node)
        return outbreaks_nodes

    @staticmethod
    def get_unrelated_edges_weights(graph, edges):
        unrelated_edges_weights = list()
        for e in edges:
            if e[0][:2] != e[1][:2] or e[0][:2] == GraphAnalyzer.UNRELATED:
                unrelated_edges_weights.append(graph[e[0]][e[1]]['weight'])
        return unrelated_edges_weights

    @staticmethod
    def get_related_edges_weights(graph, edges):
        related_edges_weights = list()
        for e in edges:
            if e[0][:2] == e[1][:2] and e[0][:2] != GraphAnalyzer.UNRELATED:
                related_edges_weights.append(graph[e[0]][e[1]]['weight'])
        return related_edges_weights

    def extract_verified_sources(self, sources_file_name):
        with open(sources_file_name) as f:
            sources = [l.strip('\n') for l in f.readlines()]
        nodes = list()
        for s in sources:
            for node in self.outbreaks_nodes[s[:2]]:
                if node[:len(s)] == s:
                    nodes.append(node)
                    continue
        return nodes

    def how_many_edges_directions_is_correct(self, source):
        true_positive = 0
        for d in self.outbreaks_nodes[source[:2]]:
            if d == source:
                continue
            if self.graph[source][d]['weight'] < self.graph[d][source]['weight']:
                true_positive += 1
        return true_positive

    def check_sources(self, verified_sources_file_name):
        true_positive = list()
        total = list()
        verified_sources = self.extract_verified_sources(verified_sources_file_name)
        for s in verified_sources:
            true_positive.append(self.how_many_edges_directions_is_correct(s))
            total.append(len(self.outbreaks_nodes[s[:2]]))
        return [float(x[0])/x[1] for x in zip(true_positive, total)]


class DirectedGraphAnalyzer(GraphAnalyzer):
    def __init__(self, edges_list_file_name):
        super(DirectedGraphAnalyzer, self).__init__(edges_list_file_name,
                                                    nx.read_edgelist(edges_list_file_name, create_using=nx.DiGraph()))
        self.transmission_edges = self.get_transmission_edges(self.graph)
        self.dist_diff = self.get_dist_diff(self.graph)
        self.unrelated_edges_weights = self.get_unrelated_edges_weights(self.graph, self.transmission_edges)
        self.related_edges_weights = self.get_related_edges_weights(self.graph, self.transmission_edges)

    @staticmethod
    def get_transmission_edges(graph):
        transmission_edges = list()
        nodes = graph.nodes()
        for i in range(len(nodes)-1):
            for j in range(i+1, len(nodes)):
                transmission_edges.append((nodes[i], nodes[j])
                                          if graph[nodes[i]][nodes[j]]['weight'] < graph[nodes[j]][nodes[i]]['weight']
                                          else (nodes[j], nodes[i]))
        return transmission_edges

    @staticmethod
    def get_dist_diff(graph):
        dist_diff = list()
        nodes = graph.nodes()
        for i in range(len(nodes) - 1):
            for j in range(i + 1, len(nodes)):
                diff = graph[nodes[i]][nodes[j]]['weight'] - graph[nodes[j]][nodes[i]]['weight']
                dist_diff.append((nodes[i], nodes[j], -diff)
                                 if diff < 0
                                 else (nodes[j], nodes[i], diff))
        return dist_diff


class UndirectedGraphAnalyzer(GraphAnalyzer):
    def __init__(self, edges_list_file_name):
        super(UndirectedGraphAnalyzer, self).__init__(edges_list_file_name,
                                                      nx.read_edgelist(edges_list_file_name))
        self.unrelated_edges_weights = self.get_unrelated_edges_weights(self.graph, self.graph.edges())
        self.related_edges_weights = self.get_related_edges_weights(self.graph, self.graph.edges())


#def parse_arguments():
#    arguments_parser = argparse.ArgumentParser('Outbreak graph analyzer.')
#    arguments_parser.add_argument('edges_list_file_name', help='file with edges of an outbreak graph')
#    return arguments_parser.parse_args()


def main():
#    args = parse_arguments()
    simulation_analyzer = DirectedGraphAnalyzer(SIMULATION_EDGE_LIST)
    min_dist_analyzer = UndirectedGraphAnalyzer(MIN_DIST_EDGE_LIST)

    y1, binEdges1 = np.histogram(simulation_analyzer.unrelated_edges_weights, 100, normed=False)
    y2, binEdges2 = np.histogram(simulation_analyzer.related_edges_weights, 40, normed=False)

    bincenters1 = 0.5*(binEdges1[1:]+binEdges1[:-1])
    bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])

    plt.plot(bincenters1, y1, 'g')
    plt.plot(bincenters2, y2, 'r')
    plt.show()

    y1, binEdges1 = np.histogram(min_dist_analyzer.unrelated_edges_weights, 100, normed=False)
    y2, binEdges2 = np.histogram(min_dist_analyzer.related_edges_weights, 40, normed=False)

    bincenters1 = 0.5 * (binEdges1[1:] + binEdges1[:-1])
    bincenters2 = 0.5 * (binEdges2[1:] + binEdges2[:-1])

    plt.plot(bincenters1, y1, 'g')
    plt.plot(bincenters2, y2, 'r')
    plt.show()


#    plt.plot(sorted([d[2] for d in outbreak_analyzer.dist_diff]))
#    plt.ylabel('diff')
#    plt.show()
#    plt.plot(sorted(outbreak_analyzer.unrelated_edges_weights))
#    plt.ylabel('unrelated_weights')
#    plt.show()

#    plt.plot(sorted(outbreak_analyzer.related_edges_weights))
#    plt.ylabel('related_weights')
#    plt.show()


if __name__ == '__main__':
    main()
