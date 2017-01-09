#!/usr/bin/env python3
"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/1/17
"""

#import argparse
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


SOURCES_FILE_NAME = 'data/sources.txt'
CENTRALITY_GRAPH = 'data/anti_centrality_graph.txt'
ANTI_CENTRALITY_GRAPH = 'data/anti_centrality_graph.txt'
MIN_DIST_EDGE_LIST = 'data/all_clipped_min_dist.txt'
SIMULATION_EDGE_LIST_FULL = 'data/all_clipped_full_graph_simulation.txt'
SIMULATION_EDGE_LIST_FILTERES = 'data/all_clipped_filtered_simulation.txt'
SIMULATION_EDGE_LIST = 'data/all_clipped_simulation.txt'
VIRIFIED_OUTBREAKS = ['AA', 'AC', 'AI', 'AJ', 'AW', 'BA', 'BB', 'BC', 'BJ']
GRAPH = SIMULATION_EDGE_LIST


class SourceFinder(object):
    def __init__(self, outbreak_graph):
        self.outbreak_graph = outbreak_graph

    def get_sum_of_weights_of_outgoing_edges(self, vertex):
        out_edges_weight_sum = 0
        for e in self.outbreak_graph.out_edges(vertex):
            out_edges_weight_sum = self.outbreak_graph[e[0]][e[1]]['weight']
        return out_edges_weight_sum

    def get_sum_of_outgoing_edges(self, vertex):
        DIRECTION_COEF = 1.5
        out_edges_sum = 0
        for e in self.outbreak_graph.out_edges(vertex):
            if self.outbreak_graph[e[0]][e[1]]['weight'] * DIRECTION_COEF < self.outbreak_graph[e[1]][e[0]]['weight']:
                out_edges_sum -= 1
            elif self.outbreak_graph[e[0]][e[1]]['weight'] > DIRECTION_COEF * self.outbreak_graph[e[1]][e[0]]['weight']:
                out_edges_sum += 1
        return out_edges_sum

    def get_sum_of_outgoing_edges1(self, vertex):
        out_edges_sum = 0
        for e in self.outbreak_graph.out_edges(vertex):
            out_edges_sum += self.outbreak_graph[e[0]][e[1]]['weight'] - self.outbreak_graph[e[1]][e[0]]['weight']
        return out_edges_sum

    def get_shortest_path_tree(self, vertex):
        edges_dict = dict()
        dijkstra_paths = nx.single_source_dijkstra_path(self.outbreak_graph, vertex)
        for dest in dijkstra_paths:
            v = vertex
            for w in dijkstra_paths[dest]:
                if v != w and (v, w) not in edges_dict:
                    edges_dict[(v, w)] = self.outbreak_graph[v][w]['weight']
                v = w
        shortest_path_tree = nx.DiGraph()
        shortest_path_tree.add_weighted_edges_from([(x[0][0], x[0][1], x[1]) for x in edges_dict.items()])
        return shortest_path_tree

    def get_weight_of_shortest_path_tree(self, vertex):
        shortes_paths_tree_cost = 0
        shortes_paths_tree = self.get_shortest_path_tree(vertex)
        for e in shortes_paths_tree.edges():
            shortes_paths_tree_cost += shortes_paths_tree[e[0]][e[1]]['weight']
        return shortes_paths_tree_cost

    def find_source(self, f):
        source_cost = float('inf')
        source = None
        for v in self.outbreak_graph.nodes():
            new_source_cost = f(v)
            print('{0} cost: {1}'.format(v, new_source_cost))
            if new_source_cost < source_cost:
                source = v
                source_cost = new_source_cost
        return source

    def find_source_by_star_sum_of_edges(self):
        return self.find_source(self.get_sum_of_outgoing_edges1)

    def find_source_by_star_total_cost(self):
        return self.find_source(self.get_sum_of_weights_of_outgoing_edges)

    def find_sorce_by_shortest_path_tree(self):
        return self.find_source(self.get_weight_of_shortest_path_tree)

    def get_source_true_positive(self, outbreak_verified_source):
        verified_source_node = list(
            filter(lambda x: x.split('_')[0] == outbreak_verified_source, self.outbreak_graph.nodes()))[0]
        #verified_source_node = verified_source_node[0]
        out_edges = self.outbreak_graph.out_edges(verified_source_node)
        true_determined = len(list(
            filter(lambda v: self.outbreak_graph[v[0]][v[1]]['weight'] < self.outbreak_graph[v[1]][v[0]]['weight'],
                   out_edges)))
        return true_determined, len(out_edges)


class GraphAnalyzer(object):
    UNRELATED = 'XX'

    def __init__(self, edges_list_file_name, graph):
        self.edges_list_file_name = edges_list_file_name
        self.graph = graph
        self.outbreaks_nodes_dict = self.infer_outbreaks_nodes(self.graph)
        self.unrelated_edges = self.get_unrelated_edges(self.graph, self.graph.edges())
        self.related_edges = self.get_related_edges(self.graph, self.graph.edges())
#        self.thresholds = self.get_thresholds(self.graph)

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
    def get_unrelated_edges(graph, edges):
        unrelated_edges = list()
        for e in edges:
            if e[0][:2] != e[1][:2] or e[0][:2] == GraphAnalyzer.UNRELATED:
                unrelated_edges.append((e[0], e[1], graph[e[0]][e[1]]['weight']))
        return unrelated_edges

    @staticmethod
    def get_related_edges(graph, edges):
        related_edges = list()
        for e in edges:
            if e[0][:2] == e[1][:2] and e[0][:2] != GraphAnalyzer.UNRELATED:
                related_edges.append((e[0], e[1], graph[e[0]][e[1]]['weight']))
        return related_edges

    def extract_verified_sources(self, sources_file_name):
        with open(sources_file_name) as f:
            sources = [l.strip('\n') for l in f.readlines()]
        nodes = list()
        for s in sources:
            for node in self.outbreaks_nodes_dict[s[:2]]:
                if node[:len(s)] == s:
                    nodes.append(node)
                    continue
        return nodes

    def how_many_edges_directions_is_correct(self, source):
        true_positive = 0
        for d in self.outbreaks_nodes_dict[source[:2]]:
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
            total.append(len(self.outbreaks_nodes_dict[s[:2]]))
        return [float(x[0])/x[1] for x in zip(true_positive, total)]

    def get_zero_type_I_error_thresholds(self):
        '''
        min_unrelated_edge_weight = min([e[2] for e in self.unrelated_edges])
        max_related_less_min_unrelated = 0
        for e in self.related_edges:
            if e[2] <= min_unrelated_edge_weight:
                if e[2] > max_related_less_min_unrelated:
                    max_related_less_min_unrelated = e[2]
        return np.mean([min_unrelated_edge_weight, max_related_less_min_unrelated])
        '''
        return min([e[2] for e in self.unrelated_edges])

    @staticmethod
    def get_outbreak_edges(graph, outbreak_name):
        edges = list()
        for e in graph.edges():
            if e[0][:2] == outbreak_name and e[1][:2] == outbreak_name:
                edges.append((e[0], e[1], graph[e[0]][e[1]]['weight']))
        return edges

    def find_bridge_length(self, outbreak_name):
        graph = nx.Graph()
        graph.add_nodes_from(self.outbreaks_nodes_dict[outbreak_name])
        sorted_outbreak_edges = sorted(self.get_outbreak_edges(self.graph, outbreak_name), key=lambda x: x[2])
        for e in sorted_outbreak_edges:
            graph.add_weighted_edges_from([e])
            if nx.is_connected(graph):
                return e[2]
        return None

    def get_zero_type_II_error_threshold(self):
        threshold = 0
        for k in self.outbreaks_nodes_dict.keys():
            if k[:2] != self.UNRELATED:
                new_threshold = self.find_bridge_length(k)
                if new_threshold > threshold:
                    threshold = new_threshold
        return threshold


class DirectedGraphAnalyzer(GraphAnalyzer):
    def __init__(self, edges_list_file_name):
        self.di_graph = nx.read_edgelist(edges_list_file_name, create_using=nx.DiGraph())
        graph = nx.Graph()
        graph.add_weighted_edges_from(self.get_transmission_edges(self.di_graph))
        super(DirectedGraphAnalyzer, self).__init__(edges_list_file_name, graph)
        self.dist_diff = self.get_dist_diff(self.graph)

    def get_outbreak_graph(self, outbreak_name):
        graph = nx.DiGraph()
        graph.add_weighted_edges_from(self.get_outbreak_edges(self.di_graph, outbreak_name))
        return graph

    @staticmethod
    def get_transmission_edges(graph):
        transmission_edges = list()
        nodes = graph.nodes()
        for i in range(len(nodes)-1):
            for j in range(i+1, len(nodes)):
                transmission_edges.append((nodes[i], nodes[j], graph[nodes[i]][nodes[j]]['weight'])
                                          if graph[nodes[i]][nodes[j]]['weight'] < graph[nodes[j]][nodes[i]]['weight']
                                          else (nodes[j], nodes[i], graph[nodes[j]][nodes[i]]['weight']))
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


def get_outbreak_verified_sources(file_name):
    with open(file_name) as f:
        sources = [l.strip('\n') for l in f.readlines()]
    return dict(zip([x[:2] for x in sources], sources))


#def parse_arguments():
#    arguments_parser = argparse.ArgumentParser('Outbreak graph analyzer.')
#    arguments_parser.add_argument('edges_list_file_name', help='file with edges of an outbreak graph')
#    return arguments_parser.parse_args()


def main():
#    args = parse_arguments()
    simulation_analyzer = DirectedGraphAnalyzer(GRAPH)
    min_dist_analyzer = UndirectedGraphAnalyzer(MIN_DIST_EDGE_LIST)
    print(simulation_analyzer.get_zero_type_I_error_thresholds())
    print(simulation_analyzer.get_zero_type_II_error_threshold())

    print(min_dist_analyzer.get_zero_type_I_error_thresholds())
    print(min_dist_analyzer.get_zero_type_II_error_threshold())

    y1, binEdges1 = np.histogram([e[2] for e in simulation_analyzer.unrelated_edges], 100, normed=False)
    y2, binEdges2 = np.histogram([e[2] for e in simulation_analyzer.related_edges], 40, normed=False)

    bincenters1 = 0.5*(binEdges1[1:]+binEdges1[:-1])
    bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])

    plt.plot(bincenters1, y1, 'g')
    plt.plot(bincenters2, y2, 'r')
    plt.show()

    y1, binEdges1 = np.histogram([e[2] for e in min_dist_analyzer.unrelated_edges], 100, normed=False)
    y2, binEdges2 = np.histogram([e[2] for e in min_dist_analyzer.related_edges], 40, normed=False)

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


def main1():
    simulation_analyzer = DirectedGraphAnalyzer(GRAPH)
    outbreak_verified_sources = get_outbreak_verified_sources(SOURCES_FILE_NAME)

#    print("Source by cost")
#    for o in VIRIFIED_OUTBREAKS:
#        outbreak_graph = SourceFinder(simulation_analyzer.get_outbreak_graph(o))
#        print(outbreak_graph.find_source_by_star_total_cost())
    print("Source by edges number")
    total_found_sources = 0
    total_nodes_count = 0
    for o in VIRIFIED_OUTBREAKS:
        print('{0}:'.format(o))
        outbreak_graph = SourceFinder(simulation_analyzer.get_outbreak_graph(o))
        outbreak_source = outbreak_graph.find_sorce_by_shortest_path_tree()
        found_sources, nodes_count = outbreak_graph.get_source_true_positive(outbreak_verified_sources[o])
        print(outbreak_source)
        total_found_sources += found_sources
        total_nodes_count += nodes_count

    print(float(total_found_sources/total_nodes_count))


def report_wrong_directions(simulation_analyzer, outbreak_verified_sources):
    print("Errors in source detection:")
    print("1) Before error correction:")

    total_found_sources = 0
    total_nodes_count = 0
    for o in VIRIFIED_OUTBREAKS:
        print('{0}:'.format(o))
        outbreak_graph = SourceFinder(simulation_analyzer.get_outbreak_graph(o))
        outbreak_source = outbreak_graph.find_sorce_by_shortest_path_tree()
        found_sources, nodes_count = outbreak_graph.get_source_true_positive(outbreak_verified_sources[o])
        print(outbreak_source)
        total_found_sources += found_sources
        total_nodes_count += nodes_count
    print("True positive: {0}".format(float(total_found_sources / total_nodes_count)))
    print("Wrong directions:")


def main2():
    simulation_analyzer = DirectedGraphAnalyzer(GRAPH)
    outbreak_verified_sources = get_outbreak_verified_sources(SOURCES_FILE_NAME)

    report_wrong_directions(simulation_analyzer, outbreak_verified_sources)


if __name__ == '__main__':
#    main()
#    main1()
    main2()
