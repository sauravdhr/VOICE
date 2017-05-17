#!/usr/bin/env python3
"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/1/17
"""

#import argparse
import statistics
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import utils.roc2img
import sklearn.metrics
import csv

COLORS = ['red', 'green', 'blue']
SOURCES_FILE_NAME = 'data/sources.txt'
CENTRALITY_GRAPH = 'data/anti_centrality_graph.txt'
ANTI_CENTRALITY_GRAPH = 'data/anti_centrality_graph.txt'
#MIN_DIST_EDGE_LIST = 'data/all_clipped_min_dist.txt'
#MIN_DIST_EDGE_LIST = 'data/all_clipped_min_dist_vicinity3.txt'
MIN_DIST_PLUS_BORDER_EDGE_LIST = 'data/all_clipped_min_dist_v_4_log_3.txt'
MIN_DIST_EDGE_LIST = 'data/all_clipped_min_dist.txt'
SIMULATION_EDGE_LIST_FULL = 'data/all_clipped_full_graph_simulation.txt'
MONTECARLO_10 = 'results/montecarlo_10.txt'
MONTECARLO_4 = 'results/montecarlo_4.txt'
SIMULATION_EDGE_LIST_FILTERES = 'data/all_clipped_filtered_simulation.txt'
SIMULATION_EDGE_LIST_BORDER_CONSENSUS = 'data/all_clipped_border_consensus_simulation_graph.txt'
SIMULATION_EDGE_LIST = 'data/all_clipped_simulation.txt'
AW = 'results/AW.txt'
AQ = 'results/AQ.txt'
#VERIFIED_OUTBREAKS = ['AA', 'AC', 'AI', 'AJ', 'AW', 'BA', 'BB', 'BC', 'BJ', 'AQ']
#VERIFIED_OUTBREAKS = ['AA', 'AC', 'AI', 'AJ', 'AW', 'BA', 'BB', 'BC', 'BJ']
VERIFIED_OUTBREAKS = ['AW']
#VERIFIED_OUTBREAKS = ['AQ']
#GRAPH = MONTECARLO_10
GRAPH = SIMULATION_EDGE_LIST
#GRAPH = AW


class SourceFinder(object):
    def __init__(self, outbreak_graph):
        self.outbreak_graph = outbreak_graph
        self.centrality = None

    def get_sum_of_weights_of_outgoing_edges(self, vertex):
        out_edges_weight_sum = 0
        for e in self.outbreak_graph.out_edges(vertex):
            out_edges_weight_sum = self.outbreak_graph[e[0]][e[1]]['weight']
        return out_edges_weight_sum

    def get_sum_of_weights_of_edges(self, vertex):
        out_edges_weight_sum = 0
        for e in self.outbreak_graph.edges(vertex):
            out_edges_weight_sum = self.outbreak_graph[e[0]][e[1]]['weight']
        return out_edges_weight_sum

    def get_median_of_weights_of_edges(self, vertex):
        out_edges_weight_sum = 0
        weights = list()
        for e in self.outbreak_graph.edges(vertex):
            weights.append(self.outbreak_graph[e[0]][e[1]]['weight'])
        return statistics.median(weights)

    def get_number_of_recipients(self, vertex):
        number_of_recipients = 0
        for e in self.outbreak_graph.edges(vertex):
            if self.outbreak_graph[e[0]][e[1]]['weight'] < self.outbreak_graph[e[1]][e[0]]['weight']:
                number_of_recipients += 1
        return number_of_recipients

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
        print("Vertex: {0}, spt: {1}".format(vertex, shortes_paths_tree_cost))
        return shortes_paths_tree_cost

    def find_source(self, f):
        source_cost = float('inf')
        source = None
        for v in self.outbreak_graph.nodes():
            new_source_cost = f(v)
#            print('{0} cost: {1}'.format(v, new_source_cost))
            if new_source_cost < source_cost:
                source = v
                source_cost = new_source_cost
        return source

    def find_source_by_star_total_cost(self):
        return self.find_source(self.get_sum_of_weights_of_edges)

    def find_source_by_star_median(self):
        return self.find_source(self.get_median_of_weights_of_edges)

    def find_source_by_shortest_path_tree(self):
        return self.find_source(self.get_weight_of_shortest_path_tree)

    def find_source_by_number_of_recipients(self):
        return self.find_source(self.get_number_of_recipients)

    def find_source_by_centrality(self):
#        self.centrality = nx.eigenvector_centrality_numpy(self.outbreak_graph)
        self.centrality = nx.closeness_centrality(self.outbreak_graph)
#        self.centrality = nx.betweenness_centrality(self.outbreak_graph)
#        self.centrality = nx.current_flow_closeness_centrality(self.outbreak_graph)
        return self.find_source(self.get_node_centrality)

    def get_node_centrality(self, v):
        return self.centrality[v]
#        return -self.centrality[v]

    def get_true_directions(self, outbreak_verified_source):
        verified_source_node = list(
            filter(lambda x: x.split('_')[0] == outbreak_verified_source, self.outbreak_graph.nodes()))[0]
        #verified_source_node = verified_source_node[0]
        out_edges = self.outbreak_graph.out_edges(verified_source_node)
        true_determined = len(list(
            filter(lambda v: self.outbreak_graph[v[0]][v[1]]['weight'] < self.outbreak_graph[v[1]][v[0]]['weight'],
                   out_edges)))
        wrong_determined = list(
            filter(lambda v: self.outbreak_graph[v[0]][v[1]]['weight'] >= self.outbreak_graph[v[1]][v[0]]['weight'],
                   out_edges))
        for i in wrong_determined:
            print("{0}-{1}".format(i[0].split("_")[0], i[1].split("_")[0]))
        return true_determined, len(out_edges)


class GraphAnalyzer(object):
    UNRELATED = 'XX'

    def __init__(self, edges_list_file_name, graph):
        self.edges_list_file_name = edges_list_file_name
        self.graph = graph
        self.outbreaks_nodes_dict = self.infer_outbreaks_nodes(self.graph)
        self.unrelated_edges = self.get_unrelated_edges(self.graph, self.graph.edges())
        self.related_edges = self.get_related_edges(self.graph, self.graph.edges())

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

    def get_zero_type_II_error_threshold_for_pairs(self):
        threshold = 0
        for e in self.related_edges:
            w = min(self.graph[e[0]][e[1]]['weight'], self.graph[e[1]][e[0]]['weight'])
            if w > threshold:
                threshold = w
        return threshold

    def get_false_positive_rate_pairs(self, thr):
        return len(self.get_false_related_edges(thr))/len(self.unrelated_edges)

    def get_true_positive_rate_pairs(self, thr):
        return len(self.get_true_related_edges(thr))/len(self.related_edges)

    def get_number_of_clusters(self, threshold):
        clusters_count = 0
        for k in self.outbreaks_nodes_dict.keys():
            if k[:2] != self.UNRELATED:
                clusters_count += self.get_number_of_components(k, threshold)
        return clusters_count

    def get_number_of_components(self, outbreak_name, threshold):
        graph = nx.Graph()
        graph.add_nodes_from(self.outbreaks_nodes_dict[outbreak_name])
        outbreak_edges = list(filter(lambda x: x[2] < threshold, self.get_outbreak_edges(self.graph, outbreak_name)))
        graph.add_weighted_edges_from(outbreak_edges)
        return nx.number_connected_components(graph)

#    def get_true_related_number(self, threshold):
#        return sum(map(lambda x: x[2] >= threshold, self.unrelated_edges)), len(self.unrelated_edges)

    def get_false_related_edges(self, threshold):
        return list(filter(lambda x: x[2] <= threshold, self.unrelated_edges))

    def get_true_related_edges(self, threshold):
        return list(filter(lambda x: x[2] <= threshold, self.related_edges))

    def get_related_error(self, threshold):
        return sum(map(lambda x: x[2] > threshold, self.related_edges)), len(self.related_edges)


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

    def get_outbreak_graph(self, outbreak_name):
        graph = nx.Graph()
        graph.add_weighted_edges_from(self.get_outbreak_edges(self.graph, outbreak_name))
        return graph


def get_outbreak_verified_sources(file_name):
    with open(file_name) as f:
        sources = [l.strip('\n') for l in f.readlines()]
    return dict(zip([x[:2] for x in sources], sources))


def depict_ROC_curve(rocs, colors, labels, fname, randomline):
    plt.figure(figsize=(4, 4), dpi=80)
    utils.roc2img.SetupROCCurvePlot(plt)
    for i in range(len(rocs)):
        plt.plot(rocs[i][0], rocs[i][1], color=colors[i], linewidth=2, label=labels[i])
    utils.roc2img.SaveROCCurvePlot(plt, fname, randomline)


def get_rocs(analyzers):
    roc = [[] for _ in range(len(analyzers))]
    for i in range(len(analyzers)):
        x = [0.0]
        y = [0.0]
        thr = 0
        tpr = 0
        while tpr != 1:
            tpr = analyzers[i].get_true_positive_rate_pairs(thr)
            fpr = analyzers[i].get_false_positive_rate_pairs(thr)
            x.append(fpr)
            y.append(tpr)
            thr += 1
        x.append(1.0)
        y.append(1.0)
        roc[i] = (x, y)
    return roc


def get_aucs(rocs):
    aucs = []
    for roc in rocs:
        aucs.append(sklearn.metrics.auc(roc[0], roc[1]))
    return aucs


def draw_ROC(analyzers, colors, labels, fname, randomline):
    rocs = get_rocs(analyzers)
    aucs = get_aucs(rocs)
    for i in range(len(labels)):
        labels[i] += ' (auc = %0.3f)' % aucs[i]
    depict_ROC_curve(rocs, colors, labels, fname, randomline)


def main0():
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

def report_source_finding_quality(graph_analyzer, outbreak_verified_sources, source_finding_method):
    found_sources_count = 0
    for o in VERIFIED_OUTBREAKS:
        outbreak_source = None
        outbreak_graph = SourceFinder(graph_analyzer.get_outbreak_graph(o))
        if source_finding_method == 'spt':
            outbreak_source = outbreak_graph.find_source_by_shortest_path_tree()
        elif source_finding_method == 'centrality':
            outbreak_source = outbreak_graph.find_source_by_centrality()
        elif source_finding_method == 'star':
            outbreak_source = outbreak_graph.find_source_by_star_median()
        elif source_finding_method == 'recipients':
            outbreak_source = outbreak_graph.find_source_by_number_of_recipients()
        print("Source: found - {0}, real - {1}".format(outbreak_source.split('_')[0], outbreak_verified_sources[o]))
        if outbreak_source.split('_')[0] == outbreak_verified_sources[o]:
            found_sources_count += 1
    print("#######")
    print("Found true sources: {0} of {1}".format(found_sources_count, len(VERIFIED_OUTBREAKS)))
    print("True rate: {0}".format(float(found_sources_count) / len(VERIFIED_OUTBREAKS)))
    print("#######")


def report_direction_finding_quality(simulation_analyzer, outbreak_verified_sources):
    total_true_directions = 0
    total_all_directions = 0
    for o in VERIFIED_OUTBREAKS:
        outbreak_graph = SourceFinder(simulation_analyzer.get_outbreak_graph(o))
        true_directions, all_directions = outbreak_graph.get_true_directions(outbreak_verified_sources[o])
        total_true_directions += true_directions
        total_all_directions += all_directions
    print("#######")
    print("Found true directions: {0} of {1}".format(total_true_directions, total_all_directions))
    print("True rate: {0}".format(float(total_true_directions) / total_all_directions))
    print("#######")


def report_relatedness(simulation_analyzer, min_dist_analyzer):
    zero_unrelated_thr_simulation = simulation_analyzer.get_zero_type_I_error_thresholds()
    zero_broken_outbreaks_thr_simulation = simulation_analyzer.get_zero_type_II_error_threshold()

    zero_unrelated_thr_min_dist = min_dist_analyzer.get_zero_type_I_error_thresholds()
    zero_broken_outbreaks_thr_min_dist = min_dist_analyzer.get_zero_type_II_error_threshold()

    print("-------")
    print("Thresholds for simulation:")
    print("No unrelated edges: {0}".format(zero_unrelated_thr_simulation))
    print("No outbreak breakage: {0}".format(zero_broken_outbreaks_thr_simulation))
    print("-------")
    print("Thresholds for min dist method:")
    print("No unrelated edges: {0}".format(zero_unrelated_thr_min_dist))
    print("No outbreak breakage: {0}".format(zero_broken_outbreaks_thr_min_dist))
    print("-------")

    number_of_clusters_simulation = simulation_analyzer.get_number_of_clusters(zero_unrelated_thr_simulation)
    true_related_number_simulation = simulation_analyzer.get_true_related_number(zero_broken_outbreaks_thr_simulation)
    related_error_simulation = simulation_analyzer.get_related_error(zero_unrelated_thr_simulation)

    number_of_clusters_min_dist = min_dist_analyzer.get_number_of_clusters(zero_unrelated_thr_min_dist)
    true_related_number_min_dist = min_dist_analyzer.get_true_related_number(zero_broken_outbreaks_thr_min_dist)
    related_error_min_dist = min_dist_analyzer.get_related_error(zero_unrelated_thr_min_dist)

    print("#######")
    print('Relatedness correctness for simulations:')
    print("#######")
    print('Related error with thr {0}: {1} of {2}'.format(zero_unrelated_thr_simulation,
                                                             related_error_simulation[0],
                                                             related_error_simulation[1]))
    print("Sensitivity: {0}".format(1 - float(related_error_simulation[0]) / related_error_simulation[1]))
    print("-------")
    print('Number of clusters with thr {0}: {1}'.format(zero_unrelated_thr_simulation,
                                                        number_of_clusters_simulation))
    print("-------")
    print('Number of true positive related pairs with thr {0}: {1} of {2}'.format(zero_broken_outbreaks_thr_simulation,
                                                                                  true_related_number_simulation[0],
                                                                                  true_related_number_simulation[1]))
    print("True rate: {0}".format(float(true_related_number_simulation[0]) / true_related_number_simulation[1]))
    print("#######")
    print('Relatedness correctness for min dist method:')
    print("#######")
    print('Related error with thr {0}: {1} of {2}'.format(zero_unrelated_thr_min_dist,
                                                          related_error_min_dist[0],
                                                          related_error_min_dist[1]))
    print("Sensitivity: {0}".format(1 - float(related_error_min_dist[0]) / related_error_min_dist[1]))
    print("-------")
    print('Number of clusters with thr {0}: {1}'.format(zero_unrelated_thr_min_dist,
                                                        number_of_clusters_min_dist))
    print("-------")
    print('Number of true positive related pairs with thr {0}: {1} of {2}'.format(zero_broken_outbreaks_thr_min_dist,
                                                                                  true_related_number_min_dist[0],
                                                                                  true_related_number_min_dist[1]))

    print("True rate: {0}".format(float(true_related_number_min_dist[0]) / true_related_number_min_dist[1]))
    print("-------")

    print("Number of unrelated samples: {0}".format(len(simulation_analyzer.outbreaks_nodes_dict['XX'])))


def export_classifiers(analyzers, file_names):
    for i in range(len(analyzers)):
        with open(file_names[i], 'w') as csv_f:
            w = csv.writer(csv_f)
            w.writerow(['l', 'related'])
            for e in analyzers[i].related_edges:
                w.writerow([e[2], 1])
            for e in analyzers[i].unrelated_edges:
                w.writerow([e[2], 0])


def main():
    simulation_analyzer = DirectedGraphAnalyzer(GRAPH)
    min_dist_analyzer = UndirectedGraphAnalyzer(MIN_DIST_EDGE_LIST)
    min_dist_plus_border_analyzer = UndirectedGraphAnalyzer(MIN_DIST_PLUS_BORDER_EDGE_LIST)
    outbreak_verified_sources = get_outbreak_verified_sources(SOURCES_FILE_NAME)

#    report_direction_finding_quality(simulation_analyzer, outbreak_verified_sources)
#    print('Simulations')
#    report_source_finding_quality(simulation_analyzer, outbreak_verified_sources, 'recipients')
#    report_source_finding_quality(simulation_analyzer, outbreak_verified_sources, 'spt')
#    report_source_finding_quality(simulation_analyzer, outbreak_verified_sources, 'centrality')
#    report_source_finding_quality(simulation_analyzer, outbreak_verified_sources, 'star')
#    report_relatedness(simulation_analyzer, min_dist_analyzer)

#    draw_ROC([min_dist_analyzer, min_dist_plus_border_analyzer, simulation_analyzer],
#             ['red', 'green', 'blue'], ['min dist', 'min dist + border', 'simulation'], 'roc.png', True)

    export_classifiers([min_dist_analyzer, min_dist_plus_border_analyzer, simulation_analyzer],
                       ['min_dist.csv', 'min_dist_plus_border.csv', 'simulation.csv'])

#    thrII = min_dist_analyzer.get_zero_type_II_error_threshold_for_pairs()
#    print(thrII)
#    print(min_dist_analyzer.get_false_positive_rate_pairs(10))
#    print(min_dist_analyzer.get_true_positive_rate_pairs(10))

#    print('Min dist')
#    report_source_finding_quality(min_dist_analyzer, outbreak_verified_sources, 'star')


if __name__ == '__main__':
    main()
