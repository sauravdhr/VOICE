#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 01/03/2017
"""
import argparse
import networkx as nx
import os
import statistics
import glob
import random


K_MAX_DEFAULT = 4


class ResultsParser(object):
    def __init__(self, dir, pair_name):
        self.dir = dir
        self.pair_name = pair_name
        self.nodes = pair_name.split('_vs_')
        self.distances = {self.nodes[0]: [], self.nodes[1]: []}

    def parse_edges(self):
        pass


class SimulationResultsParser(ResultsParser):
    def __init__(self, dir, pair_name):
        super(SimulationResultsParser, self).__init__(dir, pair_name)
        self.edges = self.parse_edges()

    def parse_edges(self):
        with open(os.path.join(self.dir, self.pair_name, 'results.log')) as f:
            f.readline()
            for l in f.readlines():
                c = l.split()
                self.distances[c[0]].append(int(c[2]))
        return [(self.nodes[0], self.nodes[1], statistics.median(self.distances[self.nodes[1]])),
                (self.nodes[1], self.nodes[0], statistics.median(self.distances[self.nodes[0]]))]


class SimulationResultsParserMonteCarlo(ResultsParser):
    MONTE_CARLO_RUNS = 501

    def __init__(self, dir, pair_name, k_max):
        super(SimulationResultsParserMonteCarlo, self).__init__(dir, pair_name)
        self.k_max = k_max
        self.edges = self.parse_edges()

    def parse_edges(self):
        return [(self.nodes[0], self.nodes[1], statistics.median(self.get_distances(self.nodes[1]))),
                (self.nodes[1], self.nodes[0], statistics.median(self.get_distances(self.nodes[0])))]

    def get_distances(self, recipient):
        simulation_logs = glob.glob(os.path.join(self.dir, self.pair_name, 'simulation', recipient + '*'))
        distances = []
        for l in simulation_logs:
            distances.append(self.get_distance_monte_carlo(l))
        return distances

    def get_distance_monte_carlo(self, recipient_log):
        nodes_infection_times_sorted = self.get_node_infection_times_sorted(recipient_log)
        return statistics.median([self.run_montecarlo_tact(nodes_infection_times_sorted)
                                  for _ in range(self.MONTE_CARLO_RUNS)])

    def run_montecarlo_tact(self, nodes_infection_times_sorted):
        if not nodes_infection_times_sorted:
            return 0
        if len(nodes_infection_times_sorted) <= self.k_max:
            return nodes_infection_times_sorted[-1]
        ids_shuffled = list(range(len(nodes_infection_times_sorted)))
        random.shuffle(ids_shuffled)
        return nodes_infection_times_sorted[max(ids_shuffled[:self.k_max])]

    @staticmethod
    def get_node_infection_times_sorted(recipient_log):
        nodes_infection_times = []
        with open(recipient_log) as f:
            for l in f.readlines():
                t = int(l.split(', ')[1])
                if t != 0:
                    nodes_infection_times.append(t)
        return nodes_infection_times


def build_graph(outbreak_dir_name, k_max):
    graph = nx.DiGraph()
    out_dir = os.path.join(outbreak_dir_name, 'out')
    pairs = os.listdir(out_dir)
    edges_list = list()
    for pair in pairs:
        if not k_max:
            pair_parser = SimulationResultsParser(out_dir, pair)
        else:
            pair_parser = SimulationResultsParserMonteCarlo(out_dir, pair, k_max)
        edges_list.extend(pair_parser.edges)
    graph.add_weighted_edges_from(edges_list)
    return graph


def parse_arguments():
    parser = argparse.ArgumentParser("Outbreak graph builder.")
    parser.add_argument('-i', dest='outbreak_dir_name', type=str, required=True,
                        help='Folder with simulation results.')
    parser.add_argument('-o', dest='graph_base_name', type=str, default=None,
                        help='Name of base of out file where json graph will be written.')
    parser.add_argument('-m', dest='k_max', type=int, default=K_MAX_DEFAULT,
                        help='Count of sequences for Monte-Carlo normalization.')
    return parser.parse_args()


def determine_graph_file_name(base_file_name):
    return base_file_name if base_file_name.split('.')[-1] == 'gz' else base_file_name + '.gz'


def main():
    args = parse_arguments()
#    out_graph_file_name = determine_graph_file_name(args.graph_base_name)
    graph = build_graph(args.outbreak_dir_name, args.k_max)
    nx.write_edgelist(graph, args.graph_base_name)

if __name__ == '__main__':
    main()
