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


class SimulationResultsParser(object):
    def __init__(self, dir, pair_name):
        self.dir = dir
        self.pair_name = pair_name
        self.nodes = pair_name.split('_vs_')
        self.edges = self.parse_edges(dir, pair_name)

    def parse_edges(self, dir, pair_name):
        distances = {self.nodes[0]: [], self.nodes[1]: []}
        with open(os.path.join(dir, pair_name, 'results.log')) as f:
            f.readline()
            for l in f.readlines():
                c = l.split()
                distances[c[0]].append(int(c[2]))
        return [(self.nodes[0], self.nodes[1], statistics.median(distances[self.nodes[1]])),
                (self.nodes[1], self.nodes[0], statistics.median(distances[self.nodes[0]]))]


def build_graph(outbreak_dir_name):
    graph = nx.DiGraph()
    out_dir = os.path.join(outbreak_dir_name, 'out')
    pairs = os.listdir(out_dir)
    edges_list = list()
    for pair in pairs:
        pair_parser = SimulationResultsParser(out_dir, pair)
        edges_list.extend(pair_parser.edges)
    graph.add_weighted_edges_from(edges_list)
    return graph


def parse_arguments():
    parser = argparse.ArgumentParser("Outbreak graph builder.")
    parser.add_argument('-i', dest='outbreak_dir_name', type=str, required=True,
                        help='Folder with simulation results.')
    parser.add_argument('-o', dest='graph_base_name', type=str, default=None,
                        help='Name of base of out file where json graph will be written.')
    return parser.parse_args()


def determine_graph_file_name(base_file_name):
    return base_file_name if base_file_name.split('.')[-1] == 'gz' else base_file_name + '.gz'


def main():
    args = parse_arguments()
#    out_graph_file_name = determine_graph_file_name(args.graph_base_name)
    graph = build_graph(args.outbreak_dir_name)
    nx.write_edgelist(graph, args.graph_base_name)

if __name__ == '__main__':
    main()
