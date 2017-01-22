#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/4/17
"""
import argparse
import networkx as nx
import os

import hamming_dist_graph
import network_creator
import math


INF = 10e6


class MinDistanceGraph(object):
    def __init__(self, outbreak_folder, vicinity, log_scale):
        self.outbreak_folder = outbreak_folder
        self.vicinity = vicinity
        self.log_scale = log_scale
        self.edges = self.get_edges(outbreak_folder, self.vicinity, self.log_scale)
        self.graph = nx.Graph()
        self.graph.add_weighted_edges_from(self.edges)

    @staticmethod
    def get_edges(outbreak_folder, vicinity, log_scale):
        edges = list()
        fasta_files_names = os.listdir(outbreak_folder)
        for i in range(len(fasta_files_names)-1):
            for j in range(i+1, len(fasta_files_names)):
                min_dist = MinDistanceGraph.get_min_dist(os.path.join(outbreak_folder, fasta_files_names[i]),
                                                         os.path.join(outbreak_folder, fasta_files_names[j]))
                vicinity = MinDistanceGraph.get_vicinity(os.path.join(outbreak_folder, fasta_files_names[i]),
                                                         os.path.join(outbreak_folder, fasta_files_names[j]),
                                                         min_dist - vicinity)
                edges.append((fasta_files_names[i], fasta_files_names[j], min_dist + log_scale * math.log(vicinity)))
            print("{0} is done out of {1}".format(i+1, len(fasta_files_names)-1))
        return edges

    @staticmethod
    def get_min_dist(fasta1_file_name, fasta2_file_name):
        seqs1 = network_creator.parse_fasta(fasta1_file_name)
        seqs2 = network_creator.parse_fasta(fasta2_file_name)
        min_dist = INF
        for i in range(len(seqs1)):
            for j in range(len(seqs2)):
                dist = hamming_dist_graph.hamming_distance(seqs1[i], seqs2[j])
                if dist < min_dist:
                    min_dist = dist
        return min_dist

    @staticmethod
    def get_vicinity(fasta1_file_name, fasta2_file_name, max_edge_length):
        seqs1 = network_creator.parse_fasta(fasta1_file_name)
        seqs2 = network_creator.parse_fasta(fasta2_file_name)
        vicinity = 0
        for i in range(len(seqs1)):
            for j in range(len(seqs2)):
                dist = hamming_dist_graph.hamming_distance(seqs1[i], seqs2[j])
                if dist >= max_edge_length:
                    vicinity += 1
        return vicinity


def parse_arguments():
    arguments_parser = argparse.ArgumentParser('Generate min dist graph.')
    arguments_parser.add_argument('-i', dest='in_dir', type=str, required=True,
                                  help='Directory with outbreak fasta files')
    arguments_parser.add_argument('-o', dest='out_file', type=str, required=True,
                                  help='Name of output file')
    arguments_parser.add_argument('-v', dest='vicinity', type=int, default=1,
                                  help='border_vicinity')
    arguments_parser.add_argument('-l', dest='log_scale', type=int, default=3,
                                  help='log scale')
    return arguments_parser.parse_args()


def main():
    args = parse_arguments()
    graph = MinDistanceGraph(args.in_dir, args.vicinity, args.log_scale)
    nx.write_edgelist(graph.graph, args.out_file)


if __name__ == '__main__':
    main()
