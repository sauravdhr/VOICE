#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/4/17
"""
import argparse
import networkx as nx
import os
import network_creator


INF = 10e6


class MinDistanceGraph(object):
    def __init__(self, outbreak_folder):
        self.outbreak_folder = outbreak_folder
        self.edges = self.get_edges(outbreak_folder)
        self.graph = nx.Graph()
        self.graph.add_weighted_edges_from(self.edges)

    @staticmethod
    def get_edges(outbreak_folder):
        edges = list()
        fasta_files_names = os.listdir(outbreak_folder)
        for i in range(len(fasta_files_names)-1):
            for j in range(i+1, len(fasta_files_names)):
                weight = MinDistanceGraph.get_min_dist(os.path.join(outbreak_folder, fasta_files_names[i]),
                                                       os.path.join(outbreak_folder,fasta_files_names[j]))
                edges.append((fasta_files_names[i], fasta_files_names[j], weight))
            print("{0} is done out of {1}".format(i+1, len(fasta_files_names)-1))
        return edges

    @staticmethod
    def get_min_dist(fasta1_file_name, fasta2_file_name):
        seqs1 = network_creator.parse_fasta(fasta1_file_name)
        seqs2 = network_creator.parse_fasta(fasta2_file_name)
        min_dist = INF
        for i in range(len(seqs1)):
            for j in range(len(seqs2)):
                dist = network_creator.hamming_distance(seqs1[i], seqs2[j])
                if dist < min_dist:
                    min_dist = dist
        return min_dist



def parse_arguments():
    arguments_parser = argparse.ArgumentParser('Generate min dist graph.')
    arguments_parser.add_argument('-i', dest='in_dir', type=str, required=True,
                                  help='Directory with outbreak fasta files')
    arguments_parser.add_argument('-o', dest='out_file', type=str, required=True,
                                  help='Name of output file')
    return arguments_parser.parse_args()


def main():
    args = parse_arguments()
    graph = MinDistanceGraph(args.in_dir)
    nx.write_edgelist(graph.graph, args.out_file)


if __name__ == '__main__':
    main()
