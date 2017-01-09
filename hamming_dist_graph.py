#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/8/17
"""
import networkx as nx
import graph_utils
import operator


def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def find_consensus(seqs):
    consensus = list()
    for i in range(len(seqs[0])):
        counter = dict(zip('ACGT', [0] * 4))
        for j in range(len(seqs)):
            counter[seqs[j][i]] += 1
        consensus.append(max(counter.items(), key=operator.itemgetter(1))[0])
    return "".join(consensus)


def find_closest_sequences(seqs1, seqs2, vicinity=0):
    min_dist = float('inf')
    seqs_pairs = []
    for s1 in seqs1:
        for s2 in seqs2:
            new_dist = hamming_distance(s1, s2)
            if new_dist < min_dist:
                min_dist = new_dist
                seqs_pairs = [(s1, s2)]
            elif new_dist == min_dist:
                seqs_pairs.append((s1, s2))
    if vicinity:
        seqs_pairs = []
        for s1 in seqs1:
            for s2 in seqs2:
                if min_dist <= hamming_distance(s1, s2) + vicinity:
                    seqs_pairs.append((s1, s2))
    return seqs_pairs


class HammingDistGraph(graph_utils.Graph):
    def __init__(self, sequences):
        super(HammingDistGraph, self).__init__(nx.Graph)
        self.graph.add_nodes_from(sequences)
        self.build_edges()
        self.shortest_path_trees = self.find_shortest_path_trees()

    def build_edges(self):
        nodes = self.graph.nodes()
        for i in range(len(nodes)-1):
            for j in range(i+1, len(nodes)):
                self.graph.add_edge(nodes[i], nodes[j])
                self.graph[nodes[i]][nodes[j]]['weight'] = hamming_distance(nodes[i], nodes[j])
