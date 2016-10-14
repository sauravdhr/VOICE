#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

import itertools
import math
import copy
from graphviz import Digraph
import networkx as nx
from networkx.readwrite import json_graph
import json

OUT_DIR = "out/graphs"
#FASTA = "data/AC122_unique_1a_48.fas"
FASTA = "data/AA45_unique_1b_161.fas"
#FASTA = "data/examples.fas"
MAX_DIST = 4


class SymMatrixWithoutDiagonal(object):
# it is upper-packed matrix (see http://www.ibm.com/support/knowledgecenter/SSFHY8_5.3.0/com.ibm.cluster.essl.v5r3.essl100.doc/am5gr_upsm.htm#am5gr_upsm).
    def __init__(self, vector):
        self.vector = vector
        self.n = int((1+math.sqrt(1+8*len(self.vector)))/2)

    def __setitem__(self, index, value):
        i = self.get_index(index[0], index[1])
        self.vector[i] = value

    def __getitem__(self, index):
        if len(index) == 1:
            return self.get_row(index)
        i = self.get_index(index[0], index[1])
        if i is not None:
            return self.vector[i]
        return 0

    @staticmethod
    def get_index(i, j):
        if i == j:
            return None
        if i < j:
            return SymMatrixWithoutDiagonal.convert_index(i, j-1)
        else:
            return SymMatrixWithoutDiagonal.convert_index(j, i-1)

    @staticmethod
    def convert_index(i, j):
        return int(i+j*(j+1)/2)

    def __iter__(self):
        return MatrixIterator(self)

    def get_row(self, i):
        return [self[i, j] for j in range(0, self.n)]

    def __len__(self):
        return self.n


class MatrixIterator(object):
    def __init__(self, matrix):
        self.i = 0
        self.matrix = matrix

    def __iter__(self):
        return self

    def __next__(self):
        if self.i < self.matrix.n:
            self.i += 1
            return self.matrix.get_row(self.i-1)
        else:
            raise StopIteration()


def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def filter_repeated_sequences(sequences):
    if not sequences:
        return sequences
    return list(sequences[i] for i in
                [0] + list(
                    filter(
                        lambda i: not(True in map(lambda s2: sequences[i] == s2, sequences[:i])),
                        list(range(1, len(sequences))))))


def get_sequences_distance_matrix(sequences):
    return SymMatrixWithoutDiagonal(list(itertools.chain.from_iterable(map(
        lambda i: list(map(
            lambda s2: hamming_distance(sequences[i], s2),
            sequences[:i])),
        list(range(1, len(sequences)))))))


def get_median(vector1, vector2, vector3):
    median = ""
    for ch1, ch2, ch3 in zip(vector1, vector2, vector3):
        if ch1 == ch2 or ch1 == ch3:
            median += ch1
        elif ch2 == ch3:
            median += ch2
        else:
            return None
    return median


def find_all_medians(sequences):
        medians = []
        n = len(sequences)
        for i in range(n-2):
            for j in range(i+1, n-1):
                for k in range(j+1, n):
                    median = get_median(sequences[i], sequences[j], sequences[k])
                    if median:
                        medians.append(median)
        return filter_repeated_sequences(medians)


def find_medians_for_similar_sequences(sequences, max_dist):
    distance_matrix = get_sequences_distance_matrix(sequences)

    similar_sequences = []
    for i in range(len(sequences)):
        similar_sequences_acc = []
        for j in range(i, len(sequences)):
            if distance_matrix[i, j] <= max_dist:
                similar_sequences_acc.append(sequences[j])
        if len(similar_sequences_acc) <= 2:
            continue
        similar_sequences.append(similar_sequences_acc)
    medians = []
    for seqs in similar_sequences:
        medians += find_all_medians(seqs)
    return filter_repeated_sequences(medians)


class SequencesNetworkCreator(object):
    PROBABILITY_TO_STAY = 0.38

    def __init__(self, fasta1, fasta2, max_dist_for_median_triple):
        self.sequences1 = filter_repeated_sequences(self.parse_fasta(fasta1))
        self.sequences2 = filter_repeated_sequences(self.parse_fasta(fasta2))
        self.vertices = self.infer_vertices(self.sequences1, self.sequences2, max_dist_for_median_triple)
        self.sequences_distance_matrix = get_sequences_distance_matrix(self.vertices)
        self.minimal_connected_graph_matrix = self.get_minimal_connected_graph_matrix(self.sequences_distance_matrix)
        self.network_matrix = self.get_reciprocal_distance_matrix(self.minimal_connected_graph_matrix)
        self.propagation_network = self.construct_propagation_network(self.network_matrix)

    @staticmethod
    def get_reciprocal_distance_matrix(distance_matrix):
        reciprocal_distance_matrix = copy.deepcopy(distance_matrix)
        for i in range(len(reciprocal_distance_matrix.vector)):
            if reciprocal_distance_matrix.vector[i]:
                reciprocal_distance_matrix.vector[i] = 1/reciprocal_distance_matrix.vector[i]
        return reciprocal_distance_matrix

    @staticmethod
    def parse_fasta(fasta):
        seqs = []
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seqs.append(str(seq_record.seq))
        return seqs

    @staticmethod
    def infer_vertices(sequences1, sequences2, max_dist):
        common_sequences = filter_repeated_sequences(sequences1 + sequences2)
        vertices = common_sequences[:]
        old_n = 0
        new_n = len(vertices)
        while True:
            vertices += find_medians_for_similar_sequences(vertices[old_n:new_n], max_dist)
            vertices = filter_repeated_sequences(vertices)
            old_n = new_n
            new_n = len(vertices)
            if old_n == new_n:
                break
            print(new_n - old_n)
        return vertices

    @staticmethod
    def get_minimal_connected_graph_matrix(distance_matrix):
        mst = minimum_spanning_tree(csr_matrix(distance_matrix))
        min_length = max([max(e) for e in mst.toarray().astype(int)])
        print(min_length)
        minimal_connected_graph_matrix = copy.deepcopy(distance_matrix)
#        min_edges = []
#        for row in minimal_connected_graph_matrix:
#            min_edges.append(min([e for e in row if e != 0]))
#        min_length = max(min_edges)
        for i, row in enumerate(minimal_connected_graph_matrix):
            for j, length in enumerate(row):
                if length and length > min_length:
                    minimal_connected_graph_matrix[i, j] = None
        return minimal_connected_graph_matrix

    def construct_propagation_network(self, distance_matrix):
        network = list()
        for vertex_ind in range(len(self.vertices)):
            adj_vertices_count, sum_weight_of_edges = self.get_info_about_adj_edges(distance_matrix, vertex_ind)
            edges = list()
            edges.append((vertex_ind, self.PROBABILITY_TO_STAY))
            if adj_vertices_count != 0:
                prob_move_to_distant_vertex = (1 - self.PROBABILITY_TO_STAY)/sum_weight_of_edges
                for adj_vertex_ind in range(len(self.vertices)):
                    d = distance_matrix[vertex_ind, adj_vertex_ind]
                    if d is not None and d != 0:
                        edges.append((adj_vertex_ind, prob_move_to_distant_vertex * d))
            network.append(edges)
        return network

    def get_info_about_adj_edges(self, distance_matrix, vertex_ind):
        adjacent_vertices_count = 0
        sum_weight_of_edges = 0
        for j in range(len(self.vertices)):
            d = distance_matrix[vertex_ind, j]
            if d is not None:
                adjacent_vertices_count += 1
                sum_weight_of_edges += d
        return adjacent_vertices_count, sum_weight_of_edges


def export_graph_to_dot(graph, file_name):
#    for item in graph:
#        print(', '.join(map(str, item[:])))
    dot = Digraph()
    for v in range(len(graph)):
        dot.node(str(v))
    for v1 in range(len(graph)):
        for (v2, weight) in graph[v1]:
#            dot.edge(str(v1), str(v2))
            dot.edge(str(v1), str(v2), weight=str(weight))
    with open(file_name, 'w') as f:
        f.write(dot.source)


def export_graph_to_json(graph, file_name):
    g = nx.DiGraph()
    for v in range(len(graph)):
        g.add_node(v)
    for v1 in range(len(graph)):
        for (v2, weight) in graph[v1]:
            g.add_edge(v1, v2, weight=weight)
    data = json_graph.adjacency_data(g)
    with open(file_name, 'w') as f:
        json.dump(data, f)


def main(fasta1, fasta2, max_dist_for_median_triple):
    graph = SequencesNetworkCreator(fasta1, fasta2, max_dist_for_median_triple)
    f = os.path.join(OUT_DIR, os.path.splitext(os.path.basename(fasta1))[0] + '_vs_'
                     + os.path.splitext(os.path.basename(fasta1))[0] + '_' + str(max_dist_for_median_triple))
    out_file_json = f + '.json'
    out_file_dot = f + '.dot'
    export_graph_to_dot(graph.propagation_network, out_file_dot)
    export_graph_to_json(graph.propagation_network, out_file_json)

if __name__ == "__main__":
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    max_dist_for_triple_mean = int(sys.argv[3])
    main(fasta1, fasta2, max_dist_for_triple_mean)
