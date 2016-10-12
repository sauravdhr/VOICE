#!/usr/bin/env python3
from Bio import SeqIO
import itertools
import math
import copy
from graphviz import Digraph

#FASTA = "data/AC122_unique_1a_48.fas"
FASTA = "data/AA45_unique_1b_161.fas"
MAX_DIST = 5
MAX_EDGE_LENGTH = 6


class SymMatrixWithoutDiagonal(object):
# it is upper-packed matrix (see http://www.ibm.com/support/knowledgecenter/SSFHY8_5.3.0/com.ibm.cluster.essl.v5r3.essl100.doc/am5gr_upsm.htm#am5gr_upsm).
    def __init__(self, vector):
        self.vector = vector
        self.n = int((1+math.sqrt(1+8*len(self.vector)))/2)

    def __getitem__(self, index):
        if len(index) != 2:
            return self.vector[index]
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
        return i+int(j*(j+1)/2)

    def __iter__(self):
        return MatrixIterator(self)

    def get_row(self, i):
        return [self[i, j] for j in range(0, self.n)]


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
#    for item in distance_matrix:
#            print(', '.join(map(str, item[:])))

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

    def __init__(self, fasta):
        self.sequences = self.parse_fasta(fasta)
        self.vertices = self.infer_vertices(self.sequences, MAX_DIST)
        self.sequences_distance_matrix = get_sequences_distance_matrix(self.vertices)
        self.reciprocal_distance_matrix = self.get_reciprocal_distance_matrix()
        self.propagation_network = self.construct_propagation_network(self.reciprocal_distance_matrix, MAX_EDGE_LENGTH)

    def get_reciprocal_distance_matrix(self):
        reciprocal_distance_matrix = copy.deepcopy(self.sequences_distance_matrix)
        for i in range(len(reciprocal_distance_matrix.vector)):
            reciprocal_distance_matrix.vector[i] = 1/reciprocal_distance_matrix.vector[i]
        return reciprocal_distance_matrix

    @staticmethod
    def parse_fasta(fasta):
        seqs = []
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seqs.append(str(seq_record.seq))
        return seqs

    @staticmethod
    def infer_vertices(sequences, max_dist):
        medians = sequences[:]
        old_n = 0
        new_n = len(medians)
        while True:
            medians += find_medians_for_similar_sequences(medians[old_n:new_n], max_dist)
            medians = filter_repeated_sequences(medians)
            old_n = new_n
            new_n = len(medians)
            if old_n == new_n:
                break
            print(new_n - old_n)
        return medians

    def construct_propagation_network(self, distance_matrix, max_dist):
        network = list()
        for vertex_ind in range(self.vertices):
            _, sum_weight_of_edges = self.get_info_about_adj_edges(distance_matrix, vertex_ind, max_dist)
            prob_move_to_distant_vertix = (1 - self.PROBABILITY_TO_STAY)/sum_weight_of_edges
            edges = list()
            edges.append((vertex_ind, self.PROBABILITY_TO_STAY))
            for adj_vertex_ind in range(self.vertices):
                d = distance_matrix[vertex_ind, adj_vertex_ind]
                if 1 <= d <= max_dist:
                    edges.append((adj_vertex_ind, prob_move_to_distant_vertix * d))
            network.append(edges)
        return network

    def get_info_about_adj_edges(self, distance_matrix, vertex_ind, max_dist):
        adjacent_vertices_count = 0
        sum_weight_of_edges = 0
        for j in range(self.vertices):
            d = distance_matrix[vertex_ind, j]
            if 1 <= d <= max_dist:
                adjacent_vertices_count += 1
                sum_weight_of_edges += d
        return adjacent_vertices_count, sum_weight_of_edges


def export_graph_to_dot(graph, file_name):
    dot = Digraph()
    for i in range(len(graph)):
        dot.node(str(i))
    for v1, edges in enumerate(graph):
        for v2 in edges:
            dot.edge(str(v1), str(v2[0]))
    with open(file_name, 'w') as f:
        f.write(dot.source)


def main(fasta):
    graph = SequencesNetworkCreator(fasta)
    export_graph_to_dot(graph.propagation_network, "test.dot")

if __name__ == "__main__":
    main(FASTA)
