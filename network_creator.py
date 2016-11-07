#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from enum import Enum

import itertools
import math
import copy
from graphviz import Digraph
import networkx as nx
from networkx.readwrite import json_graph
import json

OUT_DIR = "out/graphs"
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


def get_sequence_sets_difference(sequences1, sequences2):
    return list(sequences1[i] for i in list(
        filter(lambda i: not(True in map(lambda s2: sequences1[i] == s2, sequences2)), list(range(len(sequences1))))))


def infer_distance_matrix(sequences):
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
    distance_matrix = infer_distance_matrix(sequences)

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


class Graph(object):
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.edges = edges


class VertexType(Enum):
    original = 1
    median = 2
    compressed_hypercube = 3


class DistanceGraphBuilder(object):
    MIN_SEQS_DIST_THRESHOLD = 6
    MAX_DIST_FOR_MEDIANS = 8

    def __init__(self, fasta, max_dist_for_medians):
        self.max_dist_for_medians = max_dist_for_medians
        self.sequences = filter_repeated_sequences(self.parse_fasta(fasta))
        self.medians = self.infer_medians(self.sequences, self.max_dist_for_medians)
        self.vertices = self.infer_vertices()
        self.distance_matrix = infer_distance_matrix(self.sequences + self.medians)
        self.minimal_connected_graph_matrix = self.construct_minimal_connected_graph_matrix(
            self.distance_matrix)

    def infer_vertices(self):
        return [{'sequence': s, 'type': VertexType.original} for s in self.sequences]\
               + [{'sequence': m, 'type': VertexType.median} for m in self.medians]

    def get_minimal_connected_graph(self):
        edges = list()
        for vertex_ind in range(len(self.sequences + self.medians)):
            adjacent_vertices = list()
            for adj_vertex_ind in range(len(self.sequences + self.medians)):
                d = self.minimal_connected_graph_matrix[vertex_ind, adj_vertex_ind]
                if d is not None and d != 0:
                    adjacent_vertices.append((adj_vertex_ind, {"weight": d}))
            edges.append(adjacent_vertices)
        return Graph(self.vertices, edges)

    @staticmethod
    def parse_fasta(fasta):
        seqs = []
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seqs.append(str(seq_record.seq))
        return seqs

    @staticmethod
    def infer_medians(sequences, max_dist_for_medians):
        filtered_sequences = filter_repeated_sequences(sequences)
        medians = []
        old_n = 0
        new_n = 0
        while True:
            medians += find_medians_for_similar_sequences(filtered_sequences + medians[old_n:new_n], max_dist_for_medians)
            medians = get_sequence_sets_difference(filter_repeated_sequences(medians), filtered_sequences)
            old_n = new_n
            new_n = len(medians)
            if old_n == new_n:
                break
            print(new_n - old_n)
        return medians

    @staticmethod
    def construct_minimal_connected_graph_matrix(distance_matrix):
        mst = minimum_spanning_tree(csr_matrix(distance_matrix))
        min_length = max([max(e) for e in mst.toarray().astype(int)] + [DistanceGraphBuilder.MIN_SEQS_DIST_THRESHOLD])
        minimal_connected_graph_matrix = copy.deepcopy(distance_matrix)
        for i, row in enumerate(minimal_connected_graph_matrix):
            for j, length in enumerate(row):
                if length and length > min_length:
                    minimal_connected_graph_matrix[i, j] = None
        return minimal_connected_graph_matrix


class ProbabilityGraphBuilder(object):
    e = 0.01
    s = e/(1-3*e)
    min_edge_probability = 10e-6

    @staticmethod
    def loop_probability(L):
        return (1 - 3 * ProbabilityGraphBuilder.e) ** L

    def edge_probability(self, m):
        return self.c*self.s ** m

    def __init__(self, distance_graph):
        self.distance_graph = distance_graph
        self.L = len(distance_graph.vertices[0]['sequence'])
        self.c = self.loop_probability(self.L)
        self.distance_graph_with_compressed_hypercubes = \
            self.infer_distance_graph_with_compressed_hypercubes(self.distance_graph)
        self.probability_graph = self.infer_probability_graph(self.distance_graph_with_compressed_hypercubes)

    @staticmethod
    def infer_distance_graph_with_compressed_hypercubes(distance_graph):
        distance_graph_with_compressed_hypercubes = copy.deepcopy(distance_graph)
        for u, adjacency_list in enumerate(distance_graph.edges):
            for v, properties in adjacency_list:
                if u < v and properties["weight"] > 1:
                    ProbabilityGraphBuilder.add_edges_of_compressed_hypercube(distance_graph_with_compressed_hypercubes,
                                                                              u, v, properties["weight"])
        return distance_graph_with_compressed_hypercubes

    @staticmethod
    def add_edges_of_compressed_hypercube(graph, u, v, weight):
        count_of_vertices = len(graph.vertices)
        for d in range(1, weight):
            graph.vertices.append({'type': VertexType.compressed_hypercube})
            graph.edges.append([(u, {'weight': d, 'multiplicity': weight}), (v, {'weight': weight-d, 'multiplicity': weight})])
            graph.edges[u].append((count_of_vertices+d-1, {'weight': d, 'multiplicity': weight}))
            graph.edges[v].append((count_of_vertices+d-1, {'weight': weight-d, 'multiplicity': weight}))
        for v in range(count_of_vertices, count_of_vertices + weight - 2):
            for u in range(v, count_of_vertices + weight - 1):
                if v != u:
                    graph.edges[v].append((u, {'weight': u-v, 'multiplicity': weight}))
                    graph.edges[u].append((v, {'weight': u-v, 'multiplicity': weight}))

    def infer_probability_graph(self, distance_graph):
        probability_graph = copy.deepcopy(distance_graph)
        for u, adjacency_list in enumerate(distance_graph.edges):
            probability_graph.edges[u] = []
            for i, (v, properties) in enumerate(adjacency_list):
                multiplicity = properties["multiplicity"] if 'multiplicity' in properties else 1
                edge_probability = self.edge_probability(properties["weight"]) * multiplicity
                if edge_probability >= self.min_edge_probability:
                    probability_graph.edges[u].append((v, {"weight": edge_probability}))
            probability_graph.edges[u].append((u, {"weight": self.c}))
        return probability_graph


class GraphExporter(object):
    @staticmethod
    def export(graph, file_name):
        pass

    @staticmethod
    def get_vertex_color(vertex_type):
        if vertex_type == VertexType.original:
            return 'red'
        if vertex_type == VertexType.median:
            return 'yellow'
        if vertex_type == VertexType.compressed_hypercube:
            return 'green'
        return None


class DotExporter(GraphExporter):
    @staticmethod
    def export(graph, file_name):
        dot = Digraph()
        for v in range(len(graph.vertices)):
            label = graph.vertices[v]['sequence'] if 'sequence' in graph.vertices[v] else 'pool'
            dot.node(str(v), label=label,
                     color=super(JsonExporter, JsonExporter).get_vertex_color(graph.vertices[v]['type']))
        for v1 in range(len(graph.edges)):
            for (v2, properties) in graph.edges[v1]:
                dot.edge(str(v1), str(v2), weight=str(properties['weight']))
        with open(file_name, 'w') as f:
            f.write(dot.source)


class JsonExporter(GraphExporter):
    @staticmethod
    def export(graph, file_name):
        g = nx.DiGraph()
        for v in range(len(graph.vertices)):
            g.add_node(v, label=graph.vertices[v]['sequence'],
                       color=super(JsonExporter, JsonExporter).get_vertex_color(graph.vertices[v]['type']))
        for v1 in range(len(graph.edges)):
            for (v2, properties) in graph.edges[v1]:
                g.add_edge(v1, v2, weight=properties['weight'])
        data = json_graph.adjacency_data(g)
        with open(file_name, 'w') as f:
            json.dump(data, f)


def main(fastas, max_dist_for_median_triple):
    graphs = [ProbabilityGraphBuilder(DistanceGraphBuilder(
        fasta, max_dist_for_median_triple).get_minimal_connected_graph()) for fasta in fastas]
    for i, graph in enumerate(graphs):
        f = os.path.join(OUT_DIR, os.path.splitext(os.path.basename(fastas[i]))[0]
                         + '_' + str(max_dist_for_median_triple))
        out_file_json = f + '.json'
        out_file_dot = f + '.dot'
        DotExporter.export(graph.probability_graph, out_file_dot)
#        DotExporter.export(graph.distance_graph_with_compressed_hypercubes, out_file_dot)
#        JsonExporter.export(graph.probability_graph, out_file_json)

if __name__ == "__main__":
    fastas = [sys.argv[1], sys.argv[2]]
    max_dist_for_triple_mean = int(sys.argv[3])
    main(fastas, max_dist_for_triple_mean)
