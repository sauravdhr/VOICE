#!/usr/bin/env python3
"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 05.30.2017
"""

import analyze_outbreak_graph
import random
from functools import reduce
import argparse
import os
import networkx as nx


class GraphType(object):
    DIRECTED = 'directed'
    UNDIRECTED = 'undirected'


class KFoldGenerator(object):
    UNRELATED_PREFIX = 'XX'

    def __init__(self, graph_file_name, graph_type=GraphType.DIRECTED):
        self.graph_file_name = graph_file_name
        self.graph_type = graph_type
        self.outbreak_graph = self.get_outbreak_graph()
        self.graph_seqs = self.get_graph_seqs()
        self.outbreak_names = list(filter(lambda x: x != self.UNRELATED_PREFIX, self.graph_seqs.keys()))

        self.node_sets = []

    def get_graph_seqs(self):
        outbreak_seqs = dict()
        cur = ''
        seqs = list()
        for n in sorted(self.outbreak_graph.graph.nodes()):
            name = n[:2]
            if not cur:
                cur = name
            if cur != name:
                outbreak_seqs[cur] = seqs
                cur = name
                seqs = list()
            seqs.append(n)
        outbreak_seqs[cur] = seqs
        return outbreak_seqs

    def get_outbreak_graph(self):
        if self.graph_type == GraphType.DIRECTED:
            return analyze_outbreak_graph.DirectedGraphAnalyzer(self.graph_file_name)
        elif self.graph_type == GraphType.UNDIRECTED:
            return analyze_outbreak_graph.UndirectedGraphAnalyzer(self.graph_file_name)
        else:
            raise NameError("Wrong graph type")

    def get_k_folds(self, k):
        random.seed()
        k_sets_of_outbreaks = self.subdivide_outbreaks(k)
        k_sets_of_nodes = self.subdivide_outbreak_seqs(k_sets_of_outbreaks)
        self.add_unrelated_seqs(k_sets_of_nodes)
        return k_sets_of_nodes

    def add_unrelated_seqs(self, k_sets_of_nodes):
        min_seqs_per_fold = int(len(self.outbreak_graph.graph.nodes())/len(k_sets_of_nodes))
        rest_seqs = len(self.outbreak_graph.graph.nodes()) % len(k_sets_of_nodes)
        unrelated_seqs_shuffled = self.graph_seqs[self.UNRELATED_PREFIX][:]
        random.shuffle(unrelated_seqs_shuffled)
        cur_ind = 0
        for i in range(len(k_sets_of_nodes)):
            diff = min_seqs_per_fold - len(k_sets_of_nodes[i])
            if i < rest_seqs:
                diff += 1
            k_sets_of_nodes[i].extend(unrelated_seqs_shuffled[cur_ind:cur_ind+diff])
            cur_ind += diff

    def subdivide_outbreak_seqs(self, sets_of_outbreaks):
        k_sets_of_seqs = list()
        for fold in sets_of_outbreaks:
            fold_seqs = list()
            for o in fold:
                fold_seqs.extend(self.graph_seqs[o])
            k_sets_of_seqs.append(fold_seqs)
        return k_sets_of_seqs

    def subdivide_outbreaks(self, k):
        shuffled_outbreaks = self.outbreak_names[:]
        random.shuffle(shuffled_outbreaks)
        seqs_count = [0] * k
        k_folds_outbreaks = [list() for _ in range(k)]
        for o in shuffled_outbreaks:
            smallest_fold_ind = reduce(lambda a, b: a if seqs_count[a] < seqs_count[b] else b,
                                       range(k), 0)
            k_folds_outbreaks[smallest_fold_ind].append(o)
            seqs_count[smallest_fold_ind] += len(self.graph_seqs[o])
        return k_folds_outbreaks

    def write_data_sets(self, folds):
        out_dir = os.path.splitext(self.graph_file_name)[0] + "_{0}_fold".format(len(folds))
        create_dir(out_dir)
        for i in range(len(folds)):
            self.write_data_set(folds, i, out_dir)

    def write_data_set(self, folds, test_fold_ind, out_dir):
        train_nodes = reduce(lambda a, i: a + folds[i],
                             filter(lambda i: i != test_fold_ind, range(len(folds))),
                             list())
        self.write_subgraph(train_nodes, os.path.join(out_dir, "train_{0}.edgelist".format(test_fold_ind)))
        self.write_subgraph(folds[test_fold_ind], os.path.join(out_dir, "test_{0}.edgelist".format(test_fold_ind)))

    def write_subgraph(self, nodes, file_name):
        nx.write_edgelist(self.outbreak_graph.graph.subgraph(nodes), file_name)


def create_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_file', type=str, required=True,
                        help='File with outbreak graph')
    parser.add_argument('-t', dest='graph_type', type=str, default=GraphType.DIRECTED,
                        help='Graph type: \'directed\' or \'undirected\'')
    parser.add_argument('-k', dest='k', type=int, default=5,
                        help='Count of folds')
    return parser.parse_args()


def main():
    args = parse_arguments()
    k_fold_generator = KFoldGenerator(args.input_file, args.graph_type)
    k_folds = k_fold_generator.get_k_folds(args.k)
    k_fold_generator.write_data_sets(k_folds)

if __name__ == '__main__':
    main()
