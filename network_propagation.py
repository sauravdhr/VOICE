#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12.01.2016
"""

import random
import os
import graph_utils

FILE_NAME = "out/graphs/AD002_unique_1a_72.json"
LOG_FILE_NAME = "out/graphs/AD002_unique_1a_72.out"
#LOGGING_PERIOD = 100


class Propagator(object):
    MAX_ORIGINAL_POPULATION = 10
    MAX_HIDDEN_POPULATION = 1

    def __init__(self, network, log_file_name):
        self.files = []
        self.network = network
        self.log_file = open(log_file_name, 'w')
        self.files.append(self.log_file)
        self.observed_seqs_count = self.get_observed_seqs_count()
        self.counter = 0
        self.population = [0] * self.network.number_of_nodes()
        self.influential_nodes = dict()
        random.seed()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        for file in self.files:
            os.unlink(file)

    def get_observed_seqs_count(self):
        observed_seqs_count = 0
        for i in self.network.nodes_iter():
            if self.network.node[i]['type'] != 'original':
                return observed_seqs_count
            observed_seqs_count += 1
        return observed_seqs_count

    def log_propagation(self, node, counter):
        self.log_file.write("{0}, {1}\n".format(node, counter))

    def propagate(self, begin_nodes):
        for begin_node in begin_nodes:
            self.log_propagation(begin_node, self.counter)
            self.population[begin_node] = 1
            self.influential_nodes[begin_node] = None
        while True:
            self.counter += 1
            visited_nodes = []
            for i, influential_node in enumerate(self.influential_nodes):
                for _ in range(self.population[influential_node]):
                    d = random.random()
                    p = 0
                    for (u, v, w) in self.network.edges(influential_node, data='weight'):
                        p += w
                        if p > d:
                            visited_nodes.append(v)
                            break

            new_influential_nodes = []
            for n in visited_nodes:
                max_population = self.MAX_HIDDEN_POPULATION if n > self.observed_seqs_count \
                    else self.MAX_ORIGINAL_POPULATION
                if self.population[n] == 0:
                    new_influential_nodes.append(n)
                if self.population[n] < max_population:
                    self.population[n] += 1

            for n in new_influential_nodes:
                if n < self.observed_seqs_count:
                    self.log_propagation(n, self.counter)
                self.influential_nodes[n] = None
            self.update_influential_nodes()

            if self.are_all_original_nodes_visited():
#                self.log_file.write(self.get_population_status_string())
                break
        # Prevent infinite simulations
        if self.counter > 3000:
            print("Exceeded the maximum number of tacts. Stopping the simulation.")
        return self.counter

    def update_influential_nodes(self):
        nodes_with_lost_influence = []
        for i, influential_node in enumerate(self.influential_nodes):
            if self.are_all_neighbours_visited(influential_node):
                nodes_with_lost_influence.append(influential_node)
        for node in nodes_with_lost_influence:
            self.influential_nodes.pop(node, None)

    def are_all_neighbours_visited(self, node):
        for (_, v) in self.network.edges(node):
            if self.population[v] == 0:
                return False
        return True

    def are_all_original_nodes_visited(self):
#        if not self.counter % LOGGING_PERIOD:
#            print(len(list(filter(lambda i: i != 0, self.population[:self.observed_seqs_count]))))
#            self.log_file.write(self.get_population_status_string())
        return min(self.population[:self.observed_seqs_count])

    def get_population_status_string(self):
        return ' '.join(str(e) for e in self.population[:self.observed_seqs_count]) + '\n'


def main():
    network = graph_utils.import_graph_from_json(FILE_NAME)
    print(Propagator(network, LOG_FILE_NAME).propagate([10]))


if __name__ == "__main__":
    main()
