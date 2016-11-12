#!/usr/bin/env python3

import random
import json

MAX_POPULATION = 1000

import networkx as nx
from networkx.readwrite import json_graph

FILE_NAME = "out/graphs/BJ30_unique_1a_6.json"


def import_graph(file_name):
    with open(file_name) as f:
        data = json.load(f)
        return json_graph.adjacency_graph(data)


def random_walk(network, begin_node):
    current_node = begin_node
    counter = 0
    visited_nodes = [False] * network.number_of_nodes()
    while True:
        counter += 1
        visited_nodes[current_node] = True
        if are_all_nodes_visited(visited_nodes):
            break
        d = random.random()
        p = 0
        for (u, v, w) in network.edges(current_node, data='weight'):
            if u == current_node:
                p += w
                if p > d:
                    current_node = v
                    break
    return counter


def propagate(network, begin_node, last_original_node):
    counter = 0
    visited_nodes = [0] * network.number_of_nodes()
    visited_nodes[begin_node] = 1
    while True:
        counter += 1
        if are_all_nodes_visited(visited_nodes, last_original_node):
            break
        for current_node, visited in enumerate(visited_nodes):
            for _ in range(visited):
                d = random.random()
                p = 0
                for (u, v, w) in network.edges(current_node, data='weight'):
                    if u == current_node:
                        p += w
                        if p > d:
                            if visited_nodes[v] < MAX_POPULATION:
                                visited_nodes[v] += 1
                            break
    return counter


def are_all_nodes_visited(nodes, last):
    print(len(list(filter(lambda i: i != 0, nodes[:last]))))
    return min(nodes[:last])


def main():
    network = import_graph(FILE_NAME)
#    print(random_walk(network, 10))
    print(propagate(network, 1, 5))


if __name__ == "__main__":
    main()
