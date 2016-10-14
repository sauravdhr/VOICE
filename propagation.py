#!/usr/bin/env python3

import random
import json

import networkx as nx
from networkx.readwrite import json_graph

FILE_NAME = "out/graphs/AA45_unique_1b_161_4.json"


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


def propagate(network, begin_node):
    counter = 0
    visited_nodes = [False] * network.number_of_nodes()
    visited_nodes[begin_node] = True
    while True:
        counter += 1
        if are_all_nodes_visited(visited_nodes):
            break
        for current_node, visited in enumerate(visited_nodes):
            if visited:
                d = random.random()
                p = 0
                for (u, v, w) in network.edges(current_node, data='weight'):
                    if u == current_node:
                        p += w
                        if p > d:
                            visited_nodes[v] = True
                            break
    return counter


def are_all_nodes_visited(nodes):
    return min(nodes)


def main():
    network = import_graph(FILE_NAME)
    print(random_walk(network, 10))
    print(propagate(network, 10))


if __name__ == "__main__":
    main()
