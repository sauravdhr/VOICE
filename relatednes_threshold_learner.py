#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12.01.2016
"""

import sys
import graph_utils
import os


def get_outbreak_breakage_preventing_threshold(full_network_edges, outbreaks_networks_edges):
    #TODO:
    return 0


def get_outbreaks_joining_preventing_threshold(outbreaks_networks_edges):
    #TODO:
    return 0


def main(full_network_filename, outbreaks_networks_dir):
    full_network_edges = graph_utils.import_graph(full_network_filename).edges()
    outbreaks_networks = [graph_utils.import_graph(f) for f in os.listdir(outbreaks_networks_dir)]
    outbreaks_breakage_preventing_threshold = get_outbreak_breakage_preventing_threshold(
        full_network_edges, [g.edges() for g in outbreaks_networks])
    outbreaks_joining_preventing_threshold = get_outbreaks_joining_preventing_threshold(outbreaks_networks)


if __name__ == '__main__':
    full_network_filename = sys.argv[1]
    outbreaks_networks_dir = sys.argv[2]
    main(full_network_filename, outbreaks_networks_dir)
