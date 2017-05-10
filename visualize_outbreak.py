#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 5/5/17
"""

import argparse
from Bio import SeqIO
import glob
import networkx as nx
import numpy as np
from itertools import chain
from networkx.drawing.nx_pydot import write_dot


#COLORS = ['grey', 'purple', 'blue', 'yellow', 'cyan', 'magenta',
#          'orange', 'green', 'red',
#          'tomato', 'lime', 'tan', 'white', 'turquoise',
#          'darkorange', 'lightskyblue', 'darkseagreen', 'hotpink', 'black']

COLORS = ['grey', 'purple', 'blue', 'yellow', 'cyan', 'magenta',
          'green', 'orange', 'turquoise',
          'tomato', 'tan', 'white', 'yellowgreen', 'seagreen',
          'darkorange', 'red', 'lightskyblue', 'darkseagreen', 'hotpink', 'black']


def parse_args():
    parser = argparse.ArgumentParser(description="Draws outbreak network")
    parser.add_argument('outbreak_dir')
    return parser.parse_args()


def hamming(a, b):
    if len(a) != len(b):
        return -1
    l = len(a)
    return sum(map(lambda i: a[i] != b[i], range(l)))


def get_dist_matr(seqs):
    dist = np.zeros((len(seqs), len(seqs)))
    for i in range(len(seqs)-1):
        for j in range(i+1, len(seqs)):
            dist[i, j] = hamming(seqs[i], seqs[j])
            dist[j, i] = dist[i, j]
    return dist


def visualize_outbreak(outbreak_dir):
    g = nx.Graph()
    files = glob.glob("%s/*.fas" % outbreak_dir)
    fastas = [list(SeqIO.parse(f, 'fasta')) for f in files]
    populations_names = [f.split('/')[-1].split('_')[0] for f in files]
    total_population = list(map(lambda x: len(x), fastas))
    seqs = sum(fastas, [])
    inner_neigh = []
    outer_neigh = []
    dist = get_dist_matr(seqs)
    low_bound = 0
    high_bound = -1
    population_ind = -1
    for i in range(len(seqs)):
        if i > high_bound:
            population_ind += 1
            low_bound = high_bound + 1
            high_bound = high_bound + total_population[population_ind]
#            print("Population: {0}".format(populations_names[population_ind]))
        if i == low_bound:
            inner = list(range(low_bound+1, high_bound+1))
        elif i == high_bound:
            inner = list(range(low_bound, high_bound))
        else:
            inner = list(chain(range(low_bound, i), range(i+1, high_bound+1)))
        if i == 0:
            outer = list(range(high_bound+1, len(dist[i])))
        elif high_bound == len(seqs) - 1:
            outer = list(range(low_bound))
        else:
            outer = list(chain(range(low_bound), range(high_bound+1, len(dist[i]))))
        inner_min_dist = min(dist[i][inner])
        outer_min_dist = min(dist[i][outer])
        inner_ind = list(filter(lambda x: dist[i][x] == inner_min_dist, inner))
        outer_ind = list(filter(lambda x: dist[i][x] == outer_min_dist, outer))
#        print("Nearest {2}: inner - {0}, outer - {1}".format(inner_min_dist, outer_min_dist, i))
#        print("Inner - {0}, outer - {1}".format(inner_ind, outer_ind))
        g.add_node(i, color=COLORS[population_ind], label=populations_names[population_ind])
        for j in sum([inner_ind, outer_ind], []):
            g.add_edge(i, j)
    write_dot(g, 'test.dot')
#        inner_neigh.append(inner_ind)
#        outer_neigh.append(outer_ind)



def main():
    args = parse_args()
    visualize_outbreak(args.outbreak_dir)


if __name__ == "__main__":
    main()
