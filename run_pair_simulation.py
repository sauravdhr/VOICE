#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 1/9/2017
"""

import os
import argparse
import networkx as nx
from networkx.drawing.nx_agraph import write_dot

import network_creator
import simulation_graph
import probability_graph
import simulator

OUT_DIR = "results"
DEFAULT_SIMULATIONS_NUMBER = 5
DEFAULT_N = 100


def main(fastas, out_dir, simulations_count, L):
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in fastas]

    if L == 0:
        L = network_creator.get_count_of_heterogeneous_positions(sequences_sets[0] + sequences_sets[1])

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = os.path.join(out_dir, fastas_basenames[0] + '_vs_' + fastas_basenames[1])
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    pair_simulation_graph = simulation_graph.SimulationDistGraph(sequences_sets[0], sequences_sets[1])
    pair_probability_graphs = [probability_graph.ProbabilityGraph(s, L)
                               for s in pair_simulation_graph.simulation_graphs]

    results_log = os.path.join(out_dir, 'results.log')
    for i in [0, 1]:
        f = os.path.join(out_dir, fastas_basenames[i])
        nx.write_gexf(pair_simulation_graph.simulation_graphs[i], f + '.gexf')

    with open(results_log, 'w') as log:
        log.write('host simulation_number cycles\n')
        for i in [0, 1]:
            for k in range(simulations_count):
                log.write(fastas_basenames[i] + ' '
                          + str(k) + ' '
                          + str(simulator.Simulator(pair_probability_graphs[i].graph).
                                simulate()) + '\n')


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta1",
                        help="Fasta file for 1st host.")
    parser.add_argument("fasta2",
                        help="Fasta file for 2nd host.")
    parser.add_argument("-o", dest='out_dir', type=str, default=OUT_DIR,
                        help="Path to an output directory relative to input dir. "
                             "By default it creates directory named \'" + OUT_DIR + "\' in the input directory")
    parser.add_argument("-n", dest='simulations_count', type=int, default=DEFAULT_SIMULATIONS_NUMBER,
                        help="Count of simulations repeats.")
    parser.add_argument("-L", dest='L', type=int, default=DEFAULT_N,
                        help="Parameter which represent a number of nucleotides prone to mutate during simulation" +
                        " (it is used in probability formula for edge weight calculation).")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main([args.fasta1, args.fasta2], args.out_dir, args.simulations_count, args.L)
