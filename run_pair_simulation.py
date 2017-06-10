#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12.01.2016
"""

import os
import network_creator
import network_propagation
import graph_utils
import argparse
import montecarlo_network

OUT_DIR = "results"
DEFAULT_SIMULATIONS_NUMBER = 5
L = 100
K_MAX_DEFAULT = 0


def determine_network_sources(sequences1, sequences2):
    VICINITY = 1
    dist_s1 = [10e6] * len(sequences1)
    dist_s2 = [10e6] * len(sequences2)
    for i, s1 in enumerate(sequences1):
        for j, s2 in enumerate(sequences2):
            d = network_creator.hamming_distance(s1, s2)
            if d < dist_s1[i]:
                dist_s1[i] = d
            if d < dist_s2[j]:
                dist_s2[j] = d
    min_dist = min(dist_s1)
    return [list(map(lambda x: x[0], filter(lambda x: 0 < x[1] <= min_dist + VICINITY, enumerate(s))))
            for s in [dist_s1, dist_s2]]


def get_sequences_intersections_indices(sequences_sets_1, sequences_sets_2):
    indices_of_copies = [[], []]
    for i, s1 in enumerate(sequences_sets_1):
        for j, s2 in enumerate(sequences_sets_2):
            if s1 == s2:
                indices_of_copies[0].append(i)
                indices_of_copies[1].append(j)
    return indices_of_copies


def get_source_nodes_indices(sequences_sets, indices_of_copies, sources):
    source_nodes_indices = []
    for i, k in [(0, 1), (1, 0)]:
        source_nodes_indices.append(
            indices_of_copies[i] + list(range(len(sequences_sets[i]), len(sequences_sets[i]) + len(sources[k]))))
    return source_nodes_indices


def run_simple_simulation(fastas, out_dir, simulations_count, L):
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in fastas]
    indices_of_copies = get_sequences_intersections_indices(sequences_sets[0], sequences_sets[1])
    sources = determine_network_sources(sequences_sets[0], sequences_sets[1])

    if L == 0:
        L = network_creator.get_count_of_heterogeneous_positions(sequences_sets[0] + sequences_sets[1])

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = os.path.join(out_dir, fastas_basenames[0] + '_vs_' + fastas_basenames[1])
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    source_nodes_indices = get_source_nodes_indices(sequences_sets, indices_of_copies, sources)

    sources_nodes_file = os.path.join(out_dir, fastas_basenames[0] + '_vs_' + fastas_basenames[1] + '.data')
    with open(sources_nodes_file, 'w') as f:
        for i, s in enumerate(source_nodes_indices):
            f.write(str(fastas_basenames[i]) + ' ' + ' '.join(str(e) for e in s) + '\n')

    graphs_sequences = [sequences_sets[0] + [sequences_sets[1][i] for i in sources[1]],
                        sequences_sets[1] + [sequences_sets[0][i] for i in sources[0]]]

    graphs = [network_creator.ProbabilityGraphBuilder(network_creator.DistanceGraphBuilder(
        sequences, False).get_minimal_connected_graph(), L) for sequences in graphs_sequences]

    dist_out_file = [None] * len(graphs)
    prob_json_out_file = [None] * len(graphs)
    prob_gexf_out_file = [None] * len(graphs)

    for i, graph in enumerate(graphs):
        f = os.path.join(out_dir, fastas_basenames[i])
        prob_json_out_file[i] = f + '_prob.json'
#        prob_gexf_out_file[i] = f + '_prob.gexf'
#        dist_out_file[i] = f + '_dist.gexf'
#        network_creator.GexfExporter.export(graph.distance_graph, network_creator.GraphType.UNDIRECTED,
#                                            dist_out_file[i])
        network_creator.JsonExporter.export(graph.probability_graph, network_creator.GraphType.DIRECTED,
                                            prob_json_out_file[i])
#        network_creator.GexfExporter.export(graph.probability_graph, network_creator.GraphType.DIRECTED,
#                                            prob_gexf_out_file[i])

    simulations_out_dir = out_dir + '/simulation'
    if not os.path.exists(simulations_out_dir):
        os.mkdir(simulations_out_dir)

    results_log = os.path.join(out_dir, 'results.log')

    with open(results_log, 'w') as log:
        log.write('host simulation_number cycles simulation_results_file\n')
        for i, j in [[0, 1], [1, 0]]:
            network = graph_utils.import_graph_from_json(prob_json_out_file[i])
            f = os.path.join(simulations_out_dir, fastas_basenames[i])
            for k in range(simulations_count):
                log_file = f + '_' + str(k) + '.out'
                log.write(fastas_basenames[i] + ' '
                          + str(k) + ' '
                          + str(network_propagation.Propagator(network, log_file).
                                propagate(source_nodes_indices[i]))
                          + ' ' + os.path.basename(log_file) + '\n')


def run_montecarlo_simulation(fastas, out_dir, simulations_count, k_max):
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in fastas]
    network = montecarlo_network.MonteCarloNetwork(sequences_sets)
    timer = network.run_simulation(simulations_count, k_max)

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = os.path.join(out_dir, fastas_basenames[0] + '_vs_' + fastas_basenames[1])
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    results_log = os.path.join(out_dir, 'results.log')

    with open(results_log, 'w') as log:
        log.write('source timer\n')
        for i in range(2):
            log.write('{0} {1}\n'.format(fastas_basenames[i], timer[i]))


def main(fastas, out_dir, simulations_count, L, k_max=0):
    if not k_max:
        run_simple_simulation(fastas, out_dir, simulations_count, L)
    else:
        run_montecarlo_simulation(fastas, out_dir, simulations_count, k_max)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta1",
                        help="Fasta file for 1st host.")
    parser.add_argument("fasta2",
                        help="Fasta file for 2nd host.")
    parser.add_argument("-o", dest='out_dir', type=str, default=OUT_DIR,
                        help="Output directory. "
                             "Default directory is \'" + OUT_DIR + "\'")
    parser.add_argument("-n", dest='simulations_count', type=int, default=DEFAULT_SIMULATIONS_NUMBER,
                        help="Count of simulations repeats.")
    parser.add_argument("-L", dest='L', type=int, default=L,
                        help="Number of variable nucleotides.")
    parser.add_argument("-k_max", dest='k_max', type=int, default=K_MAX_DEFAULT)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    print(args.k_max)
    main([args.fasta1, args.fasta2], args.out_dir, args.simulations_count, args.L, args.k_max)
