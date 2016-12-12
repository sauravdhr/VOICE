#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12.01.2016
"""

import sys
import os
import network_creator
import propagation
import graph_utils

OUT_DIR = "out/graphs"
SIMULATIONS_NUMBER = 1


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


def main(fastas, L):
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in fastas]
    indices_of_copies = get_sequences_intersections_indices(sequences_sets[0], sequences_sets[1])
    sources = determine_network_sources(sequences_sets[0], sequences_sets[1])

    if L == 0:
        L = network_creator.get_count_of_heterogeneous_positions(sequences_sets[0] + sequences_sets[1])

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = OUT_DIR + '/' + fastas_basenames[0] + '_to_' + fastas_basenames[1]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, fastas_basenames[0] + '_to_' + fastas_basenames[1] + '.data')
    l = []
    with open(out_file, 'w') as f:
        for i, k in [(0, 1), (1, 0)]:
            l = list(range(len(sequences_sets[i]), len(sequences_sets[i]) + len(sources[k])))
            f.write(str(fastas_basenames[i]) + ' ' + ' '.join(str(e) for e in l) + '\n')
            graphs_sequences = [[sequences_sets[1][i] for i in sources[1]] + sequences_sets[0],
                                [sequences_sets[0][i] for i in sources[0]] + sequences_sets[1]]

    graphs = [network_creator.ProbabilityGraphBuilder(network_creator.DistanceGraphBuilder(
        sequences, False).get_minimal_connected_graph(), L) for sequences in graphs_sequences]

    out_file_dots = [None] * len(graphs)
    out_file_jsons = [None] * len(graphs)

    for i, graph in enumerate(graphs):
        f = os.path.join(out_dir, fastas_basenames[i])
        out_file_jsons[i] = f + '.json'
        out_file_dots[i] = f + '.dot'
        network_creator.DotExporter.export(graph.distance_graph, out_file_dots[i])
        network_creator.JsonExporter.export(graph.probability_graph, out_file_jsons[i])

    simulations_out_dir = out_dir + '/simulation'
    if not os.path.exists(simulations_out_dir):
        os.mkdir(simulations_out_dir)

    results_log = os.path.join(out_dir, 'results.log')

    with open(results_log, 'w') as log:
        log.write('host simulation_number cycles simulation_results_file\n')
        for i, j in [[0, 1], [1, 0]]:
            network = graph_utils.import_graph(out_file_jsons[i])
            f = os.path.join(simulations_out_dir, fastas_basenames[i])
            for k in range(SIMULATIONS_NUMBER):
                log_file = f + '_' + str(k) + '.out'
                log.write(fastas_basenames[i] + ' '
                          + str(k) + ' '
                          + str(propagation.Propagator(network, log_file).
                                propagate(sources[i] + [ind + len(sources[i]) for ind in indices_of_copies[j]]))
                          + ' ' + os.path.basename(log_file) + '\n')


if __name__ == "__main__":
    fastas = [sys.argv[1], sys.argv[2]]
    L = float(sys.argv[4]) if len(sys.argv) > 5 else 0
    main(fastas, L)
