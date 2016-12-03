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

OUT_DIR = "out/graphs"
SIMULATIONS_NUMBER = 5


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
    return [list(map(lambda x: x[0], filter(lambda x: x[1] <= min_dist + VICINITY, enumerate(s))))
            for s in [dist_s1, dist_s2]]


def main(fastas, L, mode):
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in fastas]
    sources = determine_network_sources(sequences_sets[0], sequences_sets[1])
    if L == 0:
        L = network_creator.get_count_of_heterogeneous_positions(sequences_sets[0] + sequences_sets[1])

#    sequences = [list(filter(lambda x: x,
#                             map(lambda v: v['sequence'] if 'sequence' in v else None, g.distance_graph.vertices)))
#                 for g in graphs]

#    sources = determine_network_sources(sequences[0], sequences[1])

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = OUT_DIR + '/' + fastas_basenames[0] + '_to_' + fastas_basenames[1]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, fastas_basenames[0] + '_to_' + fastas_basenames[1] + '.data')
    l = []
    with open(out_file, 'w') as f:
        for i, j in [(0, 1), (1, 0)]:
            if mode == 's':
                l = sources[i]
            elif mode == 'a':
                l = list(range(len(sequences_sets[i]), len(sequences_sets[i]) + len(sources[j])))
            f.write(str(fastas_basenames[i]) + ' ' + ' '.join(str(e) for e in l) + '\n')

    graphs_sequences = [sequences_sets[0] + [sequences_sets[1][i] for i in sources[1]],
                        sequences_sets[1] + [sequences_sets[0][i] for i in sources[0]]] if mode == 'a' \
        else [sequences_sets[0], sequences_sets[1]]

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
        if mode == 's':
            log.write('host source_node simulation_number cycles simulation_results_file\n')
        elif mode == 'a':
            log.write('host simulation_number cycles simulation_results_file\n')
        for i, json in enumerate(out_file_jsons):
            network = propagation.import_graph(json)
            f = os.path.join(simulations_out_dir, fastas_basenames[i])
            for j in range(SIMULATIONS_NUMBER):
                if mode == 's':
                    for k in sources[i]:
                        log_file = f + '_' + str(k) + '_' + str(j) + '.out'
                        log.write(fastas_basenames[i] + ' '
                                  + str(k) + ' '
                                  + str(j) + ' '
                                  + str(propagation.Propagator(network, log_file).propagate([k])) + ' '
                                  + os.path.basename(log_file) + '\n')
                elif mode == 'a':
                    log_file = f + '_' + str(j) + '.out'
                    log.write(fastas_basenames[i] + ' '
                              + str(j) + ' '
                              + str(propagation.Propagator(network, log_file).propagate(sources[i])) + ' '
                              + os.path.basename(log_file) + '\n')


if __name__ == "__main__":
    fastas = [sys.argv[1], sys.argv[2]]
    mode = sys.argv[3] if len(sys.argv) > 3 else 'a'
    L = float(sys.argv[4]) if len(sys.argv) > 4 else 0
    main(fastas, L, mode)
