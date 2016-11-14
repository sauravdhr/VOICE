import sys
import os
import network_creator
import propagation

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
    return [list(map(lambda x: x[0], filter(lambda x: x[1] <= min_dist + VICINITY, enumerate(s))))
            for s in [dist_s1, dist_s2]]


def main(fastas):
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in fastas]
    L = network_creator.get_count_of_heterogeneous_positions(sequences_sets[0] + sequences_sets[1])
    print(L)
    graphs = [network_creator.ProbabilityGraphBuilder(network_creator.DistanceGraphBuilder(
        sequences, False).get_minimal_connected_graph(), L) for sequences in sequences_sets]

    sequences = [list(filter(lambda x: x,
                             map(lambda v: v['sequence'] if 'sequence' in v else None, g.distance_graph.vertices)))
                 for g in graphs]

    sources = determine_network_sources(sequences[0], sequences[1])

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = OUT_DIR + '/' + fastas_basenames[0] + '_to_' + fastas_basenames[1]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, fastas_basenames[0] + '_to_' + fastas_basenames[1] + '.data')
    with open(out_file, 'w') as f:
        for i in range(len(fastas_basenames)):
            f.write(str(fastas_basenames[i]) + ' '.join(str(e) for e in sources[i]) + '\n')

    out_file_dots = [None] * len(graphs)
    out_file_jsons = [None] * len(graphs)

    for i, graph in enumerate(graphs):
        f = os.path.join(out_dir, fastas_basenames[i])
        out_file_jsons[i] = f + '.json'
        out_file_dots[i] = f + '.dot'
        network_creator.DotExporter.export(graph.probability_graph, out_file_dots[i])
        network_creator.JsonExporter.export(graph.probability_graph, out_file_jsons[i])

    simulations_out_dir = out_dir + '/simulation'
    if not os.path.exists(simulations_out_dir):
        os.mkdir(simulations_out_dir)

    for i, json in enumerate(out_file_jsons):
        network = propagation.import_graph(json)
        f = os.path.join(simulations_out_dir, fastas_basenames[i])
        for j in range(SIMULATIONS_NUMBER):
            log_file = f + '_' + str(j) + '.out'
            print(propagation.Propagator(network, log_file).propagate(sources[i]))

if __name__ == "__main__":
    fastas = [sys.argv[1], sys.argv[2]]
    main(fastas)
