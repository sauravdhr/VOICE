import sys
import os
import network_creator

OUT_DIR = "out/graphs"


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


def main(fastas, max_dist_for_median_triple):
    graphs = [network_creator.ProbabilityGraphBuilder(network_creator.DistanceGraphBuilder(
        fasta, max_dist_for_median_triple).get_minimal_connected_graph()) for fasta in fastas]

    sequences = [list(filter(lambda x: x,
                             map(lambda v: v['sequence'] if 'sequence' in v else None, g.distance_graph.vertices)))
                 for g in graphs]

    sources = determine_network_sources(sequences[0], sequences[1])

    fastas_basenames = [os.path.splitext(os.path.basename(f))[0] for f in fastas]
    out_dir = OUT_DIR + '/' + fastas_basenames[0] + '_to_' + fastas_basenames[1]
    os.mkdir(out_dir)

    for i, graph in enumerate(graphs):
        f = os.path.join(out_dir, fastas_basenames[i])
        out_file_json = f + '.json'
        #out_file_dot = f + '.dot'
        #network_creator.DotExporter.export(graph.probability_graph, out_file_dot)
        network_creator.JsonExporter.export(graph.probability_graph, out_file_json)


if __name__ == "__main__":
    fastas = [sys.argv[1], sys.argv[2]]
    if len(sys.argv) <= 3:
        max_dist_for_triple_mean = -1
    else:
        max_dist_for_triple_mean = int(sys.argv[3])
    main(fastas, max_dist_for_triple_mean)
