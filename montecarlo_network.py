#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 5/9/17
"""

import network_creator
import random
import copy
import statistics
import os


MUTATION_PROBABILITY = 0.01
VICINITY = 1


class MonteCarloNetwork(object):
    e = MUTATION_PROBABILITY
    s = e / (1 - 3 * e)

    def __init__(self, seqs):
        self.seqs = seqs
        self.distance_matrix = network_creator.infer_distance_matrix(sum(self.seqs, []))

    def get_seqs_subset_inds(self, k):
        prev_ind = 0
        seqs_inds = []
        for i in range(len(self.seqs)):
            new_ind = prev_ind + len(self.seqs[i])
            seqs_inds.append(list(range(prev_ind, new_ind)))
            prev_ind = new_ind

        seqs_subset_inds = [[] for _ in range(len(self.seqs))]
        for i in range(len(self.seqs)):
            random.shuffle(seqs_inds[i])
            if len(seqs_inds[i]) >= k:
                seqs_subset_inds[i] = set(seqs_inds[i][:k])
            else:
                seqs_subset_inds[i] = set(seqs_inds[i])
        return seqs_subset_inds

    def run_simulation(self, runs, k):
        spreading_time = [[0] * runs for _ in range(len(self.seqs))]
        for i in range(runs):
            print("run#: {0}".format(i+1))
            seqs_subset_inds = self.get_seqs_subset_inds(k)
            for j in range(len(self.seqs)):
                recipient = sum([seqs_subset_inds[:j], seqs_subset_inds[j+1:]], [])
                spreading_time[j][i] = self.single_run(seqs_subset_inds[j],
                                                       set.union(*recipient))
        return list(map(lambda x: statistics.median(x), spreading_time))

    @staticmethod
    def is_mutate():
        return random.random() < MonteCarloNetwork.s

    def single_run(self, source_inds, recepient_inds):
        timer = 0
        simulation_matrix = copy.deepcopy(self.distance_matrix)
        source_inds_filtered = self.filter_sources(source_inds, recepient_inds)
#        source_inds_filtered = source_inds
        while True:
            if not recepient_inds:
                break
            new_sources = set()
            for s in source_inds_filtered:
                for r in recepient_inds:
                    if simulation_matrix[s, r] == 0:
                        new_sources.add(r)
                    else:
                        if self.is_mutate():
                            simulation_matrix[s, r] -= 1
            timer += 1
            source_inds_filtered = source_inds_filtered.union(new_sources)
            recepient_inds = recepient_inds.difference(new_sources)
        return timer

    def filter_sources(self, source_inds, recipient_inds):
        filtered_sources = set()
        dist = 10e6
        for s in source_inds:
            for r in recipient_inds:
                if self.distance_matrix[s, r] < dist:
                    dist = self.distance_matrix[s, r]
        for s in source_inds:
            for r in recipient_inds:
                if self.distance_matrix[s, r] <= dist + VICINITY:
                    filtered_sources.add(s)
        return filtered_sources


if __name__ == '__main__':
    outbreak_path = 'data/AW_clipped'
    fasta_files = [os.path.join(outbreak_path, f) for f in os.listdir(outbreak_path)
                   if os.path.isfile(os.path.join(outbreak_path, f))]
    sequences_sets = [network_creator.parse_fasta(fasta_name) for fasta_name in
                      fasta_files]
    network = MonteCarloNetwork(sequences_sets)
    timer = network.run_simulation(21, 4)
    samples = map(lambda x: x.split('/')[-1].split('_')[0], fasta_files)
    print(list(zip(samples, timer)))
