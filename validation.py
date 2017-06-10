#!/usr/bin/env python3
"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 05.30.2017
"""
import analyze_outbreak_graph
import glob
import os

SIMULATION_FOLDS = 'data/all_clipped_simulation_5_fold/'
SUBSAMPLING_FOLDS = 'data/all_clipped_subsam_mean_2_5_fold/'
MINDIST_FOLDS = 'data/all_clipped_min_dist_5_fold/'
MINDIST_PLUS_BORDER_FOLDS = 'data/all_clipped_min_dist_v_4_log_3_5_fold/'

CLASSIFIERS = [('simulation', SIMULATION_FOLDS),
               ('subsampling', SUBSAMPLING_FOLDS),
               ('MinDist', MINDIST_FOLDS),
               ('MinDistPlusBorder', MINDIST_PLUS_BORDER_FOLDS)]


class Validator(object):
    def __init__(self, train_file, test_file):
        self.train_file = train_file
        self.test_file = test_file
        self.relatedness_threshold = 0
        self.cluster_threshold = 0
        self.relatedness_sensitivity = 0
        self.relatedness_specificity = 0

        self.train_graph_analyzer = analyze_outbreak_graph.UndirectedGraphAnalyzer(train_file)
        self.test_graph_analyzer = analyze_outbreak_graph.UndirectedGraphAnalyzer(test_file)

        self.train_threshold()
        self.calc_statistics()

    def train_threshold(self):
        self.relatedness_threshold = self.train_graph_analyzer.get_threshold_with_no_false_related()
        self.cluster_threshold = self.train_graph_analyzer.get_threshold_with_no_false_outbreak_linkage()

    def calc_statistics(self):
        self.relatedness_sensitivity = self.test_graph_analyzer.get_relatedness_sensitivity(self.relatedness_threshold)
        self.relatedness_specificity = self.test_graph_analyzer.get_relatedness_specificity(self.relatedness_threshold)


class KFoldValidator(object):
    def __init__(self, data_set_dir):
        self.data_set_dir = data_set_dir
        self.data_set_file_names = self.get_data_set_file_names(self.data_set_dir)
        self.validators = [None] * len(self.data_set_file_names)
        for i in range(len(self.data_set_file_names)):
            self.validators[i] = Validator(self.data_set_file_names[i][0], self.data_set_file_names[i][1])
            self.validators[i].calc_statistics()

    @staticmethod
    def get_data_set_file_names(data_set_dir):
        k = len(glob.glob(os.path.join(data_set_dir, 'train_*')))
        return [(os.path.join(data_set_dir, 'train_{0}.edgelist'.format(i)),
                 os.path.join(data_set_dir, 'test_{0}.edgelist'.format(i)))
                for i in range(k)]


def main():
    for classifier in CLASSIFIERS:
        print(classifier[0])
        k_fold_validator = KFoldValidator(classifier[1])
        for i in range(len(k_fold_validator.validators)):
            print('Data set #{0}'.format(i))
            print('Relatedness sensitivity: {0}, relatedness specificity: {1}'.format(
                  k_fold_validator.validators[i].relatedness_sensitivity,
                  k_fold_validator.validators[i].relatedness_specificity))

if __name__ == '__main__':
    main()
