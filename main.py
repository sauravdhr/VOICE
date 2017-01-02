#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12/30/16
"""
import argparse
import simulation_tasks_manager


OUT_DIR = './results'
MINIMUM_NUMBER_OF_SEQUENCES_IN_HOST_AFTER_NORMALIZATION = 10
SIMULATIONS_NOMBER = 10


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", dest='tasks_file', type=str, required=True,
                        help="File with tasks for analyzing bunches of hosts.")
    parser.add_argument("-o", dest='out_dir', type=str, default=OUT_DIR,
                        help="Path to an output directory relative to input dir. "
                             "By default it creates directory named \'" + OUT_DIR + "\' in the input directory")
    parser.add_argument("-k", dest='k_min', type=int, default=MINIMUM_NUMBER_OF_SEQUENCES_IN_HOST_AFTER_NORMALIZATION,
                        help="Minimum number of sequences in fasta files after normalization.")
    parser.add_argument("-n", dest='simulations_count', default=SIMULATIONS_NOMBER,
                        help="Count of simulations repeats.")
    parser.add_argument("-L", dest='L',
                        help="Parameter which represent a number of nucleotides prone to mutate during simulation" +
                        " (it is used in probability formula for edge weight calculation).")
    parser.add_argument("-c", dest='cores_count',
                        help="Desired number of CPU cores to run simulations.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    tasks = simulation_tasks_manager.TasksManager(args.tasks_file, args.out_dir)
    tasks.init_tasks()
    tasks.normalize_populations(args.k_min)
    tasks.run_simulations(args.simulations_count, args.L, args.cores_count)

if __name__ == '__main__':
    main()
