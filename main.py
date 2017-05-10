#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12/30/16
"""
import argparse
import simulation_tasks_manager


OUT_DIR = './results'
SIMULATIONS_NUMBER = 51
K_MAX = 10
L = 60


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", dest='tasks_file', type=str, required=True,
                        help="File with tasks description.")
    parser.add_argument("-o", dest='out_dir', type=str, default=OUT_DIR,
                        help="Path to an output directory where results will be generated. "
                             "Default out directory is \'" + OUT_DIR + "\' ")
    parser.add_argument("-k", dest='k_max', type=int, default=K_MAX,
                        help="Maximum number of sequences in fasta files after normalization.")
    parser.add_argument("-n", dest='simulations_count', default=SIMULATIONS_NUMBER,
                        help="Number of simulations.")
    parser.add_argument("-L", dest='L', type=int, default=L,
                        help="Number of variable positions in viral sequences.")
    parser.add_argument("-c", dest='cores_count', type=int,
                        help="Desired number of CPU cores to run simulations.")
    return parser.parse_args()


#TODO:
def main():
    args = parse_arguments()
    tasks = simulation_tasks_manager.TasksManager(args.tasks_file, args.out_dir)
    tasks.init_tasks()
#    tasks.normalize_populations(args.k_max)
    tasks.run_simulations(args.simulations_count, args.L, args.k_max, args.cores_count)

if __name__ == '__main__':
    main()
