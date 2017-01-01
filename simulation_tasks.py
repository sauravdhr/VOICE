#!/usr/bin/env python3

"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12/16/16
"""
import subprocess
import sys
import os
import shutil
import run_simulation_tasks
from abc import ABCMeta, abstractmethod


NORMALIZE_POPULATION_SCRIPT_NAME = "normalize_populations.py"
NORMALIZED_POPULATION_OUT_DIR = "normalized_population"
CORES_COUNT = 3


class Task(object):
    __metaclass__ = ABCMeta

    def __init__(self, in_dir, out_dir):
        self.in_dir = in_dir
        self.normalized_in_dir = os.path.join(out_dir, NORMALIZED_POPULATION_OUT_DIR)
        self.out_dir = out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def normalize_samples(self, k_min, out_dir=None):
        if not out_dir:
            out_dir = self.out_dir
        subprocess.call([sys.executable, NORMALIZE_POPULATION_SCRIPT_NAME,
                         "-i", self.in_dir,
                         "-o", self.normalized_in_dir,
                         "-k", str(k_min)])


class TwoHostsTask(Task):
    IN_DIR = 'in'

    def __init__(self, fasta1_file_name, fasta2_file_name, out_dir):
        out_path = os.path.join(out_dir,
                                os.path.basename(fasta1_file_name).split('.')[0] +
                                '_vs_' +
                                os.path.basename(fasta2_file_name).split('.')[0])
        in_path = os.path.join(out_path, self.IN_DIR)
        super(TwoHostsTask, self).__init__(in_path, out_path)
        if not os.path.exists(self.in_dir):
            os.makedirs(self.in_dir)

        shutil.copy(fasta1_file_name, self.in_dir)
        shutil.copy(fasta2_file_name, self.in_dir)


class OutbreakTask(Task):
    def __init__(self, outbreak_in_dir, out_dir):
        super(OutbreakTask, self).__init__(outbreak_in_dir, out_dir)


class TasksManager(object):
    def __init__(self, tasks_file_name, out_dir):
        self.tasks_file_name = tasks_file_name
        self.tasks = None
        self.out_dir = out_dir

    def init_tasks(self):
        self.tasks = list()
        with open(self.tasks_file_name) as f:
            for l in f.readlines():
                self.tasks.append(TasksManager.parse_task(l, self.out_dir))

    @staticmethod
    def parse_task(task_string, out_dir):
        a = task_string.split()
        if a[0] == 'p':
            return TwoHostsTask(a[1], a[2], out_dir)
        if a[0] == 'd':
            return OutbreakTask(a[1], os.path.join(out_dir, os.path.basename(os.path.normpath(a[1]))))

    def normalize_populations(self, k_min):
        for t in self.tasks:
            t.normalize_samples(k_min)

    def run_simulations(self, simulations_count, L):
        run_simulation_tasks_manager = run_simulation_tasks.SimulationManager(self.tasks,
                                                                              simulations_count, L, CORES_COUNT)
        run_simulation_tasks_manager.run()
