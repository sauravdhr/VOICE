#!/usr/bin/env python3
"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 12/15/16
"""
from abc import ABCMeta, abstractmethod
import os
import subprocess
import multiprocessing
import time
import sys

NAME_OF_SIMULATION_SCRIPT = 'run_pair_simulation.py'
LAUNCH_PERIOD = 1


class TaskComposite(object):
    __metaclass__ = ABCMeta

    def __init__(self, out_dir=None):
        self.out_dir = out_dir

    @abstractmethod
    def list_tasks(self):
        pass


class PairTask(TaskComposite):
    def __init__(self, fasta1, fasta2, out_dir=None):
        super(PairTask, self).__init__(out_dir)
        self.fasta1 = fasta1
        self.fasta2 = fasta2

    def list_tasks(self):
        return [((self.fasta1, self.fasta2), self.out_dir)]


class MultiTasks(TaskComposite):
    def __init__(self, tasks, out_dir=None):
        super(MultiTasks, self).__init__(out_dir)
        self.tasks = tasks

    def list_tasks(self):
        tasks_list = list()
        for t in self.tasks:
            tasks_list.extend(t.list_tasks())
        return tasks_list


class TasksContainer(MultiTasks):
    def __init__(self, out_dir=None):
        super(TasksContainer, self).__init__(list(), out_dir)

    def add_task(self, task):
        self.tasks.append(task)


class OutbreakTask(MultiTasks):
    def __init__(self, outbreak_dir, out_dir=None):
        super(OutbreakTask, self).__init__(self.get_outbreak_tasks(outbreak_dir, out_dir), out_dir)

    @staticmethod
    def get_outbreak_tasks(outbreak_dir, out_dir):
        outbreak_tasks = list()
        files = os.listdir(outbreak_dir)
        for i in range(len(files) - 1):
            for j in range(i+1, len(files)):
                outbreak_tasks.append(PairTask(
                    os.path.join(outbreak_dir, files[i]), os.path.join(outbreak_dir, files[j]), out_dir))
        return outbreak_tasks


class ProcessPool(object):
    def __init__(self, cores_to_use, simulations_count, L):
        self.cores_to_use = cores_to_use
        self.simulations_count = simulations_count
        self.L = L
        self.subprocesses = [None] * cores_to_use

    def launch_new_task(self, task, out_dir):
        for i in range(len(self.subprocesses)):
            if not self.subprocesses[i]:
                print('Launching simulation for pair: ' + task[0] + ' ' + task[1])
                cmd = [sys.executable, NAME_OF_SIMULATION_SCRIPT, task[0], task[1],
                       '-o', out_dir, '-n', str(self.simulations_count)]
                if self.L:
                    cmd.extend(['-L', str(self.L)])
                print(cmd)
                self.subprocesses[i] = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                return True
        return False

    def check_results(self):
        for i in range(len(self.subprocesses)):
            if self.subprocesses[i]:
                if self.subprocesses[i].poll() is not None:
                    self.subprocesses[i] = None

    def wait(self):
        for p in self.subprocesses:
            if p:
                p.wait()


class SimulationManager(object):
    CORES_RESERVE = 4

    def __init__(self, tasks_list, simulations_count, L, cores_number=0):
        self.simulations_count = simulations_count
        self.L = L
        self.cores_number = self.get_number_of_available_cores(cores_number)
        self.outbreaks_tasks = self.determine_outbreaks_tasks(tasks_list)

    #TODO:
    @staticmethod
    def determine_outbreaks_tasks(tasks_list):
        outbreak_tasks = list()
        for t in tasks_list:
#            outbreak_tasks.append(OutbreakTask(t.normalized_in_dir, t.out_dir))
            outbreak_tasks.append(OutbreakTask(t.in_dir, t.out_dir))
        return outbreak_tasks

    @staticmethod
    def get_number_of_available_cores(desired_number_of_processes):
        available_number_of_cores = multiprocessing.cpu_count()
        number_of_processes = desired_number_of_processes
        if number_of_processes == 0:
            number_of_processes = available_number_of_cores - SimulationManager.CORES_RESERVE
        if number_of_processes > available_number_of_cores:
            number_of_processes = available_number_of_cores
        if number_of_processes < 1:
            number_of_processes = 1
        return number_of_processes

    def run_simulations(self):
        pool = ProcessPool(self.cores_number, self.simulations_count, self.L)
        i = 0
        tasks_list = list()
        for t in self.outbreaks_tasks:
            tasks_list.extend(t.list_tasks())
        while i != len(tasks_list):
            if pool.launch_new_task(tasks_list[i][0], tasks_list[i][1]):
                print('Launched {0} tasks of {1}'.format(i, len(tasks_list)))
                i += 1
                continue
            time.sleep(LAUNCH_PERIOD)
            pool.check_results()
        pool.wait()

    def run(self):
        self.run_simulations()
