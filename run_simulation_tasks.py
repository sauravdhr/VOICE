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


class TaskComposite(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def list_tasks(self):
        pass


class PairTask(TaskComposite):
    def __init__(self, fasta1, fasta2):
        self.fasta1 = fasta1
        self.fasta2 = fasta2

    def list_tasks(self):
        return [(self.fasta1, self.fasta2)]


class MultiTasks(TaskComposite):
    def __init__(self, tasks):
        self.tasks = tasks

    def list_tasks(self):
        tasks_list = list()
        for t in self.tasks:
            tasks_list += t.list_tasks()
        return tasks_list


class TasksContainer(MultiTasks):
    def __init__(self):
        super(TasksContainer, self).__init__(list())

    def add_task(self, task):
        self.tasks.append(task)


class OutbreakTask(MultiTasks):
    def __init__(self, outbreak_dir):
        super(OutbreakTask, self).__init__(self.get_outbreak_tasks(outbreak_dir))

    @staticmethod
    def get_outbreak_tasks(outbreak_dir):
        outbreak_tasks = list()
        files = os.listdir(outbreak_dir)
        for i in range(len(files) - 1):
            for j in range(i+1, len(files)):
                outbreak_tasks.append(PairTask(
                    os.path.join(outbreak_dir, files[i]), os.path.join(outbreak_dir, files[j])))
        return outbreak_tasks


class ProcessPool(object):
    def __init__(self, cores_to_use):
        self.cores_to_use = cores_to_use
        self.subprocesses = [None] * cores_to_use

    def add_new_task(self, task):
        for i in range(len(self.subprocesses)):
            if not self.subprocesses[i]:
                print('Launching simulation for pair: ' + task[0] + ' ' + task[1])
                self.subprocesses[i] = subprocess.Popen(
                    [sys.executable, 'run_pair_simulation.py', task[0], task[1], 'a', '100'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

    def __init__(self, tasks_file_name, desired_number_of_processes=0):
        self.tasks_file_name = tasks_file_name
        self.number_of_processes = self.calculate_number_of_processes(desired_number_of_processes)
        self.tasks = self.parse_tasks_file(tasks_file_name)

    @staticmethod
    def calculate_number_of_processes(desired_number_of_processes):
        available_number_of_cores = multiprocessing.cpu_count()
        number_of_processes = desired_number_of_processes
        if number_of_processes == 0:
            number_of_processes = available_number_of_cores - SimulationManager.CORES_RESERVE
        if number_of_processes > available_number_of_cores:
            number_of_processes = available_number_of_cores
        if number_of_processes < 1:
            number_of_processes = 1
        return number_of_processes

    @staticmethod
    def parse_tasks_file(tasks_file_name):
        tasks = TasksContainer()
        with open(tasks_file_name) as f:
            for l in f.readlines():
                tasks.add_task(SimulationManager.parse_task(l))
        return tasks

    @staticmethod
    def parse_task(task_string):
        a = task_string.split()
        if a[0] == 'p':
            return PairTask(a[1], a[2])
        if a[1] == 'd':
            return OutbreakTask(a[1])

    def run_simulations(self):
        pool = ProcessPool(self.number_of_processes)
        i = 0
        tasks_list = self.tasks.list_tasks()
        while i != len(tasks_list):
            print(tasks_list[i])
            if pool.add_new_task(tasks_list[i]):
                i += 1
                continue
            time.sleep(3)
            pool.check_results()
        pool.wait()

    def run(self):
        self.run_simulations()


def main(tasks_file_name, desired_number_of_processes):
    simulation_manager = SimulationManager(tasks_file_name, desired_number_of_processes)
    simulation_manager.run()


if __name__ == "__main__":
    tasks_file_name = sys.argv[1]
    desired_number_of_processes = sys.argv[2] if len(sys.argv) > 2 else 0
    main(tasks_file_name, desired_number_of_processes)
