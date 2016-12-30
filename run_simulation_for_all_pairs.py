# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 21:27:46 2016

@author: melni
"""

import sys
from run_pair_simulation import main
import math
import os
from os import listdir, path, walk
from os.path import isfile, join
from itertools import combinations
import fnmatch
import itertools
from collections import defaultdict
import multiprocessing
import subprocess
from subprocess import Popen, STDOUT
from threading import Timer
from time import sleep
import ntpath

class RepeatedTimer(object):
    def __init__(self, interval, function, *args, **kwargs):
        self._timer     = None
        self.interval   = interval
        self.function   = function
        self.args       = args
        self.kwargs     = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self):
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False

def replace_completed_subprocesses(all_pairs, simulations_subprocesses, simulations_subprocesses_pairs, completed_pairs, init_range):
    for i in range (0, len(all_pairs)):
        if simulations_subprocesses[i] and simulations_subprocesses[i].poll() == 0:
            # This simulation is completed, run another one
            # Before this, remove completed simulation from running pairs array and add it to completed ones
            pair_to_launch = []
            completed_pairs.append(simulations_subprocesses_pairs[i])
            for pair in all_pairs:
                if list(pair) not in completed_pairs:
                    pair_to_launch = list(pair)
                    break
                if not pair_to_launch:
                    sys.exit()
                simulations_subprocesses_pairs[i] = pair_to_launch
                simulations_subprocesses[i] = subprocess.Popen(['python3.4', 'run_pair_simulation.py', pair[0], pair[1], 'a', '100'])
        
        
# directory with all sequences
in_folder = sys.argv[1]
# All in one directory
#files = [f for f in listdir(in_folder) if isfile(join(in_folder, f))]
# Recursive
files = []
for root, dirnames, filenames in walk(in_folder):
    for filename in fnmatch.filter(filenames, '*.fas'):
        files.append(path.join(root, filename))          

pairs = list(itertools.combinations(files, 2))       

pairs_num = len(pairs)

cores_num = multiprocessing.cpu_count()

s = 'Total input pairs: ' + str(pairs_num)
s = s + '\nCores: ' + str(cores_num)
with open("all_pairs_simulations_log.txt","a+") as f:
    f.write(s)
    
print('Total input pairs: ' + str(pairs_num))

print('Number of available cores : ' + str(cores_num))

if (cores_num - 7 > 0):
    cores_to_use = cores_num - 7
else:
    cores_to_use = cores_num    

print('Number of cores to be used: ' + str(cores_to_use))

#subprocesses = [None] * len(pairs)
#subprocesses_pairs = [None] * len(pairs)
#completed_pairs = [None] * len(pairs)

print("Number of pairs: " + str(len(pairs)))

class ProcessPool(object):
    def __init__(self, cores_to_use):
        self.cores_to_use = cores_to_use
        self.subprocesses = [None] * cores_to_use

    def add_new_task(self, task):
        for i in range(len(self.subprocesses)):
            if not self.subprocesses[i]:
                print('Launching simulation for pair: ' + task[0] + ' ' + task[1])
                self.subprocesses[i] = Popen(['python3.4', 'run_pair_simulation.py', task[0], task[1], 'a', '100'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                return True
        return False

    def check_results(self):
        for i in range(len(self.subprocesses)):
            if self.subprocesses[i]:
                if self.subprocesses[i].poll() != None:
                    #self.subprocesses[i].kill()
                    self.subprocesses[i] = None
                    #print('Finished simulation for pair: ' + task[0] + ' ' + task[1])

    def wait(self):
        for p in self.subprocesses:
            if p:
                p.wait()

out = sys.argv[2]
completed_pairs_folders = [dI for dI in os.listdir(out) if os.path.isdir(os.path.join(out,dI))]
completed_pairs = []

for folder in completed_pairs_folders:
      completed_pairs.append(folder.split('_to_'))

print('Already completed pairs: ' + str(len(completed_pairs)))

pairs_to_complete = []
for pair in pairs:
    already_completed = False
    for completed_pair in completed_pairs:
        if (os.path.splitext(ntpath.basename(list(pair)[0]))[0] == completed_pair[0] and os.path.splitext(ntpath.basename(list(pair)[1]))[0] == completed_pair[1]):
            already_completed = True
    if not already_completed:
        pairs_to_complete.append(pair)
    else:
        print("Simulation for pair " + list(pair)[0]  + '_to_' + list(pair)[1] + ' has already been completed. Skipping.')

print('Pairs to be completed: ' + str(len(pairs_to_complete)))

pool = ProcessPool(cores_to_use)
i = 0
while i != len(pairs_to_complete):
    print(pairs_to_complete[i])
    if pool.add_new_task(pairs_to_complete[i]):
        i += 1
        continue
    sleep(3)
    pool.check_results()
pool.wait()


