# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 21:27:46 2016

@author: melni
"""

import sys
from main import main
import math
from os import listdir, path, walk
from os.path import isfile, join
from itertools import combinations
import fnmatch
import itertools
from collections import defaultdict
import multiprocessing
import subprocess
from subprocess import STDOUT
from threading import Timer
from time import sleep

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
    for i in range (0, init_range):
        if simulations_subprocesses[i].poll() == 0:
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
            simulations_subprocesses[i] = subprocess.Popen(['D:\Program Files\WinPython-64bit-3.4.4.5Qt5\python-3.4.4.amd64\python.exe', 'D:\bioinfmatics\from_git\voicerep\main.py', pair[0], pair[1], 'a', '100'])
        
        
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

print(cores_num)

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

subprocesses = []
subprocesses_pairs = []
completed_pairs = [None] * len(pairs)

# Launch first simulations
init_range = 0
if (pairs_num > cores_to_use):
    init_range = cores_to_use
else:
    init_range = pairs_num
    
for i in range (0, init_range):
    #for pair in pairs:    
        current_pair = list(pairs[i])
        #print('Launching simulation for pair: ' + current_pair[0] + ' ' + current_pair[1])
        subprocesses_pairs.append(pairs[i])
        p = subprocess.Popen(['D:\Program Files\WinPython-64bit-3.4.4.5Qt5\python-3.4.4.amd64\python.exe', 'D:\\bioinfmatics\\from_git\\voicerep\\main.py', pairs[i][0], pairs[i][1], 'a', '100'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocesses.append(p)
        for s in subprocesses:
            for line in s.stderr:
                print(line, end='')
            for line in s.stdout:
                print(line, end='')

        #p =  subprocess.Popen(['D:\Program Files\WinPython-64bit-3.4.4.5Qt5\python-3.4.4.amd64\python.exe', 'D:\bioinfmatics\from_git\voicerep\main.py', pairs[i][0], pairs[i][1], 'a', '100'], stdout=subprocess.PIPE)

# Check for completed simulations every 3 seconds
#rt = RepeatedTimer(3, replace_completed_subprocesses, pairs, subprocesses, subprocesses_pairs, completed_pairs, init_range)       
        

#for pair in pairs:
 #   main(pair, 100, 'a')