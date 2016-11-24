# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:43:57 2016

@author: melni
input parameters:
sys.argv[1] - graphs folder (e.g, out/graphs/)
"""
import os
import sys
import statistics
from itertools import islice

graphs_folder = sys.argv[1]

target_subfolders = [name for name in os.listdir(graphs_folder)
            if os.path.isdir(os.path.join(graphs_folder, name))]

#print(target_subfolders)

target_result_files = []            
for folder in target_subfolders:
    for name in os.listdir(os.path.join(graphs_folder, folder)):
        if (name =="results.log"):
            target_result_files.append(os.path.join(graphs_folder, folder, name))
    
    
    #target_results_files.append([name for name in os.listdir(os.path.join(graphs_folder, folder))
    #        if os.path.isfile(os.path.join(graphs_folder, name))])
first_source_name = ''
second_source_name = ''

for infile in target_result_files:
    current_pair_name = os.path.basename(os.path.normpath(os.path.dirname(infile)))
    #print (current_pair_name)
    first_source_cycles = []
    second_source_cycles = []
    
    with open(infile) as f_tmp:
      # get first source name from second line (first line of file contains parameters' names)
      for line in islice(f_tmp, 1, 2):
        #print (line)
        first_source_name = line.split()[0]
        #print(first_source_name)
    
    with open(infile) as f:                
      for line in f:
         #print (line)
         if line.split()[0] == 'host':
             continue
         else:
             current_source_name = line.split()[0]
             current_cycles_num = int(line.split()[3])
             #print (current_cycles_num)
             
             if current_source_name == first_source_name:
                 first_source_cycles.append(current_cycles_num)
             else:
                 second_source_name = current_source_name
                 second_source_cycles.append(current_cycles_num)
                 
             #elif current_source_name == first_source_name:
                 #first_source_name = current_source_name
                 #first_source_cycles.append(current_cycles_num)
             #elif current_source_name == second_source_name:
                 #first_source_name = current_source_name
                 #first_source_cycles.append(current_cycles_num)
             #print (first_source_name)
             #print (second_source_name)
         #debug empty replacement
         #first_source_cycles = [0]        
         #second_source_cycles = [0] 
    
    print('debug1')
    print(first_source_name)
    print(first_source_cycles) 
    print('debug2')
    print(second_source_name)
    print(second_source_cycles)
     
    if first_source_cycles:
     first_source_median_cycles = statistics.median(first_source_cycles)
     first_source_mean_cycles = statistics.mean(first_source_cycles)
     first_source_total_cycles = sum(first_source_cycles)
              
    if second_source_cycles:
     second_source_median_cycles = statistics.median(second_source_cycles)
     second_source_mean_cycles = statistics.mean(second_source_cycles)
     second_source_total_cycles = sum(second_source_cycles)
             
    #print(first_source_cycles)
    #print(second_source_cycles)
             
    stats_dir = 'stats'
    if not os.path.exists(stats_dir):
     os.makedirs(stats_dir)
             
    #print (first_source_name)
    #print (second_source_name)
    #stats = open(stats_dir +'/'+first_source_name + '_to_' + second_source_name + '_stats.txt', 'w+')
    stats = open(os.path.join(stats_dir, current_pair_name +'_stats.txt'), 'w+')
    
    s = first_source_name + ' cycles stats. Mean: ' + str(first_source_mean_cycles) + '; median: ' + str(first_source_median_cycles)         
    #s = s + '\n'+ first_source_name + ' total_cycles: ' + str(first_source_total_cycles)
    #print(s)
    s = s + '\n' + second_source_name + ' cycles stats. Mean: ' + str(second_source_mean_cycles) + '; median: ' + str(second_source_median_cycles)         
    #s = s + '\n'+ second_source_name + ' total_cycles ' + str(second_source_total_cycles)
    #print(s)
    if (first_source_mean_cycles > second_source_mean_cycles):
        source = first_source_name
    else:
        source = second_source_name
    
    s = s + '\nSource: ' + source
    stats.write(s)
         
 
        
#print(target_result_files)
