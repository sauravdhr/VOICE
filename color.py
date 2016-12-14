# -*- coding: utf-8 -*-
"""
input parameters:
sys.argv[1] - simulation folder (e.g, out/graphs/BB44_unique_1a_12_to_BB45_unique_1a_177/simulation/)
sys.argv[2] - first graph template (e.g, out/graphs/BB44_unique_1a_12_to_BB45_unique_1a_177/BB44_unique_1a_12.dot)
sys.argv[3] - second graph template (e.g. out/graphs/BB44_unique_1a_12_to_BB45_unique_1a_177/BB45_unique_1a_177.dot)
"""
import os
from os import listdir
from os.path import isfile, join
import sys
#import pydot
#import pyparsing
#import json
import shutil
import pydotplus
import pygraphviz as pgv
from networkx.drawing import nx_pydot
from networkx.drawing import nx_pylab as nxp
import networkx as nx
import matplotlib.pyplot as plt
#import glob
#import collections
from collections import defaultdict

infolder = sys.argv[1]
first_graph_template = sys.argv[2]
second_graph_template = sys.argv[3]

#files = glob.glob(infolder)

files = [f for f in listdir(infolder) if isfile(join(infolder,f)) and f.endswith('.out')]

print(files)

infiles = []

for f in files:
    infiles.append(os.path.join(infolder, f))

print(infiles)

if infiles:
    first_result_graph = os.path.splitext(infiles[0])[0]+'.dot'
    shutil.copyfile(first_graph_template, first_result_graph)    
    G = nx.Graph(nx_pydot.read_dot(first_result_graph))
    graph = pydotplus.graphviz.graph_from_dot_file(first_result_graph)
    
    print(os.path.basename(infiles[0]))    
    print(os.path.basename(second_graph_template))
    
    for infile in infiles:
     if os.path.basename(infile)[:4] == os.path.basename(second_graph_template)[:4]:
         second_result_graph = os.path.splitext(infile)[0]+'.dot'
         shutil.copyfile(second_graph_template, second_result_graph)
         G = nx.Graph(nx_pydot.read_dot(second_result_graph))
         graph = pydotplus.graphviz.graph_from_dot_file(second_result_graph)
     idx = []
     idx_in_line = 0
     lines_num = 0
     
     print(infile)
     
     with open(infile) as f:
      for line in f: # read rest of lines
          idx.append([int(x) for x in line.split()])
                  
     print(idx)     
          
     nodes = graph.get_node_list()
     full_labels = []
     for node in nodes:
         node.set('color', 'green')
         current_label = node.get('label')
         full_labels.append(current_label)
        
     step = 1
     i = 0
     if len(idx)/10 > 0:
        step = len(idx)/10
     while (i < len(idx)):
      colors = idx[i]
      print(i)
      print(len(idx[i]))
      print(len(G.nodes()))
      print(colors)
      i+=step
      labels_keys = []
      labels_values = []
      
      #dict_labels = defaultdict(dict)
      dict_labels = {}
      for c in range (0, len(colors)):
          #dict_labels[c] = str(c)
          dict_labels[str(c)] = str(c)
          #labels_values.append(str(c))    
      
      print(G.nodes())
      print(dict_labels)
      print(colors)
      color_val_map = {}
      curr_idx = 0;
      for node in G.nodes():
          if curr_idx < len(colors):
            color_val_map[node] = colors[curr_idx]
          else:
            break;
          curr_idx += 1   
      color_values = [color_val_map.get(node) for node in G.nodes()]   
      
      #drawing
      #nx.draw_graphviz(G, cmap=plt.get_cmap('Reds'), node_color=color_values, labels=dict_labels, with_labels = True)
      ordered_nodes = G.nodes()
      #print("before")
      #print (ordered_nodes)      
      ordered_nodes = sorted(ordered_nodes, key=lambda x: int(x))
      #print("after")
      #print (ordered_nodes)
      nx.draw_graphviz(G, cmap=plt.get_cmap('Reds'), nodelist=ordered_nodes, labels = dict_labels, node_color=color_values, with_labels = True)
      #nx.draw_networkx_labels(G, dict_labels, font_size=16)
      plt.savefig(os.path.splitext(infile)[0]+'_cycle_'+str(i)+'.png', bbox_inches='tight')
      plt.clf()
      
      Matrix = nx.to_numpy_matrix(G)
      
      print(Matrix)
