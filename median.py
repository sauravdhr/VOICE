# -*- coding: utf-8 -*-
"""
Created on Thu Oct 01 10:24:10 2015

@author: Olga
"""

'''
# prints a list of all fileNAMES (NOT contents) in the folder
import glob, os
os.chdir("data")
combined = ""
for file in glob.glob("*.fas"):
    print(file) # prints filename!
'''



'''
#------------------------------------------------------------------------------------------------------------------------- 
from os import walk # makes a list of all fileNAMES (NOT contents)
f = []
for (dirpath, dirnames, filenames) in walk("data"):
    f.extend(filenames)
    break
print f # f= ['file1.fas', 'file2.fas', ...],=


with open('combined', 'w') as outfile: # creates a file, combining all the contents of the filename array created earlier
    for fname in f:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

#------------------------------------------------------------------------------------------------------------------------- 
'''


'''
# reading the files                
myfile1 = open('AA8_unique_1b_88.fas', 'r')
myfile2 = open('AA20_unique_1b_43.fas', 'r')
file_contents1 = myfile1.read()
file_contents2 = myfile1.read()
#print (file_contents)
myfile1.close()
'''
#-------------------------------------------------------------------------------------------------------------------------

'''
#saving header info in a list, and sequence info into the other list
# something wrong with sys.argv[1]


from Bio import SeqIO
import sys

headerList = []
seqList = []

inFile = open(sys.argv[1],'r')
for record in SeqIO.parse('AA8_unique_1b_88.fas','fasta'):
   headerList.append(record.id)
   seqList.append(str(record.seq))

print headerList
print "OMG"
print seqList

'''

#-------------------------------------------------------------------------------------------------------------------------

#import os
#myfilename = raw_input('Enter a file name with extension: ')

'''
INPUT: fasta file
OUTPUT: list of seqs_names; dictionary of pairs "seq_name: seqs" 

This procedure creates a list 'order' with all names of fasta sequences from a file myfilename in an order they were incountered;
it also creates a dictionary with pairs of names and respective sequences from a file myfilename. 
NOTE that the order of this pairs in a dictionary is NOT the order in which it was placed in the original file!
NOTE that the command sequences.keys() will produce the names of all fasta files, but NOT in the original order. 
So we keep a list 'order' - just in case,not sure if we need it that crucially :)
#print len(myfilename) # print the length of a file name; not the length of the contents
'''
def FASTA(filename):           
    try:
        f = file(filename)
    except IOError:
        print "The file, %s, does not exist" % filename
        return
    order = [] # keeps identifiers-names-tags in a list, length = # sequences in file
    sequences = {} # so far dictionary - so it keeps pairs of "identifier-name-tag: sequences"
    seqs_values = []

    for line in f:
        if line.startswith('>'):
            name = line[1:].rstrip('\n') # getting rid of sequence identifiers
            #name = name.replace('_', ' ') # underscores ('_') are replaced by spaces
            order.append(name)
            sequences[name] = '' #remembers sequence identifier w/o '>'
        else:
            sequences[name] += line.rstrip('\n').rstrip('*')
    seqs_values = sequences.values()
    
    #print "length of vector:", len(value)  
    #print "%d sequences found" % len(order) # prints number of sequences in fasta file
    #print "# of sequences in each file:", len (seqs_values)     ## same thing: number of sequences in each file 
    
    return order, seqs_values, sequences # 1 keys, 2 values, 3 dictionary of pairs of keys and values        

#-------------------------------------------------------------------------------------------------------------------------
   
    
'''

# looping through all fils in a directory
import glob
path = "*.fas"  # "path/to/dir/*.fas"
#i=0
for fname in glob.glob(path):
    tags, vectors, pairs = FASTA(fname)

    #print "a", a
    #a=len(sequences.FASTA(fname))
    #i=i+1
    #print(fname)
#print "files found:", i

'''
#-------------------------------------------------------------------------------------------------------------------------
    
#header = FASTA.order
#vectors = FASTA.sequences
#
#------------------------------------------------------------- 
#myfilename='BB45_unique_1a_177.fas' # 'BB45_unique_1a_177.fas' has 92 sequences - the biggest; 'BJ30_unique_1a_6' has 5 seqs
#header, vectors = FASTA (myfilename)

#------------------------------------------------------------- 
#print header
#print vectors
#fasta3 = "CACCTACAGCAGCCCTAGTGGTATCGCAGTTACTCCGGATCCCACAAGCTGTCGTGGATATGGTAGCGGGGGCCCACTGGGGAGTCCTGGCGGGCCTCGCATACTATTCCATGGTGGGGAACTGGGCTAAGGTTTTGATTGTGATGCTGCTCTTTGCCGGCGTTGACGGGGTGACCTACACGACGGGGGGGGCGACGGCCCGTAATACTCACAGGCTGACGTCCTTCTTATCGACTGGGTCGGCTCAGAACATCCAGCTTATAA"
#fasta4 = "CCCCTACGACGGCGCTGGTAGTAGCTCAGCTGCTCCGGATCCCACAAGCCATCTTGGACATGATCGCTGGTGCCCACTGGGGAGTCCTAGCGGGCATGGCGTATTTCTCCATGGTGGGGAACTGGGCGAAGGTCGTGGTAGTGCTGCTGCTATTCGCCGGCGTCGACGCGGAAACCCGCGTCTCCGGGGGAACTGCTGGCCGCACTACGGCTGGATTTGTCGGGTTCCTCACACAAGGCGCCAAGCAGAACATCCAGCTGATCA"
#print "length vector", len(fasta4)


#-------------------------------------------------------------------------------------------------------------------------

'''
looping through all subsets of order m=3 from a set of vectors S                                                                                                                                                                                                                                                                                                 
'''
import itertools
def findsubsets(S,m):
    return set(itertools.combinations(S, m))
#-------------------------------------------------------------------------------------------------------------------------


'''
finding median out of 3 strings (consensus, majority vote)
'''
#fasta1 = "CACCTACAGCAGCCCTAGTGGTATCGCAGTTACTCCGGATCCCACAAGCTGTCGTGGATATGGTAGCGGGGGCCCACTGGGGAGTCCTGGCGGGCCTTGCATACTATTCCATGGTGGGGAACTGGGCTAAGGTTTTGATTGTGATGCTGCTCTTTGCCGGCGTTGACGGGGTGACCTACACGACGGGGGGGGCGACGGCCCGTAATACTCACAGGCTGACGTCCTTCTTATCGACTGGGTCGGCTCAGAACATCCAGCTTATAA"
#fasta2 = "CACCTACAGCAGCCCTAGTGGTATCGCAGTTACTCCGGATCCCACAAGCTGTCGTGGATATGGTAGCGGGGGCCCACTGGGGAGTCCTGGCGGGCCTTGCATACTATTCCATGGTGGGGAACTGGGCTAAGGTTTTGATTGTGATGCTGCTCTTTGCCGGCGTTGACGGGGTGACCTACACGACGGGGGGGGCGACGGCCCGTAATACTCACAGGCTGACGTCCTTCTTATCGACTGGGTCGGCTCAGACCATCCAGCTTATAA"
#fasta3 = "CACCTACAGCAGCCCTAGTGGTATCGCAGTTACTCCGGATCCCACAAGCTGTCGTGGATATGGTAGCGGGGGCCCACTGGGGAGTCCTGGCGGGCCTCGCATACTATTCCATGGTGGGGAACTGGGCTAAGGTTTTGATTGTGATGCTGCTCTTTGCCGGCGTTGACGGGGTGACCTACACGACGGGGGGGGCGACGGCCCGTAATACTCACAGGCTGACGTCCTTCTTATCGACTGGGTCGGCTCAGAACATCCAGCTTATAA"


#fasta1 ="CCCCTACGGCGGCGTTGGTAATGGCTCAGCTGCTCCGGATCCCGCAAGCCATCGTGGACATGATCGCTGGTGCTCACTGGGGAGTCCTAGCGGGCATAGCGTATTTCTCCATGGTGGGGAACTGGGCGAAGGTCCTGGTAGTGCTGCTGCTATTTGCTGGCGTCGACGCGGAAACCCAGGTCACCGGGGGAAGTGCCGCTCGCGCCGCGTCTGGGCTTGCAAGCTTTTTCTCACCAGGCGCCAAGCAGAACATCCAGCTGATTA"
#fasta2 ="CCCCTACAGCGGCGTTGGCAATGGCTCAGCTGCTCCGGATCCCACAAGCCATCGTGGACATGATCGCTGGTGCTCACTGGGGAGTCCTAGCGGGCATAGCGTATTTCTCCATGGTGGGGAACTGGGCGAAGGTCCTGGTAGTGCTGCTGCTATTTGCTGGCGTCGACGCGGAAACCCAGGTCACCGGGGGAAGTGCCGCTCGCGCCGCGTCTGGGCTTGCAAGCTTTTTCTCACCAGGCGCCAAGCAGAACATCCAGCTGATTA"
#fasta3 ="CGCCCACAGCGGCGCTGGCAATGGCTCAGCTGCTCCGGATCCCACAAGCCATCGTGGACATGATCGCTGGTGCTCACTGGGGAGTCCTAGCGGGCATAGCGTATTTCTCCATGGTGGGGAACTGGGCGAAGGTCCTGGTAGTGCTGCTGCTATTTGCTGGCGTCGACGCGGAAACCCAGGTCACCGGGGGAAGTGCCGCTCGCGCCGCGTCTGGGCTTGCAAGCTTTTTCTCACCAGGCGCCAAGCAGAACATCCAGCTGATTA"


def findMedian(vector1, vector2, vector3):
    no_median_counter = 0
    median = ""
    if len(vector1)==len(vector2) and len(vector1)==len(vector3):
        #print "lengthes match"
        for i in range (0, len(vector1)):
            #print "i:", i
            if vector1[i]==vector2[i] or vector1[i]==vector3[i]:
                median+=vector1[i]
            elif vector2[i]==vector3[i]:
                median+=vector2[i]
                #print "median so far: ", median
            else:
                #print "no median possible for given vectors"
                no_median_counter = no_median_counter +1
    #print "no medians were constructed", no_median_counter, "times"        
    return median, no_median_counter


'''
#finding ???

def difference(vector, median):
    sum=0
    if len(vector)==len(median):
        for i in range (0, len(vector)):
            if vector[i]!=median[i]:
                sum=sum+1
        return sum
'''                
#def Unique_FASTA(filename):
#    tags, vectors, pairs = FASTA(filename)
    
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------    
'''
# looping through all fils in a directory
import glob
path = "*.fas"  # "path/to/dir/*.fas"
file_counter=0
for fname in glob.glob(path):
    #print "filename: ", path
    print "file considered:", file_counter+1 #i+1
    tags, vectors, pairs = FASTA(fname)   
    triplets_list = findsubsets(set(vectors),3) # for each three sequences, find a median; output is a set([ (), ... ,() ])

    medians_list = []
    no_median_sum=0
    for triplet in triplets_list: ## triplet is a list of 3 strings
        a,b=findMedian(triplet[0],triplet[1],triplet[2])
        medians_list.append(a)
        no_median_sum=no_median_sum+b
    file_counter=file_counter+1
    print "no medians were constructed", no_median_sum, "times"    
#print "original sequences list", vectors 
    print "length of vectors list:", len(vectors)

#print "SET of original sequences list", set(vectors) 
    print "length of SET vectors list [=]:", len(set(vectors))
#-----------------------------  
#print "original medians list", medians_list
    print "length of medians list:", len(medians_list)  
  
#print "SET of original medians list", set(medians_list) 
    print "length of SET medians list:", len(set(medians_list))

#print medians_list      #medians, some coinciding!!!   
#print set(medians_list) #unique medians from the list   

    medians_list.extend(vectors)
#print "combining initial sequences with medians:", medians_list
#print "length of combined medians list:", len(medians_list)

#print "SET of all combined:", set(medians_list)
    print "length of combined SET(medians_list) list:", len(set(medians_list))

#-------------------------------------------------------------------------------------
'''


import networkx as nx
G=nx.Graph()
G.add_node("spam")
G.add_edge(1,2)
print(list(G.nodes()))



from networkx import *

G = lollipop_graph(4,6)

pathlengths=[]

print("source vertex {target:length, }")
for v in G.nodes():
    spl=single_source_shortest_path_length(G,v)
    print('%s %s' % (v,spl))
    for p in spl.values():
        pathlengths.append(p)

print('')
print("average shortest path length %s" % (sum(pathlengths)/len(pathlengths)))

dist={}
for p in pathlengths:
    if p in dist:
        dist[p]+=1
    else:
        dist[p]=1

print('')
print("length #paths")
verts=dist.keys()
for d in sorted(verts):
    print('%s %d' % (d,dist[d]))

print("radius: %d" % radius(G))
print("diameter: %d" % diameter(G))
print("eccentricity: %s" % eccentricity(G))
print("center: %s" % center(G))
print("periphery: %s" % periphery(G))
print("density: %s" % density(G))

print("====================================================")
#-------------------------------------------------------------------------------------


import matplotlib.pyplot as plt
#import networkx as nx

G = nx.gnp_random_graph(100,0.02)

degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
#print "Degree sequence", degree_sequence
dmax=max(degree_sequence)

plt.loglog(degree_sequence,'b-',marker='o')
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")

# draw graph in inset
plt.axes([0.45,0.45,0.45,0.45])
Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
pos=nx.spring_layout(Gcc)
plt.axis('off')
nx.draw_networkx_nodes(Gcc,pos,node_size=20)
nx.draw_networkx_edges(Gcc,pos,alpha=0.4)

plt.savefig("degree_histogram.png")
plt.show()
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
'''
#myfilename='AI004_unique_2a_39.fas' # 'BB45_unique_1a_177.fas' has 92 sequences - the biggest; 'BJ30_unique_1a_6' has 5 seqs
#header, vectors, pairs = FASTA (myfilename)

# NOTE: it is important to save time doing SET operation to exclude repetitive sequences!
triplets_list = findsubsets(set(vectors),3) # for each three sequences, find a median; output is a set([ (), ... ,() ])
#print triplets_list, "\nhow many triplets? -- ", len(triplets_list)

medians_list = []
no_median_sum=0
for triplet in triplets_list: ## triplet is a list of 3 strings
    a,b=findMedian(triplet[0],triplet[1],triplet[2])
    #print len(triplets_list)
    #print triplet[1]
    medians_list.append(a)
    no_median_sum=no_median_sum+b
    #print a
print "no medians were constructed", no_median_sum, "times"    
#print "original sequences list", vectors 
print "length of vectors list:", len(vectors)

#print "SET of original sequences list", set(vectors) 
print "length of SET vectors list [=]:", len(set(vectors))
#-----------------------------  
#print "original medians list", medians_list
print "length of medians list:", len(medians_list)  
  
#print "SET of original medians list", set(medians_list) 
print "length of SET medians list:", len(set(medians_list))

#print medians_list      #medians, some coinciding!!!   
#print set(medians_list) #unique medians from the list   

medians_list.extend(vectors)
#print "combining initial sequences with medians:", medians_list
#print "length of combined medians list:", len(medians_list)

#print "SET of all combined:", set(medians_list)
print "length of combined SET(medians_list) list:", len(set(medians_list))


'''
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#print set(mycombined)      
#ymedian = 
#print findMedian(fasta1,fasta2, fasta3)                                                                                                                   
#difference1= difference(fasta1, findMedian(fasta1,fasta2, fasta3))
#difference2= difference(fasta2, findMedian(fasta1,fasta2, fasta3))
#difference3= difference(fasta3, findMedian(fasta1,fasta2, fasta3))

#print "3diffs = ", difference1, difference2, difference3

#print "length of vector is:", len(fasta3)


#for i in range (0,len(vectors[i])):
#    print "step #", i
#    if vectors[i]!=vectors[i+1]:
#        dist+=1
 
'''
dist=0       
vector_len = len(vectors[0])
        
for vector in vectors:
    for i in range (0,vector_len):
        print "length of vector is:", vector_len
        if vectors[i]!=vectors[i+1]:
            dist+=1
            '''
            