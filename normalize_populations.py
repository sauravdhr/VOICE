"""
Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
Created: 11/28/2016
"""

from Bio import SeqIO
from Bio.Seq import Seq
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import (fcluster, linkage)
from collections import defaultdict
import numpy as np
import argparse
import glob
import sys
import os
import re

nucls='ACTG'
method='average'
criterion='maxclust'
count_re_pattern='' #'[0-9]*$'

def distance(a, b):
    if len(a) != len(b):
        return -1
    l = len(a)
    return sum(map(lambda i: a[i]!=b[i], range(l)))

def to_consensus(fseqs):
    consen = fseqs[0]
    profile = [defaultdict(float) for _ in  range(len(consen))]
    for s in fseqs:
        cnt = get_count(s, count_re_pattern)
        for i in range(len(consen)):
            profile[i][s.seq[i]] += cnt

    constr = ''.join(map(lambda x: max(x, key=x.get), profile))
    cc = Seq(constr, consen.seq.alphabet)
    consen.seq = cc
    return consen

def process_file(fasta, k):
    vectors = list(map(lambda x: list(map(lambda y: nucls.index(y), x.seq)), fasta))
    vs = np.array(vectors)

    Y = pdist(vs, distance)
    Z = linkage(Y, method=method)

    labels = fcluster(Z, k, criterion=criterion)

    bag = defaultdict(lambda: [])

    for i in range(len(labels)):
        bag[labels[i]].append(fasta[i])

    return [to_consensus(bag[key]) for key in bag]

def get_count(fas, patt):
    if patt == "": return 1.0
    m = re.search(patt, fas.description)
    if m is not None:
        return float(m.group(0))
    return 1.0

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input_dir', type=str, required=True,
                        help="Path to a directory with fasta files of samples to be normalized.")
    parser.add_argument("-o", dest='out_dir', type=str, default='./out',
                        help="Path to an output directory relative to input dir. "
                             "By default it creates directory named \'out\' in the input directory")
    parser.add_argument("-c", dest='sample_code', type=str, default='',
                        help="Sample code (first 2 letters case sensitive in filename)")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        sys.stderr.write("Input directory path does not exist or unavailable!")
        sys.exit(-1)

    output = "%s/%s" % (args.input_dir, args.out_dir)

    if os.path.isdir(output):
        os.system('rm -r %s' % output)
    if os.path.exists(output):
        os.remove(output)
    os.makedirs(output)

    fastas = {os.path.basename(f):list(SeqIO.parse(f, 'fasta'))
              for f in glob.glob("%s/%s*.fas" % (args.input_dir, args.sample_code))}

    k = min(map(len, fastas.values()))

    for fname in fastas:
        fasta = fastas[fname]
        res = process_file(fasta, k) if len(fasta) > k+1 else fasta
        SeqIO.write(res, open("%s/%s" % (output, fname), 'w+'), 'fasta')