"""
Author: Alexander Artyomenko <aartyomenko@cs.gsu.edu>
Created: 12/2/2016
"""

import argparse
import glob
import sys
import os

nucls='ACTG'
unrelated_code='XX'


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input_dir', type=str, required=True,
                        help="Path to a directory with fasta files of samples to be normalized.")
    parser.add_argument("-k", dest='k', type=int, default=6,
                        help="Amount of folds (k) for k-fold cross validation. Default : 6")
    parser.add_argument("-o", dest='out_dir', type=str, default='./out',
                        help="Path to an output directory relative to input dir. "
                             "By default it creates directory named \'out\' in the input directory")

    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        sys.stderr.write("Input directory path does not exist or unavailable!")
        sys.exit(-1)

    output = "%s/%s" % (args.input_dir, args.out_dir)

    k = args.k

    if os.path.isdir(output):
        os.system('rm -r %s' % output)
    if os.path.exists(output):
        os.remove(output)

    paths = [f for f in glob.glob("%s/*" % args.input_dir)]

    os.makedirs(output)

    for i in range(k):
        tmp_out = "%s/%d" % (output, i)
        os.makedirs(tmp_out)
        tmp_out_unr = "%s/%s" % (tmp_out, unrelated_code)
        os.makedirs(tmp_out_unr)
        for j, p in enumerate(paths):
            if not unrelated_code in p:
                if j % k != i:
                    os.system("cp -r %s %s" % (p, tmp_out))
            else:
                unrelated_samples = [f for f in glob.glob("%s/*.fas" % p)]
                for sj, sp in enumerate(unrelated_samples):
                    if sj % k != i:
                        os.system("cp %s %s" % (sp, tmp_out_unr))
