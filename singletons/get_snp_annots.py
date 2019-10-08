"""
Last Edits: Aug 11th, 2019

Example command:

python get_snp_annots.py --fweight /storage/pandaman/project/AgentBind-GM12878-DanQ-unfixed-rnn-trans/storage/AgentBind-GM12878-DanQ/tmp/CTCF+GM12878/seqs_one_hot_c/vis-weights-total/weight.txt --resultdir /storage/mgymrek/agent-bind/singletons --TFname CTCF
"""

import warnings
warnings.filterwarnings("ignore")

from os import listdir
from os.path import join
import os
import subprocess
import random
from copy import deepcopy
from scipy.special import comb
import scipy.stats
from math import pow
from math import log10 as log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
import numpy as np
import sys

def print_score_map(weight_file, outfile, igvfile):
    '''read predicted scores
    output:
        file with: chrom, pos, score, norm score, and other annotations
    '''
    outf = open(outfile, "w")
    igvf = open(igvfile, "w")
    outf.write("\t".join(["chrom","start","score.raw","score.SNR","score.rank","is.core"])+"\n")
    igvf.write("track type=bedGraph\n")
    with open(weight_file) as f:
        while True:
            position_info = f.readline()
            weight_info = f.readline()

            if (position_info != "") and (weight_info != ""):
                chromID, motif_start, motif_end, seq_start, seq_end, strand =\
                                (position_info.strip()).split(";")
                motif_start = int(motif_start)
                motif_end = int(motif_end)
                seq_start = int(seq_start)
                seq_end = int(seq_end)

                weights_arr = [max([float(wt), 0]) for wt in (weight_info.strip()).split(";")] # set neg numbers to 0
                if strand == "-": # TODO is this correct?
                    weights_arr = weights_arr[::-1]
                ranks = scipy.stats.rankdata(weights_arr)
                for i in range(seq_end-seq_start):
                    is_core = int(seq_start+i>= motif_start and seq_start+i<= motif_end)
                    outf.write("\t".join(map(str, [chromID, seq_start+i, weights_arr[i], weights_arr[i]/np.mean(weights_arr), ranks[i], is_core]))+"\n")
                    igvf.write("\t".join(map(str, [chromID, seq_start+i, seq_start+i+1, weights_arr[i]]))+"\n")
            else:
                break
    outf.close()

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--fweight', dest='weight_file', default=None,
            help='The file which stores weights of positions of test samples [Default: %default]')
    parser.add_option('--resultdir', dest='result_dir', default=None,
            help='The file where to save the residue sequences [Default: %default]')
    parser.add_option('--TFname', dest='TF_name', default=None,
            help='The name of the core TF [Default: %default]')
    (options, args) = parser.parse_args()

    result_dir_for_TF = "%s/%s" %(options.result_dir, options.TF_name)
    if not os.path.exists(result_dir_for_TF):
        os.makedirs(result_dir_for_TF)

    #####
    # Calculate weight scores for each TF
    #####
    print_score_map(options.weight_file, \
                    os.path.join(result_dir_for_TF, "scores.tab"),
                    os.path.join(result_dir_for_TF, "%s-scores.bedGraph"%options.TF_name))

if __name__ == "__main__":
    main()

##### END OF FILE #####################################
