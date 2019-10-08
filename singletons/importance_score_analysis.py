"""
Last Edits: Aug 11th, 2019

Example command:

python importance_score_analysis.py --fweight /storage/pandaman/project/AgentBind-GM12878-DanQ-unfixed-rnn-trans/storage/AgentBind-GM12878-DanQ/tmp/FOS+GM12878/seqs_one_hot_c/vis-weights-total/weight.txt --resultdir /storage/mgymrek/agent-bind/singletons --datadir /storage/pandaman/project/AgentBind-GM12878-analysis/data/ --TFname FOS
"""


from os import listdir
from os.path import join
import os
import subprocess
import random
from copy import deepcopy
from scipy.special import comb
from math import pow
from math import log10 as log
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
import numpy as np
#from converter import converter_template as converter

def read_score_map(weight_file):
    '''read predicted scores
    output:
        score_map: which records the predicted scores
    '''
    score_map = {}
    weights_ttl = []
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

                weights_arr = [abs(float(wt)) for wt in (weight_info.strip()).split(";")]
                if strand == "-":
                    weights_arr = weights_arr[::-1]
                weights_ttl += [wt for wt in weights_arr]
                weights_arr_sort = [wt for wt in weights_arr]
                weights_arr_sort.sort(reverse=True)

                # weight_thres_upper and weight_thres_lower are now deprecated
                weight_thres_upper = weights_arr_sort[len(weights_arr_sort)//20] #depr
                weight_thres_lower = weights_arr_sort[-len(weights_arr_sort)//20] # depr

                if chromID in score_map:
                    score_map[chromID].append((seq_start, seq_end, weights_arr, weight_thres_upper, weight_thres_lower, motif_start, motif_end, strand))
                else:
                    score_map[chromID] = [(seq_start, seq_end, weights_arr, weight_thres_upper, weight_thres_lower, motif_start, motif_end, strand)]
            else:
                break

    # sort based on seq_start
    for chromID in score_map:
        score_map[chromID].sort(key=lambda x:x[0])

    # threshold for weight significance
    threshold_control = 200 # 200 -> top 0.005; 100 -> top 0.01; 20 -> top 0.05
    weights_ttl.sort(reverse=True)
    weight_thres_ttl = weights_ttl[len(weights_ttl)//threshold_control]
    print("Weight threshold: %s"%(weight_thres_ttl))
    return score_map, weight_thres_ttl

def read_allele(allele_dir, figure_dir):
    allele_dict = {}
    freq_dict = {}
    for chrom_idx in range(1,23):
        chrom_ID = "chr%d" %(chrom_idx)
        allele_path = "%s/allele_%s.txt" %(allele_dir, chrom_ID)
        allele_dict[chrom_ID] = []
        for line in open(allele_path):
            # 1 10177   A   AC  0.425319
            elems = (line.strip()).split()
            nuc_list = ['A', 'C', 'G', 'T']
            #if (len(elems[2]) >= 4) and (len(elems[2]) <= 9)  and (len(elems[3]) == 1) and \
            if (len(elems[2]) == 1) and (len(elems[3]) == 1) and \
                    (np.all([nt in nuc_list for nt in elems[2]])) and\
                    (np.all([nt in nuc_list for nt in elems[3]])):
                allele_pos = int(elems[1])
                allele_freq = float(elems[4])
                len_ref = len(elems[2])
                # for score analysis
                #if (allele_freq > 0.001) and (allele_freq < 0.1):
                allele_dict[chrom_ID].append((allele_pos, allele_freq, len_ref))

                # frequency distribution stats
                if allele_freq not in freq_dict:
                    freq_dict[allele_freq] = 1
                else:
                    freq_dict[allele_freq] += 1
        allele_dict[chrom_ID].sort(key=lambda x:x[0])

    freq_arr = []
    for allele_freq in freq_dict:
        freq_arr.append((allele_freq, log(freq_dict[allele_freq])))
    freq_arr.sort(key=lambda x:x[0])
#    draw_bar_graph(freq_arr, figure_dir, figure_name="frequency_distribution_all.png")
    return allele_dict

def search_valid_point(score_map, allele_dict, score_thres_ttl):
    valid_points = []
    invalid_points = []
    core_points = []
    for chromID in allele_dict:
        pts_arr = allele_dict[chromID]
        scores_arr = score_map[chromID]
        read_idx = 0
        for (allele_pos, allele_freq, len_ref) in pts_arr:
            # jump to the overlapped read
            while read_idx < (len(scores_arr)-1):
                # break only if the allele pos is smaller than the read end
                if allele_pos < scores_arr[read_idx][1]:
                    break
                read_idx += 1

            # compare
            if allele_pos < scores_arr[read_idx][0]:
                # no read is overlapped with this SNP
                continue
            elif allele_pos+len_ref <= scores_arr[read_idx][1]:
                allele_relative_pos = allele_pos - scores_arr[read_idx][0]
                importance_score = scores_arr[read_idx][2][allele_relative_pos]
                print("debug %s:%s:%s"%(allele_pos, allele_freq, importance_score))
                if (scores_arr[read_idx][5] <= allele_pos) and (allele_pos+len_ref <= scores_arr[read_idx][6]):
                    # skip the core motifs
                    core_points.append((allele_freq, importance_score))
                else:
                    if importance_score >= score_thres_ttl:
                        valid_points.append((allele_freq, importance_score))
                    elif importance_score < score_thres_ttl:
                        invalid_points.append((allele_freq, importance_score))
            else:
                # no read is left; have reached the end of this chromosome
                break
  
    return valid_points, invalid_points, core_points

def cmp_valid_vs_invalid(valid_points, invalid_points, figure_dir):
    valid_freq = [allele_freq for (allele_freq, importance_score) in valid_points]
    invalid_freq = [allele_freq for (allele_freq, importance_score) in invalid_points]
    x_ticks_arr = ['valid', 'other']
    draw_boxplot([valid_freq, invalid_freq], x_ticks_arr,\
                        figure_dir, "valid-invalid-freq-distribution")
    return

def draw_point_boxplot(points, figure_dir, figure_name):
    point_dict = {}
    xticks_arr = []
    for (allele_freq, importance_score) in points:
        if allele_freq not in point_dict:
            point_dict[allele_freq] = [importance_score]
            xticks_arr.append(allele_freq)
        else:
            point_dict[allele_freq].append(importance_score)

    xticks_arr.sort()
    point_arr = []
    for xtick in xticks_arr:
        point_arr.append(point_dict[xtick])

    draw_boxplot(point_arr, xticks_arr, figure_dir, figure_name)
    return

def draw_boxplot(data, x_ticks_arr, figure_dir, figure_name):
    plt.figure()
    plt.boxplot(data, showfliers=False)
    plt.xticks([idx+1 for idx in range(len(x_ticks_arr))], x_ticks_arr)
    plt.ylabel('importance score')
    plt.title(figure_name)
    plt.savefig('%s/%s.png' %(figure_dir, figure_name))
    plt.close()
    return

def draw_scatter(points, figure_dir, figure_name):
    plt.figure(figsize=(11,11))
    plt.scatter([pt[0] for pt in points],\
                    [pt[1] for pt in points],\
                    marker=".")
    plt.title(figure_name)
    plt.xlabel('frequency')
    plt.ylabel('importance score')
    plt.savefig('%s/%s.png' %(figure_dir, figure_name))
    plt.close()
    return

def draw_bar_graph(data_list, figure_dir, figure_name):
    x = []
    y = []
    for (allele_freq, n_samples) in data_list:
        x.append(allele_freq)
        y.append(n_samples)
    fig, ax = plt.subplots()
    index = np.arange(0.2*len(x), step=0.2)
    bar_width = 0.2
    opacity = 0.4
    
    rects1 = ax.bar(index, y, bar_width,
                        alpha=opacity, color='r')
    #rects2 = ax.bar(index + bar_width, data_list['ab-unfixed-1'], bar_width,
    #                    alpha=opacity, color='g',
    #                    label='AgentBind-unfixed-1')
    #rects3 = ax.bar(index + 2*bar_width, data_list['ab-unfixed-2'], bar_width,
    #                    alpha=opacity, color='g',
    #                    label='AgentBind-unfixed-2')
    ax.set_xlabel('frequency')
    ax.set_ylabel('number of samples')
    ax.set_title(figure_name)
    ax.set_xticks(index + bar_width/2.0)
    ax.set_xticklabels(x, rotation=90)
    fig.tight_layout()
    plt.savefig('%s/%s.png' %(figure_dir, figure_name))
    plt.close()
    return

def draw_histogram(points, figure_dir, figure_name):
    freq_arr = [allele_freq for (allele_freq, importance_score) in points]
    plt.figure(figsize=(11,11))
    plt.hist(freq_arr, bins = 50)
    plt.title(figure_name)
    plt.xlabel('frequency')
    plt.ylabel('#samples')
    plt.savefig('%s/%s.png' %(figure_dir, figure_name))
    plt.close()
    return

def calculate_singleton_rate(points):
    n_ttl = len(points)
    n_singletons = 0
    for (allele_freq, importance_score) in points:
        if allele_freq < 0.0002:
            n_singletons += 1
    sgt_rate = float(n_singletons) / float(n_ttl)
    return (sgt_rate, n_singletons)

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--datadir', dest='dir_data', default=None,
            help='The directory where all input and dependent data are stored [Default: %default]')
    parser.add_option('--fweight', dest='weight_file', default=None,
            help='The file where stores weights of positions of test samples [Default: %default]')
    parser.add_option('--resultdir', dest='result_dir', default=None,
            help='The file where to save the residue sequences [Default: %default]')
    parser.add_option('--TFname', dest='TF_name', default=None,
            help='The name of the core TF [Default: %default]')
    (options, args) = parser.parse_args()

    genome_limit_path = "%s/genomes/hg19/hg19.fa.fai" %(options.dir_data)
    ref_genomes_dir = "%s/genomes/hg19/" %(options.dir_data)
    allele_dir = "%s/1000genomes/" %(options.dir_data)
    result_dir_for_TF = "%s/%s/1000genomes/" %(options.result_dir, options.TF_name)
    if not os.path.exists(result_dir_for_TF):
        os.makedirs(result_dir_for_TF)

    #####
    # Calculate weight scores for each TF
    #####
    # load data
    print("read score map")
    score_map, weight_thres_ttl = read_score_map(options.weight_file)
    print("get allele dict")
    allele_dict = read_allele(allele_dir, result_dir_for_TF)
    print("search valid point")
    valid_points, invalid_points, core_motif_points = search_valid_point(score_map, allele_dict, weight_thres_ttl)
    print("get singleton rate valid")
    sgt_rate_valid, n_singletons_valid = calculate_singleton_rate(valid_points)
    print("get singleton rate invalid")
    sgt_rate_invalid, n_singletons_invalid = calculate_singleton_rate(invalid_points)
    print("get singleton rate core")
    sgt_rate_core, n_singletons_core = calculate_singleton_rate(core_motif_points)
    
    #draw_scatter(valid_points, result_dir_for_TF, figure_name="freq_vs_IS")
    #draw_point_boxplot(valid_points, result_dir_for_TF, figure_name="freq_vs_IS_valid_points_boxplot")
    #draw_point_boxplot(invalid_points, result_dir_for_TF, figure_name="freq_vs_IS_invalid_points_boxplot")
    #cmp_valid_vs_invalid(valid_points, invalid_points, result_dir_for_TF)
    #draw_histogram(valid_points, result_dir_for_TF, figure_name="freq_hist_valid.png")
    #draw_histogram(invalid_points, result_dir_for_TF, figure_name="freq_hist_invalid.png")

    data_size_file = "%s/data_size.txt" %(result_dir_for_TF)
    with open(data_size_file, 'w') as ofile:
        line = "size of valid points: %d, number of singletons: %d, singleton rate: %f\n" %(len(valid_points), n_singletons_valid, sgt_rate_valid)
        ofile.write(line)
        line = "size of other points: %d, number of singletons: %d, singleton rate: %f\n" %(len(invalid_points), n_singletons_invalid, sgt_rate_invalid)
        ofile.write(line)
        line = "size of core motif points: %d, number of singletons: %d, singleton rate: %f\n" %(len(core_motif_points), n_singletons_core, sgt_rate_core)
        ofile.write(line)
    return

if __name__ == "__main__":
    main()

##### END OF FILE #####################################
