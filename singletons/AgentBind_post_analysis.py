'''
    -> AgentBind <-
    
'''

from optparse import OptionParser
import os
import subprocess
import time

def main():
    ############
    # user-defined options
    ############
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--motif', dest='f_core_motif',
            default="/storage/pandaman/project/AgentBind-GM12878-analysis/data/table_matrix/table_core_motifs.txt",
            help='The file address where you save your motifs of interst [Default: %default]')
    parser.add_option('--scope', dest='len_seq', default=1000, type='int',
            help='The scope of contexts to be examined for each core motif [Default: %default]')
    parser.add_option('--workdir', dest='dir_work',
            default='/storage/pandaman/project/AgentBind-GM12878-analysis/tmp/',
            help='The directory where all temporary and final results are saved [Default: %default]')
    parser.add_option('--datadir', dest='dir_data',
            default='/storage/pandaman/project/AgentBind-GM12878-analysis/data/',
            help='The directory where all input and dependent data are stored [Default: %default]')
    parser.add_option('--resultdir', dest='dir_result',
            default='/storage/pandaman/project/AgentBind-GM12878-analysis/results-singlton-0-005-global-threshold/',
            help='The directory where all input and dependent data are stored [Default: %default]')
    (options, args) = parser.parse_args()

    ############
    # analysis of genomic contexts around motifs
    #############
    if options.f_core_motif != None and\
            os.path.isfile(options.f_core_motif):
        TF_list = []
        for line in open(options.f_core_motif):
            (TF_name, TF_ID, TF_chipseq) = line.strip().split(";")
            TF_list.append((TF_name, TF_ID, TF_chipseq))
    else:
        exit("Cannot find the given motif file address: %s" %(options.f_core_motif))

    ########
    # step 1: identify all matched locations across the genome
    #         for each motif
    
    ref_genomes = "%s/genomes/hg19/hg19.fa" %(options.dir_data)
    bgfile = "%s/genomes/hg19/hg19.fna.bfile" %(options.dir_data)
    dir_jaspar_motif = "%s/Jaspar-motifs/JASPAR2018_CORE_vertebrates_redundant_pfms_meme/"\
                            %(options.dir_data)
    dir_fimo_out = "%s/fimo_out/" %(options.dir_work)
    if not os.path.exists(dir_fimo_out):
        os.makedirs(dir_fimo_out)

    ######
    for (TF_name, TF_ID, TF_chipseq) in TF_list:
        start_time = time.time()

        dir_work_TF = "%s/%s/" %(options.dir_work, TF_name)
        f_fimo = os.path.join(dir_fimo_out, TF_ID)
        ################
        ## Model training
        ################

        for suffix in ['c']: #'b' stand for block, 'c' stand for core
            #######
            # step 2: prepare data for model training
            dir_seqs = "%s/seqs_one_hot_%s/" %(dir_work_TF, suffix)
            f_chipseq = "%s/TFBS_ENCODE/%s" %(options.dir_data, TF_chipseq)
            dir_result_with_suffix = "%s/%s/" %(options.dir_result, suffix)
            f_weight = "%s/vis-weights-total/weight.txt" %(dir_seqs)

            # singleton analysis
            IS_analysis_cmd = "python importance_score_analysis.py "\
                    " --datadir %s --fweight %s --resultdir %s --TFname %s " %(options.dir_data,\
                            f_weight, dir_result_with_suffix, TF_name)
            subprocess.check_call(IS_analysis_cmd, shell=True)

            # kmer analysis
            #f_neg_coord = "%s/vis-weights-total/coordinates_neg.txt" %(dir_seqs)
            #f_vis_samples = "%s/vis-samples/data.txt" %(dir_seqs)
            #dir_analysis_tmp_data = "%s/data_analysis/" %(dir_work_TF)
            #kmer_analysis_cmd = "python ./data_analysis/data_analysis.py "\
            #        " --datadir %s --storage %s --fweight %s --fseq %s --negcoord %s "\
            #        " --fimodir %s --resultdir %s "\
            #        " --TFname %s " %(options.dir_data, dir_analysis_tmp_data, f_weight, f_vis_samples, f_neg_coord,
            #                dir_fimo_out, dir_result_with_suffix, TF_name)
            #subprocess.check_call(kmer_analysis_cmd, shell=True)

            # GC content analysis
            #f_pos_coord = "%s/vis-weights-total/coordinates_pos.txt" %(dir_seqs)
            #GC_analysis_cmd = "python ./data_analysis/CG_content_ratio.py "\
            #        " --datadir %s --negcoord %s --poscoord %s "\
            #        " --resultdir %s --TFname %s " %(options.dir_data, f_neg_coord, f_pos_coord,
            #                dir_result_with_suffix, TF_name)
            #subprocess.check_call(GC_analysis_cmd, shell=True)

            
        end_time = time.time()
        print ("Analysis of %s used %f seconds" %(TF_name, end_time-start_time))
        start_time = end_time

        #################
if __name__ == "__main__":
    main()
