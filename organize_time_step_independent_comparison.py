import numpy as np
import os
import sys
import pdb
import math
import random

def get_variant_gene_pairs(dynamic_qtl_all_hits_file):
    f = open(dynamic_qtl_all_hits_file)
    dicti = {}
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        pvalue = float(data[-3])
        rs_id = data[2]
        ensamble_id = data[5]
        dicti[ensamble_id + '_' + rs_id] = np.zeros(17)
        dicti[ensamble_id + '_' + rs_id][-1] = pvalue
    return dicti

def fill_in_dictionary_with_time_step_independent_file(variant_gene_pairs, time_step_independent_file, time_step):
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        rs_id = data[3]
        pvalue = float(data[-1])
        if ensamble_id + '_' + rs_id in variant_gene_pairs:
            variant_gene_pairs[ensamble_id + '_' + rs_id][time_step] = pvalue
    return variant_gene_pairs

def make_all_tests_comparison_file(dynamic_qtl_all_hits_file, time_step_independent_stem, all_tests_comparison_file):
    variant_gene_pairs = get_variant_gene_pairs(dynamic_qtl_all_hits_file)
    names = []
    names.append('test_name')
    for time_step in range(16):
        print(time_step)
        time_step_independent_file =  time_step_independent_stem + str(time_step) + '_eqtl_results.txt'
        variant_gene_pairs = fill_in_dictionary_with_time_step_independent_file(variant_gene_pairs, time_step_independent_file, time_step)
        names.append('time_step_' + str(time_step))
    names.append('dynamic')
    t = open(all_tests_comparison_file, 'w')
    t.write('\t'.join(names) + '\n')
    for test_name in variant_gene_pairs.keys():
        t.write(test_name + '\t')
        pvalues = variant_gene_pairs[test_name].astype(str)
        t.write('\t'.join(pvalues) + '\n')
    t.close()


# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def get_background_variant_gene_pairs(input_dicti, input_file):
    distance_bin_size = 10000
    maf_bin_size = .05
    eqtl_distance = 50000
    ####################
    # Initialize object
    ####################
    background_qtls = []
    # number of bins needed for maf and distance
    num_distance_bins = int(math.ceil(eqtl_distance/distance_bin_size + 1))
    num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))
    maf_mapping = {}
    distance_mapping = {}
    # Add each possible bin
    for distance_bin in range(num_distance_bins):
        background_qtls.append([])
        for maf_bin in range(num_maf_bins):
            background_qtls[distance_bin].append([])
    output_dicti = {}
    f = open(input_file)
    head_count = 0
    background = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        rs_id = data[3]
        test_name = ensamble_id + '_' + rs_id
        distance = abs(float(data[2]) - float(data[4]))
        maf = float(data[5])
        if test_name in input_dicti:
            maf_mapping[test_name] = maf
            distance_mapping[test_name] = distance
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        background_qtls[distance_bin][maf_bin].append(test_name)
    f.close()

    output_dicti = {}
    for test_name in input_dicti.keys():
        arr = input_dicti[test_name]
        distance = distance_mapping[test_name]
        maf = maf_mapping[test_name]
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)

        tests_to_choose_from = background_qtls[distance_bin][maf_bin]
        if len(tests_to_choose_from) < 10:
            print('short supply')
            pdb.set_trace()
        working = False
        while working == False:
            new_test = random.choice(tests_to_choose_from)
            if new_test not in output_dicti:
                output_dicti[new_test] = arr
                working = True
    return output_dicti
        

def make_background_tests_comparison_file(dynamic_qtl_all_hits_file, time_step_independent_stem, all_tests_comparison_file):
    variant_gene_pairs_real = get_variant_gene_pairs(dynamic_qtl_all_hits_file)
    variant_gene_pairs = get_background_variant_gene_pairs(variant_gene_pairs_real, time_step_independent_stem + '0_eqtl_results.txt' )
    names = []
    names.append('test_name')
    for time_step in range(16):
        print(time_step)
        time_step_independent_file = time_step_independent_stem + str(time_step) + '_eqtl_results.txt'
        variant_gene_pairs = fill_in_dictionary_with_time_step_independent_file(variant_gene_pairs, time_step_independent_file, time_step)
        names.append('time_step_' + str(time_step))
    names.append('dynamic')
    t = open(all_tests_comparison_file, 'w')
    t.write('\t'.join(names) + '\n')
    for test_name in variant_gene_pairs.keys():
        t.write(test_name + '\t')
        pvalues = variant_gene_pairs[test_name].astype(str)
        t.write('\t'.join(pvalues) + '\n')
    t.close()




# Output files
dynamic_egenes_comparison_file = sys.argv[1]
dynamic_egenes_background_comparison_file = sys.argv[2]
# Input files
time_step_independent_stem = sys.argv[3]
egenes_file = sys.argv[4]



make_all_tests_comparison_file(egenes_file, time_step_independent_stem, dynamic_egenes_comparison_file)


make_background_tests_comparison_file(egenes_file, time_step_independent_stem, dynamic_egenes_background_comparison_file)
