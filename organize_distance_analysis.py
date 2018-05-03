import numpy as np 
import os
import sys
import pdb
import math
import random

# First create dictionary list of the significant variant gene pairs where each key is of form $variantID"_"$geneID
def extract_significant_variant_gene_pairs(file_name):
    f = open(file_name)
    dicti = {}  # dictionary to keep variant gene pairs
    head_count = 0  # Used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        #  A significant variant-gene pairs
        rs_id = data[2]
        ensamble_id = data[5]
        pvalue = float(data[-3])
        # Simple check
        if rs_id + '_' + ensamble_id in dicti:
            print('fundamental assumption error')
            pdb.set_trace()
        # Add variant gene pair to dictionary
        dicti[rs_id + '_' + ensamble_id] = 1
    return dicti


def extract_significant_variant_gene_pairs_from_time_ind(file_name):
    f = open(file_name)
    dicti = {}  # dictionary to keep variant gene pairs
    head_count = 0  # Used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        #  A significant variant-gene pairs
        rs_id = data[3]
        ensamble_id = data[1]
        # Simple check
        if rs_id + '_' + ensamble_id in dicti:
            print('fundamental assumption error')
            pdb.set_trace()
        # Add variant gene pair to dictionary
        dicti[rs_id + '_' + ensamble_id] = 1
    return dicti


# Input list of variant gene pairs, return list of distance to tss's 
def get_distance_to_tss_for_variant_gene_pair_list(sig_variant_gene_pairs, time_step_independent_file):
    # Array to keep track of distances
    distances = []
    # Used to skip header
    head_count = 0
    # Open time_step_independent file that contains distance to tss knowledge for each variant gene pair
    f = open(time_step_independent_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Extract info
        ensamble_id = data[1]
        rs_id = data[3]
        pair_name = rs_id + '_' + ensamble_id
        # Throw out variant gene pairs that are not in input dictionary
        if pair_name not in sig_variant_gene_pairs:
            continue
        # compute distance to tss
        gene_tss_pos = float(data[2])
        rs_pos = float(data[4])
        disty = abs(rs_pos - gene_tss_pos)
        distances.append(disty)

    # Simple check
    if len(distances) != len(sig_variant_gene_pairs):
        print('Fundamental assumption errro')
        pdb.set_trace()

    return np.asarray(distances)

# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))


# Make object that takes as input a MAF and returns a list of variant-gene pairs in this MAF
def make_background_pairs_object(time_step_independent_file):
    maf_bin_size = .05
    ####################
    # Initialize object
    ####################
    background_qtls = []
    # number of bins needed for maf and distance
    num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))
    maf_mapping = {}
    # Add each possible bin
    for maf_bin in range(num_maf_bins):
        background_qtls.append([])

    # Fill in object
    f = open(time_step_independent_file)
    head_count = 0
    background = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        maf = float(data[5])
        ensamble_id = data[1]
        rs_id = data[3]
        test_name = rs_id + '_' + ensamble_id
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        background_qtls[maf_bin].append(test_name)
    return background_qtls

# Extract dictionary list of len(sig_variant_gene_pairs) that is matched for MAF
def sample_background_variant_gene_pairs(time_step_independent_file, background_pairs_object, sig_variant_gene_pairs):
    maf_bin_size = .05
    # output dictionary 
    background_variant_gene_pairs = {}
    # open file with MAF info
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        # extract relevent info
        maf = float(data[5])
        ensamble_id = data[1]
        rs_id = data[3]
        test_name = rs_id + '_' + ensamble_id
        # ignore lines (variant-gene pairs) that are not in sig_variant_gene_pairs
        if test_name not in sig_variant_gene_pairs:
            continue
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        # Randomly select a variant gene pair
        converged = False
        while converged == False:
            randomly_selected_pair = random.choice(background_pairs_object[maf_bin])
            if randomly_selected_pair not in background_variant_gene_pairs:
                background_variant_gene_pairs[randomly_selected_pair] = 1
                converged = True 
    f.close()
    # simple check
    if len(background_variant_gene_pairs) != len(sig_variant_gene_pairs):
        print('FUNDAMENTAL ASSUMPTIONER ERROR')
        pdb.set_trace()
    return background_variant_gene_pairs




# Return list of length num_permutations where each element of the list is a vector of length real_distances.
def get_distance_to_tss_for_background_variant_gene_pair_list(sig_variant_gene_pairs, time_step_independent_file, num_permutations):
    # Make object that takes as input a MAF and returns a list of variant-gene pairs in this MAF
    background_pairs_object = make_background_pairs_object(time_step_independent_file)
    distances = []
    # Loop through all permutations
    for perm_number in range(num_permutations):
        # Extract dictionary list of len(sig_variant_gene_pairs) that is matched for MAF
        background_variant_gene_pairs = sample_background_variant_gene_pairs(time_step_independent_file, background_pairs_object, sig_variant_gene_pairs)
        # Input list of variant gene pairs, return list of distance to tss's 
        perm_distances = get_distance_to_tss_for_variant_gene_pair_list(background_variant_gene_pairs, time_step_independent_file)

        distances.append(perm_distances)
    return distances    

significant_variant_gene_pairs_file = sys.argv[1]
time_step_independent_stem = sys.argv[2]
distance_results_file = sys.argv[3]
num_permutations = int(sys.argv[4])


# Results file for time step 0 (could have been any time step)
# Only use this to extract maf and distance to tss info
time_step_independent_file = time_step_independent_stem + '0_eqtl_results.txt'


# First create dictionary list of the significant variant gene pairs where each key is of form $variantID"_"$geneID
sig_variant_gene_pairs = extract_significant_variant_gene_pairs(significant_variant_gene_pairs_file)
#sig_variant_gene_pairs = extract_significant_variant_gene_pairs_from_time_ind(significant_variant_gene_pairs_file)

# Input list of variant gene pairs, return list of distance to tss's 
real_distances = get_distance_to_tss_for_variant_gene_pair_list(sig_variant_gene_pairs, time_step_independent_file)
# Return list of length num_permutations where each element of the list is a vector of length real_distances.
perm_distances = get_distance_to_tss_for_background_variant_gene_pair_list(sig_variant_gene_pairs, time_step_independent_file, num_permutations)

# Save results to output file
output_mat = np.transpose(np.vstack((real_distances,np.asmatrix(perm_distances))))
np.savetxt(distance_results_file, output_mat, delimiter='\t',fmt="%s")
