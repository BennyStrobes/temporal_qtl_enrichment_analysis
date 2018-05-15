import numpy as np 
import os
import sys
import pdb
import math
import gzip
import random
import scipy.stats


# First create dictionary list of the significant variant gene pairs where each key is of form $variantID"_"$geneID
def extract_significant_variant_gene_pairs(file_name, variant_gene_pair_time_step_info, hits_version):
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
        if hits_version == 'all_hits':
            dicti[rs_id + '_' + ensamble_id] = 1
        elif hits_version == 'early_time_step_hits' and variant_gene_pair_time_step_info[rs_id + '_' + ensamble_id] == 'early_time_step_hits':
            dicti[rs_id + '_' + ensamble_id] = 1
        elif hits_version == 'late_time_step_hits' and variant_gene_pair_time_step_info[rs_id + '_' + ensamble_id] == 'late_time_step_hits':
            dicti[rs_id + '_' + ensamble_id] = 1
        elif hits_version == 'change_in_sign_hits' and variant_gene_pair_time_step_info[rs_id + '_' + ensamble_id] == 'change_in_sign_hits':
            dicti[rs_id + '_' + ensamble_id] = 1
    return dicti



# Create Mapping from variant-gene pair to quartuple (chrom_num, variant_position, distToTss, MAF)
def extract_variant_gene_pair_info(time_step_independent_file):
    head_count = 0
    # Open time_step_independent file that contains distance to tss knowledge for each variant gene pair
    f = open(time_step_independent_file)
    # Dictionary for each variant-gene pair
    variant_gene_pair_info = {}
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
        # compute distance to tss
        gene_tss_pos = float(data[2])
        rs_pos = float(data[4])
        disty = abs(rs_pos - gene_tss_pos)
        maf = float(data[5])
        # Add relevent info to small dictionary
        mini_dicti = {}
        mini_dicti['chrom_num'] = int(data[0])
        mini_dicti['dist_to_tss'] = disty
        mini_dicti['maf'] = maf
        mini_dicti['variant_position'] = rs_pos
        # Add variant-gene pair to dictionary
        if pair_name in variant_gene_pair_info:
            print('assumption errror')
            pdb.set_trace()
        variant_gene_pair_info[pair_name] = mini_dicti
    return variant_gene_pair_info


# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def make_background_object(time_step_independent_file):
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
    # Add each possible bin
    for distance_bin in range(num_distance_bins):
        background_qtls.append([])
        for maf_bin in range(num_maf_bins):
            background_qtls[distance_bin].append([])
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
        test_name = rs_id + '_' + ensamble_id
        distance = abs(float(data[2]) - float(data[4]))
        maf = float(data[5])
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        background_qtls[distance_bin][maf_bin].append(test_name)
    f.close()
    return background_qtls

def sample_background_variant_gene_pairs(background_object, variant_gene_pair_info, sig_variant_gene_pairs):
    distance_bin_size = 10000
    maf_bin_size = .05
    eqtl_distance = 50000
    background_variant_gene_pairs = {}
    for test_name in sig_variant_gene_pairs.keys():
        dist_to_tss = variant_gene_pair_info[test_name]['dist_to_tss']
        maf = variant_gene_pair_info[test_name]['maf']
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(dist_to_tss, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        converged = False
        while converged == False:
            if len(background_object[distance_bin][maf_bin]) < 10:
                print('small backgroudn: should investigate')
                pdb.set_trace()
            randomly_selected_pair = random.choice(background_object[distance_bin][maf_bin])
            if randomly_selected_pair not in background_variant_gene_pairs:
                background_variant_gene_pairs[randomly_selected_pair] = 1
                converged = True
    if len(background_variant_gene_pairs) != len(sig_variant_gene_pairs):
        print('ASSUMPTION EROROR')
    return background_variant_gene_pairs


# Return list of length num_permutations where each element of the dictionary list of variant-gene pairs (of len(sig_variant_gene_pairs)) matched for dist_to_tss and maf
def extract_perm_variant_gene_pairs(sig_variant_gene_pairs, time_step_independent_file, num_permutations, variant_gene_pair_info):
    wrapper_list = []
    background_object = make_background_object(time_step_independent_file)
    for perm_num in range(num_permutations):
        background_variant_gene_pairs = sample_background_variant_gene_pairs(background_object, variant_gene_pair_info, sig_variant_gene_pairs)
        wrapper_list.append(background_variant_gene_pairs)
    return wrapper_list

# Extract list of cell line ids used for this cell_line_version
def get_cell_line_ids(cell_line_version, chrom_hmm_input_dir):
    cell_line_ids = []
    if cell_line_version == 'heart_cell_lines':
        cell_line_ids.append('E095')
        #cell_line_ids.append('E104')
        cell_line_ids.append('E105')
        cell_line_ids.append('E083')
    elif cell_line_version == 'ipsc_cell_lines':
        cell_line_ids.append('E018')
        cell_line_ids.append('E019')
        cell_line_ids.append('E020')
        cell_line_ids.append('E021')
        cell_line_ids.append('E022')

    elif cell_line_version == 'all_cell_lines':
        for file_name in os.listdir(chrom_hmm_input_dir):
            if file_name.endswith('mnemonics.bed.gz') == False:
                continue
            cell_line_id = file_name.split('_')[0]
            cell_line_ids.append(cell_line_id)
    return np.asarray(cell_line_ids)

def get_valid_markers(marker_type):
    valid_markers = {}
    if marker_type == 'promotor':
        valid_markers['1_TssA'] = 1
        valid_markers['2_TssAFlnk'] = 1
        valid_markers['10_TssBiv'] = 1
        valid_markers['11_BivFlnk'] = 1
    elif marker_type == 'enhancer':
        valid_markers['7_Enh'] = 1
        valid_markers['6_EnhG'] = 1
        valid_markers['12_EnhBiv'] = 1
        valid_markers['11_BivFlnk'] = 1
    return valid_markers


# Make binary array length of a chromosome. If array == 0, no marker there. If array == 1, there is a marker there
def make_binary_chromosome(chrom_num, chrom_hmm_input_dir, cell_line_ids, marker_type):
    chrom_num_string = 'chr' + str(chrom_num)
    chromosome = np.zeros(259250621)
    valid_markers = get_valid_markers(marker_type)
    for cell_line_id in cell_line_ids:
        chrom_hmm_file = chrom_hmm_input_dir + cell_line_id + '_15_coreMarks_mnemonics.bed.gz'
        f = gzip.open(chrom_hmm_file)
        for line in f:
            line = line.rstrip()
            data = line.split()
            chromer = data[0]
            # Only consider marks on this chromsome
            if chromer != chrom_num_string:
                continue
            line_marker = data[3]
            if line_marker in valid_markers:
                start = int(data[1])
                end = int(data[2])
                if start > end:
                    print('assumption errororo')
                chromosome[start:end] = np.ones(end-start)
    return chromosome


# Count number of variants that overlap a marker on this chromosome
def count_variant_overlap(chrom_num, chromosome, sig_variant_gene_pairs, variant_gene_pair_info):
    county = 0
    for variant_gene in sig_variant_gene_pairs.keys():
        variant_gene_pair_dict = variant_gene_pair_info[variant_gene]
        if variant_gene_pair_dict['chrom_num'] != chrom_num:
            continue
        variant_position = int(variant_gene_pair_dict['variant_position'])
        county = county + chromosome[variant_position]
    return county

def count_variant_overlap_specificity(chrom_num, heart_chromosome, ipsc_chromosome, sig_variant_gene_pairs, variant_gene_pair_info, cell_line_version):
    county = 0
    for variant_gene in sig_variant_gene_pairs.keys():
        variant_gene_pair_dict = variant_gene_pair_info[variant_gene]
        if variant_gene_pair_dict['chrom_num'] != chrom_num:
            continue
        variant_position = int(variant_gene_pair_dict['variant_position'])
        if cell_line_version == 'heart_and_ipsc_cell_lines':
            if heart_chromosome[variant_position] == 1.0 and ipsc_chromosome[variant_position] == 1.0:
                county = county + 1
        elif cell_line_version == 'heart_only_cell_lines':
            if heart_chromosome[variant_position] == 1.0 and ipsc_chromosome[variant_position] == 0.0:
                county = county + 1
        elif cell_line_version == 'ipsc_only_cell_lines':
            if heart_chromosome[variant_position] == 0.0 and ipsc_chromosome[variant_position] == 1.0:
                county = county + 1
    return county

def in_same_direction(estimated_t0, estimated_t15):
    if estimated_t0 <= .5 and estimated_t15 <= .5:
        return True
    elif estimated_t0 > .5 and estimated_t15 > .5:
        return True
    else:
        return False

def extract_variant_gene_pair_time_step_info(time_step_independent_stem):
    # Dictionary to keep track of mapping from variant-gene pairs to vector of length time-steps taht contains pvalues
    dicti = {}
    dicti_pval = {}
    # loop through all time steps
    for time_step in range(16):
        file_name = time_step_independent_stem + str(time_step) + '_eqtl_results.txt'
        head_count = 0  # Used to skip header
        f = open(file_name)
        for line in f:
            line = line.rstrip()
            data = line.split()
            if head_count == 0:
                head_count = head_count + 1
                continue
            # Extract relevent info
            rs_id = data[3]
            ensamble_id = data[1]
            pvalue = float(data[-1])
            alpha = float(data[-3])
            beta = float(data[-2])
            allelic_fraction = alpha/(alpha+beta)
            test_name = rs_id + '_' + ensamble_id
            # If we've never seen this test before
            if test_name not in dicti:
                if time_step != 0:
                    print('assumpriton erroro')
                    pdb.set_trace()
                dicti[test_name] = np.zeros(16)
            dicti[test_name][time_step] = allelic_fraction
            # If we've never seen this test before
            if test_name not in dicti_pval:
                if time_step != 0:
                    print('assumpriton erroro')
                    pdb.set_trace()
                dicti_pval[test_name] = np.zeros(16)
            dicti_pval[test_name][time_step] = pvalue
        f.close()
    dicti2 = {}
    counter = 0
    for test_name in dicti.keys():
        counter = counter + 1
        allelic_vec = dicti[test_name]
        pvalue_vec = dicti_pval[test_name]
        time_steps = np.arange(16)
        result = scipy.stats.linregress(time_steps,allelic_vec)
        slope = result[0]
        intercept = result[1]
        estimated_t0 = intercept
        estimated_t15 = intercept + (slope*15.0)
        pvalue = result[3]
        '''
        if pvalue < .2:
            if abs(estimated_t15 - .5) > abs(estimated_t0 - .5):
                dicti2[test_name] = 'late_time_step_hits'
            else:
                dicti2[test_name] = 'early_time_step_hits'
        else:
            dicti2[test_name] = 'neither'
        '''
        if abs(estimated_t15 - .5) > abs(estimated_t0 - .5) and in_same_direction(estimated_t0, estimated_t15):
            dicti2[test_name] = 'late_time_step_hits'
        elif abs(estimated_t15 - .5) <= abs(estimated_t0 - .5) and in_same_direction(estimated_t0, estimated_t15):
            dicti2[test_name] = 'early_time_step_hits'
        else:
            dicti2[test_name] = 'change_in_sign_hits'
    return dicti2
############################################
# Command Line Args!
############################################
marker_type = sys.argv[1]
cell_line_version = sys.argv[2]
num_permutations = int(sys.argv[3])
chrom_hmm_input_dir = sys.argv[4]
significant_variant_gene_pairs_file = sys.argv[5]
time_step_independent_stem = sys.argv[6]
output_root = sys.argv[7]
hits_version = sys.argv[8]

# Extract list of cell line ids used for this cell_line_version
cell_line_ids = get_cell_line_ids(cell_line_version, chrom_hmm_input_dir)


# Results file for time step 0 (could have been any time step)
# Only use this to extract maf and distance to tss info
time_step_independent_file = time_step_independent_stem + '0_eqtl_results.txt'

# Creating mapping from variant-gene pairs to a vector of length num_time_steps where each element in vector is a pvalue corresponding to that time stpe in the time step independent analysis
variant_gene_pair_time_step_info = extract_variant_gene_pair_time_step_info(time_step_independent_stem)

# Create Mapping from variant-gene pair to quartuple (chrom_num, variant_position, distToTss, MAF)
variant_gene_pair_info = extract_variant_gene_pair_info(time_step_independent_file)

# First create dictionary list of the significant variant gene pairs where each key is of form $variantID"_"$geneID
sig_variant_gene_pairs = extract_significant_variant_gene_pairs(significant_variant_gene_pairs_file, variant_gene_pair_time_step_info, hits_version)
print(len(sig_variant_gene_pairs))

# Return list of length num_permutations where each element of the dictionary list of variant-gene pairs (of len(sig_variant_gene_pairs)) matched for dist_to_tss and maf
perm_variant_gene_pairs = extract_perm_variant_gene_pairs(sig_variant_gene_pairs, time_step_independent_file, num_permutations, variant_gene_pair_info)
#############################
# Count up number of times a variant overlaps one of our markers
###############################
real_overlaps = 0  # Initialze real counts
perm_overlaps = np.zeros(num_permutations) # initialize vector of perm counts

# Do seperately for each chromosome
for chrom_num in range(1,23):
    print(chrom_num)
    if cell_line_version == 'heart_cell_lines' or cell_line_version == 'ipsc_cell_lines' or cell_line_version == 'all_cell_lines':
        # Make binary array length of a chromosome. If array == 0, no marker there. If array == 1, there is a marker there
        chromosome = make_binary_chromosome(chrom_num, chrom_hmm_input_dir, cell_line_ids, marker_type)

        # Count number of variants that overlap a marker on this chromosome
        real_overlaps = real_overlaps + count_variant_overlap(chrom_num, chromosome, sig_variant_gene_pairs, variant_gene_pair_info)
     
        # Count number of variants that overlap a marker on this chromosome
        for perm_num in range(num_permutations):
            perm_overlaps[perm_num] = perm_overlaps[perm_num] + count_variant_overlap(chrom_num, chromosome, perm_variant_gene_pairs[perm_num], variant_gene_pair_info)
    elif cell_line_version == 'heart_and_ipsc_cell_lines' or cell_line_version == 'heart_only_cell_lines' or cell_line_version == 'ipsc_only_cell_lines':
        heart_cell_line_ids = get_cell_line_ids('heart_cell_lines', chrom_hmm_input_dir)
        ipsc_cell_line_ids = get_cell_line_ids('ipsc_cell_lines', chrom_hmm_input_dir)
        # Make binary array length of a chromosome. If array == 0, no marker there. If array == 1, there is a marker there
        heart_chromosome = make_binary_chromosome(chrom_num, chrom_hmm_input_dir, heart_cell_line_ids, marker_type)
        ipsc_chromosome = make_binary_chromosome(chrom_num, chrom_hmm_input_dir, ipsc_cell_line_ids, marker_type)


        # Count number of variants that overlap a marker on this chromosome
        real_overlaps = real_overlaps + count_variant_overlap_specificity(chrom_num, heart_chromosome, ipsc_chromosome, sig_variant_gene_pairs, variant_gene_pair_info, cell_line_version)

        # Count number of variants that overlap a marker on this chromosome
        for perm_num in range(num_permutations):
            perm_overlaps[perm_num] = perm_overlaps[perm_num] + count_variant_overlap_specificity(chrom_num, heart_chromosome, ipsc_chromosome, perm_variant_gene_pairs[perm_num], variant_gene_pair_info, cell_line_version)

# For each permutation run, compute fischers exact test and print to output file
t = open(output_root + '_enrich.txt', 'w')
num_samp = len(sig_variant_gene_pairs)
t.write('real_overlaps\treal_misses\tperm_overlaps\tperm_misses\todds_ratio\tpvalue\n')
for perm_num in range(num_permutations):
    aa = real_overlaps
    bb = num_samp - real_overlaps
    cc = perm_overlaps[perm_num]
    dd = num_samp - cc
    odds_ratio, pvalue = scipy.stats.fisher_exact([[aa,bb],[cc,dd]])
    t.write(str(aa) + '\t' + str(bb) + '\t' + str(cc) + '\t' + str(dd) + '\t' + str(odds_ratio) + '\t' + str(pvalue) + '\n')
t.close()

