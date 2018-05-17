import numpy as np 
import os
import sys
import pdb
import gzip
import scipy.stats

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




def extract_ensamble_ids(time_step_independent_file, significant_variant_gene_pairs_file, variant_gene_pair_time_step_info, hits_version):
    ########################
    # Extract dictionary containing names of all tested genes
    all_genes = {}
    f = open(time_step_independent_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[1]
        all_genes[ensamble_id] = 1
    f.close()
    ########################
    # Extract dictionary list of genes that are signficant in this hits_version
    significant_genes = {}
    f = open(significant_variant_gene_pairs_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        ensamble_id = data[5]
        rs_id = data[2]
        if hits_version == 'all_hits':
            significant_genes[ensamble_id] = 1
        elif hits_version == 'early_time_step_hits' and variant_gene_pair_time_step_info[rs_id + '_' + ensamble_id] == 'early_time_step_hits':
            significant_genes[ensamble_id] = 1
        elif hits_version == 'late_time_step_hits' and variant_gene_pair_time_step_info[rs_id + '_' + ensamble_id] == 'late_time_step_hits':
            significant_genes[ensamble_id] = 1
        elif hits_version == 'change_in_sign_hits' and variant_gene_pair_time_step_info[rs_id + '_' + ensamble_id] == 'change_in_sign_hits':
            significant_genes[ensamble_id] = 1
    f.close()
    return significant_genes, all_genes



def convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file):
    f = gzip.open(gencode_file)
    gene_symbol_hits = []
    for line in f:
        line = line.decode('utf-8').rstrip()
        data = line.split()
        if line.startswith('#'):
            continue
        line_ensamble_id = data[9].split('"')[1].split('.')[0]
        line_gene_symbol = data[17].split('"')[1]
        if line_ensamble_id in ensamble_hits:
            gene_symbol_hits.append(line_gene_symbol)
    return np.unique(gene_symbol_hits)


def print_array(file_name, array):
    t = open(file_name,'w')
    for ele in array:
        t.write(ele + '\n')
    t.close()

def sort_gsea(save_file, new_save_file):
    f = open(save_file)
    t = open(new_save_file,'w')
    pairs = []
    for i,line in enumerate(f):
        line = line.rstrip()
        data = line.split()
        if i < 4:
            continue
        pvalue = float(data[6])
        pairs.append((pvalue, line))
    sorted_pairs = sorted(pairs, key=lambda x: x[0])
    for pair in sorted_pairs:
        liner = pair[1]
        t.write(liner + '\n')
    t.close()


#########################
# Command Line args
#########################
parameter_string = sys.argv[1]  # String that describes current version of the data
hits_version = sys.argv[2]  # what type of dynmaic qtl hits to subset to
significant_variant_gene_pairs_file = sys.argv[3]  # File containing significant dynamic qtl hits
time_step_independent_stem = sys.argv[4]  # Stem to time step independent files
gencode_file = sys.argv[5]  # Gencode gene annotation file to be used to convert from ensamble id to gene-symbol id
gene_set_enrichment_directory = sys.argv[6]  # Output directory
gsea_data_dir = sys.argv[7]  # Directory containing all necessary gsea dat

print(hits_version)

parameter_string = parameter_string + hits_version

# Results file for time step 0 (could have been any time step)
# Only use this to extract maf and distance to tss info
time_step_independent_file = time_step_independent_stem + '0_eqtl_results.txt'

# Creating mapping from variant-gene pairs to a vector of length num_time_steps where each element in vector is a pvalue corresponding to that time stpe in the time step independent analysis
variant_gene_pair_time_step_info = extract_variant_gene_pair_time_step_info(time_step_independent_stem)

ensamble_hits, ensamble_background = extract_ensamble_ids(time_step_independent_file, significant_variant_gene_pairs_file, variant_gene_pair_time_step_info, hits_version)


gene_symbol_hits = convert_from_ensamble_to_gene_symbol(ensamble_hits, gencode_file)

gene_symbol_background = convert_from_ensamble_to_gene_symbol(ensamble_background, gencode_file)

print(len(gene_symbol_hits))
print(len(gene_symbol_background))


hits_file = gene_set_enrichment_directory + parameter_string + '_hit_genes.txt'

background_file = gene_set_enrichment_directory + parameter_string + '_background_genes.txt'


print_array(hits_file, gene_symbol_hits)
print_array(background_file, gene_symbol_background)
#np.savetxt(hits_file, gene_symbol_hits,fmt="%s",delimiter="\n")
#np.savetxt(background_file, gene_symbol_background,fmt="%s",delimiter="\n")


genesets = ['h.all.v5.1.symbols.gmt.txt', 'c2.cp.biocarta.v5.1.symbols.gmt.txt', 'c2.cp.kegg.v5.1.symbols.gmt.txt']
names = ['hallmark', 'biocarta', 'kegg']
for i, val in enumerate(genesets):
    name = names[i]
    geneset_file = gsea_data_dir + val
    save_file = gene_set_enrichment_directory + parameter_string + '_' + name + '_gsea_output.txt'
    #geneset_file = '/project2/gilad/bstrober/tools/tools/gsea/data/' + 'c2.cp.biocarta.v5.1.symbols.gmt.txt'
    os.system('gsea ' + hits_file + ' ' + background_file + ' ' + geneset_file + ' ' + save_file)


    new_save_file = gene_set_enrichment_directory + parameter_string + '_' + name + '_gsea_sorted_output.txt'
    sort_gsea(save_file, new_save_file)
    # Remove un-sorted file
    os.system('rm ' + save_file)

# Remove some unnecessary files
os.system('rm ' + hits_file)
os.system('rm ' + background_file)

