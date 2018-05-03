#####################################################################
# Perform downstream/enrichment analysis on temporal qtl analyses
# Ben Strober (bstrober3@gmail.com)
#####################################################################


#####################################################################
#INPUT DATA
#####################################################################
# All output files will have this parameter string
parameter_string="te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.01_"

#File containing list of variant gene pairs that are significant
significant_variant_gene_pairs_file="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.01_significant_egenes.txt"

all_significant_variant_gene_pairs_file="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.01_significant.txt"



# File containing results for all dynamic qtls tested
real_dynamic_qtl_results_file="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_none_permute_False_merged_dynamic_qtl_results.txt"

# File containing permuted results for all dynamic qtls tested
perm_dynamic_qtl_results_file="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_permute_True_merged_dynamic_qtl_results.txt"



# Time step independent data
# Suffixes are:
#### 1. "$time_step"_eqtl_results.txt, containing pvalues for each time step, along with MAF, and dist-toTSS for our variant-gene pairs
#### 2. "$time_step"_efdr_thresh_.1_significant_egenes.txt giving list of significant variant gene pairs at this time step
time_step_independent_stem="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/cht_output/cht_results_cis_distance_50000_maf_cutoff_0.1_min_reads_90_min_as_reads_20_min_het_counts_4_num_pc_1_time_"

# Directory containing chromHMM results
# Each cell line has its own file with suffix $cell_line_identifier'_15_coreMarks_mnemonics.bed.gz'
chrom_hmm_input_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/chrom_hmm/"

# File containing dosage based genotypes
genotype_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess/genotype/YRI_genotype.vcf"


#####################################################################
#OUTPUT DATA
#####################################################################

# ROOT DIRECTORY FOR ALL OUTPUT FILES
output_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/temporal_qtl_enrichment_analysis/"


# Directory containing results for distance to tss analysis
distance_to_tss_directory=$output_root"distance_to_tss/"

# Directory containing results for chromHmm enrichment analysis
chrom_hmm_enrichment_directory=$output_root"chrom_hmm_enrichment/"

# Directory containing results for distance to tss analysis
time_step_independent_comparison_directory=$output_root"time_step_independent_comparison/"

# Directory containing results for cell line overlap analysis
cell_line_overlap_directory=$output_root"cell_line_overlap/"


# Directory containing visualizations from all analyses
visualization_directory=$output_root"visualization/"





#########################################################
# PART 1: Distance Analysis
# Compute distance to TSS for all hits.
# Compute distance to TSS for matched background set.
#########################################################

# Number of background sets/permutations to run
num_permutations="100"
if false; then
sh distance_analysis.sh $parameter_string $significant_variant_gene_pairs_file $time_step_independent_stem $distance_to_tss_directory $visualization_directory $num_permutations
fi





#########################################################
# PART 2: ChromHmm enrichment analysis
# Compute enrichment of:
### A. 'promotor'
### B. 'enhancer'
# based on Chromhmm cell types for hits vs randomly matched (for maf and disttoTss) background
# Perform analysis based on:
### 'all_cell_lines'
### 'heart_cell_lines'
### 'ipsc_cell_lines'
#########################################################
num_permutations="100"
sh chrom_hmm_enrichment_analysis.sh $parameter_string $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory




#########################################################
# PART 3: Cell Line overlap analysis
# Compute how often a given pair of cell lines overlap in their genotype in hits compared to background
#########################################################

if false; then
sh cell_line_overlap_analysis.sh $all_significant_variant_gene_pairs_file $real_dynamic_qtl_results_file $genotype_file $parameter_string $cell_line_overlap_directory $visualization_directory
fi




#########################################################
# PART 4: time step independent comparison
# Compare dynamic qtl results with time step independent analysis
#########################################################
if false; then
sh time_step_independent_comparison.sh $parameter_string $significant_variant_gene_pairs_file $time_step_independent_stem $time_step_independent_comparison_directory $visualization_directory
fi

