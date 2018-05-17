#####################################################################
# Perform downstream/enrichment analysis on temporal qtl analyses
# Ben Strober (bstrober3@gmail.com)
#####################################################################


#####################################################################
#INPUT DATA
#####################################################################
#### FDR <= .01
# All output files will have this parameter string
parameter_string="te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.01_"

#File containing list of variant gene pairs that are significant
significant_variant_gene_pairs_file="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.01_significant_egenes.txt"

all_significant_variant_gene_pairs_file="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.01_significant.txt"


#### FDR <= .01
# All output files will have this parameter string
parameter_string_05="te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.05_"

#File containing list of variant gene pairs that are significant
significant_variant_gene_pairs_file_05="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.05_significant_egenes.txt"

all_significant_variant_gene_pairs_file_05="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/ipsc_data_te/qtl_results/te_log_linear_environmental_variable_time_steps_optimizer_LBFGS_genotype_round_covariate_method_cell_line_pc1Xtime_permutation_scheme_sample_null_efdr_.05_significant.txt"



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

# File containing conversions from ensamble ids to gene symbols
gencode_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/gencode.v19.annotation.gtf.gz"


# Directory containing gsea data
gsea_data_dir="/project2/gilad/bstrober/tools/tools/gsea/data/"

#####################################################################
#OUTPUT DATA
#####################################################################

# ROOT DIRECTORY FOR ALL OUTPUT FILES
output_root="/project2/gilad/bstrober/ipsc_differentiation/dynamic_qtl_pipelines/temporal_qtl_enrichment_analysis/"


# Directory containing results for chromHmm enrichment analysis
chrom_hmm_enrichment_directory=$output_root"chrom_hmm_enrichment/"

# Directory containing results for distance to tss analysis
time_step_independent_comparison_directory=$output_root"time_step_independent_comparison/"

# Directory containing results for cell line overlap analysis
cell_line_overlap_directory=$output_root"cell_line_overlap/"

# Directory containing gene set enrichment analysis results
gene_set_enrichment_directory=$output_root"gene_set_enrichment/"


# Directory containing visualizations from all analyses
visualization_directory=$output_root"visualization/"




#########################################################
# PART 1: Run gene set enrichment analysis
# Run gene set enrichment on dynamic qtl result
#########################################################
if false; then
sbatch gene_set_enrichment_analysis.sh $parameter_string $significant_variant_gene_pairs_file $time_step_independent_stem $gencode_file $gene_set_enrichment_directory $gsea_data_dir
fi


#########################################################
# PART 2: ChromHmm enrichment analysis
# Compute enrichment of:
### A. 'promotor'
### B. 'enhancer'
# based on Chromhmm cell types for hits vs randomly matched (for maf and disttoTss) background
# Perform analysis based on:
### 'heart_cell_lines'
### 'ipsc_cell_lines'
# Also, do seperate analysis investigating:
### 'change_in_sign_hits'
### 'early_time_step_hits'
### 'late_time_step_hits'
## Define 'change_in_sign_hits' using error_bound

##########################
# We use hits defined at eFDR <= .05 to increase our sample size in the enrichment analysis
#########################################################
num_permutations="100"

if false; then
# Run for a range of error bounds
error_bound=".0"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound

error_bound=".005"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound

error_bound=".01"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound

error_bound=".015"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound

error_bound=".02"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound

error_bound=".025"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound

error_bound=".03"
sbatch chrom_hmm_enrichment_analysis.sh $parameter_string_05 $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file_05 $time_step_independent_stem $chrom_hmm_enrichment_directory $visualization_directory $error_bound
fi


if false; then

Rscript global_odds_ratio_plotter.R $parameter_string $num_permutations $chrom_hmm_enrichment_directory $visualization_directory

fi

#########################################################
# PART 3: Cell Line overlap analysis
# Compute how often a given pair of cell lines overlap in their genotype in hits compared to background
#########################################################
if false; then
sh cell_line_overlap_analysis.sh $significant_variant_gene_pairs_file $real_dynamic_qtl_results_file $genotype_file $parameter_string $cell_line_overlap_directory $visualization_directory $time_step_independent_stem
fi




#########################################################
# PART 4: time step independent comparison
# Compare dynamic qtl results with time step independent analysis
#########################################################
if false; then
sh time_step_independent_comparison.sh $parameter_string $significant_variant_gene_pairs_file $time_step_independent_stem $time_step_independent_comparison_directory $visualization_directory
fi

