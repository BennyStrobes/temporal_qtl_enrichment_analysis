#!/bin/bash
#SBATCH --time=9:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
num_permutations="$2"
chrom_hmm_input_dir="$3"
significant_variant_gene_pairs_file="$4"
time_step_independent_stem="$5"
chrom_hmm_enrichment_directory="$6"
visualization_directory="$7"



marker_type="enhancer"
marker_type="promotor"
cell_line_version="heart_cell_lines"
cell_line_version="ipsc_cell_lines"
cell_line_version="all_cell_lines"



#######################
marker_type="enhancer"
#######################

cell_line_version="heart_cell_lines"
echo $marker_type"  "$cell_line_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string


cell_line_version="ipsc_cell_lines"
echo $marker_type"  "$cell_line_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string



cell_line_version="all_cell_lines"
echo $marker_type"  "$cell_line_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string




#######################
marker_type="promotor"
#######################

cell_line_version="heart_cell_lines"
echo $marker_type"  "$cell_line_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string


cell_line_version="ipsc_cell_lines"
echo $marker_type"  "$cell_line_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string



cell_line_version="all_cell_lines"
echo $marker_type"  "$cell_line_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string

