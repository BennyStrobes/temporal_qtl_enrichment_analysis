#!/bin/bash
#SBATCH --time=9:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
num_permutations="$2"
chrom_hmm_input_dir="$3"
significant_variant_gene_pairs_file="$4"
time_step_independent_stem="$5"
chrom_hmm_enrichment_directory="$6"
visualization_directory="$7"
hits_version="$8"


#######################
marker_type="enhancer"
#######################

cell_line_version="heart_cell_lines"
echo $marker_type"  "$cell_line_version" "$hits_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_hits_version_"$hits_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string $hits_version




cell_line_version="ipsc_cell_lines"
echo $marker_type"  "$cell_line_version" "$hits_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_hits_version_"$hits_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string $hits_version



cell_line_version="all_cell_lines"
echo $marker_type"  "$cell_line_version" "$hits_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_hits_version_"$hits_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string $hits_version




#######################
marker_type="promotor"
#######################

cell_line_version="heart_cell_lines"
echo $marker_type"  "$cell_line_version" "$hits_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_hits_version_"$hits_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string $hits_version


cell_line_version="ipsc_cell_lines"
echo $marker_type"  "$cell_line_version" "$hits_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_hits_version_"$hits_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string $hits_version



cell_line_version="all_cell_lines"
echo $marker_type"  "$cell_line_version" "$hits_version
mod_parameter_string=$parameter_string"marker_type_"$marker_type"_cell_line_version_"$cell_line_version"_hits_version_"$hits_version"_num_perm_"$num_permutations
python organize_chrom_hmm_enrichment_analysis.py $marker_type $cell_line_version $num_permutations $chrom_hmm_input_dir $significant_variant_gene_pairs_file $time_step_independent_stem $chrom_hmm_enrichment_directory$mod_parameter_string $hits_version



if false; then
Rscript visualize_chrom_hmm_enrichment_analysis.R $chrom_hmm_enrichment_directory$parameter_string $visualization_directory$parameter_string
fi
