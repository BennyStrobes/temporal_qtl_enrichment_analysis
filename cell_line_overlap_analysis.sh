#!/bin/bash
#SBATCH --time=9:00:00 --partition=broadwl --mem=5GB


significant_variant_gene_pairs_file="$1"
real_dynamic_qtl_results_file="$2"
genotype_file="$3"
parameter_string="$4"
cell_line_overlap_directory="$5"
visualization_directory="$6"


real_overlap_matrix=$cell_line_overlap_directory$parameter_string"real_overlap_matrix.txt"
perm_overlap_matrix=$cell_line_overlap_directory$parameter_string"perm_overlap_matrix.txt"

python organize_cell_line_overlap_analysis.py $significant_variant_gene_pairs_file $real_dynamic_qtl_results_file $genotype_file $real_overlap_matrix $perm_overlap_matrix


overlap_output_root=$visualization_directory$parameter_string
Rscript visualize_cell_line_overlap.R $real_overlap_matrix $perm_overlap_matrix $overlap_output_root