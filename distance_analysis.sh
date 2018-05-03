#!/bin/bash
#SBATCH --time=7:00:00 --partition=broadwl --mem=5GB

parameter_string="$1"
significant_variant_gene_pairs_file="$2"
time_step_independent_stem="$3"
distance_to_tss_directory="$4"
visualization_directory="$5"
num_permutations="$6"


distance_results_file=$distance_to_tss_directory$parameter_string"distance_to_tss_"$num_permutations"_perms.txt"

python organize_distance_analysis.py $significant_variant_gene_pairs_file $time_step_independent_stem $distance_results_file $num_permutations



Rscript visualize_distance_analysis.R $distance_results_file $visualization_directory $parameter_string"_"$num_permutations