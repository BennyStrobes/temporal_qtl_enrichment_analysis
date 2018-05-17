#!/bin/bash
#SBATCH --time=9:00:00 --partition=broadwl --mem=5GB

parameter_string="$1"
significant_variant_gene_pairs_file="$2"
time_step_independent_stem="$3"
gencode_file="$4"
gene_set_enrichment_directory="$5"
gsea_data_dir="$6"



hits_version="all_hits"
#python gene_set_enrichment_analysis.py $parameter_string $hits_version $significant_variant_gene_pairs_file $time_step_independent_stem $gencode_file $gene_set_enrichment_directory $gsea_data_dir


hits_version="early_time_step_hits"
# python gene_set_enrichment_analysis.py $parameter_string $hits_version $significant_variant_gene_pairs_file $time_step_independent_stem $gencode_file $gene_set_enrichment_directory $gsea_data_dir

hits_version="late_time_step_hits"
#  python gene_set_enrichment_analysis.py $parameter_string $hits_version $significant_variant_gene_pairs_file $time_step_independent_stem $gencode_file $gene_set_enrichment_directory $gsea_data_dir


hits_version="change_in_sign_hits"
python gene_set_enrichment_analysis.py $parameter_string $hits_version $significant_variant_gene_pairs_file $time_step_independent_stem $gencode_file $gene_set_enrichment_directory $gsea_data_dir