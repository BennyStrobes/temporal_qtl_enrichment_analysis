#!/bin/bash
#SBATCH --time=9:00:00 --partition=broadwl --mem=5GB


parameter_string="$1"
egenes_file="$2"
time_step_independent_stem="$3"
time_step_independent_comparison_directory="$4"
visualization_directory="$5"


dynamic_egenes_comparison_file=$time_step_independent_comparison_directory$parameter_string"egenes_comparison.txt"
dynamic_egenes_background_comparison_file=$time_step_independent_comparison_directory$parameter_string"background_egenes_comparison.txt"

python organize_time_step_independent_comparison.py $dynamic_egenes_comparison_file $dynamic_egenes_background_comparison_file $time_step_independent_stem $egenes_file


Rscript visualize_time_step_independent_comparison.R $dynamic_egenes_comparison_file $dynamic_egenes_background_comparison_file $visualization_directory$parameter_string