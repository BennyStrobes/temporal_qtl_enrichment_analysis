args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)




load_in_odds_ratios <- function(file_name, adding_constant) {
    aa <- read.table(file_name,header=TRUE)

    real_overlaps <- as.numeric(aa$real_overlaps) + adding_constant
    real_misses <- as.numeric(aa$real_misses) + adding_constant
    perm_overlaps <- as.numeric(aa$perm_overlaps) + adding_constant
    perm_misses <- as.numeric(aa$perm_misses) + adding_constant
    odds_ratios <- (real_overlaps/real_misses)/(perm_overlaps/perm_misses)
    return(odds_ratios)
}


odds_ratio_cell_line_specific_boxplot <- function(ipsc_early_or, ipsc_late_or, cardio_early_or, cardio_late_or, ipsc_change_or, cardio_change_or, output_file, marker_type) {
    odds_ratios <- c()
    roadmap_cell_types <- c()
    dynamic_qtl_versions <- c()


    odds_ratios <- c(odds_ratios, ipsc_early_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("ipsc", length(ipsc_early_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("early_qtl", length(ipsc_early_or)))


    odds_ratios <- c(odds_ratios, ipsc_late_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("ipsc", length(ipsc_late_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("late_qtl", length(ipsc_late_or)))

    odds_ratios <- c(odds_ratios, cardio_early_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("heart", length(cardio_early_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("early_qtl", length(cardio_early_or)))

    odds_ratios <- c(odds_ratios, cardio_late_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("heart", length(cardio_early_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("late_qtl", length(cardio_late_or)))


    odds_ratios <- c(odds_ratios, ipsc_change_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("ipsc", length(ipsc_change_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("change_qtl", length(ipsc_change_or)))

    odds_ratios <- c(odds_ratios, cardio_change_or)
    roadmap_cell_types <- c(roadmap_cell_types, rep("heart", length(cardio_change_or)))
    dynamic_qtl_versions <- c(dynamic_qtl_versions, rep("change_qtl", length(cardio_change_or)))

    df <- data.frame(odds_ratio=odds_ratios,roadmap_cell_type=factor(roadmap_cell_types,levels=c("ipsc","heart")), qtl_version=factor(dynamic_qtl_versions))
    print(df)
    # PLOT
    boxplot <- ggplot(df, aes(x=roadmap_cell_type,y=odds_ratio,fill=qtl_version)) + geom_boxplot() + labs(x = "Roadmap Cell Type", y = "Odds Ratio", title= marker_type)
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot 
    boxplot <- boxplot
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")


}


input_root = args[1]
output_root = args[2]
error_bound = args[3]
num_permutations = args[4]


adding_constant <- 0

########################
# Promoter
########################
promoter_ipsc_file_early <- paste0(input_root, "mark_promotor_cl_ipsc_cell_lines_hits_early_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
promoter_ipsc_file_late <- paste0(input_root, "mark_promotor_cl_ipsc_cell_lines_hits_late_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
promoter_ipsc_file_change <- paste0(input_root, "mark_promotor_cl_ipsc_cell_lines_hits_change_in_sign_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )


promoter_cardio_file_early <- paste0(input_root, "mark_promotor_cl_heart_cell_lines_hits_early_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
promoter_cardio_file_late <- paste0(input_root, "mark_promotor_cl_heart_cell_lines_hits_late_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
promoter_cardio_file_change <- paste0(input_root, "mark_promotor_cl_heart_cell_lines_hits_change_in_sign_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )






promoter_ipsc_early_or <- load_in_odds_ratios(promoter_ipsc_file_early, adding_constant)
promoter_ipsc_late_or <- load_in_odds_ratios(promoter_ipsc_file_late, adding_constant)
promoter_ipsc_change_or <- load_in_odds_ratios(promoter_ipsc_file_change, adding_constant)

promoter_cardio_early_or <- load_in_odds_ratios(promoter_cardio_file_early, adding_constant)
promoter_cardio_late_or <- load_in_odds_ratios(promoter_cardio_file_late, adding_constant)
promoter_cardio_change_or <- load_in_odds_ratios(promoter_cardio_file_change, adding_constant)




output_file <- paste0(output_root, "promoter_cell_line_specific_num_perm_", num_permutations,"_error_bound_", error_bound, "_odds_ratios.png")

odds_ratio_cell_line_specific_boxplot(promoter_ipsc_early_or, promoter_ipsc_late_or, promoter_cardio_early_or, promoter_cardio_late_or, promoter_ipsc_change_or, promoter_cardio_change_or, output_file, "promoter")



########################
# Enhancer
########################



enhancer_ipsc_file_early <- paste0(input_root, "mark_enhancer_cl_ipsc_cell_lines_hits_early_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
enhancer_ipsc_file_late <- paste0(input_root, "mark_enhancer_cl_ipsc_cell_lines_hits_late_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
enhancer_ipsc_file_change <- paste0(input_root, "mark_enhancer_cl_ipsc_cell_lines_hits_change_in_sign_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )


enhancer_cardio_file_early <- paste0(input_root, "mark_enhancer_cl_heart_cell_lines_hits_early_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
enhancer_cardio_file_late <- paste0(input_root, "mark_enhancer_cl_heart_cell_lines_hits_late_time_step_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )
enhancer_cardio_file_change <- paste0(input_root, "mark_enhancer_cl_heart_cell_lines_hits_change_in_sign_hits_num_", num_permutations, "_bound_",error_bound, "_enrich.txt" )






enhancer_ipsc_early_or <- load_in_odds_ratios(enhancer_ipsc_file_early, adding_constant)
enhancer_ipsc_late_or <- load_in_odds_ratios(enhancer_ipsc_file_late, adding_constant)
enhancer_ipsc_change_or <- load_in_odds_ratios(enhancer_ipsc_file_change, adding_constant)

enhancer_cardio_early_or <- load_in_odds_ratios(enhancer_cardio_file_early, adding_constant)
enhancer_cardio_late_or <- load_in_odds_ratios(enhancer_cardio_file_late, adding_constant)
enhancer_cardio_change_or <- load_in_odds_ratios(enhancer_cardio_file_change, adding_constant)




output_file <- paste0(output_root, "enhancer_cell_line_specific_num_perm_", num_permutations,"_error_bound_", error_bound, "_odds_ratios.png")

odds_ratio_cell_line_specific_boxplot(enhancer_ipsc_early_or, enhancer_ipsc_late_or, enhancer_cardio_early_or, enhancer_cardio_late_or, enhancer_ipsc_change_or, enhancer_cardio_change_or, output_file, "enhancer")

