args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)








odds_ratio_boxplot <- function(promotor_ipsc, promotor_cardio, promotor_all, output_file, marker_type, hits_version) {
    odds_ratios <- c()
    cell_lines <- c()

    odds_ratios <- c(odds_ratios, promotor_ipsc$odds_ratio)
    cell_lines <- c(cell_lines, rep("ipsc_lines", length(promotor_ipsc$odds_ratio)))


    odds_ratios <- c(odds_ratios, promotor_cardio$odds_ratio)
    cell_lines <- c(cell_lines, rep("heart_lines", length(promotor_cardio$odds_ratio)))

    odds_ratios <- c(odds_ratios, promotor_all$odds_ratio)
    cell_lines <- c(cell_lines, rep("all_lines", length(promotor_all$odds_ratio)))

    df <- data.frame(odds_ratio=odds_ratios,cell_lines=factor(cell_lines))
        # PLOT
    boxplot <- ggplot(df, aes(x=cell_lines,y=odds_ratios,fill=cell_lines)) + geom_boxplot() + labs(x = "Cell Line", y = "Odds Ratio")
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot + theme(legend.position="none")
    boxplot <- boxplot + labs(x = "chromHMM cell lines", y = "Odds Ratio", title = paste0("hits = ", hits_version, " / marker = ",marker_type))
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}


odds_ratio_boxplot_v2 <- function(promotor_heart_ipsc, promotor_heart_only, promotor_ipsc_only, output_file, marker_type, hits_version) {
    odds_ratios <- c()
    cell_lines <- c()

    odds_ratios <- c(odds_ratios, promotor_heart_ipsc)
    cell_lines <- c(cell_lines, rep("heart_ipsc_lines", length(promotor_heart_ipsc)))


    odds_ratios <- c(odds_ratios, promotor_heart_only)
    cell_lines <- c(cell_lines, rep("heart_only_lines", length(promotor_heart_only)))

    odds_ratios <- c(odds_ratios, promotor_ipsc_only)
    cell_lines <- c(cell_lines, rep("ipsc_only_lines", length(promotor_ipsc_only)))

    df <- data.frame(odds_ratio=odds_ratios,cell_lines=factor(cell_lines))
        # PLOT
    boxplot <- ggplot(df, aes(x=cell_lines,y=odds_ratios,fill=cell_lines)) + geom_boxplot() + labs(x = "Cell Line", y = "Odds Ratio")
    boxplot <- boxplot + theme(text = element_text(size=18))
    boxplot <- boxplot + geom_hline(yintercept = 1.0) 
    boxplot <- boxplot + theme(legend.position="none")
    boxplot <- boxplot + labs(x = "chromHMM cell lines", y = "Odds Ratio", title = paste0("hits = ", hits_version, " / marker = ",marker_type))
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")


}






load_in_odds_ratios <- function(file_name, adding_constant) {
    aa <- read.table(file_name,header=TRUE)

    real_overlaps <- as.numeric(aa$real_overlaps) + adding_constant
    real_misses <- as.numeric(aa$real_misses) + adding_constant
    perm_overlaps <- as.numeric(aa$perm_overlaps) + adding_constant
    perm_misses <- as.numeric(aa$perm_misses) + adding_constant
    odds_ratios <- (real_overlaps/real_misses)/(perm_overlaps/perm_misses)
    return(odds_ratios)
}




input_root = args[1]
output_root = args[2]
hits_version = args[3]

promotor_ipsc_file <- paste0(input_root, "mark_promotor_cl_ipsc_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
promotor_cardio_file <- paste0(input_root, "mark_promotor_cl_heart_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
promotor_all_file <- paste0(input_root, "mark_promotor_cl_all_cell_lines_hits_", hits_version,"_num_100_enrich.txt")

enhancer_ipsc_file <- paste0(input_root, "mark_enhancer_cl_ipsc_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
enhancer_cardio_file <- paste0(input_root, "mark_enhancer_cl_heart_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
enhancer_all_file <- paste0(input_root, "mark_enhancer_cl_all_cell_lines_hits_", hits_version,"_num_100_enrich.txt")


promotor_heart_ipsc_file <- paste0(input_root, "mark_promotor_cl_heart_and_ipsc_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
promotor_heart_only_file <- paste0(input_root, "mark_promotor_cl_heart_only_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
promotor_ipsc_only_file <- paste0(input_root, "mark_promotor_cl_ipsc_only_cell_lines_hits_", hits_version,"_num_100_enrich.txt")


enhancer_heart_ipsc_file <- paste0(input_root, "mark_enhancer_cl_heart_and_ipsc_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
enhancer_heart_only_file <- paste0(input_root, "mark_enhancer_cl_heart_only_cell_lines_hits_", hits_version,"_num_100_enrich.txt")
enhancer_ipsc_only_file <- paste0(input_root, "mark_enhancer_cl_ipsc_only_cell_lines_hits_", hits_version,"_num_100_enrich.txt")


adding_constant <- 1

promotor_heart_ipsc <- load_in_odds_ratios(promotor_heart_ipsc_file, adding_constant)
promotor_heart_only <- load_in_odds_ratios(promotor_heart_only_file, adding_constant)
promotor_ipsc_only <- load_in_odds_ratios(promotor_ipsc_only_file, adding_constant)

enhancer_heart_ipsc <- load_in_odds_ratios(enhancer_heart_ipsc_file, adding_constant)
enhancer_heart_only <- load_in_odds_ratios(enhancer_heart_only_file, adding_constant)
enhancer_ipsc_only <- load_in_odds_ratios(enhancer_ipsc_only_file, adding_constant)

#promotor_ipsc <- read.table(promotor_ipsc_file, header=TRUE)
#promotor_cardio <- read.table(promotor_cardio_file, header=TRUE)
#promotor_all <- read.table(promotor_all_file, header=TRUE)

#enhancer_ipsc <- read.table(enhancer_ipsc_file, header=TRUE)
#enhancer_cardio <- read.table(enhancer_cardio_file, header=TRUE)
#enhancer_all <- read.table(enhancer_all_file, header=TRUE)



output_file <- paste0(output_root, hits_version, "_marker_type_promotor_odds_ratio_enrichment_v2.png")
odds_ratio_boxplot_v2(promotor_heart_ipsc, promotor_heart_only, promotor_ipsc_only, output_file, "promotor", hits_version)



output_file <- paste0(output_root, hits_version, "_marker_type_enhancer_odds_ratio_enrichment_v2.png")
odds_ratio_boxplot_v2(enhancer_heart_ipsc, enhancer_heart_only, enhancer_ipsc_only, output_file, "enhancer", hits_version)














output_file <- paste0(output_root, hits_version, "_marker_type_promotor_odds_ratio_enrichment.png")
#odds_ratio_boxplot(promotor_ipsc, promotor_cardio, promotor_all, output_file, "promotor", hits_version)



output_file <- paste0(output_root, hits_version, "_marker_type_enhancer_odds_ratio_enrichment.png")
#odds_ratio_boxplot(enhancer_ipsc, enhancer_cardio, enhancer_all, output_file, "enhancer", hits_version)