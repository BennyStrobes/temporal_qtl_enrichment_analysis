args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)








odds_ratio_boxplot <- function(promotor_ipsc, promotor_cardio, promotor_all, output_file, marker_type) {
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
    boxplot <- boxplot + labs(x = "chromHMM cell lines", y = "Odds Ratio", title = paste0("marker = ",marker_type))
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}






input_root = args[1]
output_root = args[2]

promotor_ipsc_file <- paste0(input_root, "marker_type_promotor_cell_line_version_ipsc_cell_lines_num_perm_100_chrom_hmm_enrichments.txt")
promotor_cardio_file <- paste0(input_root, "marker_type_promotor_cell_line_version_heart_cell_lines_num_perm_100_chrom_hmm_enrichments.txt")
promotor_all_file <- paste0(input_root, "marker_type_promotor_cell_line_version_all_cell_lines_num_perm_100_chrom_hmm_enrichments.txt")


enhancer_ipsc_file <- paste0(input_root, "marker_type_enhancer_cell_line_version_ipsc_cell_lines_num_perm_100_chrom_hmm_enrichments.txt")
enhancer_cardio_file <- paste0(input_root, "marker_type_enhancer_cell_line_version_heart_cell_lines_num_perm_100_chrom_hmm_enrichments.txt")
enhancer_all_file <- paste0(input_root, "marker_type_enhancer_cell_line_version_all_cell_lines_num_perm_100_chrom_hmm_enrichments.txt")


promotor_ipsc <- read.table(promotor_ipsc_file, header=TRUE)
promotor_cardio <- read.table(promotor_cardio_file, header=TRUE)
promotor_all <- read.table(promotor_all_file, header=TRUE)

enhancer_ipsc <- read.table(enhancer_ipsc_file, header=TRUE)
enhancer_cardio <- read.table(enhancer_cardio_file, header=TRUE)
enhancer_all <- read.table(enhancer_all_file, header=TRUE)



output_file <- paste0(output_root, "marker_type_promotor_odds_ratio_enrichment.png")
odds_ratio_boxplot(promotor_ipsc, promotor_cardio, promotor_all, output_file, "promotor")



output_file <- paste0(output_root, "marker_type_enhancer_odds_ratio_enrichment.png")
odds_ratio_boxplot(enhancer_ipsc, enhancer_cardio, enhancer_all, output_file, "enhancer")
