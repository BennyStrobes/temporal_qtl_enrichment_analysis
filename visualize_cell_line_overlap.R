args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)




load_in_data <- function(matrix) {
    aa <- read.table(matrix,header=TRUE)
    rownamers <- aa[,1]
    bb <- aa[,2:(dim(aa)[2])]
    rownames(bb) <- rownamers
    colnames(bb) <- rownamers
    num_rows <- dim(bb)[1]
    for (row_num in 1:num_rows) {
        bb[row_num, row_num] = NA
    }
    return(as.matrix(bb))
}

plot_overlap_heatmap <- function(overlap_mat, output_file, version, minimum, maximum) {
    melted_corr <- melt(overlap_mat)
    print(colnames(overlap_mat))

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    ord <- hclust( dist(scale(overlap_mat), method = "euclidean"), method = "ward.D" )$order

    print(ord)

    #  Use factors to represent covariate and pc name
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(overlap_mat)[ord])
    melted_corr$X2 <- factor(melted_corr$X2, levels = rownames(overlap_mat)[ord])


    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", limits=c(minimum,maximum))
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + labs(x = "Cell Lines", y = "Cell Lines", title=paste0(version))

    ggsave(heatmap, file=output_file,width = 18,height=13.5,units="cm")
}

plot_overlap_heatmap_t15_troponin <- function(overlap_mat, output_file, version, minimum, maximum) {
    melted_corr <- melt(overlap_mat)
    print(colnames(overlap_mat))

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    ord <- hclust( dist(scale(overlap_mat), method = "euclidean"), method = "ward.D" )$order

    ord <- c(8,2,6,13,7,12,14,10,11,4,1,5,9,3)
    print(length(unique(ord)))

    #  Use factors to represent covariate and pc name
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(overlap_mat)[ord])
    melted_corr$X2 <- factor(melted_corr$X2, levels = rownames(overlap_mat)[ord])


    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", limits=c(minimum,maximum))
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + labs(x = "Cell Lines", y = "Cell Lines", title=paste0(version))

    ggsave(heatmap, file=output_file,width = 18,height=13.5,units="cm")
}


plot_overlap_heatmap_cell_line_pc1 <- function(overlap_mat, output_file, version, minimum, maximum) {
    melted_corr <- melt(overlap_mat)
    print(colnames(overlap_mat))

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1)
    melted_corr$X2 <- factor(melted_corr$X2)

    ord <- hclust( dist(scale(overlap_mat), method = "euclidean"), method = "ward.D" )$order

    ord <- c(12, 11, 4, 1, 3, 10, 5, 9, 14, 7, 6, 2, 13, 8)
    print(length(unique(ord)))

    #  Use factors to represent covariate and pc name
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(overlap_mat)[ord])
    melted_corr$X2 <- factor(melted_corr$X2, levels = rownames(overlap_mat)[ord])


    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", limits=c(minimum,maximum))
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + labs(x = "Cell Lines", y = "Cell Lines", title=paste0(version))

    ggsave(heatmap, file=output_file,width = 18,height=13.5,units="cm")
}



###################
# Command line args
###################
real_overlap_matrix <- args[1]
perm_overlap_matrix <- args[2]
time_step_independent_overlap_matrix <- args[3]
time_step_independent_perm_overlap_matrix <- args[4]
overlap_output_root <- args[5]



real_overlap <- load_in_data(real_overlap_matrix)

perm_overlap <- load_in_data(perm_overlap_matrix)

time_overlap <- load_in_data(time_step_independent_overlap_matrix)

time_perm_overlap <- load_in_data(time_step_independent_perm_overlap_matrix)

minimum <- min(min(real_overlap, na.rm=TRUE), min(perm_overlap, na.rm=TRUE))
maximum <- max(max(real_overlap, na.rm=TRUE), max(perm_overlap, na.rm=TRUE))

minimum_time <- min(min(time_overlap, na.rm=TRUE), min(time_perm_overlap, na.rm=TRUE))
maximum_time <- max(max(time_overlap, na.rm=TRUE), max(time_perm_overlap, na.rm=TRUE))


real_heatmap_output <- paste0(overlap_output_root, "real_heatmap_cell_line_pc1.png")
plot_overlap_heatmap_cell_line_pc1(real_overlap, real_heatmap_output, "real_cell_line_pc1", minimum, maximum)


# Make heatmap for real data
real_heatmap_output <- paste0(overlap_output_root, "real_heatmap.png")
plot_overlap_heatmap(real_overlap, real_heatmap_output, "real", minimum, maximum)


# Make heatmap for real data
perm_heatmap_output <- paste0(overlap_output_root, "perm_heatmap.png")
plot_overlap_heatmap(perm_overlap, perm_heatmap_output, "perm", minimum, maximum)

# Make heatmap for real data
time_heatmap_output <- paste0(overlap_output_root, "time_step_independent_real_heatmap.png")
plot_overlap_heatmap(time_overlap, time_heatmap_output, "time_ind_real", minimum_time, maximum_time)

# Make heatmap for real data
time_heatmap_output <- paste0(overlap_output_root, "time_step_independent_perm_heatmap.png")
plot_overlap_heatmap(time_perm_overlap, time_heatmap_output, "time_ind_perm", minimum_time, maximum_time)




real_heatmap_output <- paste0(overlap_output_root, "real_heatmap_t15_trop.png")
plot_overlap_heatmap_t15_troponin(real_overlap, real_heatmap_output, "real_t15_troponin", minimum, maximum)