args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)






distance_histogram <- function(real_data, output_file) {
    df <- data.frame(distance_to_tss=real_data)
    histo <- ggplot(data=df, aes(df$distance_to_tss)) +
    geom_histogram(breaks=seq(0, 50000, by = 2000),col="red", fill="green",alpha = .2) +
            labs(x="Distance to TSS (BP)", y="Count")
    ggsave(histo, file=output_file, width=20, height=10.5, units="cm")
}








distance_results_file = args[1]
visualization_directory = args[2]
parameter_string = args[3]


data <- read.table(distance_results_file, header=FALSE)

real_data <- data$V1


# histogram showing distance to tss distribution
output_file <- paste0(visualization_directory, parameter_string, "_real_distance_to_tss_histogram.png")
print(output_file)
distance_histogram(real_data, output_file)