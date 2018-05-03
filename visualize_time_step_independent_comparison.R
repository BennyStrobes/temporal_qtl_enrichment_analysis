args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(reshape)




scatter_colored_by_time <- function(all_tests, output_file) {
    nnn <- 80000
    if (nrow(all_tests) < nnn) {
        all_tests_sub <- all_tests
    } else {
        all_tests_sub <- all_tests[sample(nrow(all_tests),nnn),]
    }
    dynamic_pvalz <- c()
    time_step_pvalz <- c()
    time_step_arr <- c()
    for (time_step in 0:15){
        dynamic_pvalz <- c(dynamic_pvalz, -log10(all_tests_sub[,18] + .000000000001))
        time_step_pvalz <- c(time_step_pvalz, -log10(all_tests_sub[,(time_step + 2)] + .000000000001))
        time_step_arr <- c(time_step_arr, rep(time_step, length(all_tests_sub[,18])))
    }
    df <- data.frame(dynamic_pvalue=dynamic_pvalz, time_step_pvalue=time_step_pvalz, time_step=(time_step_arr))
    max_val <-max(max(-log10(dynamic_pvalz + .000000000001)), max(-log10(time_step_pvalz + .000000000001)))
    # PLOT
    scatter <- ggplot(df, aes(x = dynamic_pvalue, y = time_step_pvalue,colour=time_step)) + geom_point(size=.000001) 
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + labs(x = "Dynamic", y = "Time step ind.")
    scatter <- scatter + geom_abline() + scale_color_gradient(low="pink",high="blue")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))

    ggsave(scatter, file=output_file, width=20, height=10.5, units="cm")


}

boxplot_comparing_time_steps <- function(data, output_file) {
    time_step_pvalz <- c()
    time_step_arr <- c()
    for (time_step in 0:15) {
        time_step_pvalz <- c(time_step_pvalz, -log10(data[,(time_step+2)]+ .000000000001))
        time_step_arr <- c(time_step_arr, rep(time_step, length(data[,18])))
    }
    df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr)
        # PLOT
    boxplot <- ggplot(df, aes(x=time_step,y=pvalue,fill=int_time_step)) + geom_boxplot() + labs(x = "Time Step", y = "-log10(pvalue)") + scale_x_discrete(name="Time Step") 
    boxplot <- boxplot + theme(text = element_text(size=18)) + scale_fill_gradient(low="pink",high="blue") + theme(legend.position="none")
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}

boxplot_comparing_time_steps_with_background <- function(data, data_background, output_file) {
    time_step_pvalz <- c()
    time_step_arr <- c()
    real_v_background <- c()
    for (time_step in 0:15) {
        time_step_pvalz <- c(time_step_pvalz, -log10(data[,(time_step+2)]+ .000000000001))
        time_step_pvalz <- c(time_step_pvalz, -log10(data_background[,(time_step+2)]+ .000000000001))
        
        real_v_background <- c(real_v_background, rep("real",length(data[,18])))
        real_v_background <- c(real_v_background, rep("background",length(data[,18])))


        time_step_arr <- c(time_step_arr, rep(time_step, 2*length(data[,18])))
    }
    df <- data.frame(pvalue=time_step_pvalz,time_step=factor(time_step_arr),int_time_step=time_step_arr, data_version=factor(real_v_background))
        # PLOT
    boxplot <- ggplot(df, aes(x=time_step,y=pvalue,fill=data_version)) + geom_boxplot() + labs(x = "Time Step", y = "-log10(pvalue)") + scale_x_discrete(name="Time Step") 
    boxplot <- boxplot + theme(text = element_text(size=18))
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}

boxplot_with_best_time_step <- function(dynamic_egenes, dynamic_egenes_background, output_file) {
    pvalz <- c()
    time_steps <- c()
    real_v_background <- c()
    for (i in 1:nrow(dynamic_egenes)){
        vec <- dynamic_egenes[i,2:17]
        index <- which(vec==min(vec))
        pvalz <- c(pvalz, -log10(min(vec) + .000000000001))
        real_v_background <- c(real_v_background, "real")
        time_step <- index - 1
        time_steps <- c(time_steps, time_step[1])

        vec <- dynamic_egenes_background[i,2:17]
        index <- which(vec==min(vec))
        pvalz <- c(pvalz, -log10(min(vec) + .000000000001))
        real_v_background <- c(real_v_background, "background")
        time_step <- index - 1
        time_steps <- c(time_steps, time_step[1])
    }

    print(length(time_steps))
    print(length(real_v_background))
    print(length(pvalz))

    df <- data.frame(pvalue=pvalz, time_step=time_steps,data_version=factor(real_v_background))
    boxplot <- ggplot(df, aes(x=data_version,y=pvalue,fill=data_version)) + geom_violin() + labs(x = "Data Version", y = "-log10(pvalue)") + scale_x_discrete(name="Data Version") 
    #boxplot <- boxplot + geom_dotplot(binaxis='y', stackdir='center', dotsize=.00000001) + geom_jitter(aes(colour=time_step),shape=16,size=.00000001, position=position_jitter(0.4))
    #boxplot <- boxplot + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")
    boxplot <- boxplot + theme(text = element_text(size=18)) + theme(legend.position="none")
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}

histogram_showing_nominal_time_steps <- function(dynamic_egenes, output_file, pvalue_threshold) {
    time_steps <- c()
    num_genes <- c()
    for (time_step in 0:15) {
        time_steps <- c(time_steps, time_step)
        num_hits <- sum(dynamic_egenes[,(time_step+2)] <= pvalue_threshold)
        num_genes <- c(num_genes, num_hits)
    }
    df <- data.frame(time_steps=time_steps, num_dynamic_qtls=num_genes)
    barplot <- ggplot(df, aes(time_steps, num_dynamic_qtls)) + geom_bar(stat = "identity",aes(fill=time_steps))
    barplot <- barplot + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    barplot <- barplot + labs(x = "time step", y = paste0("Dynamic qtls in this time step"),title = paste0("pvalue <= ",pvalue_threshold))
    barplot <- barplot + scale_fill_gradient(low="pink",high="blue")
    ggsave(barplot, file=output_file,width = 20,height=10.5,units="cm")

}

histogram_showing_number_nominal_time_steps <- function(dynamic_egenes, output_file,pvalue_threshold) {
    time_step_pvalz <- dynamic_egenes[,2:17]
    num_time_steps <- rowSums(time_step_pvalz <= pvalue_threshold)
    df <- data.frame(num_time_steps=num_time_steps)
    histo <- ggplot(data=df, aes(df$num_time_steps)) +
    geom_histogram(breaks=seq(-.5, 16.5, by = 1),col="red", fill="green",alpha = .2) +
            labs(x=paste0("Number of time steps dynamic eqtl is significant (<=", pvalue_threshold,") in"), y="Count")
    ggsave(histo, file=output_file, width=20, height=10.5, units="cm")
}




dynamic_egenes_comparison_file = args[1]
dynamic_egenes_background_comparison_file = args[2]
output_directory = args[3]


dynamic_egenes <- read.table(dynamic_egenes_comparison_file,header=TRUE)
dynamic_egenes_background <- read.table(dynamic_egenes_background_comparison_file, header=TRUE)




#################################################################################
# Histogram showing which time-steps have pvalue < $pvalue_threshold for our dynamic qtls 
##################################################################################
pvalue_threshold <- .05
output_file <- paste0(output_directory, "histogram_showing_time_step_independent_pvalue_",pvalue_threshold,"_hits_distribution_for_dynamic_qtls.png")
histogram_showing_nominal_time_steps(dynamic_egenes, output_file,pvalue_threshold)


#################################################################################
# Histogram showing number of time-steps have pvalue < $pvalue_threshold for a given dynamic qtl
##################################################################################

output_file <- paste0(output_directory, "histogram_number_of_time_steps_per_dynamic_egenes_",pvalue_threshold,".png")
histogram_showing_number_nominal_time_steps(dynamic_egenes, output_file,pvalue_threshold)




output_file <- paste0(output_directory, "dynamic_egenes_boxplot_comparing_time_step_independent_with_background.png")
boxplot_comparing_time_steps_with_background(dynamic_egenes, dynamic_egenes_background, output_file)


output_file <- paste0(output_directory, "dynamic_egenes_best_time_step_boxplot_cmp_to_background.png")
boxplot_with_best_time_step(dynamic_egenes, dynamic_egenes_background, output_file)


















####################
# OLD
####################





#################################################################################
# Scatter plot showing -log10 pvalues in dynamic and time step (colored by time)
##################################################################################


output_file <- paste0(output_directory, "all_tests_scatter_colored_by_time.png")
#scatter_colored_by_time(all_tests, output_file)

output_file <- paste0(output_directory, "dynamic_eqtl_tests_scatter_colored_by_time.png")
#scatter_colored_by_time(dynamic_eqtls, output_file)

output_file <- paste0(output_directory, "dynamic_eqtls_boxplot_comparing_time_step_independent.png")
#boxplot_comparing_time_steps(dynamic_eqtls, output_file)


output_file <- paste0(output_directory, "dynamic_egenes_boxplot_comparing_time_step_independent.png")
#boxplot_comparing_time_steps(dynamic_egenes, output_file)
