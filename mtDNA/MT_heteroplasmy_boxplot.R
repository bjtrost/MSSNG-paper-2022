#!/usr/bin/env Rscript

#####################################################
# Generate boxplots for mtDNA heteroplasmy analysis #
#####################################################

suppressMessages(library(ggplot2))

args = commandArgs(TRUE)
MSSNG_filename = args[1]
SSC_filename = args[2]
output_filename = args[3]

MSSNG_data = read.table(MSSNG_filename, header=TRUE, sep=" ", check.names=FALSE)
SSC_data = read.table(SSC_filename, header=TRUE, sep=" ", check.names=FALSE)
MSSNG_data$dataset = "MSSNG"
SSC_data$dataset = "SSC"

MSSNG_data = MSSNG_data[c("delta", "status", "dataset")]
SSC_data = SSC_data[c("delta", "status", "dataset")]
all_data = rbind(MSSNG_data, SSC_data)
all_data[all_data == "proband"] = "Child with ASD"
all_data[all_data == "Proband"] = "Child with ASD"
all_data[all_data == "unaffected sibling"] = "Sibling without ASD"
all_data[all_data == "UnaffectedSibling"] = "Sibling without ASD"

theme = theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10)) +
        theme(axis.text.y=element_text(size=10)) +
        theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_line(colour = "lightgrey", size=0.15)) +
        theme(text=element_text(size=20)) +
        theme(axis.title=element_text(size=14)) +
        theme(legend.position = "none") +
        theme(plot.title = element_text(hjust = 0.5,size=16))

p = ggplot(all_data, aes(x=status, y=delta)) +
    geom_jitter(aes(color=status), alpha=0.2) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap("dataset") +
    theme +
    xlab("") +
    ylab("Delta")

suppressMessages(ggsave(output_filename, plot=p))
