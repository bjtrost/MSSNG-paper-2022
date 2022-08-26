#!/usr/bin/env Rscript

#################################################
# Generate plot for variant count distributions #
#################################################

suppressMessages(library(BTlib))
suppressMessages(library(extrafont))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
require(scales)

args = commandArgs(TRUE)
file = args[1]
outfile = args[2]

data = read.table(file, header=TRUE, check.names=FALSE, sep="\t")
data = subset(data, select=-c(`Damaging missense:LoF ratio`))

names(data)[names(data) == "Number of variants"] = "Variants"
names(data)[names(data) == "Number of rare variants"] = "Rare variants"
names(data)[names(data) == "Number of rare exonic variants"] = "Rare exonic variants"
names(data)[names(data) == "Number of LoF variants"] = "Loss of function variants"
names(data)[names(data) == "Number of de novo variants"] = "De novo variants"

if (grepl("indel", file, fixed=TRUE)) {
    data = subset(data, select=-c(`Number of damaging missense variants`))
} else {
    names(data)[names(data) == "Number of damaging missense variants"] = "Damaging missense variants"
}

data = melt(data, id.vars=c("Sample", "Category"))

if (grepl("indel", file, fixed=TRUE)) {
    data$variable = factor(data$variable, levels=c("Variants", "Rare variants", "Rare exonic variants", "Loss of function variants", "De novo variants"))
} else {
    data$variable = factor(data$variable, levels=c("Variants", "Rare variants", "Rare exonic variants", "Damaging missense variants", "Loss of function variants", "De novo variants"))
}

data$Category = factor(data$Category, levels=c("MSSNG.Complete Genomics.-", "MSSNG.Illumina HiSeq X.PCR-free", "MSSNG.Illumina HiSeq X.PCR-based", "MSSNG.Illumina HiSeq 2000.PCR-based", "MSSNG.Illumina HiSeq 2500.PCR-based", "SSC.Illumina HiSeq X.PCR-free", "SSC.Illumina HiSeq 2500.PCR-free", "1000G.Illumina NovaSeq 6000.PCR-free"))

text_size=25
x_text_size=11

to = 250
by = 50
plot_boundary = 3

p =    ggplot(data, aes_string(x="`Category`", y="`value`")) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=x_text_size)) +
        theme(plot.margin=unit(c(plot_boundary,plot_boundary,plot_boundary,plot_boundary),"cm")) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_line(colour = "lightgrey", size=0.15)) +
        facet_wrap("variable", scales="free") +
        theme(strip.text = element_text(size = 10)) +
        scale_y_continuous(labels = comma) +
        theme(text=element_text(size=text_size))

size=12
suppressMessages(ggsave(paste(outfile, ".pdf", sep=""), plot=p, height=size, width=size))
