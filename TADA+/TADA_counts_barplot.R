#!/usr/bin/env Rscript

#####################################################################
# Generate the barplot showing the evidence for the new TADA+ genes #
#####################################################################

suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(openxlsx))

args = commandArgs(TRUE)
counts_filename = args[1]
outfile_new = args[2]
outfile_ASC102 = args[3]

theme = theme(axis.text.x=element_text(vjust=1, size=16)) +
        theme(axis.text.y=element_text(size=10)) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_blank()) +
        theme(panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_line(colour = "lightgrey", size=0.15)) +
        theme(axis.title=element_text(size=16)) +
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))


counts = read.xlsx(counts_filename, check.names=TRUE)

counts$cc_diff.ptv = counts$case.ptv - counts$control.ptv
counts$total = counts$dn.ptv + counts$dn.misb + counts$dn.misa + counts$cc_diff.ptv
colnames(counts)[which(names(counts) == "dn.ptv")] <- "De novo PTV"
colnames(counts)[which(names(counts) == "dn.misb")] <- "De novo misB"
colnames(counts)[which(names(counts) == "dn.misa")] <- "De novo misA"
colnames(counts)[which(names(counts) == "cc_diff.ptv")] <- "Case-control difference"
counts_melted = melt(counts, id.vars=c("total", "gene", "FDR", "case.ptv", "control.ptv", "ASC102", "SFARI", "EAGLE", "pLI", "X..of.clinical.labs.offering.ASD.panel.with.this.gene", "ClinGen.Disease", "ClinGen.Gene.Disease.Classification", "Comments.MR", "OMIM.P"))
counts_melted$variable = factor(counts_melted$variable, levels=c("De novo PTV", "De novo misB", "De novo misA", "Case-control difference"))


ptv = "firebrick1"
misb = "darkblue"
misa = "lightskyblue"
cc = "forestgreen"

dn.ptv = expression(paste(italic("De novo "), "PTV", sep=""))
dn.misA = expression(paste(italic("De novo "), "MisA", sep=""))
dn.misB = expression(paste(italic("De novo "), "MisB", sep=""))
cc_diff = "Case-control difference"

p = ggplot(counts_melted[!counts_melted$ASC102,], aes(x=reorder(gene, -total), y=value)) + theme +
    geom_col(aes(fill=variable), position=position_stack(reverse = TRUE)) +
    coord_cartesian(ylim=c(-3,15)) +
    xlab("") + ylab("Number of variants") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(face="italic", size=12, angle=90, vjust=0.5, hjust=1)) +
    scale_fill_manual(values=c(ptv, misb, misa, cc), guide=guide_legend(reverse=TRUE), labels=c(dn.ptv, dn.misB, dn.misA, cc_diff)) + theme(legend.text.align=0)

suppressMessages(ggsave(outfile_new, plot=p, width=16, height=6))

p2 = ggplot(counts_melted[counts_melted$ASC102,], aes(x=reorder(gene, -total), y=value)) + theme +
    geom_col(aes(fill=variable), position=position_stack(reverse = TRUE)) +
    #coord_cartesian(ylim=c(-3,15)) +
    xlab("") + ylab("Number of variants") + theme(legend.title=element_blank()) + theme(axis.text.x=element_text(face="italic", size=12, angle=90, vjust=0.5, hjust=1)) +
    scale_fill_manual(values=c(ptv, misb, misa, cc), guide=guide_legend(reverse=TRUE))

suppressMessages(ggsave(outfile_ASC102, plot=p2, width=16, height=6))
