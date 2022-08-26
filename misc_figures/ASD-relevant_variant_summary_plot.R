#!/usr/bin/env Rscript

#########################################################################################################
# Generate the barplot figuring giving the detailed breakdown of the ASD-associated variants discovered #
#########################################################################################################

suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(ggtext))

options(warn=-1)

args = commandArgs(TRUE)
SNV_indel_counts_filename = args[1]
recessive_counts_filename = args[2]
SV_counts_filename = args[3]
UPD_counts_filename = args[4]
TRE_counts_filename = args[5]
MT_counts_filename = args[6]

MSSNG_border_size=0.3

theme = theme(axis.text.x=element_text(vjust=1, size=6)) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_blank()) +
        theme(panel.grid.major.x = element_line(colour = "lightgrey", size=0.15)) +
        theme(legend.title=element_blank()) +
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))

LoF_color = rgb(red=245, green=194, blue=66, alpha=255, maxColorValue=255)
missense_color = rgb(red=125, green=146, blue=72, alpha=255, maxColorValue=255)
deletion_color = rgb(red=179, green=87, blue=81, maxColorValue=255)
duplication_color = rgb(red=90, green=128, blue=184, maxColorValue=255)
insertion_color = "green"
inversion_color = rgb(red=78, green=157, blue=158, maxColorValue=255) # Same as IGV (light blue)
complex_color = rgb(red=124, green=102, blue=158, maxColorValue=255)

alpha_vals=c(1,0.7)
plot_width=2
legend_box_margin=margin(-8,-8,-8,-8)
legend_justification="center"
legend_direction = "vertical"
legend_spacing_x = unit(0.1, "cm")

################ SNVs+indels plot ################
cat("Now generating SNV/indel plot...\n")
counts = read.table(SNV_indel_counts_filename, sep="\t", header=TRUE)
#counts = counts[counts$inheritance_category == "de novo",]
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts = counts[order(counts$Dataset),]
counts$variant_category = factor(counts$variant_category, levels=c("LoF", "missense"))
SNV_indel_plot = ggplot(counts, aes(x=reorder(gene, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=variant_category, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_fill_manual(values=c(LoF_color, missense_color)) +
    scale_y_continuous(name=NULL, breaks=seq(0,20,by=2)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(face="italic", size=4, margin=margin(-5,-5,-5,-5))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("SNVs+indels.pdf", plot=SNV_indel_plot, width=plot_width, height=6))

################ Recessive plot ################
cat("Now generating recessive plot...\n")
counts = read.table(recessive_counts_filename, sep="\t", header=TRUE)
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts = counts[order(counts$Dataset),]
counts$variant_category = factor(counts$variant_category, levels=c("LoF", "missense"))
recessive_plot = ggplot(counts, aes(x=reorder(gene, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=variant_category, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_fill_manual(values=c(LoF_color, missense_color)) +
    scale_y_continuous(name=NULL, breaks=seq(0,20,by=2)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(face="italic", size=4, margin=margin(-5,-5,-5,-5))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("recessive.pdf", plot=recessive_plot, width=plot_width, height=6))

################ Chromosomal abnormalities plot ################
cat("Now generating chromosomal abnormalities plot...\n")
counts = read.table(SV_counts_filename, sep="\t", header=TRUE)
counts = counts[counts$Class == "Chromosomal abnormality",]
counts$Type = factor(counts$Type, levels=c("DEL", "DUP", "CPX"))
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts = counts[order(counts$Dataset),]
chromosomal_abnormalities_plot = ggplot(counts, aes(x=reorder(Description, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=Type, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_fill_manual(values=c(deletion_color, duplication_color, complex_color)) +
    scale_y_continuous(name=NULL, breaks=seq(0,20,by=2)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(size=4, margin=margin(-3,-3,-3,-3))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("chromosomal_abnormalities.pdf", plot=chromosomal_abnormalities_plot, width=plot_width, height=1.5))

################ Genomic disorders plot ################
cat("Now generating genomic disorders plot...\n")
counts = read.table(SV_counts_filename, sep="\t", header=TRUE)
counts = counts[counts$Class == "Genomic disorder",]
counts$Type = factor(counts$Type, levels=c("DEL", "DUP"))
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts = counts[order(counts$Dataset),]
genomic_disorders_plot = ggplot(counts, aes(x=reorder(Description, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=Type, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_fill_manual(values=c(deletion_color, duplication_color)) +
    scale_y_continuous(name=NULL, limits=c(0,35), breaks=seq(0,35,by=5)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(size=4, margin=margin(-3,-3,-3,-3))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("genomic_disorders.pdf", plot=genomic_disorders_plot, width=plot_width, height=2.2))

################ Large or gene-rich CNVs plot ################
cat("Now generating large/gene-rich CNVs plot...\n")
counts = read.table(SV_counts_filename, sep="\t", header=TRUE)
counts = counts[counts$Class == "Large or gene-rich CNV" & counts$Description != "",]
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts$Type = factor(counts$Type, levels=c("DEL", "DUP"))

large_gene_rich_plot = ggplot(counts, aes(x=reorder(Description, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=Type, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    theme(legend.title=element_blank()) +
    scale_fill_manual(values=c(deletion_color, duplication_color)) +
    scale_y_continuous(name=NULL, breaks=seq(0,10)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(size=4, margin=margin(-3,-3,-3,-3))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("large_gene_rich.pdf", plot=large_gene_rich_plot, width=plot_width, height=2.2))

################ ASD-risk SV plot ################
cat("Now generating ASD-risk SVs plot...\n")
counts = read.table(SV_counts_filename, sep="\t", header=TRUE)
counts = counts[counts$Class == "SV disrupting ASD-associated gene",]
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts$Type[counts$Type != "DEL" & counts$Type != "DUP" & counts$Type != "INS" & counts$Type != "INV"] = "CPX"
counts$Type = factor(counts$Type, levels=c("DEL", "DUP", "INS", "INV", "CPX"))

ASD_risk_SV_plot = ggplot(counts, aes(x=reorder(Description, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=Type, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    theme(legend.title=element_blank()) +
    #theme(axis.text.y = element_markdown(size=4, margin=margin(-3,-3,-3,-3))) +
    theme(axis.text.y = element_text(face="italic", size=4, margin=margin(-3,-3,-3,-3))) +
    scale_fill_manual(values=c(deletion_color, duplication_color, insertion_color, inversion_color, complex_color)) +
    scale_y_continuous(name=NULL, breaks=seq(0,25,by=5)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("ASD_risk_SVs.pdf", plot=ASD_risk_SV_plot, width=plot_width, height=6))

################ UPD plot ################
cat("Now generating UPD plot...\n")
counts = read.table(UPD_counts_filename, sep="\t", header=TRUE)
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts$Type = factor(counts$Type, levels=rev(c("chr2 (complete)", "chr4 (complete)", "chr8 (complete)", "chr1 (partial)", "chr8 (partial)")))

UPD_plot = ggplot(counts, aes(x=Type, y=Count)) + theme +
    geom_col(aes(fill=Type, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_fill_manual(values=c(complex_color, complex_color, complex_color, complex_color, complex_color)) +
    theme(legend.title=element_blank()) + theme(axis.text.y = element_markdown(size=4, margin=margin(-3,-3,-3,-3))) +
    scale_y_continuous(name=NULL, breaks=seq(0,5,by=1), limits=c(0,2)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=FALSE, alpha=guide_legend(order=0))
suppressMessages(ggsave("UPDs.pdf", plot=UPD_plot, width=plot_width, height=1))

################ TRE plot ################
cat("Now generating TREs plot...\n")
counts = read.table(TRE_counts_filename, sep="\t", header=TRUE)
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts = counts[order(counts$Dataset),]
counts$gene_set = factor(counts$gene_set, levels=c("Nervous system development", "Cardiovascular system or muscle", "Both"))
TRE_plot = ggplot(counts, aes(x=reorder(gene, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=gene_set, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_y_continuous(name=NULL, breaks=seq(0,20,by=2)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(face="italic", size=4, margin=margin(-5,-5,-5,-5))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("TREs.pdf", plot=TRE_plot, width=plot_width, height=6))

################ MT plot ################
cat("Now generating MT plot...\n")
counts = read.table(MT_counts_filename, sep="\t", header=TRUE)
counts$Dataset = factor(counts$Dataset, levels=c("MSSNG", "SSC"))
counts = counts[order(counts$Dataset),]
MT_plot = ggplot(counts, aes(x=reorder(gene, total_for_category), y=Count)) + theme +
    geom_col(aes(fill=Disorder, alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    #geom_col(aes(alpha=Dataset), position=position_stack(reverse=TRUE), width=0.65) +
    scale_y_continuous(name=NULL, breaks=seq(0,20,by=2)) +
    theme(legend.position="top", legend.key.size=unit(0.15, "cm"), legend.text=element_text(size=5.), legend.justification=legend_justification, legend.direction=legend_direction, legend.margin=margin(0,0,0,0), legend.box.margin=legend_box_margin, legend.spacing.x=legend_spacing_x, axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(face="italic", size=4, margin=margin(-5,-5,-5,-5))) +
    scale_alpha_discrete(range=alpha_vals) + scale_x_discrete(name=NULL) + coord_flip() +
    guides(fill=guide_legend(order=1), alpha=guide_legend(order=0))
suppressMessages(ggsave("MT.pdf", plot=MT_plot, width=plot_width, height=1))


#blank = ggplot() + theme_void()
#A = ggarrange(SNV_indel_plot, labels="A")
#BCD = ggarrange(chromosomal_abnormalities_plot, genomic_disorders_plot, large_gene_rich_plot, nrow=3, labels=c("B","C","D"), heights=c(3.5,4.5,4))
#EF = ggarrange(UPD_plot, ASD_risk_SV_plot, nrow=2, labels=c("E", "F"), heights=c(1.2,10.5), label.y=c(1.05, 1))
#GH = ggarrange(TRE_plot, MT_plot, blank, nrow=3, heights=c(5,1.3,2), labels=c("G","H",""), label.x=c(-0.05, 0))
#all_plot = ggarrange(A, BCD, EF, GH, ncol=4)
#suppressMessages(ggsave("ASD_relevant_variants.pdf", plot=all_plot, width=7, height=8))

blank = ggplot() + theme_void()
A = ggarrange(SNV_indel_plot, labels="A")
BCD = ggarrange(recessive_plot, chromosomal_abnormalities_plot, genomic_disorders_plot, nrow=3, labels=c("B","C","D"), heights=c(4.5,2.5,2.7))
EF = ggarrange(large_gene_rich_plot, ASD_risk_SV_plot, nrow=2, labels=c("E","F"), heights=c(3.7,6.5))
GHI = ggarrange(UPD_plot, TRE_plot, MT_plot, blank, nrow=4, heights=c(1,5.1,1.6,1.7), labels=c("G","H","I"))
all_plot = ggarrange(A, BCD, EF, GHI, ncol=4)

#BCD = ggarrange(recessive_plot, chromosomal_abnormalities_plot, genomic_disorders_plot, nrow=3, labels=c("B","C","D"), heights=c(3.5,4.5,4))
#EF = ggarrange(UPD_plot, ASD_risk_SV_plot, nrow=2, labels=c("E", "F"), heights=c(1.2,10.5), label.y=c(1.05, 1))
#GH = ggarrange(TRE_plot, MT_plot, blank, nrow=3, heights=c(5,1.3,2), labels=c("G","H",""), label.x=c(-0.05, 0))
#all_plot = ggarrange(A, BCD, EF, GH, ncol=4)
suppressMessages(ggsave("ASD_relevant_variants.pdf", plot=all_plot, width=8.5, height=8))
