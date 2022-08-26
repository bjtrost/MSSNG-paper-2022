#!/usr/bin/env Rscript

#############################################################
# Generate the pLI distribution plot of the new TADA+ genes #
#############################################################

suppressMessages(library(tidyverse))

args = commandArgs(TRUE)
filename = args[1]
output_filename = args[2]

theme = theme(axis.text.x=element_text(vjust=1, size=20)) + ### Define ggplot theme###
        theme(axis.text.y=element_text(size=20)) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_line(colour = "lightgrey", size=0.15)) +
        theme(axis.title=element_text(size=20)) +
        theme(legend.position = "none") +
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))

###########################################
### Plot overall distribution of points ###
###########################################


mat = read.table(filename, sep="\t", header=TRUE)

bins <- 5
cols <- c("forestgreen","red1")
colGradient <- colorRampPalette(cols)
cut.cols <- colGradient(bins)
cuts <- cut(mat$pLI,bins)
names(cuts) <- sapply(cuts,function(t) cut.cols[which(as.character(t) == levels(cuts))])

density = ggplot(mat, aes(x=pLI, fill=cut(pLI, bins))) +
          geom_histogram() + theme + ylab("Number of genes") +
          scale_color_manual(values=cut.cols,labels=levels(cuts)) +
          scale_fill_manual(values=cut.cols,labels=levels(cuts))

suppressMessages(ggsave(output_filename, plot=density))
