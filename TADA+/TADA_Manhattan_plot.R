#!/usr/bin/env Rscript

###################################################
# Generate the Manhattan plot for the TADA+ genes #
###################################################

suppressMessages(library(qqman))

args = commandArgs(TRUE)
TADA_results_filename = args[1]
new_genes = args[2]
outfile = args[3]

new_genes = as.vector(read.table(new_genes, header=FALSE)[,1])

TADA_results = read.table(TADA_results_filename, header=TRUE, sep="\t", check.names=FALSE)
TADA_results$chr = as.integer(gsub("chr", "", TADA_results$chr))

pdf(outfile, width=24, height=9)
margin = 10
par(mar = c(margin,margin,margin,margin))
manhattan(TADA_results, chr="chr", snp="gene", bp="start_hg38", p="qval_dnccPTV", ylab=expression(paste("-lo", g[10], "(Q)")), ylim=c(0,20), chrlabs=c("1", "", "3", "", "5", "", "7", "", "9", "", "11", "", "13", "", "15", "", "17", "", "19", "", "21", ""), annotatePval=0.1, annotateTop=FALSE, col=c("dodgerblue", "dodgerblue4"), suggestiveline=-log10(0.1), genomewideline=FALSE, cex.lab=1.5, cex.axis=1.3, highlight=new_genes)
abline(h=c(5,10,15,20), col="gray", lty=3)
invisible(dev.off())
