#!/usr/bin/env Rscript

###############################################
# Wrapper script to perform the TADA analysis #
# PLEASE NOTE that the required libraries     #
# TADA_v2.R and gamma_pli.R are not provided, #
# as they were not written by us. We kindly   #
# received them from Dr. Kyle Satterstrom     #
###############################################

source(paste(Sys.getenv("SCRIPTS"), "TADA_v2.R", sep="/"))
source(paste(Sys.getenv("SCRIPTS"), "gamma_pli.R", sep="/"))
suppressMessages(library(investr))
options(stringsAsFactors=F)

args = commandArgs(TRUE)
variant_data_filename = args[1]
gene_data_filename = args[2]
n.trio = as.numeric(args[3])
n.case = as.numeric(args[4])

n.cc.dbs = data.frame(ca = n.case, cn = 5214)

variant_data = read.table(variant_data_filename, sep="\t", header=TRUE)
gene_data = read.table(gene_data_filename, sep="\t", header=TRUE)

TADA_results = tada2(variant_data, gene_data, n.trio=n.trio, n.cc.dbs=n.cc.dbs)
rownames(TADA_results) = gene_data$gene
TADA_results = TADA_results[!is.infinite(TADA_results$BF_dn_misa),]
TADA_results = TADA_results[!is.infinite(TADA_results$BF_dn_misb),]
TADA_results$qval_dnccPTV = Bayesian.FDR(apply(TADA_results[,1:4], 1, prod), pi0=0.95)
TADA_results = cbind(gene = rownames(TADA_results), TADA_results)
write.table(TADA_results, stdout(), sep="\t", quote=FALSE, row.names=FALSE)
invisible(file.remove("Rplots.pdf"))
