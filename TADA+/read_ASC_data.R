#!/usr/bin/env Rscript

##########################################################
# Read ASC data from Excel file, write data to TSV files #
##########################################################

suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))

args = commandArgs(TRUE)
ASC.variant_data.filename = args[1]
ASC.gene_data.filename = args[2]

ASC.sample_data = read.xlsx(ASC.variant_data.filename, sheet="Ped for all TADA samples") %>% select(Sample_ID, Sex)
ASC.sample_data$Sex = tolower(ASC.sample_data$Sex)
write.table(ASC.sample_data, file="ASC/ASC.sample_data.tsv", sep="\t", row.names=FALSE, quote=FALSE)

ASC.variant_data = read.xlsx(ASC.variant_data.filename, sheet="De novo variants") %>%
    select(Variant, Child_ID, Mom_ID, Dad_ID, Affected_Status, GENE_NAME, pLI, VEP_functional_class_canonical_simplified, MPC)
ASC.variant_data["Dataset"] = "ASC"
write.table(ASC.variant_data, file="ASC/ASC.variant_data.tsv", sep="\t", row.names=FALSE, quote=FALSE)

ASC.gene_data = read.xlsx(ASC.gene_data.filename, sheet="Autosomal", cols=1:27)
ASC.gene_data = ASC.gene_data[ASC.gene_data$entrez_id != ".",] # Remove genes with entrez_id = "."
ASC.gene_data = filter(ASC.gene_data, !duplicated(ASC.gene_data[,4])) # For some reason ASC has some genes with duplicate entrez_id, so remove all but the first one. There are only 12 and I manually checked them all, none are even close to being "ASD genes"
write.table(ASC.gene_data, file="ASC/ASC.gene_data.tsv", sep="\t", row.names=FALSE, quote=FALSE)
