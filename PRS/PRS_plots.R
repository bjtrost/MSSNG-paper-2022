#!/usr/bin/env Rscript

###############################################
# Perform the PRS analsyis and generate plots #
###############################################

suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(glue))
suppressMessages(library(grid))
suppressMessages(library(ggplotify))
suppressMessages(library(stringr))

args = commandArgs(TRUE)
PRS_data_filename = args[1]
metadata_filename = args[2]
clinically_significant_CNVs_filename = args[3]
LoF_filename = args[4]
variant_phenotype_table = args[5]

#######################################################
### Create directories for PDF and TSV output files ###
#######################################################
pdf_dir = "pdf"
tsv_dir = "tsv"
dirs = c(pdf_dir, tsv_dir)
for (d in dirs) {
    if (!dir.exists(d)) {
        dir.create(d)
    }
}

########################
### Define constants ###
########################
tip.length = 0.015
stat.label = "p = {sprintf('%.1e', p.adj)}"
stat.label.noadj = "p = {sprintf('%.1e', p)}"
ylab_PRS = "PRS"
PRS_breaks = c(-20, -10, 0, 10, 20)
# Colors: MSSNG = green, SSC = blue, 1000G = red
# Affected = dark, Unaffected = light
MSSNG_affected = rgb(83, 182, 75, maxColorValue=255)
MSSNG_unaffected = "darkseagreen3"
SSC_affected = rgb(110, 157, 248, maxColorValue=255)
SSC_unaffected = rgb(85, 188, 195, maxColorValue=255)
OneKG = rgb(232, 125, 114, maxColorValue=255)
colors = c(
"MSSNG" = MSSNG_affected,
"SSC" = SSC_affected,
"MSSNG\n(with ASD)" = MSSNG_affected,
"MSSNG\n(without ASD)" = MSSNG_unaffected,
"MSSNG parent" = MSSNG_unaffected,
"MSSNG father" = MSSNG_unaffected,
"MSSNG mother" = MSSNG_unaffected,
"MSSNG male\n(with ASD)" = MSSNG_affected,
"MSSNG female\n(with ASD)" = MSSNG_affected,
"MSSNG male\n(without ASD)" = MSSNG_unaffected,
"MSSNG female\n(without ASD)" = MSSNG_unaffected,
"SSC\n(with ASD)" = SSC_affected,
"SSC\n(without ASD)" = SSC_unaffected,
"SSC parent" = SSC_unaffected,
"SSC father" = SSC_unaffected,
"SSC mother" = SSC_unaffected,
"SSC male\n(with ASD)" = SSC_affected,
"SSC female\n(with ASD)" = SSC_affected,
"SSC male\n(without ASD)" = SSC_unaffected,
"SSC female\n(without ASD)" = SSC_unaffected,
"female (with ASD)\nwith female\nsibling (without ASD)" = SSC_affected,
"female (with ASD)\nwith male\nsibling (without ASD)" = SSC_affected,
"male (with ASD)\nwith female\nsibling (without ASD)" = SSC_affected,
"male (with ASD)\nwith male\nsibling (without ASD)" = SSC_affected,
"1000G" = OneKG,
"1000G " = OneKG,
"Multiplex" = MSSNG_affected,
"Simplex" = SSC_affected
)

theme = theme(axis.text.x=element_text(vjust=1, size=10)) + ### Define ggplot theme###
        theme(axis.text.y=element_text(size=10)) +
        theme(panel.background=element_rect(fill = "white")) +
        theme(panel.border=element_rect(color="black", fill=NA)) +
        theme(panel.grid.major.y = element_line(colour = "lightgrey", size=0.15)) +
        theme(axis.title=element_text(size=14)) +
        theme(legend.position = "none") +
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))

#######################################################
### Read in PRS data and merge with sample metadata ###
#######################################################
PRS_data = read.xlsx(PRS_data_filename, rowNames=TRUE)
PRS_data[,"Sample ID"] = rownames(PRS_data)
metadata = read.table(metadata_filename, header=TRUE, sep="\t", check.names=FALSE)
rownames(metadata) = metadata[,"Sample ID"]
merged = merge(PRS_data, metadata, by="Sample ID")
families_with_non_sibling_affected = unique(merged[merged$Relation != "proband" & merged$Relation != "affected sibling" & merged$Relation != "unaffected sibling" & merged$Affection == 2,"Family ID"])
merged = merged[!merged[,"Family ID"] %in% families_with_non_sibling_affected,]

colnames(merged) = make.names(colnames(merged)) # Replace spaces with periods in column names
merged$Dataset = factor(merged$Dataset, levels=c("MSSNG", "SSC", "1000G"))

clinically_significant_CNVs = read.xlsx(clinically_significant_CNVs_filename)
merged$clinically_significant_CNV_bool = merged$Sample.ID %in% clinically_significant_CNVs$sample
merged$clinically_significant_CNV = ifelse(merged$Sample.ID %in% clinically_significant_CNVs$sample, "CNV=yes", "CNV=no")

LoF = read.table(LoF_filename, header=TRUE, sep="\t", check.names=FALSE, comment.char="")
merged$LoF = ifelse(merged$Sample.ID %in% LoF$Sample, "LoF=yes", "LoF=no")

summary(merged$PRS_1)
sd(merged$PRS_1)

######################################################
### Generate derived columns for merged data frame ###
######################################################
affection_map = function(affection) {
    if (affection == "0") {
        return("affected - other")
    } else if (affection == "1") {
        return("(without ASD)")
    } else if (affection == "2") {
        return("(with ASD)")
    }
}

## Generate derived columns ##
parent_map = function(relation) {
    if (relation == "mother" || relation == "father") {
        return("parent")
    } else if (relation == "-") {
        return("")
    } else {
        return("non-parent")
    }
}

merged$Affection_text                                   = sapply(merged$Affection, affection_map)
merged$Parent                                           = sapply(merged$Relation, parent_map)
merged$Dataset.affection                                = paste(merged$Dataset, merged$Affection_text, sep="\n") # Make new column containing dataset and affection status; e.g. "MSSNG affected"
merged$Dataset.affection.sex                            = paste(paste(merged$Dataset, merged$Sex, sep=" "), merged$Affection_text, sep="\n") # Make new column containing dataset and affection status and sex; e.g. "MSSNG affected male"
merged$Dataset.relation                                 = paste(merged$Dataset, merged$Relation, sep=" ") # Make new column containing dataset and relation; e.g. MSSNG father
merged$Dataset.parent                                   = str_trim(paste(merged$Dataset, merged$Parent, sep=" "))
merged$is_affected                                      = merged$Affection == "2"
merged$is_1000G                                         = merged$Dataset == "1000G"
merged$Dataset.CNV                                      = paste(merged$Dataset, merged$clinically_significant_CNV, sep=" ")
merged$Dataset.LoF                                      = paste(merged$Dataset, merged$LoF, sep=" ")

merged$PRS_rank[order(merged$PRS_1, decreasing=TRUE)]   = 1:nrow(merged)

##########################################################
### Compare affected versus unaffected siblings in SSC ###
##########################################################
SSC_affected = merged[merged$Affection == "2" & merged$Dataset == "SSC",]
SSC_affected$sibling_PRS = NA
for (i in 1:nrow(SSC_affected)) {
    unaffected_siblings = merged[merged$Family.ID == SSC_affected[i,]$Family.ID & merged$Relation == "unaffected sibling",]
    if (nrow(unaffected_siblings) == 1) {
        SSC_affected[i,]$sibling_PRS = unaffected_siblings[1,]$PRS_1
    }
}
SSC_affected = SSC_affected[!is.na(SSC_affected$sibling_PRS),] %>% select(Sample.ID, PRS_1, sibling_PRS)
for (sd_mult in c(1,2)) {
    num_affected = nrow(SSC_affected[SSC_affected$PRS_1 > SSC_affected$sibling_PRS + sd_mult * sd(merged$PRS_1),])
    num_unaffected = nrow(SSC_affected[SSC_affected$sibling_PRS > SSC_affected$PRS_1 + sd_mult * sd(merged$PRS_1),])
    binom.test(c(num_affected, num_unaffected), alternative="greater")
}

###################################################################################################################################
########################################################### Begin plots ###########################################################
###################################################################################################################################

###########################################
### Plot overall distribution of points ###
###########################################
filename = "PRS.rank_dist"
cat("Generating plot", filename, "...\n")
PRS.rank_dist = ggplot(merged, aes(x=PRS_rank, y=PRS_1)) +
                geom_point() + theme + xlab("Rank") + ylab(ylab_PRS)
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.rank_dist))

####################
### Density plot ###
####################
filename = "PRS.density"
cat("Generating plot", filename, "...\n")
PRS.density = ggplot(merged, aes(x=PRS_1, color=Dataset)) +
              geom_density() + theme + xlab(ylab_PRS) + ylab("Density") + scale_colour_manual(values=colors) +
              theme(legend.position = "top", legend.text = element_text(size=8), legend.title=element_blank())
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.density))

##########################################################
### MSSNG affected/unaffected, SSC affected/unaffected ###
##########################################################
filename = "PRS.affected_versus_unaffected_children"
cat("Generating plot", filename, "...\n")
children = merged[merged$Affection != "0" & merged$Relation %in% c("proband", "affected sibling", "unaffected sibling"),]

MSSNG_affected_versus_MSSNG_unaffected_pval = coef(summary(glm(is_affected ~ Sex + PRS_1, data=children[children$Dataset == "MSSNG",])))[3,4]
SSC_affected_versus_SSC_unaffected_pval = coef(summary(glm(is_affected ~ Sex + PRS_1, data=children[children$Dataset == "SSC",])))[3,4]
MSSNG_affected_versus_SSC_unaffected_pval = coef(summary(glm(is_affected ~ Sex + PRS_1, data=children[(children$Dataset == "MSSNG" & children$Affection == "2") | (children$Dataset == "SSC" & children$Affection == "1"),])))[3,4]
SSC_affected_versus_MSSNG_unaffected_pval = coef(summary(glm(is_affected ~ Sex + PRS_1, data=children[(children$Dataset == "SSC" & children$Affection == "2") | (children$Dataset == "MSSNG" & children$Affection == "1"),])))[3,4]
MSSNG_affected_versus_SSC_affected_pval = coef(summary(glm(is_affected ~ Sex + PRS_1, data=children[children$Affection == "2",])))[3,4]
p = c(MSSNG_affected_versus_MSSNG_unaffected_pval, SSC_affected_versus_SSC_unaffected_pval, MSSNG_affected_versus_SSC_unaffected_pval, SSC_affected_versus_MSSNG_unaffected_pval, MSSNG_affected_versus_SSC_affected_pval)
p.adj = p.adjust(p, method="holm")
stat.test = data.frame(group1=c("MSSNG\n(with ASD)", "SSC\n(with ASD)", "MSSNG\n(with ASD)", "SSC\n(with ASD)", "MSSNG\n(with ASD)"),
                       group2=c("MSSNG\n(without ASD)", "SSC\n(without ASD)", "SSC\n(without ASD)", "MSSNG\n(without ASD)", "SSC\n(with ASD)"),
                       p=p,
                       p.adj=p.adj,
                       y.position=c(NA,18,22,NA,NA))

PRS.affected_versus_unaffected_children = ggplot(children, aes(x=Dataset.affection, y=PRS_1)) +
                                          geom_jitter(aes(color=Dataset.affection, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(ylab_PRS) + theme(axis.text.x=element_text(vjust=1, size=9)) +
                                          scale_colour_manual(values=colors) + scale_y_continuous(limits=c(NA, 27), breaks=PRS_breaks) +
                                          stat_pvalue_manual(stat.test, label=stat.label, hide.ns=TRUE, tip.length=tip.length)
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.affected_versus_unaffected_children))

##############################
### Mothers versus fathers ###
##############################
filename = "PRS.mothers_versus_fathers"
cat("Generating plot", filename, "...\n")
parents = merged[merged$Affection == "1" & merged$Relation %in% c("father", "mother"),] # Get only unaffected parents
stat.test = parents %>% t_test(PRS_1 ~ Dataset.relation) %>% add_xy_position(x="Dataset.relation")
write.table(stat.test[,!names(stat.test) == "groups"], file=paste(tsv_dir, "/", filename, ".tsv", sep=""))

PRS.mothers_versus_fathers = ggplot(parents, aes(x=str_wrap(Dataset.relation,7), y=PRS_1)) +
                             geom_jitter(aes(color=Dataset.relation, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(ylab_PRS) +
                             scale_colour_manual(values=colors) + stat_pvalue_manual(stat.test, label=stat.label, hide.ns=TRUE, tip.length=tip.length)
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.mothers_versus_fathers))


###################################
### Clinically significant CNVs ###
###################################
filename = "PRS.clinically_significant_CNVs"
cat("Generating plot", filename, "...\n")
affected = merged[merged$Affection == "2",]
MSSNG_clinsig_CNV_pval = coef(summary(glm(clinically_significant_CNV_bool ~ Sex + PRS_1, data=affected[affected$Dataset == "MSSNG",])))[3,4]
SSC_clinsig_CNV_pval = coef(summary(glm(clinically_significant_CNV_bool ~ Sex + PRS_1, data=affected[affected$Dataset == "SSC",])))[3,4]
MSSNG_clinsig_CNV_pval
SSC_clinsig_CNV_pval

p = c(MSSNG_clinsig_CNV_pval, SSC_clinsig_CNV_pval)
p.adj = p.adjust(p, method="holm")
stat.test = data.frame(group1=c("MSSNG\nCNV=no", "SSC\nCNV=no"),
                       group2=c("MSSNG\nCNV=yes", "SSC\nCNV=yes"),
                       p=p,
                       p.adj=p.adj,
                       y.position=c(18,22))

PRS.clinically_significant_CNVs = ggplot(affected, aes(x=str_wrap(Dataset.CNV, 8), y=PRS_1)) +
                                  geom_jitter(aes(color=Dataset.CNV, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(ylab_PRS) +
                                  stat_pvalue_manual(stat.test, label=stat.label, hide.ns=TRUE, tip.length=tip.length)
                                  #scale_colour_manual(values=colors) +
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.clinically_significant_CNVs))


###################################
### LoF ####
###################################
filename = "PRS.LoF"
cat("Generating plot", filename, "...\n")
affected = merged[merged$Affection == "2",]
stat.test = affected %>% t_test(PRS_1 ~ Dataset.LoF) %>% add_xy_position(x="Dataset.LoF")
write.table(stat.test[,!names(stat.test) == "groups"], file=paste(tsv_dir, "/", filename, ".tsv", sep=""))

PRS.LoF = ggplot(affected, aes(x=str_wrap(Dataset.LoF, 8), y=PRS_1)) +
                                  geom_jitter(aes(color=Dataset.LoF, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(ylab_PRS) +
                                  stat_pvalue_manual(stat.test, label=stat.label, hide.ns=TRUE, tip.length=tip.length)
                                  #scale_colour_manual(values=colors) +
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.LoF))



#####################################################
### MSSNG parents versus SSC parents versus 1000G ###
#####################################################
filename = "PRS.1000G_versus_parents"
cat("Generating plot", filename, "...\n")
parents_plus_1000G = merged[merged$Affection == "1" & (merged$Parent == "parent" | merged$Relation == "-"),] # Get only unaffected parents plus 1000G
parents_plus_1000G$Dataset.parent = factor(parents_plus_1000G$Dataset.parent, levels=c("MSSNG parent", "SSC parent", "1000G"))


OneKG_versus_MSSNG_parent_pval = coef(summary(glm(is_1000G ~ Sex + PRS_1, data=parents_plus_1000G[parents_plus_1000G$Dataset == "1000G" | parents_plus_1000G$Dataset == "MSSNG",])))[3,4]
OneKG_versus_SSC_parent_pval = coef(summary(glm(is_1000G ~ Sex + PRS_1, data=parents_plus_1000G[parents_plus_1000G$Dataset == "1000G" | parents_plus_1000G$Dataset == "SSC",])))[3,4]
p = c(OneKG_versus_MSSNG_parent_pval, OneKG_versus_SSC_parent_pval)
p.adj = p.adjust(p, method="holm")
stat.test = data.frame(group1=c("1000G", "1000G"),
                       group2=c("MSSNG\nparent", "SSC\nparent"),
                       p=p,
                       p.adj=p.adj,
                       y.position=c(18,22))

PRS.parents_plus_1000G = ggplot(parents_plus_1000G, aes(x=str_wrap(Dataset.parent,7), y=PRS_1)) +
                         geom_jitter(aes(color=Dataset.parent, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(ylab_PRS) +
                         scale_colour_manual(values=colors) + stat_pvalue_manual(stat.test, label=stat.label, hide.ns=TRUE, tip.length=tip.length)
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.parents_plus_1000G))


############
### pTDT ###
############
children$Parent_mean = NA
for (i in 1:nrow(children)) {
    parents <- merged[merged$Family.ID == children[i,]$Family.ID & merged$Relation %in% c("father", "mother"),] # Get the parents
    if (nrow(parents) == 2) { # If we have both parents...
        children$Parent_mean[i] = mean(parents$PRS_1)
    }
}
children_pTDT = children[!is.na(children$Parent_mean),]
children_pTDT$Dataset.affection= factor(children_pTDT$Dataset.affection, levels=c("MSSNG\n(with ASD)", "MSSNG\n(without ASD)", "SSC\n(with ASD)", "SSC\n(without ASD)"))
children_pTDT$overtransmission = children_pTDT$PRS_1 - children_pTDT$Parent_mean

t.MSSNG.affected <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "MSSNG" & children_pTDT$Affection == "2"], alternative = "greater", mu = 0)
t.MSSNG.unaffected <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "MSSNG" & children_pTDT$Affection == "1"], alternative = "greater", mu = 0)
t.SSC.affected <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "SSC" & children_pTDT$Affection == "2"], alternative = "greater", mu = 0)
t.SSC.unaffected <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "SSC" & children_pTDT$Affection == "1"], alternative = "greater", mu = 0)
adjusted_pvalues = p.adjust(c(t.MSSNG.affected$p.value, t.MSSNG.unaffected$p.value, t.SSC.affected$p.value, t.SSC.unaffected$p.value), method="holm")

filename = "PRS.pTDT"
cat("Generating plot", filename, "...\n")
PRS.pTDT = ggplot(children_pTDT, aes(x=Dataset.affection, y=overtransmission)) +
           geom_jitter(aes(color=Dataset.affection, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(paste(ylab_PRS, "overtransmission", sep=" ")) + theme(axis.text.x=element_text(vjust=1, size=9)) +
           scale_colour_manual(values=colors) + annotate("text", x=c(1,3), y=c(15), label=paste("p = ",c(signif(adjusted_pvalues[1], digits = 2), signif(adjusted_pvalues[3], digits = 2)), sep=""))
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.pTDT))


#############################
### pTDT - sex stratified ###
#############################

t.MSSNG.affected.male <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "MSSNG" & children_pTDT$Affection == "2" & children_pTDT$Sex == "male"], alternative = "greater", mu = 0)
t.SSC.affected.male <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "SSC" & children_pTDT$Affection == "2" & children_pTDT$Sex == "male"], alternative = "greater", mu = 0)
t.MSSNG.affected.female <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "MSSNG" & children_pTDT$Affection == "2" & children_pTDT$Sex == "female"], alternative = "greater", mu = 0)
t.SSC.affected.female <- t.test(children_pTDT$overtransmission[children_pTDT$Dataset == "SSC" & children_pTDT$Affection == "2" & children_pTDT$Sex == "female"], alternative = "greater", mu = 0)
adjusted_pvalues = p.adjust(c(t.MSSNG.affected.female$p.value, t.MSSNG.affected.male$p.value, t.SSC.affected.female$p.value, t.SSC.affected.male$p.value), method="holm")

filename = "PRS.pTDT.sex_stratified"
cat("Generating plot", filename, "...\n")
children_pTDT_affected = children_pTDT[children_pTDT$Affection == "2",]
children_pTDT_affected$Dataset.affection.sex = factor(children_pTDT_affected$Dataset.affection.sex, levels=c("MSSNG male\n(with ASD)", "MSSNG female\n(with ASD)", "SSC male\n(with ASD)", "SSC female\n(with ASD)"))
PRS.pTDT.sex_stratified = ggplot(children_pTDT_affected, aes(x=Dataset.affection.sex, y=overtransmission)) +
                          geom_jitter(aes(color=Dataset.affection.sex, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(paste(ylab_PRS, "overtransmission", sep=" ")) + theme(axis.text.x=element_text(vjust=1, size=9)) +
                          scale_colour_manual(values=colors) + annotate("text", x=c(2,3,4), y=c(15), label=paste("p = ",c(signif(adjusted_pvalues[2], digits = 2), signif(adjusted_pvalues[3], digits = 2), signif(adjusted_pvalues[4], digits = 2)), sep=""))
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.pTDT.sex_stratified))


##########################################
### MSSNG multiplex versus SSC simplex ###
##########################################
affected = merged[merged$Affection == "2" & merged$Relation %in% c("proband", "affected sibling"),]
for (i in 1:nrow(affected)) {
    affected_family_members = merged[merged$Family.ID == affected[i,]$Family.ID & merged$Affection == "2",]
    affected$computed_multiplex[i] = ifelse(nrow(affected_family_members) > 1, "Multiplex", "Simplex")
}
affected$computed_multiplex = factor(affected$computed_multiplex, levels=c("Simplex", "Multiplex"))
affected$is_multiplex = affected$computed_multiplex == "Multiplex"
affected = affected[!(affected$Dataset == "MSSNG" & affected$computed_multiplex == "Simplex"),] # get rid of MSSNG "simplex" families (which may not really be simplex)

filename = "PRS.multiplex"
cat("Generating plot", filename, "...\n")
multiplex_versus_simplex_pval = coef(summary(glm(is_multiplex ~ Sex + PRS_1, data=affected)))[3,4]
p = c(multiplex_versus_simplex_pval)
p.adj = p.adjust(p, method="holm")
stat.test = data.frame(group1=c("Multiplex"),
                       group2=c("Simplex"),
                       p=p,
                       p.adj=p.adj,
                       y.position=c(18))

PRS.multiplex = ggplot(affected, aes(x=computed_multiplex, y=PRS_1)) +
                geom_jitter(aes(color=computed_multiplex, alpha=0.5)) + geom_boxplot(alpha=0) + theme + xlab("") + ylab(ylab_PRS) +
                scale_colour_manual(values=colors) + stat_pvalue_manual(stat.test, label=stat.label.noadj, hide.ns=TRUE, tip.length=tip.length)
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.multiplex))

#################################################
### MSSNG MZ twin (technical reproducibility) ###
#################################################
filename = "PRS.MZ_twins"
cat("Generating plot", filename, "...\n")
MZ_twin_data = read.xlsx("MZ_twins.xlsx", sheet="Differences")
MZ_twin_data$MZ_pair = factor(MZ_twin_data$MZ_pair, levels=paste("MZ pair", 1:(length(unique(MZ_twin_data$MZ_pair)))))
PRS.MZ_twin_plot = ggplot(MZ_twin_data, aes(x=MZ_pair, y=PRS_1, color=MZ_pair)) +
                   geom_line() + geom_point(size=1.5) + theme + xlab("") + ylab(ylab_PRS) +
                   scale_y_continuous(limits=c(-20,22), breaks=seq(from=-20, to=22, by=4)) + theme(axis.text.x=element_text(angle=90, vjust=0.5))
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.MZ_twin_plot, width=3))

#################################################
### MSSNG MZ twin scatterplot (technical reproducibility) ###
#################################################
filename = "PRS.MZ_twins.scatterplot"
cat("Generating plot", filename, "...\n")
MZ_twin_data = read.xlsx("MZ_twins.xlsx", sheet="Scatterplot")
#cor(MZ_twin_data$Twin_PRS_1, MZ_twin_data$Twin_PRS_2)
PRS.MZ_twin_plot.scatter = ggplot(MZ_twin_data, aes(x=Twin_1_PRS, y=Twin_2_PRS)) +
                           scale_x_continuous(breaks=seq(-20,20,2)) +
                           scale_y_continuous(breaks=seq(-20,20,2)) +
                           theme(panel.grid.major.x = element_line(colour = "lightgrey", size=0.15)) +
                           geom_point(size=1.5) + theme + xlab("Twin 1 PRS") + ylab("Twin 2 PRS") + geom_smooth(method='lm', formula=y~x) + stat_cor(aes(label=..r.label..), r.digits=2, label.x=0, label.y=2)

suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.MZ_twin_plot.scatter))

########################################################
### Decile analysis (ratio of affected/unaffected) #####
########################################################
children = merged[merged$Affection != "0" & merged$Relation %in% c("proband", "affected sibling", "unaffected sibling"),] # Get only children
#children_1000G = merged[merged$Affection != "0" & merged$Relation %in% c("proband", "affected sibling", "unaffected sibling", "-"),] # Get only children+1000G
children_affected_1000G = merged[merged$Affection != "0" & merged$Relation %in% c("proband", "affected sibling", "-"),] # Get only affected children+1000G

probs = seq(0, 1, 0.1)
decile_counts = data.frame()
filename = "PRS.decile"
for (d in list(list(children, "Siblings without ASD"), list(children_affected_1000G, "1000G"))) {
    df = d[[1]]
    label = d[[2]]
    quant = quantile(df$PRS_1, probs=probs)

    num_affected_decile1 = nrow(df[df$PRS_1 >= quant[1] & df$PRS_1 < quant[2] & df$Affection == "2",])
    num_unaffected_decile1 = nrow(df[df$PRS_1 >= quant[1] & df$PRS_1 < quant[2] & df$Affection == "1",])
    decile1_ratio = num_affected_decile1 / num_unaffected_decile1

    for (i in 2:length(probs)) {
        num_affected = nrow(df[df$PRS_1 >= quant[i-1] & df$PRS_1 < quant[i] & df$Affection == "2",])
        num_unaffected = nrow(df[df$PRS_1 >= quant[i-1] & df$PRS_1 < quant[i] & df$Affection == "1",])
        OR = (num_affected / num_unaffected) / (num_affected_decile1 / num_unaffected_decile1)
        decile = 10*probs[i]
        decile_counts = rbind(decile_counts, list(Decile=decile, PRS_cutoff=quant[i], num_affected=num_affected, num_unaffected=num_unaffected, total_in_decile=num_affected+num_unaffected, OR=OR, label=label))
    }
}
decile_counts$label = factor(decile_counts$label, levels=c("1000G", "Siblings without ASD"))
write.table(decile_counts, file=paste(tsv_dir, "/", filename, ".tsv", sep=""), sep="\t", quote=FALSE)

cat("Generating plot", filename, "...\n")
PRS.decile = ggplot(decile_counts, aes(x=Decile, y=OR, color=str_wrap(label, 18))) +
             geom_point(size=1.5) + geom_line() + theme + xlab("Decile") + ylab("OR") +
             scale_y_continuous(breaks=seq(0, 10, by=0.1)) + scale_x_continuous(breaks=10*probs) + ylab("Odds ratio") +
             theme(legend.position = "top", legend.text = element_text(size=8), legend.title = element_text(size=10)) + labs(color="Individuals without ASD")
suppressMessages(ggsave(paste(pdf_dir, "/", filename, ".pdf", sep=""), plot=PRS.decile))





########################################################
################ PRS_phenotype plots ###################
########################################################


dt <- read.delim(variant_phenotype_table, stringsAsFactors=F)
dt$variant_harvest = ifelse(dt$variant_category == "None", FALSE, TRUE)
dt = dt[!is.na(dt$PRS), ]

dt.out <- data.frame()

do_plot = function(feat) {
  dt.tmp <- dt
  dt.tmp$feature <- dt.tmp[, feat]
  dt.tmp <- dt.tmp[which(dt.tmp$feature != 0), ]

  lm.ref <- lm("PRS ~ Sex + variant_harvest + K1 + K2 + K3", dt.tmp)
  lm.add <- lm("PRS ~ Sex + variant_harvest + K1 + K2 + K3 + feature", dt.tmp)
  conf <- confint(lm.add)
  CI <- paste(signif(conf["feature", 1], digits = 2), signif(conf["feature", 2], digits = 2), sep=",")

  p <- signif(anova(lm.ref, lm.add, test = "Chisq")[2, "Pr(>Chi)"], digits = 1)
  coeff <- signif(lm.add$coefficients["feature"], digits = 1)

  dt.out <- rbind(dt.out, data.frame(feat, coeff, CI, p))
  max <- max(dt.tmp$feature)
  p <- ggplot(dt.tmp, aes(x = PRS, y = feature)) + geom_point() +
    ylab(gsub("_", "\\ ", feat)) + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size = 11)) + geom_smooth(method = "lm") +
    annotate("text", x = -14, y = max-5, label = sprintf("italic(B) == %s", coeff), parse = T) +
    annotate("text", x = -14, y = max-13, label = sprintf("italic(p) == %s", p), parse = T) + xlim(-20, 20)
    return(p)
}

p1 = do_plot("Adaptive_behaviour_standard_score")
p2 = do_plot("Full_scale_IQ")
p3 = do_plot("Global_ability_composite_estimate")
p4 = do_plot("Socialization_standard_score")

PRS_phenotype_plot = ggarrange(p1, p2, p3, p4, ncol=4, nrow=1)
write.table(dt.out, "Sup.Fig.16.underlying.tsv", sep="\t", row.names=F, quote=F, col.names=T)


#####################
#####################
### Combine plots ###
#####################
#####################
empty_plot = ggplot() + theme_void()
filename = "PRS_plot_main_text"
cat("Generating plot", filename, "...\n")
PRS_plot_main_text_p1 = ggarrange(PRS.MZ_twin_plot, PRS.density, PRS.affected_versus_unaffected_children, PRS.pTDT, PRS.decile, empty_plot, labels=LETTERS, ncol=3, nrow=2)
PRS_plot_main_text_p2 = ggarrange(PRS_plot_main_text_p1, PRS_phenotype_plot, ncol=1, nrow=2, heights=c(20,7), labels=c("", "G"))
suppressMessages(ggsave(paste(filename, ".pdf", sep=""), plot=PRS_plot_main_text_p2, width=12, height=14))

filename = "PRS_plot_supp"
cat("Generating plot", filename, "...\n")
#PRS_plot_supp = ggarrange(PRS.rank_dist, PRS.pTDT.sex_stratified, PRS.multiplex, PRS.mothers_versus_fathers, PRS.parents_plus_1000G, labels=LETTERS, ncol=2, nrow=3)
PRS_plot_supp = ggarrange(PRS.rank_dist, PRS.pTDT.sex_stratified, PRS.multiplex, PRS.mothers_versus_fathers, PRS.parents_plus_1000G, labels=LETTERS, ncol=2, nrow=3)
suppressMessages(ggsave(paste(filename, ".pdf", sep=""), plot=PRS_plot_supp, width=10, height=19))

if (file.exists("Rplots.pdf")) {
    invisible(file.remove("Rplots.pdf"))
}
system("crop_PDFs.sh .")
