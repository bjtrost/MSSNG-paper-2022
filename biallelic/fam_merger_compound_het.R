#!/usr/bin/env Rscript

library(optparse)
library(tidyr)
library(dplyr)

`%!in%` <- Negate(`%in%`)

option_list = list(
  make_option(c("-n", "--indexname"), type="character", default=NULL, 
              help="list of index", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

index_list <-opt$indexname
path<-opt$path
output<-opt$output
MSSNG_metadata <- read.delim(opt$metadata)
# 
#index_list <- list("SSC04996")
# 
#path<-"C:/Users/ealiyev/Desktop//"
#output<-"C:/Users/ealiyev/Desktop/"
# 
#MSSNG_metadata <- read.delim("C:/Users/ealiyev/Desktop/SSC_metadata.tsv")

for (variable in index_list) {
  index_file<-read.delim(paste(path,variable,".tsv.gz",sep = ""),header = TRUE,comment.char = "",row.names = NULL)
  colnames(index_file)
  
  family_id<-filter(MSSNG_metadata,Sample.ID==variable)$Family.ID
  
  fam_members<-filter(MSSNG_metadata,Family.ID==family_id,Relation!="proband")$Individual.ID
  
  
  
  for (fam_member in fam_members) {
    fam_member_file<-read.delim(paste(path,fam_member,".tsv.gz",sep = ""),header = TRUE,comment.char = "",row.names = NULL)
    if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "proband") {
      suffix<-"proband"
    } else if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "mother") {
      suffix<-"mother"
    } else if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "affected sibling") {
      suffix<-paste("aff_sibling_",fam_member,sep = "")
    } else if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "father") {
      suffix<-"father"
    } else if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "unaffected sibling") {
      suffix<-paste("unaff_sibling_",fam_member,sep = "")
    }
    
    id<-paste(suffix,"id",sep="")
    gt<-paste(suffix,"_GT",sep="")
    fam_member_file <- fam_member_file %>%
      rename(!!gt := GT,
             !!id := X.id) 
    
    
    
    index_file<-merge(x = index_file, 
                      y = fam_member_file[ , c(id,gt)], 
                      by.x="X.id", by.y=id,all.x = TRUE)
    index_file[which(is.na(index_file[[gt]])),][[gt]]<-"0/0"
  }
}



filtered<-index_file %>%
  filter(FILTER == "PASS") %>%
  filter(CHROM != "chrX") %>%
  filter(CHROM != "chrY") %>%
  filter(freq_max < 0.01) %>%
  filter(gnomAD_exome_freq_max < 0.01) %>%
  filter(gnomAD_genome_freq_max < 0.01) %>%
  filter(A1000g_freq_max < 0.01) %>%
  filter(cg_freq_max < 0.01) %>%
  filter(inheritance == "paternal" | inheritance == "maternal") %>%
  filter(GT == "0/1")%>%
  filter(high_quality == "true") %>%
  #filter(if_all(starts_with("aff_sibling_"), ~ .x == "0/1")) %>%
  filter(if_all(starts_with("unaff_sibling_"), ~ .x != "1/1")) %>%
  filter(substr(effect,1,4) == "nons" & MPC_score >= 1 | 
           substr(effect,1,3) == "sto" |
           substr(effect,1,3) == "fra" |
           substr(typeseq,1,3) == "int" )


gene_count<-as.data.frame(table(filtered$gene_symbol))

filtered_gene_list<-gene_count %>%
  filter(Freq == 2)

final_tsv<-filter(filtered,filtered$gene_symbol %in% filtered_gene_list$Var1)

optimus<-final_tsv %>%
  group_by(gene_symbol,inheritance) %>%
  tally()

gene_count<-as.data.frame(table(optimus$gene_symbol))

filtered_gene_list<-gene_count %>%
  filter(Freq > 1)

final_final_tsv<-filter(final_tsv,final_tsv$gene_symbol %in% filtered_gene_list$Var1)

bad_genes_total = data.frame(matrix(ncol = 1, nrow = 0))
x <- c("gene_symbol")
colnames(bad_genes_total) <- x

for (fam_member in fam_members) {
  
  
  fam_member_file<-read.delim(paste(path,fam_member,".tsv.gz",sep = ""),header = TRUE,comment.char = "",row.names = NULL)
  if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "unaffected sibling") {
      suffix<-paste("unaff_sibling_",fam_member,"_GT",sep = "")
      grp <- c(suffix)
      gt_count<-final_final_tsv %>%
      group_by(gene_symbol,across(all_of(grp))) %>%
      tally()
      
      bad_genes<-gt_count %>%
        filter((!!sym(suffix)) == "0/1") %>%
        filter(n == 2) 
      
      filtered_gene_list<-bad_genes$gene_symbol
      
      bad_genes_total <- rbind(bad_genes_total,filtered_gene_list)
  }
  
}

names(bad_genes_total) <- c('Var1')

final_final_final_tsv<-filter(final_final_tsv,final_final_tsv$gene_symbol %!in% bad_genes_total$Var1)

write.table(final_final_final_tsv,paste(output,"FAM-proband-",family_id,".tsv",sep=""),quote = FALSE, sep='\t',row.names = FALSE)

affected_siblings<-filter(MSSNG_metadata,Family.ID==family_id,Relation=="affected sibling")$Individual.ID

for (affected_sibling in affected_siblings) {
  suffix<-paste("aff_sibling_",affected_sibling,sep = "")
  gt_aff_sibling<-paste(suffix,"_GT",sep="")
  
  filtered<-index_file %>%
    filter(FILTER == "PASS") %>%
    filter(CHROM != "chrX") %>%
    filter(CHROM != "chrY") %>%
    filter(freq_max < 0.01) %>%
    filter(gnomAD_exome_freq_max < 0.01) %>%
    filter(gnomAD_genome_freq_max < 0.01) %>%
    filter(A1000g_freq_max < 0.01) %>%
    filter(cg_freq_max < 0.01) %>%
    filter(inheritance == "paternal" | inheritance == "maternal") %>%
    #filter(GT == "0/1")%>%
    filter(get(gt_aff_sibling)  == "0/1") %>%
    #filter(if_all(starts_with("aff_sibling_"), ~ .x == "0/1")) %>%
    filter(if_all(starts_with("unaff_sibling_"), ~ .x != "1/1")) %>%
    filter(high_quality == "true") %>%
    filter(substr(effect,1,4) == "nons" & MPC_score >= 1 | 
             substr(effect,1,3) == "sto" |
             substr(effect,1,3) == "fra" |
             substr(typeseq,1,3) == "int" )
  
  gene_count<-as.data.frame(table(filtered$gene_symbol))
  
  filtered_gene_list<-gene_count %>%
    filter(Freq == 2)
  
  final_tsv<-filter(filtered,filtered$gene_symbol %in% filtered_gene_list$Var1)
  
  optimus<-final_tsv %>%
    group_by(gene_symbol,inheritance) %>%
    tally()
  
  gene_count<-as.data.frame(table(optimus$gene_symbol))
  
  filtered_gene_list<-gene_count %>%
    filter(Freq > 1)
  
  final_final_tsv<-filter(final_tsv,final_tsv$gene_symbol %in% filtered_gene_list$Var1)
  
  bad_genes_total = data.frame(matrix(ncol = 1, nrow = 0))
  x <- c("gene_symbol")
  colnames(bad_genes_total) <- x
  
  filtered$X.Sample<-affected_sibling
  
  for (fam_member in fam_members) {
    
    
    fam_member_file<-read.delim(paste(path,fam_member,".tsv.gz",sep = ""),header = TRUE,comment.char = "",row.names = NULL)
    if (MSSNG_metadata[which(MSSNG_metadata$Individual.ID == fam_member),]$Relation == "unaffected sibling") {
      suffix<-paste("unaff_sibling_",fam_member,"_GT",sep = "")
      grp <- c(suffix)
      gt_count<-final_final_tsv %>%
        group_by(gene_symbol,across(all_of(grp))) %>%
        tally()
      
      bad_genes<-gt_count %>%
        filter((!!sym(suffix)) == "0/1") %>%
        filter(n == 2) 
      
      filtered_gene_list<-bad_genes$gene_symbol
      
      bad_genes_total <- rbind(bad_genes_total,filtered_gene_list)
    }
    
  }
  
  names(bad_genes_total) <- c('Var1')
  
  final_final_final_tsv<-filter(final_final_tsv,final_final_tsv$gene_symbol %!in% bad_genes_total$Var1)
  
  write.table(final_final_final_tsv,paste(output,"FAM-",family_id,"affected_sib",suffix,".tsv",sep=""),quote = FALSE, sep='\t',row.names = FALSE)  
}



