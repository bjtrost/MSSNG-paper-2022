#!/usr/bin/env Rscript

library(optparse)
library(tidyr)
library(dplyr)

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

#index_list <- list("1-0206-004")

#path<-"C:/Users/ealiyev/Desktop//Work/MSSNG/"
#output<-"C:/Users/ealiyev/Desktop/Work/MSSNG/"

#MSSNG_metadata <- read.delim("C:/Users/ealiyev/Desktop/Work/MSSNG/MSSNG_metadata_rev.tsv")

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
                      by.x="X.id", by.y=id)
  }
}

if(!"mother_GT" %in% colnames(index_file))
{
  index_file$mother_GT<-NA
}

if(!"father_GT" %in% colnames(index_file))
{
  index_file$father_GT<-NA
}

filtered<-index_file %>%
  filter(FILTER == "PASS") %>%
  filter(CHROM != "chrX") %>%
  filter(CHROM != "chrY") %>%
  filter(mother_GT == "0/1" | is.na(mother_GT)) %>%
  filter(father_GT == "0/1" | is.na(father_GT)) %>%
  filter(freq_max < 0.01) %>%
  filter(gnomAD_exome_freq_max < 0.01) %>%
  filter(gnomAD_genome_freq_max < 0.01) %>%
  filter(A1000g_freq_max < 0.01) %>%
  filter(cg_freq_max < 0.01) %>%
  filter(GT == "1/1") %>%
  filter(high_quality == "true") %>%
  filter(substr(effect,1,4) == "nons" & MPC_score >= 1 | 
           substr(effect,1,3) == "sto" |
           substr(effect,1,3) == "fra" |
           substr(typeseq,1,3) == "int" ) %>%
  filter(if_all(starts_with("unaff_sibling_"), ~ .x != "1/1"))

write.table(filtered,paste(output,"FAM-proband-",family_id,".tsv",sep=""),quote = FALSE, sep='\t',row.names = FALSE)

#affected sibling part

affected_siblings<-filter(MSSNG_metadata,Family.ID==family_id,Relation=="affected sibling")$Individual.ID

for (affected_sibling in affected_siblings) {
  suffix<-paste("aff_sibling_",affected_sibling,sep = "")
  gt_aff_sibling<-paste(suffix,"_GT",sep="")
  
  
  filtered_aff_sibling<-index_file %>%
    filter(FILTER == "PASS") %>%
    filter(CHROM != "chrX") %>%
    filter(CHROM != "chrY") %>%
    filter(mother_GT == "0/1" | is.na(mother_GT)) %>%
    filter(father_GT == "0/1" | is.na(father_GT)) %>%
    filter(freq_max < 0.01) %>%
    filter(gnomAD_exome_freq_max < 0.01) %>%
    filter(gnomAD_genome_freq_max < 0.01) %>%
    filter(A1000g_freq_max < 0.01) %>%
    filter(cg_freq_max < 0.01) %>%
    filter(get(gt_aff_sibling)  == "1/1") %>%
    filter(high_quality == "true") %>%
    filter(substr(effect,1,4) == "nons" & MPC_score >= 1 | 
             substr(effect,1,3) == "sto" |
             substr(effect,1,3) == "fra" |
             substr(typeseq,1,3) == "int" ) %>%
    filter(if_all(starts_with("unaff_sibling_"), ~ .x != "1/1")) 
  
  filtered_aff_sibling$X.Sample<-affected_sibling
  
write.table(filtered_aff_sibling,paste(output,"FAM-",family_id,"affected_sib",suffix,".tsv",sep=""),quote = FALSE, sep='\t',row.names = FALSE)
  
}