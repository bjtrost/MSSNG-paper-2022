# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(GenomicRanges)
files <- list.files("denovo", full.names = T)
message(Sys.time())
manifest <- read.delim("cases.eur.info.tsv", stringsAsFactors = F)

####annotation files
tscript <- read.delim("../../ReferenceData/geneInfo2019/hg38_refGene_20200708.transcript.txt", stringsAsFactors = F, header = F,
                      col.names = c("chr", "start", "end", "isoid", "gsymbol", "enzid"))


lof <- read.delim("../../ReferenceData/constraint/gnomAD_oe_lof0.35_v2.1_Feb14_b38_coordinate.txt", stringsAsFactors = F)
ndd <- read.delim("../../ReferenceData/TADs/NDD.Genes.20200916.Shared.WithWorawatEngchuan.txt", stringsAsFactors = F)
asd102 <- data.frame(readxl::read_excel("../../ReferenceData/TADs/satterstrom_102ASD.xlsx", sheet = 3), stringsAsFactors = F)[1:102, ]
tad <- read.delim("../NoncodingAnnotation/TAD_hg19_annotation_track.tsv", stringsAsFactors = F) ## only use TAD ID
tad.hg38 <- read.delim("../NoncodingAnnotation/TAD_hg38.bed", stringsAsFactors = F, header = F)
tad.hg38$hg38ID <- paste(tad.hg38$V1, tad.hg38$V2, tad.hg38$V3, sep="#")
tad <- merge(tad, tad.hg38, by.x = "Dixon2012_TAD_boundary_ID", by.y = "V4", all.x = T)
tad <- tad[!is.na(tad$hg38ID), ]
tad$Dixon2012_TAD_boundary_ID <- tad$hg38ID
tad <- tad[, c("Dixon2012_TAD_boundary_ID", "Dixon2012_TAD_boundary_ASD135")]

ultra.conserved <- read.delim("ultra_conserved_elements_hg38.bed", stringsAsFactors = F, header = F)
ucg <- GRanges(ultra.conserved$V1, IRanges(ultra.conserved$V2, ultra.conserved$V3), "*")
deepsea <- fread("DeepSEA/Denovo_DeepSEA_reduced_binding_sites.tsv", data.table = F)

#### get promoters
asd135 <- read.delim("../../ReferenceData/TADs/TADA MSSNG+SPARK+ASC ASD gene list - TADA MSSNG+SPARK+ASC ASD gene list.tsv", stringsAsFactors = F)
asd135$gene[asd135$gene == "SRPR"] <- "SRPRA"
asd135$gene[asd135$gene == "SUV420H1"] <- "KMT5B"
asd135.promoter <- tscript[which(tscript$gsymbol %in% asd135$gene), ]
asd102.promoter <- tscript[which(tscript$enzid %in% asd102$entrez_id), ]
ndd.promoter <- tscript[which(tscript$enzid %in% ndd$egID), ]
lof.promoter <- tscript[which(tscript$gsymbol %in% lof$gene), ]

asd102.promoter$promoterstart <- ifelse(asd102.promoter$strand == "+", asd102.promoter$start-10000, asd102.promoter$end)
asd102.promoter$promoterend <- ifelse(asd102.promoter$strand == "+", asd102.promoter$start, asd102.promoter$end+10000)
asd102.promoter <- unique(asd102.promoter[, c("gsymbol", "enzid", "chr", "promoterstart", "promoterend")])

asd135.promoter$promoterstart <- ifelse(asd135.promoter$strand == "+", asd135.promoter$start-10000, asd135.promoter$end)
asd135.promoter$promoterend <- ifelse(asd135.promoter$strand == "+", asd135.promoter$start, asd135.promoter$end+10000)
asd135.promoter <- unique(asd135.promoter[, c("gsymbol", "enzid", "chr", "promoterstart", "promoterend")])

ndd.promoter$promoterstart <- ifelse(ndd.promoter$strand == "+", ndd.promoter$start-10000, ndd.promoter$end)
ndd.promoter$promoterend <- ifelse(ndd.promoter$strand == "+", ndd.promoter$start, ndd.promoter$end+10000)
ndd.promoter <- unique(ndd.promoter[, c("gsymbol", "enzid", "chr", "promoterstart", "promoterend")])

lof.promoter$promoterstart <- ifelse(lof.promoter$strand == "+", lof.promoter$start-10000, lof.promoter$end)
lof.promoter$promoterend <- ifelse(lof.promoter$strand == "+", lof.promoter$start, lof.promoter$end+10000)
lof.promoter <- unique(lof.promoter[, c("gsymbol", "enzid", "chr", "promoterstart", "promoterend")])

overlapPromoter <- function(gr.g, promoter){
  promoter.g <- GRanges(promoter$chr, IRanges(promoter$promoterstart, promoter$promoterend), "*")
  olap <- findOverlaps(gr.g, promoter.g)
  return((1:length(gr.g)) %in% unique(olap@from))
}


dt.out <- data.frame()
for(i in 1:nrow(manifest)){
  message(sprintf("%s/%s families processed %s", i, nrow(manifest), Sys.time()))
  proband.id <- manifest$Sample.ID[i]
  proband.found <- grep(sprintf("/%s", proband.id), files)[1]
  
  if(length(na.omit(c(proband.found))) > 0){
    message(sprintf("%s - proband=%s ", i, proband.found))
    proband <- data.table::fread(files[proband.found], data.table = F)
    proband <- proband[which(proband$phastCons_placental > 0), ]
    
    transmitted <- proband
    ## add tad 135
    transmitted <- merge(transmitted, tad, by = "Dixon2012_TAD_boundary_ID", all.x = T)
    
    if(nrow(transmitted) > 0){
      
      transmitted.g <- GRanges(transmitted$CHROM, IRanges(transmitted$POS, transmitted$POS+1), "*")
      
      olap <- findOverlaps(transmitted.g, ucg)
      transmitted$ultra_conserved <- 1:nrow(transmitted) %in% olap@from
      transmitted$ASD102_promoter <- overlapPromoter(transmitted.g, asd102.promoter)
      transmitted$ASD135_promoter <- overlapPromoter(transmitted.g, asd135.promoter)
      transmitted$NDD1250_promoter <- overlapPromoter(transmitted.g, ndd.promoter)
      transmitted$LOFIntorelance_promoter <- overlapPromoter(transmitted.g, lof.promoter)
      
      tmp <- data.frame("Sample" = proband.id, "All_Denovo" = nrow(transmitted), stringsAsFactors = F)
      
      
      tmp[, c("Downstream")] <- c(sum(transmitted$typeseq_priority == "downstream", na.rm = T))
      tmp[, c("Intronic")] <- c(sum(transmitted$typeseq_priority == "intronic", na.rm = T))
      tmp[, c("Upstream")] <- c(sum(transmitted$typeseq_priority == "upstream", na.rm = T))
      tmp[, c("UTR3")] <- c(sum(transmitted$typeseq_priority == "UTR3", na.rm = T))
      tmp[, c("UTR5")] <- c(sum(transmitted$typeseq_priority == "UTR5", na.rm = T)) 
      
      tmp[, c("Ultra_conserved")] <- sum(transmitted$ultra_conserved) 
      ### GeneHancer
      tmp[, c("GeneHancer")] <- 
        c(sum(!is.na(transmitted$GeneHancer_ID)))
      
      transmitted.gh.genes <- sapply(na.omit(transmitted$GeneHancer_Target), strsplit, ",")

      tmp[, c("GeneHancer_ASD102")] <- 
        sum(sapply(sapply(transmitted.gh.genes, "%in%", asd102$gene), sum) > 0)
      
      tmp[, c("GeneHancer_ASD135")] <- 
        sum(sapply(sapply(transmitted.gh.genes, "%in%", asd135$gene), sum) > 0)
      
      tmp[, c("GeneHancer_NDD")] <- 
        sum(sapply(sapply(transmitted.gh.genes, "%in%", ndd$gsymbol), sum) > 0)
      
      tmp[, c("GeneHancer_LOFIntolerance")] <- 
        sum(sapply(sapply(transmitted.gh.genes, "%in%", lof$gene), sum) > 0)
          
      ### Remap
      tmp[, c("Remap")] <- 
        sum(transmitted$Remap != "")
      
      tmp[, c("Remap_CTCF")] <- 
        length(grep("CTCF", transmitted$Remap))
      
      tmp[, c("Remap_RAD21")] <- 
        length(grep("RAD21", transmitted$Remap))
      
      ### dashr
      # tmp[, c("Dashr_Transcript")] <- 
      #   sum(!is.na(transmitted$dashr_transcript))
      # tmp[, c("Dashr_Exon")] <- 
      #   sum(!is.na(transmitted$dashr_exon))
      tmp[, c("LNCipedia_lncRNA")] <-
        sum(!is.na(transmitted$LNCipedia))
      
      ### TAD
      tmp[, c("TAD_Boundary")] <- 
        sum(!is.na(transmitted$Dixon2012_TAD_boundary_ID))
      
      tmp[, c("TAD_HighBrainExpressionDifference")] <-
        sum(transmitted$Dixon2012_TAD_boundary_diff_score > 0.75, na.rm = T)
      
      tmp[, c("TAD_ASD102_one_domain")] <-
        sum(transmitted$Dixon2012_TAD_boundary_ASD102 == 1, na.rm = T)
      
      tmp[, c("TAD_ASD102_both_domains")] <-
        sum(transmitted$Dixon2012_TAD_boundary_ASD102 == 2, na.rm = T)
      
      tmp[, c("TAD_ASD135_one_domain")] <-
        sum(transmitted$Dixon2012_TAD_boundary_ASD135 == 1, na.rm = T)
      
      tmp[, c("TAD_ASD135_both_domains")] <-
        sum(transmitted$Dixon2012_TAD_boundary_ASD135 == 2, na.rm = T)
      
      tmp[, c("TAD_NDD_one_domain")] <-
        sum(transmitted$Dixon2012_TAD_boundary_NDD1250 == 1, na.rm = T)
      
      tmp[, c("TAD_NDD_both_domains")] <-
        sum(transmitted$Dixon2012_TAD_boundary_NDD1250 == 2, na.rm = T)
      
      tmp[, c("TAD_LOFInterolance_one_domain")] <-
        sum(transmitted$Dixon2012_TAD_boundary_LOFIntorelance == 1, na.rm = T)
      
      tmp[, c("TAD_LOFInterolance_both_domains")] <-
        sum(transmitted$Dixon2012_TAD_boundary_LOFIntorelance == 2, na.rm = T)
      
      ### vep_reg_BIOTYPE
      transmitted$vep_AVERAGE_motif_score <- NA
      if(nrow(transmitted >= 1))
      for(k in 1:nrow(transmitted)){
        if(!is.na(transmitted$vep_MOTIF_SCORE_CHANGE[k]))
          transmitted$vep_AVERAGE_motif_score[k] <- mean(as.numeric(strsplit(as.character(transmitted$vep_MOTIF_SCORE_CHANGE[k]), ",")[[1]]))
      }
      
      tmp[, c("vep_TF")] <- 
        sum(!is.na(transmitted$vep_reg_BIOTYPE))
      tmp[, c("vep_TF_ASD102")] <- 
        sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$ASD102_promoter)
      tmp[, c("vep_TF_ASD135")] <- 
        sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$ASD135_promoter)
      tmp[, c("vep_TF_NDD1250")] <- 
        sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$NDD1250_promoter)
      tmp[, c("vep_TF_LOFIntorelance")] <- 
        sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$LOFIntorelance_promoter)
      
      tmp[, c("vep_CTCF")] <- 
        length(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE))
      tmp[, c("vep_CTCF_ASD102")] <- 
        length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD102_promoter)))
      tmp[, c("vep_CTCF_ASD135")] <- 
        length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD135_promoter)))
      tmp[, c("vep_CTCF_NDD1250")] <- 
        length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$NDD1250_promoter)))
      tmp[, c("vep_CTCF_LOFIntorelance")] <- 
        length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$LOFIntorelance_promoter)))
      
      tmp[, c("vep_Enhancer")] <- 
        length(grep("enhancer", transmitted$vep_reg_BIOTYPE))
      tmp[, c("vep_Enhancer_ASD102")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD102_promoter)))
      tmp[, c("vep_Enhancer_ASD135")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD135_promoter)))
      tmp[, c("vep_Enhancer_NDD1250")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$NDD1250_promoter)))
      tmp[, c("vep_Enhancer_LOFIntorelance")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$LOFIntorelance_promoter)))
      
      tmp[, c("vep_reducedEnhancer")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0)))
      tmp[, c("vep_reducedEnhancer_ASD102")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD102_promoter)))
      tmp[, c("vep_reducedEnhancer_ASD135")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD135_promoter)))
      tmp[, c("vep_reducedEnhancer_NDD1250")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$NDD1250_promoter)))
      tmp[, c("vep_reducedEnhancer_LOFIntorelance")] <- 
        length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$LOFIntorelance_promoter)))
      
      tmp[, c("vep_Promoter")] <- 
        length(grep("promoter", transmitted$vep_reg_BIOTYPE))
      tmp[, c("vep_Promoter_ASD102")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD102_promoter)))
      tmp[, c("vep_Promoter_ASD135")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD135_promoter)))
      tmp[, c("vep_Promoter_NDD1250")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$NDD1250_promoter)))
      tmp[, c("vep_Promoter)LOFIntorelance")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$LOFIntorelance_promoter)))
      
      tmp[, c("vep_reducedPromote")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0)))
      tmp[, c("vep_reducedPromote_ASD102")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD102_promoter)))
      tmp[, c("vep_reducedPromote_ASD135")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD135_promoter)))
      tmp[, c("vep_reducedPromote_NDD1250")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$NDD1250_promoter)))
      tmp[, c("vep_reducedPromote_LOFIntorelance")] <- 
        length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$LOFIntorelance_promoter)))
      
      tmp[, c("DeepSEA_anyTFAnyCell", "DeepSEA_anyTFBrainCell", "DeepSEA_DNaseBrainCell", "DeepSEA_ActiveEnhancerPromoterBrainCell")] <- c(sum(transmitted$"#id" %in% deepsea$name), 
           sum(transmitted$"#id" %in% deepsea$name[deepsea$anyTFBrainCell > 0]), 
           sum(transmitted$"#id" %in% deepsea$name[deepsea$DNaseBrainCell > 0]), 
           sum(transmitted$"#id" %in% deepsea$name[deepsea$ActiveEnhancerPromoterBrainCell > 0]))
      
      tmp[, c("DeepSEA_anyTFAnyCell_LOF_promoter", "DeepSEA_anyTFBrainCell_LOF_promoter", "DeepSEA_DNaseBrainCell_LOF_promoter", "DeepSEA_ActiveEnhancerPromoterBrainCell_LOF_promoter")] <- c(sum(transmitted$"#id" %in% deepsea$name & transmitted$LOFIntorelance_promoter), 
                                                                                                                                           sum(transmitted$"#id" %in% deepsea$name[deepsea$anyTFBrainCell > 0] & transmitted$LOFIntorelance_promoter), 
                                                                                                                                           sum(transmitted$"#id" %in% deepsea$name[deepsea$DNaseBrainCell > 0] & transmitted$LOFIntorelance_promoter), 
                                                                                                                                           sum(transmitted$"#id" %in% deepsea$name[deepsea$ActiveEnhancerPromoterBrainCell > 0] & transmitted$LOFIntorelance_promoter))
      dt.out <- rbind(dt.out, tmp)
    
    }
  }
}

write.table(dt.out, sprintf("denovo_features_phastcon0_%s.tsv", Sys.Date()), sep="\t", row.names=F, quote=F, col.names=T)
