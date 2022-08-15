# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(GenomicRanges)
files <- list.files("singletons", full.names = T)
message(Sys.time())
manifest <- read.delim("cases.info.tsv", stringsAsFactors = F)

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
dt.out <- data.frame()

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

for(i in 1:nrow(manifest)){
  message(sprintf("%s/%s families processed %s", i, nrow(manifest), Sys.time()))
  proband.id <- manifest$Sample.ID[i]
  father.id <- strsplit(manifest$FAM[i], "#")[[1]][1]
  mother.id <- strsplit(manifest$FAM[i], "#")[[1]][2]
  
  proband.found <- grep(sprintf("/%s", proband.id), files)[1]
  father.found <- grep(sprintf("/%s", father.id), files)[1]
  mother.found <- grep(sprintf("/%s", mother.id), files)[1]
  
  if(length(na.omit(c(proband.found, father.found, mother.found))) == 3){
    message(sprintf("%s - proband=%s, father=%s, mother=%s", i, proband.found, father.found, mother.found))
    proband <- data.table::fread(files[proband.found], data.table = F)
    father <- data.table::fread(files[father.found], data.table = F)
    mother <- data.table::fread(files[mother.found], data.table = F)
    
    common.columns <- names(proband)[names(proband) %in% names(mother)]
    
    # transmitted <- proband[proband$inheritance %in% c("maternal", "paternal", "ambiguous"), common.columns]
    
    transmitted <- proband[proband$inheritance %in% c("maternal", "paternal", "ambiguous"), common.columns]
    transmitted <- merge(transmitted, tad, by = "Dixon2012_TAD_boundary_ID", all.x = T)
    
    if(nrow(transmitted) > 0){
      nontransmitted <- rbind(father[, common.columns], mother[, common.columns])
      nontransmitted <- merge(nontransmitted, tad, by = "Dixon2012_TAD_boundary_ID", all.x = T)
      
      # ambiguous <- nontransmitted[nontransmitted$`#id` %in% transmitted$`#id`[transmitted$inheritance == "ambiguous"], ]
      # ambiguous <- ambiguous[!duplicated(ambiguous$`#id`), ]
      nontransmitted <- nontransmitted[!nontransmitted$`#id` %in% transmitted$`#id`, ]
      
      ##remove ambiguous out from analysis
      # transmitted <- transmitted[transmitted$inheritance != "ambiguous", ]
      # nontransmitted <- rbind(nontransmitted, ambiguous)
      
      transmitted <- transmitted[which((is.na(transmitted$freq_max) | transmitted$freq_max < 0.001) &
                                 (is.na(transmitted$A1000g_freq_max) | transmitted$A1000g_freq_max < 0.001) &
                                 (is.na(transmitted$gnomAD_exome_freq_max) | transmitted$gnomAD_exome_freq_max < 0.001) &
                                 (is.na(transmitted$gnomAD_genome_freq_max) | transmitted$gnomAD_genome_freq_max < 0.001) &
                                 (is.na(transmitted$ilm_max_int_freq) | transmitted$ilm_max_int_freq < 0.001)), ]
      
      nontransmitted <- nontransmitted[which((is.na(nontransmitted$freq_max) | nontransmitted$freq_max < 0.001) &
                                         (is.na(nontransmitted$A1000g_freq_max) | nontransmitted$A1000g_freq_max < 0.001) &
                                         (is.na(nontransmitted$gnomAD_exome_freq_max) | nontransmitted$gnomAD_exome_freq_max < 0.001) &
                                         (is.na(nontransmitted$gnomAD_genome_freq_max) | nontransmitted$gnomAD_genome_freq_max < 0.001) &
                                         (is.na(nontransmitted$ilm_max_int_freq) | nontransmitted$ilm_max_int_freq < 0.001)), ]
      
      transmitted.g <- GRanges(transmitted$CHROM, IRanges(transmitted$POS, transmitted$POS+1), "*")
      nontransmitted.g <- GRanges(nontransmitted$CHROM, IRanges(nontransmitted$POS, nontransmitted$POS+1), "*")
      
      olap <- findOverlaps(transmitted.g, ucg)
      transmitted$ultra_conserved <- 1:nrow(transmitted) %in% olap@from
      
      olap <- findOverlaps(nontransmitted.g, ucg)
      nontransmitted$ultra_conserved <- 1:nrow(nontransmitted) %in% olap@from
      # transmitted <- transmitted[transmitted$`#id` %in% nontransmitted$`#id`, ]
      
      
      
      tmp <- data.frame("Sample" = proband.id, "All_Features_Transmitted" = nrow(transmitted),
                        "All_Features_Nontransmitted" = nrow(nontransmitted), stringsAsFactors = F)
      
      tmp[, c("Downstream_Transmitted", "Downstream_Nontransmitted")] <- c(sum(transmitted$typeseq_priority == "downstream", na.rm = T),
                                                                                     sum(nontransmitted$typeseq_priority == "downstream", na.rm = T))
      tmp[, c("Intronic_Transmitted", "Intronic_Nontransmitted")] <- c(sum(transmitted$typeseq_priority == "intronic", na.rm = T),
                                                                           sum(nontransmitted$typeseq_priority == "intronic", na.rm = T))
      tmp[, c("Upstream_Transmitted", "Upstream_Nontransmitted")] <- c(sum(transmitted$typeseq_priority == "upstream", na.rm = T),
                                                                           sum(nontransmitted$typeseq_priority == "upstream", na.rm = T))
      tmp[, c("UTR3_Transmitted", "UTR3_Nontransmitted")] <- c(sum(transmitted$typeseq_priority == "UTR3", na.rm = T),
                                                                           sum(nontransmitted$typeseq_priority == "UTR3", na.rm = T))
      tmp[, c("UTR5_Transmitted", "UTR5_Nontransmitted")] <- c(sum(transmitted$typeseq_priority == "UTR5", na.rm = T),
                                                                           sum(nontransmitted$typeseq_priority == "UTR5", na.rm = T)) 
      
      tmp[, c("Ultra_conserved_Transmitted", "Ultra_conserved_Nontransmitted")] <- c(sum(transmitted$ultra_conserved),
                                                                                         sum(nontransmitted$ultra_conserved)) 
      nontransmitted.g <- GRanges(nontransmitted$CHROM, IRanges(nontransmitted$POS, nontransmitted$POS+1), "*")
      
      transmitted$ASD102_promoter <- overlapPromoter(transmitted.g, asd102.promoter)
      transmitted$ASD135_promoter <- overlapPromoter(transmitted.g, asd135.promoter)
      transmitted$NDD1250_promoter <- overlapPromoter(transmitted.g, ndd.promoter)
      transmitted$LOFIntolerance_promoter <- overlapPromoter(transmitted.g, lof.promoter)
      
      nontransmitted$ASD102_promoter <- overlapPromoter(nontransmitted.g, asd102.promoter)
      nontransmitted$ASD135_promoter <- overlapPromoter(nontransmitted.g, asd135.promoter)
      nontransmitted$NDD1250_promoter <- overlapPromoter(nontransmitted.g, ndd.promoter)
      nontransmitted$LOFIntolerance_promoter <- overlapPromoter(nontransmitted.g, lof.promoter)
      
      ### GeneHancer
      tmp[, c("GeneHancer_Transmitted", "GeneHancer_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$GeneHancer_ID)), sum(!is.na(nontransmitted$GeneHancer_ID)))
      transmitted.gh.genes <- sapply(na.omit(transmitted$GeneHancer_Target), strsplit, ",")
      nontransmitted.gh.genes <- sapply(na.omit(nontransmitted$GeneHancer_Target), strsplit, ",")
      
      tmp[, c("GeneHancer_ASD102_Transmitted", "GeneHancer_ASD102_Nontransmitted")] <- 
        c(sum(sapply(sapply(transmitted.gh.genes, "%in%", asd102$gene), sum) > 0),
          sum(sapply(sapply(nontransmitted.gh.genes, "%in%", asd102$gene), sum) > 0))
      
      tmp[, c("GeneHancer_ASD135_Transmitted", "GeneHancer_ASD135_Nontransmitted")] <- 
        c(sum(sapply(sapply(transmitted.gh.genes, "%in%", asd135$gene), sum) > 0),
          sum(sapply(sapply(nontransmitted.gh.genes, "%in%", asd135$gene), sum) > 0))
      
      tmp[, c("GeneHancer_NDD_Transmitted", "GeneHancer_NDD_Nontransmitted")] <- 
        c(sum(sapply(sapply(transmitted.gh.genes, "%in%", ndd$gsymbol), sum) > 0),
          sum(sapply(sapply(nontransmitted.gh.genes, "%in%", ndd$gsymbol), sum) > 0))
      
      tmp[, c("GeneHancer_LOFIntolerance_Transmitted", "GeneHancer_LOFIntolerance_Nontransmitted")] <- 
        c(sum(sapply(sapply(transmitted.gh.genes, "%in%", lof$gene), sum) > 0),
          sum(sapply(sapply(nontransmitted.gh.genes, "%in%", lof$gene), sum) > 0))
      
      ### Remap
      tmp[, c("Remap_Transmitted", "Remap_Nontransmitted")] <- 
        c(sum(transmitted$Remap != ""), sum(nontransmitted$Remap != ""))
      
      tmp[, c("Remap_CTCF_Transmitted", "Remap_CTCF_Nontransmitted")] <- 
        c(length(grep("CTCF", transmitted$Remap)), length(grep("CTCF", nontransmitted$Remap)))
      
      tmp[, c("Remap_RAD21_Transmitted", "Remap_RAD21_Nontransmitted")] <- 
        c(length(grep("RAD21", transmitted$Remap)), length(grep("RAD21", nontransmitted$Remap)))
      tmp[, c("LNCipedia_lncRNA_Transmitted", "LNCipedia_lncRNA_Nontransmitted")] <-
        c(sum(!is.na(transmitted$LNCipedia)), sum(!is.na(nontransmitted$LNCipedia)))
      
      ### dashr
      # tmp[, c("Dashr_Transcript_Transmitted", "Dashr_Transcript_Nontransmitted")] <- 
      #   c(sum(!is.na(transmitted$dashr_transcript)), sum(!is.na(nontransmitted$dashr_transcript)))
      # tmp[, c("Dashr_Exon_Transmitted", "Dashr_Exon_Nontransmitted")] <- 
      #   c(sum(!is.na(transmitted$dashr_exon)), sum(!is.na(nontransmitted$dashr_exon)))
      
      ### TAD
      tmp[, c("TAD_Boundary_Transmitted", "TAD_Boundary_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$Dixon2012_TAD_boundary_ID)), sum(!is.na(nontransmitted$Dixon2012_TAD_boundary_ID)))
      
      tmp[, c("TAD_HighBrainExpressionDifference_Transmitted", "TAD_HighBrainExpressionDifference_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_diff_score > 0.75, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_diff_score > 0.75, na.rm = T))
      
      tmp[, c("TAD_ASD102_one_domain_Transmitted", "TAD_ASD102_one_domain_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_ASD102 == 1, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_ASD102 == 1, na.rm = T))
      
      tmp[, c("TAD_ASD102_both_domains_Transmitted", "TAD_ASD102_both_domains_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_ASD102 == 2, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_ASD102 == 2, na.rm = T))
      
      tmp[, c("TAD_ASD135_one_domain_Transmitted", "TAD_ASD135_one_domain_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_ASD135 == 1, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_ASD135 == 1, na.rm = T))
      
      tmp[, c("TAD_ASD135_both_domains_Transmitted", "TAD_ASD135_both_domains_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_ASD135 == 2, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_ASD135 == 2, na.rm = T))
      
      tmp[, c("TAD_NDD_one_domain_Transmitted", "TAD_NDD_one_domain_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_NDD1250 == 1, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_NDD1250 == 1, na.rm = T))
      
      tmp[, c("TAD_NDD_both_domains_Transmitted", "TAD_NDD_both_domains_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_NDD1250 == 2, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_NDD1250 == 2, na.rm = T))
      
      tmp[, c("TAD_LOFInterolance_one_domain_Transmitted", "TAD_LOFInterolance_one_domain_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_LOFIntolerance == 1, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_LOFIntolerance == 1, na.rm = T))
      
      tmp[, c("TAD_LOFInterolance_both_domains_Transmitted", "TAD_LOFInterolance_both_domains_Nontransmitted")] <-
        c(sum(transmitted$Dixon2012_TAD_boundary_LOFIntolerance == 2, na.rm = T), 
          sum(nontransmitted$Dixon2012_TAD_boundary_LOFIntolerance == 2, na.rm = T))
      
      ### vep_reg_BIOTYPE
      transmitted$vep_AVERAGE_motif_score <- NA
      if(nrow(transmitted >= 1))
        for(k in 1:nrow(transmitted)){
          if(!is.na(transmitted$vep_MOTIF_SCORE_CHANGE[k]))
            transmitted$vep_AVERAGE_motif_score[k] <- mean(as.numeric(strsplit(as.character(transmitted$vep_MOTIF_SCORE_CHANGE[k]), ",")[[1]]))
        }
      nontransmitted$vep_AVERAGE_motif_score <- NA
      if(nrow(nontransmitted >= 1))
        for(k in 1:nrow(nontransmitted)){
          if(!is.na(nontransmitted$vep_MOTIF_SCORE_CHANGE[k]))
            nontransmitted$vep_AVERAGE_motif_score[k] <- mean(as.numeric(strsplit(as.character(nontransmitted$vep_MOTIF_SCORE_CHANGE[k]), ",")[[1]]))
        }
      
      tmp[, c("vep_TF_Transmitted", "vep_TF_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$vep_reg_BIOTYPE)), sum(!is.na(nontransmitted$vep_reg_BIOTYPE)))
      tmp[, c("vep_TF_ASD102_Transmitted", "vep_TF_ASD102_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$ASD102_promoter), sum(!is.na(nontransmitted$vep_reg_BIOTYPE) & nontransmitted$ASD102_promoter))
      tmp[, c("vep_TF_ASD135_Transmitted", "vep_TF_ASD135_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$ASD135_promoter), sum(!is.na(nontransmitted$vep_reg_BIOTYPE) & nontransmitted$ASD135_promoter))
      tmp[, c("vep_TF_NDD1250_Transmitted", "vep_TF_NDD1250_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$NDD1250_promoter), sum(!is.na(nontransmitted$vep_reg_BIOTYPE) & nontransmitted$NDD1250_promoter))
      tmp[, c("vep_TF_LOFIntolerance_Transmitted", "vep_TF_LOFIntolerance_Nontransmitted")] <- 
        c(sum(!is.na(transmitted$vep_reg_BIOTYPE) & transmitted$LOFIntolerance_promoter), sum(!is.na(nontransmitted$vep_reg_BIOTYPE) & nontransmitted$LOFIntolerance_promoter))
      
      tmp[, c("vep_CTCF_Transmitted", "vep_CTCF_Nontransmitted")] <- 
        c(length(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE)), 
          length(grep("CTCF_binding_site", nontransmitted$vep_reg_BIOTYPE)) )
      tmp[, c("vep_CTCF_ASD102_Transmitted", "vep_CTCF_ASD102_Nontransmitted")] <- 
        c(length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD102_promoter))), 
          length(intersect(grep("CTCF_binding_site", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$ASD102_promoter))))
      tmp[, c("vep_CTCF_ASD135_Transmitted", "vep_CTCF_ASD135_Nontransmitted")] <- 
        c(length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD135_promoter))), 
          length(intersect(grep("CTCF_binding_site", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$ASD135_promoter))))
      tmp[, c("vep_CTCF_NDD1250_Transmitted", "vep_CTCF_NDD1250_Nontransmitted")] <- 
        c(length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$NDD1250_promoter))), 
          length(intersect(grep("CTCF_binding_site", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$NDD1250_promoter))))
      tmp[, c("vep_CTCF_LOFIntolerance_Transmitted", "vep_CTCF_LOFIntolerance_Nontransmitted")] <- 
        c(length(intersect(grep("CTCF_binding_site", transmitted$vep_reg_BIOTYPE), which(transmitted$LOFIntolerance_promoter))), 
          length(intersect(grep("CTCF_binding_site", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$LOFIntolerance_promoter))))
      
      tmp[, c("vep_Enhancer_Transmitted", "vep_Enhancer_Nontransmitted")] <- 
        c(length(grep("enhancer", transmitted$vep_reg_BIOTYPE)), 
          length(grep("enhancer", nontransmitted$vep_reg_BIOTYPE)) )
      tmp[, c("vep_Enhancer_ASD102_Transmitted", "vep_Enhancer_ASD102_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD102_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$ASD102_promoter))))
      tmp[, c("vep_Enhancer_ASD135_Transmitted", "vep_Enhancer_ASD135_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD135_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$ASD135_promoter))))
      tmp[, c("vep_Enhancer_NDD1250_Transmitted", "vep_Enhancer_NDD1250_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$NDD1250_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$NDD1250_promoter))))
      tmp[, c("vep_Enhancer_LOFIntolerance_Transmitted", "vep_Enhancer_LOFIntolerance_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$LOFIntolerance_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$LOFIntolerance_promoter))))
      
      tmp[, c("vep_reducedEnhancer_Transmitted", "vep_reducedEnhancer_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0))))
      tmp[, c("vep_reducedEnhancer_ASD102_Transmitted", "vep_reducedEnhancer_ASD102_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD102_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$ASD102_promoter))))
      tmp[, c("vep_reducedEnhancer_ASD135_Transmitted", "vep_reducedEnhancer_ASD135_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD135_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$ASD135_promoter))))
      tmp[, c("vep_reducedEnhancer_NDD1250_Transmitted", "vep_reducedEnhancer_NDD1250_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$NDD1250_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$NDD1250_promoter))))
      tmp[, c("vep_reducedEnhancer_LOFIntolerance_Transmitted", "vep_reducedEnhancer_LOFIntolerance_Nontransmitted")] <- 
        c(length(intersect(grep("enhancer", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$LOFIntolerance_promoter))), 
          length(intersect(grep("enhancer", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$LOFIntolerance_promoter))))
      
      tmp[, c("vep_Promoter_Transmitted", "vep_Promoter_Nontransmitted")] <- 
        c(length(grep("promoter", transmitted$vep_reg_BIOTYPE)), 
          length(grep("promoter", nontransmitted$vep_reg_BIOTYPE)) )
      tmp[, c("vep_Promoter_ASD102_Transmitted", "vep_Promoter_ASD102_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD102_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$ASD102_promoter))))
      tmp[, c("vep_Promoter_ASD135_Transmitted", "vep_Promoter_ASD135_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$ASD135_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$ASD135_promoter))))
      tmp[, c("vep_Promoter_NDD1250_Transmitted", "vep_Promoter_NDD1250_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$NDD1250_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$NDD1250_promoter))))
      tmp[, c("vep_Promoter_LOFIntolerance_Transmitted", "vep_Promoter_LOFIntolerance_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$LOFIntolerance_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$LOFIntolerance_promoter))))
      
      tmp[, c("vep_reducedPromote_Transmitted", "vep_reducedPromote_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0))))
      tmp[, c("vep_reducedPromote_ASD102_Transmitted", "vep_reducedPromote_ASD102_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD102_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$ASD102_promoter))))
      tmp[, c("vep_reducedPromote_ASD135_Transmitted", "vep_reducedPromote_ASD135_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$ASD135_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$ASD135_promoter))))
      tmp[, c("vep_reducedPromote_NDD1250_Transmitted", "vep_reducedPromote_NDD1250_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$NDD1250_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$NDD1250_promoter))))
      tmp[, c("vep_reducedPromote_LOFIntolerance_Transmitted", "vep_reducedPromote_LOFIntolerance_Nontransmitted")] <- 
        c(length(intersect(grep("promoter", transmitted$vep_reg_BIOTYPE), which(transmitted$vep_AVERAGE_motif_score < 0 & transmitted$LOFIntolerance_promoter))), 
          length(intersect(grep("promoter", nontransmitted$vep_reg_BIOTYPE), which(nontransmitted$vep_AVERAGE_motif_score < 0 & nontransmitted$LOFIntolerance_promoter))))
      
      dt.out <- rbind(dt.out, tmp)
    
    }
  }
}

write.table(dt.out, sprintf("transmission_features_singletons_%s.tsv", Sys.Date()), sep="\t", row.names=F, quote=F, col.names=T)
