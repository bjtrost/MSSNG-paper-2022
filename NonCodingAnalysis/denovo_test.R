setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(ggplot2)
library(survival)
message(Sys.time())
### mssng
manifest <- read.delim("cases.eur.info.tsv", stringsAsFactors = F)
dt.feat.phastcon0 <- read.delim("denovo_features_phastcon0_2021-08-13.tsv", stringsAsFactors = F)
dt.feat <- read.delim("denovo_features_2021-07-31.tsv", stringsAsFactors = F)
names(dt.feat.phastcon0)[-1] <- paste(names(dt.feat.phastcon0[-1]), "phastcon0", sep="_")
dt.feat <- merge(dt.feat, dt.feat.phastcon0, by = "Sample")

mssng.prs <- read.delim("../PRS/MSSNG_PRS_Metadata.csv", stringsAsFactors = F, sep=",")
dt.feat <- merge(dt.feat, manifest, by.x = "Sample", by.y = "Sample.ID", all.x = T)
dt.feat <- merge(dt.feat, mssng.prs[, c("SampleID", "SCORESUM")],
                 by.x = "Sample", by.y = "SampleID", all.x = T)
dt <- dt.feat

### ssc
manifest <- read.delim("SSC.eur.info.tsv", stringsAsFactors = F)
selected.fam <- names(which(table(manifest$FAM) == 2))
manifest <- manifest[manifest$FAM %in% selected.fam, ]
dt.feat.phastcon0 <- read.delim("ssc_denovo_features_phastcon0_2021-08-13.tsv", stringsAsFactors = F)
dt.feat <- read.delim("ssc_denovo_features_2021-07-31.tsv", stringsAsFactors = F)
names(dt.feat.phastcon0)[-1] <- paste(names(dt.feat.phastcon0[-1]), "phastcon0", sep="_")
dt.feat <- merge(dt.feat, dt.feat.phastcon0, by = "Sample")

dt.feat <- dt.feat[dt.feat$Sample %in% manifest$Sample.ID, ]
dt.feat <- merge(dt.feat, manifest, by.x = "Sample", by.y = "Sample.ID", all.x = T)
ssc.prs <- read.delim("../PRS/SSC_metadata+PRS_total.csv", stringsAsFactors = F, sep=",")
dt.feat <- merge(dt.feat, ssc.prs[, c("SampleID", "SCORESUM")],
                 by.x = "Sample", by.y = "SampleID", all.x = T)

dt <- rbind(dt, dt.feat)
names(dt) <- gsub("ntorelance|nterolance", "ntolerance", names(dt))
names(dt) <- gsub("Promoter\\.", "Promoter_", names(dt))
names(dt) <- gsub("reducedPromote", "reducedPromoter", names(dt))

# dt <- dt[!is.na(dt$SCORESUM), ]
dt[is.na(dt)] <- 0

outliers <- c("1-0306-004", readLines("denovo_outliers.txt"))
dt <- dt[!dt$Sample %in% outliers, ]

crv <- c(readLines("CRVs/MSSNG+SSC.ASD135_LoF.tsv"),
         readLines("CRVs/MSSNG+SSC.CNVs.tsv"))
# crv <-  readLines("CRVs/MSSNG+SSC.CNVs.tsv")

dt$CRV <- dt$Sample %in% crv
pca <- read.delim("MSSNG_SSC_Sampleinfo.tsv", stringsAsFactors = F)
dt <- merge(dt, pca[, c("Sample.ID", paste0("K", 1:5))], by.x = "Sample", by.y = "Sample.ID", all.x = T)

feats <- readLines("test.features.txt")
feats <- paste0(feats, "_phastcon0")
feats <- c(feats, names(dt)[intersect(grep("DeepSEA", names(dt)), grep("phastcon0", names(dt)))] )
# feats <- names(dt)[67:130]
for(f in feats){
  if(sd(dt[, f]) != 0)
  dt[, f] <- scale(dt[, f])
}
# feats <- names(dt)[3:61]
dt$Affection <- factor(dt$Affection)
before <- ggplot(dt, aes(x = Dataset, y = All_Denovo, fill = Affection)) + 
  geom_boxplot(outlier.alpha = 0, position="dodge", alpha = .5) + ggtitle("Before removing outliers") +
  geom_point(shape = 21, color = "black", position = position_jitterdodge()) + theme_bw() +
  theme(legend.position = "top") + scale_fill_manual(values = c("blue", "red"))

outliers.threshold <- c(quantile(dt$All_Denovo, 0.25) - 3*IQR(dt$All_Denovo),
                        quantile(dt$All_Denovo, 0.75) + 3*IQR(dt$All_Denovo))
dt <- dt[dt$All_Denovo >= outliers.threshold[1] & dt$All_Denovo <= outliers.threshold[2], ]
dt <- dt[dt$Platform == "Illumina HiSeq X", ]
after <- ggplot(dt, aes(x = Dataset, y = All_Denovo, fill = Affection)) + 
  geom_boxplot(outlier.alpha = 0, position="dodge", alpha = .5) + ggtitle("After removing outliers") +
  geom_point(shape = 21, color = "black", position = position_jitterdodge()) + theme_bw() +
  theme(legend.position = "top") + scale_fill_manual(values = c("blue", "red"))

cowplot::plot_grid(before, after, nrow = 1)
ggsave("denovo_distribution.png", width = 11, height = 7)

### start testing
#### simplex vs multiplex
dt.feat <- dt[dt$Dataset == "MSSNG" & dt$Predicted.ancestry == "EUR" & dt$Platform %in% c("Illumina HiSeq X"), ]
dt.feat$Family.type <- factor(dt.feat$Family.type, levels = c("SPX", "MPX"))
dt.out <- data.frame()
for(feat in feats){
  ref <- "Family.type ~ Sex + All_Denovo + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[feat]), 
                                     "low"= exp(conf[feat, 1]), 
                                     "up" = exp(conf[feat, 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
write.table(dt.out, "results/denovo.burden.test.MSSNG.MPX.SPX.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p1 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("SPX(0) vs MPX(1) MSSNG") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) + coord_cartesian(ylim = c(0, 2))
#ggsave("SPXvsMPX.mssng.denovo.png", height = 5, width = 7)

#### SSC proband vs sibs
dt.feat <- dt[dt$Dataset == "SSC" & dt$Predicted.ancestry == "EUR", ]
dt.feat$Affection <- factor(dt.feat$Affection, levels = c(1, 2))
# fam <- dt.feat$Family.ID[duplicated(dt.feat$Family.ID)]
# dt.feat <- dt.feat[dt.feat$Family.ID %in% fam, ]
if(sum(dt.feat$Affection == 2) > 0){
  dt.feat$Affection <- ifelse(dt.feat$Affection == 2, 1, 0)
}

dt.out <- data.frame()
for(feat in feats){
  dt.feat$test <- dt.feat[, feat]
  
  if(sd(dt.feat$test) != 0){
    lm.ref <- clogit(Affection ~ Sex + All_Denovo + CRV + K1 + K2 + K3 + strata(Family.ID), dt.feat)
    lm.alt <- clogit(Affection ~ Sex + All_Denovo + CRV + K1 + K2 + K3 + strata(Family.ID) + test, dt.feat)
    
    test <- anova(lm.ref, lm.alt, test = "Chisq")
    conf <- confint(lm.alt)
    
    dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients["test"]), 
                                       "low"= exp(conf["test", 1]), 
                                       "up" = exp(conf["test", 2]), 
                                       "p"= test$`P(>|Chi|)`[2]))
  }else{
    dt.out <- rbind(dt.out, data.frame(feat, "or" = NA, 
                                       "low"= NA, 
                                       "up" = NA, 
                                       "p"= NA))
  }
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1
dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
write.table(dt.out, "results/denovo.burden.test.SSCProband.SSCUnaffected.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p2 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Unaffected sibs(0) vs Probands(1) SSC") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim = c(0, 2))
cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 2.2))
ggsave("within.dataset.test.denovo.png", height = 8, width = 14)

#### SSC proband vs MSSNG proband
dt.feat <- dt[which((dt$Dataset == "SSC" & dt$Affection == 2) | (dt$Dataset == "MSSNG") & dt$Platform == "Illumina HiSeq X"), ]
dt.feat$Dataset <- factor(dt.feat$Dataset, levels = c("SSC", "MSSNG"))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Dataset ~ Sex + All_Denovo + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[feat]), 
                                     "low"= exp(conf[feat, 1]), 
                                     "up" = exp(conf[feat, 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
write.table(dt.out, "results/denovo.burden.test.MSSNG.MPX.SPX.SSCProband.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p0 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("SSC probands(0) vs MSSNG probands(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) + coord_cartesian(ylim = c(0, 2))

#### SSC proband vs MSSNG SPX
dt.feat <- dt[which((dt$Dataset == "SSC" & dt$Affection == 2) | (dt$Dataset == "MSSNG" & dt$Family.type == "SPX") & dt$Platform == "Illumina HiSeq X"), ]
dt.feat$Family.type <- factor(dt.feat$Family.type, levels = c("-", "SPX"))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Family.type ~ Sex + All_Denovo + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[feat]), 
                                     "low"= exp(conf[feat, 1]), 
                                     "up" = exp(conf[feat, 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
write.table(dt.out, "results/denovo.burden.test.MSSNG.SPX.SSCProband.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p1 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("SSC probands(0) vs MSSNG SPX probands(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) + coord_cartesian(ylim = c(0, 2))

#### SSC proband vs MSSNG MPX
dt.feat <- dt[which((dt$Dataset == "SSC" & dt$Affection == 2) | (dt$Dataset == "MSSNG" & dt$Family.type == "MPX") & dt$Platform == "Illumina HiSeq X"), ]
dt.feat$Family.type <- factor(dt.feat$Family.type, levels = c("-", "MPX"))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Family.type ~ Sex + All_Denovo + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[feat]), 
                                     "low"= exp(conf[feat, 1]), 
                                     "up" = exp(conf[feat, 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}

dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
write.table(dt.out, "results/denovo.burden.test.MSSNG.MPX.SSCProband.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p2 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3, na.rm = T) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("SSC probands(0) vs MSSNG MPX probands(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) + coord_cartesian(ylim = c(0, 2))

cowplot::plot_grid(p0, p1, p2, nrow = 3, rel_heights = c(1, 1, 2.5))
ggsave("mssng.mpx.spx.vs.ssc.proband.denovo.png", height = 8, width = 16)


##############MSSNG and unaffected SSC
######################################
######################################

#### MSSNG proband vs SSC unaffected
dt.feat <- dt[which((dt$Dataset == "SSC" & dt$Affection == 1) | (dt$Dataset == "MSSNG" & dt$Platform == "Illumina HiSeq X")), ]
dt.feat$Affection <- factor(dt.feat$Affection, levels = c(1, 2))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Affection ~ Sex + All_Denovo + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[feat]), 
                                     "low"= exp(conf[feat, 1]), 
                                     "up" = exp(conf[feat, 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
write.table(dt.out, "results/denovo.burden.test.MSSNG.MPX.SPX.SSCUnaffected.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p1 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("SSC unaffected(0) vs MSSNG probands(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) + coord_cartesian(ylim = c(0, 2))

#### MSSNG SSC proband vs SSC unaffected
dt.feat <- dt[dt$Platform == "Illumina HiSeq X", ]
dt.feat$Affection <- factor(dt.feat$Affection, levels = c(1, 2))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Affection ~ Sex +  All_Denovo + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[feat]), 
                                     "low"= exp(conf[feat, 1]), 
                                     "up" = exp(conf[feat, 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- -Inf
dt.out$up[is.na(dt.out$up)] <- Inf
dt.out$feat <- gsub("_phastcon0", "", dt.out$feat)
write.table(dt.out, "results/denovo.burden.test.MSSNG.MPX.SPX.SSCProband.SSCUnaffected.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p2 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("SSC unaffected(0) vs MSSNG/SSC probands(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) + coord_cartesian(ylim = c(0, 2))

cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 1.7))
ggsave("proband.vs.unaffected.denovo.png", height = 8, width = 14)
