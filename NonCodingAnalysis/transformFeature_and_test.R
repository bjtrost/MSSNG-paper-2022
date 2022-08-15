setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(ggplot2)
message(Sys.time())
manifest <- read.delim("cases.info.tsv", stringsAsFactors = F)
dt.denom <- read.delim("transmission.count.data_freq_0.1perc.2022-03-23.tsv", stringsAsFactors = F)
dt.denom <- dt.denom[dt.denom$proband.transmitted > 0, ]

outliers <- c("1-0306-004", readLines("denovo_outliers.txt"))
dt.denom <- dt.denom[!dt.denom$Sample.ID %in% outliers, ]
###plot all
dt.plot <- merge(dt.denom, manifest, by = "Sample.ID", all.x = T)
dt.plot <- na.omit(dt.plot)
dt.plot <- dt.plot[dt.plot$Platform == "Illumina HiSeq X", ]

ggplot(dt.plot, aes(x = proband.nontransmitted, y = proband.transmitted)) + 
  geom_point(shape = 21, color = "black", alpha = .9, aes(fill = Platform)) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_smooth(method = "lm", se = F, lwd=0.5, lty=2, color = "red") + xlab("#non-transmitted variants") +
  ylab("#transmitted variants") + theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave("mssng.all.transmission.png", height = 6, width = 6)


# dt.denom <- dt.denom[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Platform == "Illumina HiSeq X" &
#                                                                   manifest$Library.type == "PCR-free"], ]
dt.denom <- dt.denom[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Platform != "Complete Genomics"], ]
# ratio <- dt.denom$proband.transmitted/dt.denom$proband.nontransmitted
# dt.denom <- dt.denom[ratio <= quantile(ratio, 0.75)+1.5*IQR(ratio) &
#                        ratio >= quantile(ratio, 0.25)-1.5*IQR(ratio), ]
dt.plot <- merge(dt.denom, manifest, by = "Sample.ID", all.x = T)
dt.plot <- na.omit(dt.plot)
# dt.plot <- dt.plot[dt.plot$Platform == "Illumina HiSeq X", ]
ggplot(dt.plot, aes(x = proband.nontransmitted, y = proband.transmitted)) + 
  geom_point(shape = 21, color = "black", alpha = .9, aes(fill = Predicted.ancestry)) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_smooth(method = "lm", se = F, lwd=0.5, lty=2, color = "red") + xlab("#non-transmitted variants") +
  ylab("#transmitted variants") + theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave("mssng.ilmn.transmission.png", height = 6, width = 6)

dt.denom <- dt.denom[dt.denom$proband.transmitted < 40000, ]
names(dt.denom) <- gsub("\\.", "_", names(dt.denom))

write.table(dt.denom, "for_terra_transmission.count.data.ultrarare.2022-01-27.csv",
            sep=",", row.names=F, quote=F, col.names=T)

manifest <- manifest[manifest$Sample.ID %in% dt.denom$Sample.ID, ]
dt.feat <- read.delim("transmission_features_ultra_rare_2021-08-13.tsv", stringsAsFactors = F)
dt.feat <- dt.feat[dt.feat$Sample %in% dt.denom$Sample.ID, ]
names(dt.feat) <- gsub("reducedPromote", "reducedPromoter", names(dt.feat))
names(dt.feat) <- gsub("ntorelance|nterolance", "ntolerance", names(dt.feat))
names(dt.feat) <- gsub("Promoter\\.", "Promoter_", names(dt.feat))
####
ggplot(dt.feat, aes(x = All_Features_Nontransmitted, y = All_Features_Transmitted)) + 
  geom_point(shape = 21, color = "black", alpha = .9) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_smooth(method = "lm", se = F, lwd=0.5, lty=2, color = "red") + xlab("#non-transmitted variants") +
  ylab("#transmitted variants") + theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

feats <- names(dt.feat)
feats <- gsub("_Transmitted", "", feats[grep("Transmitted$", feats)])
feats <- feats[-1]
feats <- feats[-c(8, 18:19, 21, 23, 25:27, 31, 32, 36:40, 42, 46:50, 52)]
writeLines(feats, "test.features.txt")

feats <- readLines("test.features.txt")
feats <- c(paste(feats, c("Transmitted"), sep="_"), paste(feats, c("Nontransmitted"), sep="_"))

dt.feat <- dt.feat[, c("Sample", "All_Features_Transmitted", "All_Features_Nontransmitted", feats)]
names(dt.feat)[1] <- "Sample_ID"
write.table(dt.feat, "for_terra_transmission.features_ultrarare.2022-01-27.csv", sep=",", row.names=F, quote=F, col.names=T)

# test like expansions
### all
baseline <- fisher.test(data.frame("X1" = c(sum(dt.feat$All_Features_Transmitted), 
                                            sum(dt.feat$All_Features_Nontransmitted)),
                                   "X2" = c(sum(dt.denom$proband.transmitted)-sum(dt.feat$All_Features_Transmitted),
                                            sum(dt.denom$proband.nontransmitted) - sum(dt.feat$All_Features_Nontransmitted))))
dt.out <- data.frame()
for(feat in feats){
  feat.transmitted <- sum(dt.feat[, sprintf("%s_Transmitted", feat)])
  feat.nontransmitted <- sum(dt.feat[, sprintf("%s_Nontransmitted", feat)])

  denom.transmitted <- sum(dt.denom$proband.transmitted)-feat.transmitted
  denom.nontransmitted <- sum(dt.denom$proband.nontransmitted)-feat.nontransmitted
  
  # denom.transmitted <- sum(dt.feat$All_Features_Transmitted)-feat.transmitted
  # denom.nontransmitted <- sum(dt.feat$All_Features_Nontransmitted)-feat.nontransmitted
  
  
  test <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                "X2" = c(denom.transmitted, denom.nontransmitted)), or = baseline$estimate)
  testp <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                 "X2" = c(denom.transmitted, denom.nontransmitted)), alternative = "greater")
  dt.out <- rbind(dt.out, data.frame(feat, test$estimate, test$conf.int[1], test$conf.int[2], "p"= testp$p.value))
}

dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- round(p.adjust(dt.out$p, method = "bonferroni"), digits = 2)
dt.out$BHFDR <- round(p.adjust(dt.out$p, method = "BH"), digits = 2)
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

write.table(dt.out, "results/transmission.test.MSSNG.tsv", sep="\t", row.names=F, quote=F, col.names=T)
all <- dt.out

### mpx
baseline <- fisher.test(data.frame("X1" = c(sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "MPX"]]), 
                                            sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "MPX"]])),
                                   "X2" = c(sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "MPX"]])-sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "MPX"]]),
                                            sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "MPX"]]) - sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "MPX"]]))))

dt.out <- data.frame()
for(feat in feats){
  feat.transmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "MPX"], sprintf("%s_Transmitted", feat)])
  feat.nontransmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "MPX"], sprintf("%s_Nontransmitted", feat)])
  
  denom.transmitted <- sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "MPX"]]) - feat.transmitted
  denom.nontransmitted <- sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "MPX"]]) - feat.nontransmitted
  
  
  # denom.transmitted <- sum(dt.feat$All_Features_Transmitted)-feat.transmitted
  # denom.nontransmitted <- sum(dt.feat$All_Features_Nontransmitted)-feat.nontransmitted
  
  
  # test <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
  #                                "X2" = c(denom.transmitted, denom.nontransmitted)), or = baseline$estimate)
  testp <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                  "X2" = c(denom.transmitted, denom.nontransmitted)), alternative = "greater")
  dt.out <- rbind(dt.out, data.frame(feat, test$estimate, test$conf.int[1], test$conf.int[2], "p"= testp$p.value))
}

dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- round(p.adjust(dt.out$p, method = "bonferroni"), digits = 2)
dt.out$BHFDR <- round(p.adjust(dt.out$p, method = "BH"), digits = 2)
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1
write.table(dt.out, "results/transmission.test.MSSNG.MPX.tsv", sep="\t", row.names=F, quote=F, col.names=T)

mpx <- dt.out


### spx
baseline <- fisher.test(data.frame("X1" = c(sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "SPX"]]), 
                                            sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "SPX"]])),
                                   "X2" = c(sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "SPX"]])-sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "SPX"]]),
                                            sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "SPX"]]) - sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "SPX"]]))))

dt.out <- data.frame()
for(feat in feats){
  feat.transmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "SPX"], sprintf("%s_Transmitted", feat)])
  feat.nontransmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Family.type == "SPX"], sprintf("%s_Nontransmitted", feat)])
  
  denom.transmitted <- sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "SPX"]]) - feat.transmitted
  denom.nontransmitted <- sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Family.type == "SPX"]]) - feat.nontransmitted
  
  
  # denom.transmitted <- sum(dt.feat$All_Features_Transmitted)-feat.transmitted
  # denom.nontransmitted <- sum(dt.feat$All_Features_Nontransmitted)-feat.nontransmitted
  
  
  test <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                 "X2" = c(denom.transmitted, denom.nontransmitted)), or = baseline$estimate)
  testp <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                  "X2" = c(denom.transmitted, denom.nontransmitted)), alternative = "greater")
  dt.out <- rbind(dt.out, data.frame(feat, test$estimate, test$conf.int[1], test$conf.int[2], "p"= testp$p.value))
}

dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- round(p.adjust(dt.out$p, method = "bonferroni"), digits = 2)
dt.out$BHFDR <- round(p.adjust(dt.out$p, method = "BH"), digits = 2)
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1
write.table(dt.out, "results/transmission.test.MSSNG.SPX.tsv", sep="\t", row.names=F, quote=F, col.names=T)

spx <- dt.out

p1 <- ggplot(all, aes(x = feat, y = test.estimate, color = BHFDR <= 0.1)) +
  geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("FET MSSNG - 2,395 probands") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim=c(0.75,1.35))

p2 <- ggplot(mpx, aes(x = feat, y = test.estimate, color = BHFDR <= 0.1)) +
  geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("FET MPX - 478 probands") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim=c(0.75,1.35))

p3 <- ggplot(spx, aes(x = feat, y = test.estimate, color = BHFDR <= 0.1)) +
  geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("FET SPX - 1,917 probands") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim=c(0.75,1.35))
cowplot::plot_grid(p1, p2, p3, nrow = 3, rel_heights = c(1, 1, 1.7))
ggsave("FET.mssng.png", height = 8, width =14)

#### simplex vs multiplex
# dt.plot <- merge(dt.denom, manifest, by = "Sample.ID", all.x = T)
# ggplot(dt.plot, aes(x = proband.nontransmitted, y = proband.transmitted)) + 
#   geom_point(shape = 21, color = "black", alpha = .6, aes(fill = Family.type)) +
#   geom_abline(slope = 1, intercept = 0) + theme_bw() +
#   geom_smooth(method = "lm", se = F) + xlab("#non-transmitted variants") +
#   ylab("#transmitted variants") + theme(legend.position = "top") + 
#   guides(fill=guide_legend(nrow=2,byrow=TRUE))
# ggsave("mssng.no.outliers.transmission.spx.mpx.png", height = 6, width = 6)

dt.feat <- merge(dt.feat, dt.denom, by.x = "Sample", by.y= "Sample.ID", all.x = T)
dt.feat <- merge(dt.feat, manifest, by = "Sample", by.y= "Sample.ID", all.x = T)
dt.feat$Family.type <- factor(dt.feat$Family.type, levels = c("SPX", "MPX"))
dt.feat <- dt.feat[dt.feat$Sample %in% dt.denom$Sample.ID, ]
mssng.prs <- read.delim("../PRS/MSSNG_PRS_Metadata.csv", stringsAsFactors = F, sep=",")
dt.feat <- merge(dt.feat, mssng.prs[, c("SampleID", "SCORESUM")],
                 by.x = "Sample", by.y = "SampleID", all.x = T)
crv <- c(readLines("CRVs/MSSNG+SSC.ASD135_LoF.tsv"),
         readLines("CRVs/MSSNG+SSC.CNVs.tsv"))
dt.feat$CRV <- dt.feat$Sample %in% crv
pca <- read.delim("MSSNG_SSC_Sampleinfo.tsv", stringsAsFactors = F)
dt.feat <- merge(dt.feat, pca[, c("Sample.ID", paste0("K", 1:5))], by.x = "Sample", by.y = "Sample.ID", all.x = T)
dt.out <- data.frame()
write.table(dt.feat, "MSSNG.features.tsv", sep="\t", row.names=F, quote=F, col.names=T)

for(f in feats){
  if(sd(dt.feat[, sprintf("%s_Transmitted", f)]) != 0)
    dt.feat[, sprintf("%s_Transmitted", f)] <- scale(dt.feat[, sprintf("%s_Transmitted", f)])
}

for(feat in feats){
  ref <- "Family.type ~ Sex + All_Features_Transmitted + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s_Transmitted", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[sprintf("%s_Transmitted", feat)]), 
                                     "low"= exp(conf[sprintf("%s_Transmitted", feat), 1]), 
                                     "up" = exp(conf[sprintf("%s_Transmitted", feat), 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- dt.out$or[is.na(dt.out$low)]
dt.out$up[is.na(dt.out$up)] <- dt.out$or[is.na(dt.out$up)]
write.table(dt.out, "results/burden.test.MSSNG.MPX.SPX.tsv", sep="\t", row.names=F, quote=F, col.names=T)

ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("SPX(0) vs MPX(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim=c(0, 2))
ggsave("SPXvsMPX.mssng.png", height = 6, width = 12)

##### MSSNG vs SSC
dt <- rbind(read.delim("MSSNG.features.tsv", stringsAsFactors = F),
            read.delim("SSC.features.tsv", stringsAsFactors = F))
dt <- dt[dt$Predicted.ancestry == "EUR" & dt$Platform == "Illumina HiSeq X", ]

for(f in feats){
  if(sd(dt[, sprintf("%s_Transmitted", f)]) != 0)
    dt[, sprintf("%s_Transmitted", f)] <- scale(dt[, sprintf("%s_Transmitted", f)])
}

### MSSNG MPX vs SSC
dt.feat <- dt[(dt$Dataset == "SSC" & dt$Affection == 2) |
                (dt$Dataset == "MSSNG" & dt$Family.type == "MPX"), ]
dt.feat$Family.type <- factor(dt.feat$Family.type, levels = c("-", "MPX"))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Family.type ~ Sex + All_Features_Transmitted + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s_Transmitted", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[sprintf("%s_Transmitted", feat)]), 
                                     "low"= exp(conf[sprintf("%s_Transmitted", feat), 1]), 
                                     "up" = exp(conf[sprintf("%s_Transmitted", feat), 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}

dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- dt.out$or[is.na(dt.out$low)]
dt.out$up[is.na(dt.out$up)] <- dt.out$or[is.na(dt.out$up)]
write.table(dt.out, "results/burden.test.MSSNGMPX.SSCProband.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p1 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_blank(), legend.position = "top") + ggtitle("SSC SPX(0) vs MSSNG MPX(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim = c(0, 2))

### MSSNG SPX vs SSC
dt.feat <- dt[(dt$Dataset == "SSC" & dt$Affection == 2) |
                (dt$Dataset == "MSSNG" & dt$Family.type == "SPX" & dt$Platform != "Complete Genomics"), ]
dt.feat$Family.type <- factor(dt.feat$Family.type, levels = c("-", "SPX"))

dt.out <- data.frame()
for(feat in feats){
  ref <- "Family.type ~ Sex  + All_Features_Transmitted + CRV + K1 + K2 + K3"
  alt <- sprintf("%s + %s_Transmitted", ref, feat)
  
  lm.ref <- glm(ref, dt.feat, family = binomial())
  lm.alt <- glm(alt, dt.feat, family = binomial())
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients[sprintf("%s_Transmitted", feat)]), 
                                     "low"= exp(conf[sprintf("%s_Transmitted", feat), 1]), 
                                     "up" = exp(conf[sprintf("%s_Transmitted", feat), 2]), 
                                     "p"= test$`Pr(>Chi)`[2]))
}

dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1

dt.out$low[is.na(dt.out$low)] <- dt.out$or[is.na(dt.out$low)]
dt.out$up[is.na(dt.out$up)] <- dt.out$or[is.na(dt.out$up)]
write.table(dt.out, "results/burden.test.MSSNGSPX.SSCProband.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p2 <- ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + ggtitle("SSC SPX(0) vs MSSNG SPX(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim = c(0, 2))

cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 1.4))
ggsave("SSC.vs.MSSNG.png", height = 7, width = 7)

dt.plot <- dt[(dt$Dataset == "SSC" & dt$Affection == 2) |
                (dt$Dataset == "MSSNG"), ]
ggplot(dt.plot, aes(x = All_Features_Transmitted, y = proband.transmitted)) + 
  geom_point(shape = 21, color = "black", alpha = .9, aes(fill = Dataset)) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_smooth(aes(color = Dataset), method = "lm", se = F, lwd = 0.5, lty = 2) + xlab("#variants in noncoding features") +
  ylab("#transmitted variants") + theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave("mssng.ssc.probands.transmission.png", height = 6, width = 6)


# ###################
# ###biases test
# dt.out <- data.frame()
# for(i in 1:nrow(dt.feat)){
#   tmp.denom <- dt.denom[dt.denom$Sample.ID == dt.feat$Sample[i], ]
#   tmp.feat <- dt.feat[i, ]
# 
#   feat.out <- c()
#   for(feat in feats){
#     feat.transmitted <- tmp.feat[, sprintf("%s_Transmitted", feat)]
#     feat.nontransmitted <- tmp.feat[, sprintf("%s_Transmitted", feat)] + tmp.feat[, sprintf("%s_Nontransmitted", feat)]
#     # fisher <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
#     #                                  "X2" = c(tmp.denom$proband.transmitted - feat.transmitted,
#     #                                           tmp.denom$proband.nontransmitted - feat.nontransmitted)))
#     # transmission.bias <- signif(fisher$estimate, digits = 3)
# # 
#     transmission.bias <- signif((feat.transmitted/tmp.denom$proband.transmitted)/
#       (feat.nontransmitted/(tmp.denom$proband.nontransmitted + tmp.denom$proband.transmitted)), digits = 4)
#     # transmission.bias <- signif((feat.transmitted/tmp.feat$All_Features_Transmitted)/
#     #                               (feat.nontransmitted/(tmp.feat$All_Features_Nontransmitted + tmp.feat$All_Features_Transmitted)), digits = 4)
# 
#     feat.out <- c(feat.out, transmission.bias)
#   }
# 
# 
#   dt.tmp <-data.frame("Sample" = dt.feat$Sample[i],
#                                      "tmp" = "tmp")
#   dt.tmp[, feats] <- feat.out
#   dt.out <- rbind(dt.out, dt.tmp)
# }
# 
# dt.out <- merge(dt.out, manifest, by.x = "Sample", by.y = "Sample.ID", all.x = T)
# write.table(dt.out, "MSSNG.bias.features.tsv", sep="\t",row.names=F, quote=F, col.names=T)
# 
# ### new test on transmission bias
# test.out <- data.frame()
# for(feat in feats){
#   # test <- t.test(dt.out[, feat], mu = 1)
#   if(sum(!is.na(dt.out[, feat])) >= 3){
#     test <- t.test(dt.out[dt.out$Sample %in% manifest$Sample.ID, feat], mu = 1)
#     wilcoxtest <- t.test(dt.out[dt.out$Sample %in% manifest$Sample.ID, feat], mu = 1, alternative = "greater")
#   
#     test.out <- rbind(test.out, data.frame(feat, test$estimate, test$conf.int[1], test$conf.int[2], wilcoxtest$p.value))
#   }else{
#     tmp <- test.out[1, ]
#     tmp$feat <- feat
#     tmp[, 2:5] <- NA
#     test.out <- rbind(test.out, tmp)
#   }
# }
# 
# test.out$feat <- factor(test.out$feat, levels = test.out$feat)
# test.out$FWER <- p.adjust(test.out$wilcoxtest.p.value, method = "bonferroni")
# # dt.plot <- reshape2::melt(dt.out, id.vars = c("Sample", "tmp"))
# # dt.plot <- merge(dt.plot, test.out, by.x = "variable", by.y = "feat", all.x = T)
# #
# # ggplot(dt.plot, aes(x = variable, y = value, color = FWER < 0.01)) +
# #   geom_boxplot(outlier.alpha = 0) + geom_jitter(width = .2, alpha = .3, shape = 21) +
# #   theme_bw() + xlab("") + ylab("Average odds ratio") +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("MSSNG - 2,157 probands") +
# #   scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1)
# 
# ggplot(test.out, aes(x = feat, y = test.estimate, color = FWER < 0.01)) +
#   geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
#   geom_point() + theme_bw() + xlab("") + ylab("Average odds ratio") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("MSSNG - 2,157 probands") +
#   scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
#   coord_cartesian(ylim=c(0.9, 1.1))
# 
# ggsave(sprintf("result.%s.png", Sys.Date()), width = 10, height = 6)
