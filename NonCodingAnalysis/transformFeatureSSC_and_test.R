setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(data.table)
library(ggplot2)
library(survival)
message(Sys.time())
manifest <- read.delim("SSC.info.tsv", stringsAsFactors = F)
selected.fam <- names(which(table(manifest$FAM) == 2))
manifest <- manifest[manifest$FAM %in% selected.fam, ]

dt.denom <- read.delim("ssc.transmission.count.data_freq_0.1perc.2021-08-13.tsv", stringsAsFactors = F)
dt.denom <- dt.denom[dt.denom$proband.transmitted > 0 & dt.denom$Sample.ID %in% manifest$Sample.ID, ]

outliers <- c("1-0306-004", readLines("denovo_outliers.txt"))
dt.denom <- dt.denom[!dt.denom$Sample.ID %in% outliers, ]

# ratio <- dt.denom$proband.transmitted/dt.denom$proband.nontransmitted
# dt.denom <- dt.denom[ratio <= quantile(ratio, 0.75)+1.5*IQR(ratio) &
#                        ratio >= quantile(ratio, 0.25)-1.5*IQR(ratio), ]

# outlier <- quantile(dt.denom$proband.transmitted, 0.75)+3*IQR(dt.denom$proband.transmitted)
# dt.denom <- dt.denom[dt.denom$proband.transmitted <= outlier, ]

dt.plot <- merge(dt.denom, manifest, by = "Sample.ID", all.x = T)
ggplot(dt.plot, aes(x = proband.nontransmitted, y = proband.transmitted)) + 
  geom_point(shape = 21, color = "black", alpha = .9, aes(fill = Predicted.ancestry)) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_smooth(method = "lm", se = F, lwd=0.5, lty=2, color = "red") + xlab("#non-transmitted variants") +
  ylab("#transmitted variants") + theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave("ssc.all.transmission.png", height = 6, width = 6)

dt.denom <- dt.denom[dt.denom$proband.transmitted < 40000, ]
manifest <- manifest[manifest$Sample.ID %in% dt.denom$Sample.ID, ]

dt.feat <- read.delim("ssc.transmission_features_ultra_rare_2021-08-13.tsv", stringsAsFactors = F)
dt.feat <- dt.feat[dt.feat$Sample %in% dt.denom$Sample.ID, ]
names(dt.feat) <- gsub("reducedPromote", "reducedPromoter", names(dt.feat))
names(dt.feat) <- gsub("ntorelance|nterolance", "ntolerance", names(dt.feat))
names(dt.feat) <- gsub("Promoter\\.", "Promoter_", names(dt.feat))

###
ggplot(dt.feat, aes(x = All_Features_Nontransmitted, y = All_Features_Transmitted)) + 
  geom_point(shape = 21, color = "black", alpha = .9) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_smooth(method = "lm", se = F, lwd=0.5, lty=2, color = "red") + xlab("#non-transmitted variants") +
  ylab("#transmitted variants") + theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

feats <- names(dt.feat)
feats <- gsub("_Transmitted", "", feats[grep("Transmitted$", feats)])
feats <- feats[-1]
### test like expansions
feats <- readLines("test.features.txt")
baseline <- fisher.test(data.frame("X1" = c(sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2]]), 
                                            sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2]])),
                                   "X2" = c(sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 2]])-sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2]]),
                                            sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 2]]) - sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2]]))))

dt.out <- data.frame()
for(feat in feats){
  feat.transmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2], sprintf("%s_Transmitted", feat)])
  feat.nontransmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2], sprintf("%s_Nontransmitted", feat)])
  
  denom.transmitted <- sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 2]]) - feat.transmitted
  denom.nontransmitted <- sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 2]]) - feat.nontransmitted
  
  # denom.transmitted <- sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2]])-feat.transmitted
  # denom.nontransmitted <- sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 2]])-feat.nontransmitted
  
  test <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                 "X2" = c(denom.transmitted, denom.nontransmitted)), or = baseline$estimate)
  testp <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                  "X2" = c(denom.transmitted, denom.nontransmitted)), alternative = "greater")
  dt.out <- rbind(dt.out, data.frame(feat, test$estimate, test$conf.int[1], test$conf.int[2], "p"= testp$p.value))
}

proband <- dt.out
proband$feat <- factor(proband$feat, levels = proband$feat)
proband$FWER <- p.adjust(proband$p, method = "bonferroni")
proband$BHFDR <- p.adjust(proband$p, method = "BH")
proband$BHFDR[is.na(proband$BHFDR)] <- 1
write.table(proband, "results/transmission.test.SSCProband.tsv", sep="\t", row.names=F, quote=F, col.names=T)


baseline <- fisher.test(data.frame("X1" = c(sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1]]), 
                                            sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1]])),
                                   "X2" = c(sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 1]])-sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1]]),
                                            sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 1]]) - sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1]]))))

dt.out <- data.frame()
for(feat in feats){
  feat.transmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1], sprintf("%s_Transmitted", feat)])
  feat.nontransmitted <- sum(dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1], sprintf("%s_Nontransmitted", feat)])
  
  denom.transmitted <- sum(dt.denom$proband.transmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 1]]) - feat.transmitted
  denom.nontransmitted <- sum(dt.denom$proband.nontransmitted[dt.denom$Sample.ID %in% manifest$Sample.ID[manifest$Affection == 1]]) - feat.nontransmitted
  
  # denom.transmitted <- sum(dt.feat$All_Features_Transmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1]])-feat.transmitted
  # denom.nontransmitted <- sum(dt.feat$All_Features_Nontransmitted[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1]])-feat.nontransmitted
  test <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                 "X2" = c(denom.transmitted, denom.nontransmitted)), or = baseline$estimate)
  testp <- fisher.test(data.frame("X1" = c(feat.transmitted, feat.nontransmitted),
                                  "X2" = c(denom.transmitted, denom.nontransmitted)), alternative = "greater")
  
  ###regression
  case <- dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1], c("Sample", sprintf("%s_Transmitted", feat))]
  ctrl <- dt.feat[dt.feat$Sample %in% manifest$Sample.ID[manifest$Affection == 1], c("Sample", sprintf("%s_Nontransmitted", feat))]
  case$status <- 1
  ctrl$status <- 0
  
  names(case) <- names(ctrl) <- c("Sample", feat, "status")
  dt <- rbind(case, ctrl)
  dt <- merge(dt, dt.denom, by.x = "Sample", by.y = "Sample.ID", all.x = T)
  dt$all_count <- dt$proband.transmitted
  dt$all_count[dt$status == 0] <- dt$proband.nontransmitted[dt$status == 0]
  
  ref <- glm(status ~ all_count, dt, family = binomial())
  add <- glm(sprintf("status ~ all_count + %s", feat), dt, family = binomial())
  p <- anova(ref, add, test = "Chisq")$"Pr(>Chi)"[2]
  
  dt.out <- rbind(dt.out, data.frame(feat, test$estimate, test$conf.int[1], test$conf.int[2], "p"= testp$p.value))
}

sibs <- dt.out
sibs$feat <- factor(sibs$feat, levels = sibs$feat)
sibs$FWER <- p.adjust(sibs$p, method = "bonferroni")
sibs$BHFDR <- p.adjust(sibs$p, method = "BH")
sibs$BHFDR[is.na(sibs$BHFDR)] <- 1
write.table(sibs, "results/transmission.test.SSCUnaffected.tsv", sep="\t", row.names=F, quote=F, col.names=T)

p1 <- ggplot(proband, aes(x = feat, y = test.estimate, color = BHFDR <= 0.1)) + 
  geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Average odds ratio") +
  theme(axis.text.x = element_blank()) + ggtitle("SSC - 1,931 probands") + 
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim=c(0.75,1.35))

p2 <- ggplot(sibs, aes(x = feat, y = test.estimate, color = BHFDR <= 0.1)) + 
  geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Average odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("SSC - 1,931 unaffected siblings") + 
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim=c(0.75,1.35))

cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 1.7))
ggsave(sprintf("ssc.result.%s.png", Sys.Date()), width = 14, height = 8)

#### proband vs sibs
dt.feat <- merge(dt.feat, dt.denom, by.x = "Sample", by.y= "Sample.ID", all.x = T)
dt.feat <- merge(dt.feat, manifest, by = "Sample", by.y= "Sample.ID", all.x = T)
dt.feat$Affection <- factor(dt.feat$Affection, levels = c("1", "2"))
#### display significant result
plot(dt.feat$GeneHancer_LOFIntolerance_Transmitted[dt.feat$Affection == 1], dt.feat$GeneHancer_LOFIntolerance_Nontransmitted[dt.feat$Affection == 1])
abline(0, 1)
ssc.prs <- read.delim("../PRS/SSC_metadata+PRS_total.csv", stringsAsFactors = F, sep=",")
dt.feat <- merge(dt.feat, ssc.prs[, c("SampleID", "SCORESUM")],
                 by.x = "Sample", by.y = "SampleID", all.x = T)
crv <- c(readLines("CRVs/MSSNG+SSC.ASD135_LoF.tsv"),
         readLines("CRVs/MSSNG+SSC.CNVs.tsv"))
dt.feat$CRV <- dt.feat$Sample %in% crv

pca <- read.delim("MSSNG_SSC_Sampleinfo.tsv", stringsAsFactors = F)
dt.feat <- merge(dt.feat, pca[, c("Sample.ID", paste0("K", 1:5))], by.x = "Sample", by.y = "Sample.ID", all.x = T)
write.table(dt.feat, "SSC.features.tsv", sep="\t", row.names=F, quote=F, col.names=T)

for(f in feats){
  if(sd(dt.feat[, sprintf("%s_Transmitted", f)]) != 0)
    dt.feat[, sprintf("%s_Transmitted", f)] <- scale(dt.feat[, sprintf("%s_Transmitted", f)])
}

dt.out <- data.frame()
if(sum(dt.feat$Affection == 2) > 0){
  dt.feat$Affection <- ifelse(dt.feat$Affection == 2, 1, 0)
}

for(feat in feats){
  # ref <- "Affection ~ Sex + All_Features_Transmitted + CRV + K1 + K2 + K3"
  # alt <- sprintf("%s + %s_Transmitted", ref, feat)
  #lm.ref <- glm(ref, dt.feat, family = binomial())
  #lm.alt <- glm(alt, dt.feat, family = binomial())
  dt.feat$test <- dt.feat[, sprintf("%s_Transmitted", feat)]
  lm.ref <- clogit(Affection ~ Sex + All_Features_Transmitted + CRV + K1 + K2 + K3 + strata(Family.ID), dt.feat)
  lm.alt <- clogit(Affection ~ Sex + All_Features_Transmitted + CRV + K1 + K2 + K3 + strata(Family.ID) + test, dt.feat)
  
  test <- anova(lm.ref, lm.alt, test = "Chisq")
  conf <- confint(lm.alt)
  
  dt.out <- rbind(dt.out, data.frame(feat, "or" = exp(lm.alt$coefficients["test"]), 
                                     "low"= exp(conf["test", 1]), 
                                     "up" = exp(conf["test", 2]), 
                                     "p"= test$`P(>|Chi|)`[2]))
}


dt.out$feat <- factor(dt.out$feat, levels = dt.out$feat)
dt.out$FWER <- p.adjust(dt.out$p, method = "bonferroni")
dt.out$BHFDR <- p.adjust(dt.out$p, method = "BH")
dt.out$low[is.na(dt.out$low)] <- dt.out$or[is.na(dt.out$low)]
dt.out$up[is.na(dt.out$up)] <- dt.out$or[is.na(dt.out$up)]
dt.out$BHFDR[is.na(dt.out$BHFDR)] <- 1
write.table(dt.out, "results/burden.test.SSCProband.SSCUnaffected.tsv", sep="\t", row.names=F, quote=F, col.names=T)
dt.out

ggplot(dt.out, aes(x = feat, y = or, color = BHFDR <= 0.25)) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0, size = 3, alpha = .3) +
  geom_point() + theme_bw() + xlab("") + ylab("Odds ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("SSC unaffected sibs(0) vs probands(1)") +
  scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1) +
  coord_cartesian(ylim = c(0, 2))
ggsave("Proband.vs.Sibs.SSC.png", height = 6, width = 12)



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
# 
#     transmission.bias <- signif((feat.transmitted/tmp.denom$proband.transmitted)/
#       (feat.nontransmitted/(tmp.denom$proband.nontransmitted + tmp.denom$proband.transmitted)), digits = 4)
#     # transmission.bias <- signif((feat.transmitted/tmp.feat$All_Features_Transmitted)/
#     #   (feat.nontransmitted/(tmp.feat$All_Features_Nontransmitted + tmp.feat$All_Features_Transmitted)), digits = 4)
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
# write.table(dt.out, "SSC.bias.features.tsv", sep="\t",row.names=F, quote=F, col.names=T)
# 
# ### new test on transmission bias for probands
# test.out <- data.frame()
# for(feat in feats){
#   # test <- t.test(dt.out[, feat], mu = 1)
#   if(sum(!is.na(dt.out[dt.out$Sample %in% manifest$Sample.ID[manifest$Affection == 2], feat])) >= 3){
#     
#     test <- t.test(dt.out[dt.out$Sample %in% manifest$Sample.ID[manifest$Affection == 2], feat], mu = 1)
#     wilcoxtest <- t.test(dt.out[dt.out$Sample %in% manifest$Sample.ID[manifest$Affection == 2], feat], mu = 1, alternative = "greater")
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
# 
# p1 <- ggplot(test.out, aes(x = feat, y = test.estimate, color = FWER < 0.01)) +
#   geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
#   geom_point() + theme_bw() + xlab("") + ylab("Average odds ratio") +
#   theme(axis.text.x = element_blank()) + ggtitle("SSC - 1,499 probands") +
#   scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1)
# 
# ### new test on transmission bias for unaffected sibs
# test.out <- data.frame()
# for(feat in feats){
#   # test <- t.test(dt.out[, feat], mu = 1)
#   if(sum(!is.na(dt.out[dt.out$Sample %in% manifest$Sample.ID[manifest$Affection == 2], feat])) >= 3){
#     test <- t.test(dt.out[dt.out$Sample %in% manifest$Sample.ID[manifest$Affection == 1], feat], mu = 1)
#     wilcoxtest <- t.test(dt.out[dt.out$Sample %in% manifest$Sample.ID[manifest$Affection == 1], feat], mu = 1, alternative = "greater")
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
# 
# p2 <- ggplot(test.out, aes(x = feat, y = test.estimate, color = FWER < 0.01)) +
#   geom_errorbar(aes(ymin = test.conf.int.1., ymax = test.conf.int.2.), width = 0, size = 3, alpha = .3) +
#   geom_point() + theme_bw() + xlab("") + ylab("Average odds ratio") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("SSC - 1,499 unaffected siblings") +
#   scale_color_manual(values = c("#333333", "red")) + geom_hline(lty=2, yintercept = 1)
# 
# cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 1.7))
# ggsave(sprintf("ssc.sibs.result.%s.png", Sys.Date()), width = 10, height = 8)
