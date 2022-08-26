#!/usr/bin/env Rscript

################################################################
# Generate data/plots for rare variant phenotype distributions #
################################################################

suppressMessages(library(extrafont))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplotify))
suppressMessages(library(ggpubr))

args = commandArgs(TRUE)
file = args[1]

data = read.table(file, header=TRUE, check.names=FALSE, sep="\t")
variant_categories =  c("SNV/indel (dominant)", "SNV/indel (recessive)", "Chr abnormality", "Genomic disorder", "Large or gene-rich CNV", "Uniparental isodisomy", "SV disrupting ASD gene", "TRE", "Mitochondrial variant", "Multiple")

data$variant_category = factor(data$variant_category, levels=c(variant_categories, "None"))

do_plot = function(phenotype_measure) {
    all_pvals = c()
    d = data[data[,phenotype_measure] != 0,]

    for (category in variant_categories) {
        if (!category %in% d$variant_category) {
            next
        }
        d2 = d[d$variant_category == category | d$variant_category == "None",]
        d2$has_variant = d2$variant_category == category
        f = as.formula(paste("has_variant ~ Sex + K1 + K2 + K3 + ", phenotype_measure, sep=""))
        pval = coef(summary(glm(f, data=d2)))[6,4]
        print(summary(glm(f, data=d2)))
        print(pval)
        print("*******************************************")
        all_pvals = c(all_pvals, pval)
    }

    p.adj = p.adjust(all_pvals, method="holm")

    labels = c()
    for (i in 1:length(p.adj)) {
        if (p.adj[i] < 0.05) {
            labels = c(labels, "**")
        } else if (all_pvals[i] < 0.05) {
            labels = c(labels, "*")
        } else {
            labels = c(labels, "")
        }
    }

    x = d[d$variant_category == "None",]

    max_val = max(d[,phenotype_measure])
    min_val = min(d[,phenotype_measure])
    m = median(x[,phenotype_measure])

    p =     ggplot(d, aes_string(x="variant_category", y=phenotype_measure)) +
            geom_boxplot(aes_string(fill="variant_category"), outlier.shape=NA) +
            ylab(gsub("_", " ", phenotype_measure)) +
            theme(axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=1)) +
            theme(panel.background=element_rect(fill = "white")) +
            theme(axis.title.x=element_blank()) +
            theme(panel.border=element_rect(color="black", fill=NA)) +
            #theme(plot.margin = margin(1,1,1,2.3, "cm")) +
            theme(legend.position="none", axis.text.y=element_text(size=12)) +
            scale_y_continuous(breaks=seq(from=0, to=200, by=20), limits=c(min_val, max_val)) +
            theme(panel.grid.major.y = element_line(colour="gray", size=0.1)) +
            annotate("text", size=8, x=1:length(labels), y=c(max(d[,phenotype_measure])), label=labels) +
            geom_hline(yintercept=m, linetype="dashed", color="blue", size=0.3)

    #if (phenotype_measure == "Adaptive_behaviour_standard_score" || phenotype_measure == "Full_scale_IQ") {
    #    p = p + theme(axis.text.x = element_blank())
    #}

    suppressMessages(ggsave(paste(phenotype_measure, ".pdf", sep=""), plot=p, width=6))
    return(p)
}

p1 = do_plot("Adaptive_behaviour_standard_score")
p2 = do_plot("Full_scale_IQ")
p3 = do_plot("Global_ability_composite_estimate")
p4 = do_plot("Socialization_standard_score")

#plot_combined = ggarrange(p1, p2, p3, p4, labels=LETTERS, ncol=2, nrow=2, heights=c(7.5,15))
#plot_combined = ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, heights=c(7.5,11))
plot_combined = ggarrange(p1, p2, p3, p4, ncol=4, nrow=1)
suppressMessages(ggsave("genotype-phenotype.pdf", plot=plot_combined, width=14, height=5))
