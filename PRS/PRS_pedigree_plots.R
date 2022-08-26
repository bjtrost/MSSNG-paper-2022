#!/usr/bin/env Rscript

##############################
# Generate the PRS pedigrees #
##############################

suppressMessages(library(kinship2))

### PRS pedigree for FAM_1-0627-007 ###
mom = "10.6"
dad = "-4.7"
c1  = "4.7"
c2  = "-0.3"
c3  = "-1.7"
c4  = "9.3"
c5 = "7.7"
IDs     = c(mom, dad, c1, c2, c3, c4, c5)
dad_IDs = c(NA, NA, dad, dad, dad, dad, dad)
mom_IDs = c(NA, NA, mom, mom, mom, mom, mom)
sex     = c("female", "male", "male", "male", "male", "male", "male")
affected = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
ped     =  pedigree(IDs, dad_IDs, mom_IDs, sex, affected)

cat("Generating plot PRS.pedigree.FAM_1-0627-007.pdf...\n")
pdf("pdf/PRS.pedigree.FAM_1-0627-007.pdf", height=3, width=4)
plot(ped)
invisible(dev.off())

### PRS pedigree for FAM_AU3889305 ###
mom = "2.9"
dad = "0.9"
c1  = "4.6"
c2  = "3.5"
c3  = "1.3"
c4  = "7.7"
IDs     = c(mom, dad, c1, c2, c3, c4)
dad_IDs = c(NA, NA, dad, dad, dad, dad)
mom_IDs = c(NA, NA, mom, mom, mom, mom)
sex     = c("female", "male", "female", "female", "male", "female")
affected = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
ped     =  pedigree(IDs, dad_IDs, mom_IDs, sex, affected)

cat("Generating plot PRS.pedigree.FAM_AU3889305.pdf...\n")
pdf("pdf/PRS.pedigree.FAM_AU3889305.pdf", height=3, width=4)
plot(ped)
invisible(dev.off())
