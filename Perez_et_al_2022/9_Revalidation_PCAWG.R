#In PCAWG merge clinical info and calculate scars, then store results
library(ggplot2)
library(rpart)
library(rpart.utils)
library(rpart.plot)
library(ggpubr)
library(survival)
library(survminer)
library(coxphw)
library(mfp)
library(svglite)

####### Defining variables ###########
outputfolder = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/Cutoff-54/" #Mainly for plots

############## Loading input PCAGW data ##############
segOVA_PCAWG <- read.table(file="/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/PCAWG-OVA_segments.txt", header=T)
PCAWG_sample_info <- read.table(file="/home/fernpere/HRD/TCGA_analysis/PCAW/pcawg_Supplementary_Table1.csv", header=T, sep=",", row.names=1)
PCAWG_clinical_info <- read.csv(file="/home/fernpere/HRD/TCGA_analysis/PCAW/pcawg_donor_clinical_August2016_v9.csv", header = T, sep =",", row.names=1)
PCAWG_mutational_drivers <- read.table(file="/home/fernpere/HRD/TCGA_analysis/PCAW/TableS3_panorama_driver_mutations_ICGC_samples.public_OVA.tsv", header=T)

############## Loading input data TCGA (arrays) ##############
segOVA_TCGA <- read.table(file="/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/SegmentsOVA_TCGA.txt", header=T)
clinical_info <- read.table(file="/home/fernpere/HRD/TCGA_analysis/clinical/TCGA-CDR-SupplementalTableS1_page1-OVA.tsv", header=T, sep="\t", row.names=1)
biological_status_samples <- read.table(file="/home/fernpere/HRD/TCGA_analysis/HRD_status-samples.csv", header = T, sep =",", row.names=1)

#Convert to months
clinical_info$PFI.time <- clinical_info$PFI.time / 30.4
clinical_info$OS.time <- clinical_info$OS.time / 30.4

############## Calculating scars in input segments for TCGA (arrays) ##############
segOVA_TCGA$ploidy <- rep(2, nrow(segOVA_TCGA))
#TCGA_scars <- get_HRDs(segOVA_TCGA, chrominfo = chrominfo_grch38)
TCGA_scars <- get_HRDs(segOVA_TCGA, chrominfo = chrominfo_grch38, chrominfo = chrominfo_grch38, LST_segSize=9e6, LST_mindistance=1e6)
dfTCGA_scars <- as.data.frame(TCGA_scars)

############## Calculating scars in input segments for PCAWG ##############
PCAWG_scars <- get_HRDs(segOVA_PCAWG, chrominfo = chrominfo_grch38)
#PCAWG_scars <- get_HRDs(segOVA_PCAWG, chrominfo = chrominfo_grch38, LST_segSize=9e6, LST_mindistance=1e6)
dfPCAWG_scars <- as.data.frame(PCAWG_scars)

#Using previous criteria
PCAWG_scars_prevdef <- get_HRDs(segOVA_PCAWG, chrominfo = chrominfo_grch38, LOH_windos=c(15,1000), LST_segSize=10e6, LST_mindistance=3e6)
dfPCAWG_scars_prevdef <- as.data.frame(PCAWG_scars_prevdef)

############## Merging previous scars values with ovaHRDscar values ##############
HRDscars <- rep("Low_number", nrow(dfPCAWG_scars))

#Cutoff of 54 for new HRDscore
dfPCAWG_scars <- cbind(dfPCAWG_scars, HRDscars)
dfPCAWG_scars$HRDscars[which(dfPCAWG_scars$HRDsum >= 60)] <- "High_number"

#Cutoff of 42 for previous Telli2016 HRDscore
dfPCAWG_scars_prevdef <- cbind(dfPCAWG_scars_prevdef, HRDscars)
dfPCAWG_scars_prevdef$HRDscars[which(dfPCAWG_scars_prevdef$HRDsum >= 42)] <- "High_number"
colnames(dfPCAWG_scars_prevdef) <- c("pHRD_LOH","pLSTs","pnTAIs","pHRDsum","pHRDscars")

dfPCAWG_scars_status <- cbind(dfPCAWG_scars,dfPCAWG_scars_prevdef[,c(4,5)])

################################################ Merging scars with clinical information ##################################################
#Columns of interest from the sample info file
columns_interest <- c(1:8,10,12,34,41,45,46,51,53)
PCAWG_sample_info_subset <- PCAWG_sample_info[,c(columns_interest)]

#Merging with sample info
dfPCAWG_scars_status_info <- merge(dfPCAWG_scars_status, PCAWG_sample_info_subset, by="row.names")
dfPCAWG_scars_status_info_clin <- merge(dfPCAWG_scars_status_info, PCAWG_clinical_info, by="icgc_donor_id")

names(dfPCAWG_scars_status_info_clin)[2] <- "sample_id"

BRCAmutants <- unique(PCAWG_mutational_drivers[which(grepl("BRCA", PCAWG_mutational_drivers$gene)),"sample_id"])
BRCAdrivers <- data.frame(sample_id=BRCAmutants, BRCAness=rep(1,length(BRCAmutants)))

dfPCAWG_scars_status_info_clin_drivers <- merge(dfPCAWG_scars_status_info_clin, BRCAdrivers, by="sample_id", all.x = TRUE)
dfPCAWG_scars_status_info_clin_drivers$BRCAness[which(is.na(dfPCAWG_scars_status_info_clin_drivers$BRCAness))] <- 0

#Change to months
dfPCAWG_scars_status_info_clin_drivers$donor_survival_time <- dfPCAWG_scars_status_info_clin_drivers$donor_survival_time/30.4

write.table(dfPCAWG_scars_status_info_clin_drivers, file=paste0(outputfolder, "PCAWG_newScars_clinical-info_HGSC_LST91.csv"), sep=",", row.names = FALSE)

########################################################################################################################################################
################################################ Measure concordance between Arrays and Whole Genome ###################################################
########################################################################################################################################################

#Measuring the concordance of ovaHRDscar values in TCGA samples also used in PCAWG
library("DescTools")

#Reading input values
dfPCAWG_scars_status_info_clin <- read.table(file=paste0(outputfolder, "PCAWG_newScars_clinical-info_HGSC.csv"), sep=",", header = TRUE)
dfTCGA_scars <- read.table(file=paste0(outputfolder,"TCGA_HRD_newScars_clinical-info_HGSC.csv"), sep=",", header = TRUE, row.names = 1)

#Merging info between TCGA and PCAWG
row.names(dfPCAWG_scars_status_info_clin) <- dfPCAWG_scars_status_info_clin$submitted_donor_id
dfPCAWG_scars_status_info_clin <- dfPCAWG_scars_status_info_clin[,-2]

#41 samples merged
mergePCAWG_TCGA_scars <- merge(dfPCAWG_scars_status_info_clin, dfTCGA_scars, by="row.names")

###### Bland-Altman plot for HRDsum measured in WSG ans SNP-arrays
mean_values <- (mergePCAWG_TCGA_scars$HRDsum.x + mergePCAWG_TCGA_scars$HRDsum.y)/2
differences <- (mergePCAWG_TCGA_scars$HRDsum.x - mergePCAWG_TCGA_scars$HRDsum.y)

mean_dif <- mean(mergePCAWG_TCGA_scars$HRDsum.x - mergePCAWG_TCGA_scars$HRDsum.y)
sd_dif <- sd(abs(mergePCAWG_TCGA_scars$HRDsum.x - mergePCAWG_TCGA_scars$HRDsum.y))
CIpositive  <- mean_dif + 1.96 * sd_dif
CInegative  <- mean_dif - 1.96 * sd_dif

#Calculating the  Lin's concordance correlation coefficient for agreement
#This value is part of Supplementary Fig. 2h
HRDsum_ccc <- CCC(mergePCAWG_TCGA_scars$HRDsum.x, mergePCAWG_TCGA_scars$HRDsum.y)
HRDsum_ccc$rho.c

##Plot for article, Supplementary Fig. 2h
df <- data.frame(meanvalues = mean_values, diff = differences)
p <- ggplot(df,aes(meanvalues, diff)) + geom_point()
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)), legend.text=element_text(size=rel(1)),
               strip.text.x = element_text(size=rel(1)),
               legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(1.5)),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "white"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + labs(x="(ovaHRDscar(WGS) + ovaHRDscar(SNP-array))/2 ", y = "ovaHRDscar(WGS) - ovaHRDscar(SNP-array)")
p <- p + geom_hline(yintercept=0, color = "black")
p <- p + geom_hline(yintercept=mean_dif, color = "blue")
p <- p + geom_hline(yintercept=CIpositive, linetype="dashed", color = "red")
p <- p + geom_hline(yintercept=CInegative, linetype="dashed", color = "red")
p <- p + scale_y_continuous(breaks = round(seq(min(-20), max(10), by = 5),1), limits = c(-20, 11))
print(p)
ggsave(p, filename = paste0(outputfolder, "Bland_TCGA_PCAWG_HRDsum.svg"), width = 11, height = 11, units = "cm")


##Not in article, correlation between WGS and SNP-arrays
p <- ggplot(mergePCAWG_TCGA_scars, aes(x=HRDsum.x, y=HRDsum.y)) + geom_point()
p <- p +  geom_smooth(method = "lm") + stat_cor(method="spearman")
p <- p + ylab("HRDsum WGS") + xlab("HRDsum arrays")
p <- p + theme(axis.text=element_text(size=rel(1.5)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)), legend.text=element_text(size=rel(1.4)),
               strip.text.x = element_text(size=rel(2.8)),
               legend.title=element_text(size=rel(2)), axis.title=element_text(size=rel(2.5)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"))
p <- p + geom_abline(intercept = 0, slope = 1, color="red",  size=0.7)
p <- p + xlim(8,86) + ylim(8,86)
print(p)
ggsave(p, filename = "/home/fernpere/HRD/Figures/Correlations_TCGA_PCAWG_HRDsum.png", width = 15, height = 15, units = "cm")
