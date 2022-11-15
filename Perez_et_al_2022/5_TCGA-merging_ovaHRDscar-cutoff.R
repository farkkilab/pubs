#In this script I merge the clinical info from TCGA with the biological info and number of scars
#Then I calculate the cut-off point for establishing HRD or HRP samples
#Then also performed some exploratory survival-analysis

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


#Importing other function inside the repository
setwd("C:/Users/fernpere/ovaHRDscar_manuscript_scripts/")
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")

####### Defining variables ###########
outputfolder = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Mainly for plots
setwd(outputfolder)

############## Loading input data ##############
segOVA_TCGA <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signatures/SegmentsOVA_TCGA.txt", header=T)
clinical_info <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/clinical/TCGA-CDR-SupplementalTableS1_page1-OVA.tsv", header=T, sep="\t", row.names=1)
biological_status_samples <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/HRD_status-samples_2021.csv", header = T, sep =",", row.names=1)
TP53mut.samples <- scan(file = "C:/Users/fernpere/HRD/TCGA_analysis/TP53mutations/TP53mutdels", what = "char")

#Convert days to months
clinical_info$PFI.time <- clinical_info$PFI.time / 30.4
clinical_info$OS.time <- clinical_info$OS.time / 30.4
############################################## Calculating scars ##############################################
segOVA_TCGA$ploidy <- rep(2, nrow(segOVA_TCGA))
OVA_TCGA_scars <- get_HRDs(segOVA_TCGA, chrominfo = chrominfo_grch38)
dfOVA_TCGA_scars <- as.data.frame(OVA_TCGA_scars)

#Using previous criteria of Telli2016
OVA_TCGA_scars_previousdefinition <- get_HRDs(segOVA_TCGA, chrominfo = chrominfo_grch38, LOH_windos=c(15,1000), LST_segSize=10e6, LST_mindistance=3e6)
dfOVA_TCGA_scars_previousdefinition <- as.data.frame(OVA_TCGA_scars_previousdefinition)

############################################## Merging scars with clinical information ##############################################

#Selecting relevant columns for further analysis
relevant_columns <- c("age_at_initial_pathologic_diagnosis", "clinical_stage", "histological_grade",
                      "tumor_status", "OS", "OS.time", "PFI", "PFI.time")

clinical_info_selected <- clinical_info[,relevant_columns]

#Renaming stages
stage_III <- which(clinical_info_selected$clinical_stage == "Stage IIIA" | clinical_info_selected$clinical_stage == "Stage IIIB" | clinical_info_selected$clinical_stage == "Stage IIIC")
clinical_info_selected$clinical_stage[stage_III] <- 3
stage_IV <- which(clinical_info_selected$clinical_stage == "Stage IV")
clinical_info_selected$clinical_stage[stage_IV] <- 4

#Selecting the biological information important columns
biological_status_samples_selectedcolumns <- biological_status_samples[,-c(1:4)]

###Merging clinical information with scars detected
#For the new definition
OVA_scars_clinical <- merge(dfOVA_TCGA_scars, clinical_info_selected, by="row.names", all.x = TRUE)
row.names(OVA_scars_clinical) <- OVA_scars_clinical[,1]
OVA_scars_clinical <- OVA_scars_clinical[,-1]

#Merge for the prev definition
OVA_scars_clinical_previousdefinition <- merge(dfOVA_TCGA_scars_previousdefinition, clinical_info_selected, by="row.names", all.x = TRUE)
row.names(OVA_scars_clinical_previousdefinition) <- OVA_scars_clinical_previousdefinition[,1]
OVA_scars_clinical_previousdefinition <- OVA_scars_clinical_previousdefinition[,-1]

####Merging clinical information and scars with biological classification
OVA_scars_clinical_biology <- merge(OVA_scars_clinical, biological_status_samples_selectedcolumns, by="row.names", all.x=TRUE)
OVA_scars_clinical_previousdefinition_biology <- merge(OVA_scars_clinical_previousdefinition, biological_status_samples_selectedcolumns, by="row.names", all.x=TRUE)


####Selecting only G3 or G4, or G2 samples with TPF53mut for HGSC
high.grades <- which(OVA_scars_clinical_biology$histological_grade == "G3" | OVA_scars_clinical_biology$histological_grade == "G4")
tp53.mutants.g2 <- which(OVA_scars_clinical_biology$Row.names %in% TP53mut.samples & OVA_scars_clinical_biology$histological_grade == "G2")
selected.HGSC.samples <- sort(unique(c(high.grades, tp53.mutants.g2)))

OVA_scars_clinical_biologyHGSC <- OVA_scars_clinical_biology[selected.HGSC.samples,]
OVA_scars_clinical_previousdefinition_biologyHGSC <- OVA_scars_clinical_previousdefinition_biology[selected.HGSC.samples, ]

####Marking samples in low stages
#Those samples will be ignored for survival analysis
low.stages.samples <- which(is.na(OVA_scars_clinical_biologyHGSC$clinical_stage)
                           | OVA_scars_clinical_biologyHGSC$clinical_stage == "Stage IB"
                           | OVA_scars_clinical_biologyHGSC$clinical_stage == "Stage IA"
                           | OVA_scars_clinical_biologyHGSC$clinical_stage == "Stage IC"
                           | OVA_scars_clinical_biologyHGSC$clinical_stage == "Stage IIA"
                           | OVA_scars_clinical_biologyHGSC$clinical_stage == "Stage IIB"
                           | OVA_scars_clinical_biologyHGSC$clinical_stage == "Stage IIC")

OVA_scars_clinical_biologyHGSC[c(low.stages.samples), "clinical_stage"] <- "Low stage"
OVA_scars_clinical_previousdefinition_biologyHGSC[c(low.stages.samples), "clinical_stage"] <- "Low stage"

########################################################################################################################################
### Sample TCGA-13-1511 (HRP), is an HRP outlier, we have to tag it as "Undefined" ###
OVA_scars_clinical_biologyHGSC[OVA_scars_clinical_biologyHGSC$Row.names == "TCGA-13-1511","HRDstatus"] <- "Undefined"
OVA_scars_clinical_previousdefinition_biologyHGSC[OVA_scars_clinical_previousdefinition_biologyHGSC$Row.names == "TCGA-13-1511","HRDstatus"] <- "Undefined"

#There are samples with not HRDstatus value, those correspond to samples with SNP-array info, but not mutation records
OVA_scars_clinical_biologyHGSC[which(is.na(OVA_scars_clinical_biologyHGSC$HRDstatus)),"HRDstatus"] <- "Undefined"
OVA_scars_clinical_previousdefinition_biologyHGSC[which(is.na(OVA_scars_clinical_previousdefinition_biologyHGSC$HRDstatus)),"HRDstatus"] <- "Undefined"
#Not_clear tag correspond to undefined

#Storing tables to avoid calculating again the HRD scars, this tables are used also in the script `9_Survival_analysis.R`
write.table(OVA_scars_clinical_biologyHGSC, file=paste0(outputfolder, "TCGA_HRD_newScars_clinical-info_HGSC.csv"), sep=",", row.names = FALSE)
write.table(OVA_scars_clinical_previousdefinition_biologyHGSC, file=paste0(outputfolder, "TCGA_HRD_prevScars_clinical-info_HGSC.csv"), sep=",", row.names = FALSE)

#########################################################################################################################
########################################## Establishing cut-off for ovaHRDscar ##########################################
#########################################################################################################################

#Reading merge of clinical info and ovaHRDscar, add tag of "Not_clear" to the samples with no biological info
HGSC.OVA.scars.clinical.biology <- read.table(paste0(outputfolder,"TCGA_HRD_newScars_clinical-info_HGSC.csv"), sep=",", header=TRUE, row.names = 1)

###Reading previous scars (Telli2016) for plotting distributions and cut-offs
prev.scars.clinical.biology <- read.table(paste0(outputfolder,"TCGA_HRD_prevScars_clinical-info_HGSC.csv"), sep=",", header=TRUE, row.names = 1)


#Jackknife/Bootstraping estimation of intersection using 1000 iterations
intersecting.points <- NULL

#Get as input pre-annotated samples HRD/HRP and the HRDscore values
#Then, samples as labeled as TP, TN using the input cut-off
#values1 are the corresponding for the pre-annotated HRD
#values2 are the corresponding for the pre-annotated HRP
get_localdispertion_acc <- function(values1, values2, cutoff, status1="HRD", status2="HRP"){
  TP <- length(which(values1 >= cutoff))
  FN <- length(which(values1 < cutoff))
  FP <- length(which(values2 >= cutoff))
  TN <- length(which(values2 < cutoff))
  TPR <- TP/(TP+FN)
  TNR <- TN/(TN+FP)
  BA <- (TPR + TNR)/2
  return(BA)
}

cut.off.values <- c(30:80) #Cutoff values to explore
Accuracies <- sapply(cut.off.values, function(x){
  newHRDscores_HRDs <- HGSC.OVA.scars.clinical.biology[HGSC.OVA.scars.clinical.biology$HRDstatus %in% "HRD","HRDsum"]
  newHRDscores_HRPs <- HGSC.OVA.scars.clinical.biology[HGSC.OVA.scars.clinical.biology$HRDstatus %in% "HRP","HRDsum"]
  BA.values  <- NULL
  for (i in 1:10000){
    HRD.HRDsum.sampled <- sort(sample(newHRDscores_HRDs, 29, replace = TRUE))
    HRP.HRDsum.sampled <- sort(sample(newHRDscores_HRPs, 29, replace = TRUE))
    BA <- round(get_localdispertion_acc(HRD.HRDsum.sampled, HRP.HRDsum.sampled, status1="HRD", status2="HRP", cutoff=x),2)
    BA.values <- c(BA.values,BA)
  }
  mean.BA <- round(mean(BA.values),3)
  CI.BA1 <- round(quantile(BA.values, probs = c(0.025)),4)
  CI.BA2 <- round(quantile(BA.values, probs = c(0.975)),4)
  df.BA <- data.frame(mean= mean.BA, CI1=CI.BA1, CI2=CI.BA2)
  return(df.BA)
})

#Generating Fig.2d
BA.values <- data.frame(Number.scars =cut.off.values,
                        Accuracy= unlist(Accuracies[1,]), CI1=unlist(Accuracies[2,]),
                        CI2=unlist(Accuracies[3,]))
p <- ggplot(BA.values, aes(x=cut.off.values, y=Accuracy))
p <- p + geom_errorbar(aes(ymin=CI1, ymax=CI2), width=.1, colour="grey50")
p <- p + geom_line() + geom_point()
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)),
               legend.text=element_text(size=rel(1.1)), legend.title=element_text(size=rel(1.7)),
               axis.title=element_text(size=rel(1.7)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p +  guides(size = guide_legend(keyheight  = 3))
p <- p + ylim(0.5,1)
p <- p + geom_vline(xintercept = 54, colour="red", linetype = "longdash", size=1)
p <- p + xlab("\nHR status cut-off point")
print(p)
ggsave(p, filename = paste0(outputfolder,"HRstatus_cut-off-point.svg"), width = 12, height = 12, units = "cm")


#Accuracy values for Figure 2d
cutoff.ovaHRDscar=54
cutoff.telli2016=42

#This is the balanced accuracy of Telli2016 reported in Fig.2e left panel
prevHRDscores_HRDs <- prev.scars.clinical.biology[prev.scars.clinical.biology$HRDstatus %in% "HRD","HRDsum"]
prevHRDscores_HRPs <- prev.scars.clinical.biology[prev.scars.clinical.biology$HRDstatus %in% "HRP","HRDsum"]
round(get_localdispertion_acc(prevHRDscores_HRDs, prevHRDscores_HRPs, status1="HRD", status2="HRP", cutoff=cutoff.telli2016),2)

#This is the balanced accuracy of ovaHRDscar reported in Fig.2e right panel
prevHRDscores_HRDs <- HGSC.OVA.scars.clinical.biology[HGSC.OVA.scars.clinical.biology$HRDstatus %in% "HRD","HRDsum"]
prevHRDscores_HRPs <- HGSC.OVA.scars.clinical.biology[HGSC.OVA.scars.clinical.biology$HRDstatus %in% "HRP","HRDsum"]
round(get_localdispertion_acc(prevHRDscores_HRDs, prevHRDscores_HRPs, status1="HRD", status2="HRP", cutoff=cutoff.ovaHRDscar),2)

#This is the balanced accuracy of Telli2016 but using a cut-off of 54
prevHRDscores_HRDs <- prev.scars.clinical.biology[prev.scars.clinical.biology$HRDstatus %in% "HRD","HRDsum"]
prevHRDscores_HRPs <- prev.scars.clinical.biology[prev.scars.clinical.biology$HRDstatus %in% "HRP","HRDsum"]
round(get_localdispertion_acc(prevHRDscores_HRDs, prevHRDscores_HRPs, status1="HRD", status2="HRP", cutoff=cutoff.ovaHRDscar),2)


#Plotting right panel for Figure 2e
#Plotting the distribution of ovaHRDscar levels between HRD and HRP samples
#Marking a value of 54
p <- ggplot(HGSC.OVA.scars.clinical.biology, aes(x=HRDsum, fill=HRDstatus)) + geom_histogram(alpha=1, position="dodge", color="black", aes(y = ..density..), binwidth=5.5)
p <- p + geom_density(alpha=0.3)
p <- p + geom_vline(xintercept = 54, colour="red", linetype = "longdash", size=1)
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)),
               legend.text=element_text(size=rel(1.1)), legend.title=element_text(size=rel(1.7)),
               axis.title=element_text(size=rel(1.7)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p +  guides(size = guide_legend(keyheight  = 3))
p <- p + ylab("Density \n") + xlab("\nHRD-AIs")
p <- p + scale_fill_manual(labels=c("HRD", "HRP", "Undefined"), values = c("#FA5A41", "#34A0D3", "#00b8384d"))
p <- p + scale_x_continuous(limits=c(3,193))
p <- p + ylim(0,0.039)
print(p)
ggsave(p, filename = paste0(outputfolder,"Cutoff-54/","HRD-HRP_histogram-cutoff_ovaHRDscar_HGSC.svg"), width = 17, height = 12, units = "cm")


#Plotting left panel for Figure 2e
#Plotting the distribution of HRD-AIs levels (Telli 2016) between HRD and HRP samples
#A value of HRDsum at 42 is the Telli2016 cutoff
p <- ggplot(prev.scars.clinical.biology, aes(x=HRDsum, fill=HRDstatus)) + geom_histogram(alpha=1, position="dodge", color="black", aes(y = ..density..), binwidth=6)
p <- p + geom_density(alpha=0.3)
p <- p + geom_vline(xintercept = 42, colour="red", linetype = "longdash", size=1)
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.5)), axis.text.y=element_text(size=rel(1.5)),
               legend.text=element_text(size=rel(1.1)), legend.title=element_text(size=rel(1.7)),
               axis.title=element_text(size=rel(1.7)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p +  guides(size = guide_legend(keyheight  = 3))
p <- p + ylab("Density \n") + xlab("\nHRD-AIs")
p <- p + scale_fill_manual(labels=c("HRD", "HRP", "Undefined"), values = c("#FA5A41", "#34A0D3", "#00b8384d") )
p <- p + scale_x_continuous(limits=c(3,193))
p <- p + ylim(0,0.039)
print(p)
ggsave(p, filename = paste0(outputfolder,"HRD-HRP_histogram-cutoff_Telli2016_HGSC.svg"), width = 17, height = 12, units = "cm")
