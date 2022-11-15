#Validation in Hercules
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
library(stringr)

setwd("C:/Users/fernpere/ovaHRDscar_manuscript_scripts/")

#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")

####### Defining variables ###########
outputfolder = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Mainly for plots

############## Loading input data ##############
Segments <- read.table(file="Z:/Documents/HRD_Hercules/input/combinedAscatSegments-out_NoCLbigseq_allmatched.tsv", header=T)

Hercules_sample_info <- read.table(file="Z:/Documents/HRD_Hercules/input/OvCa_experimentSummary_Oct2020Tis.tsv", header=T, sep="\t")
Hercules_clinical_info <- read.csv(file="Z:/Documents/HRD_Hercules/input/220129_Clinical_binary.csv", header = T, sep =",", row.names=1)
BRCAmut <- scan(file="Z:/Documents/HRD_Hercules/input/BRCAmuted_samples.txt", what = "character")

####Formatting Hercules information
#Remove NA values
Segments <- Segments[which(!is.na(Segments[,10])),]

#Formatting input data
ploy = rep(2, nrow(Segments))

Hercules_segs <- data.frame(SampleID = Segments[,1],
                           Chromosome = paste("chr",Segments[,2],sep=""),
                           Start_position = Segments[,3],
                           End_position = Segments[,4],
                           total_cn = (Segments[,9] + Segments[,10]),
                           A_cn  = Segments[,9],
                           B_cn = Segments[,10],
                           ploidy = ploy)

############## Calculating scars in input segments for Hercules ###############
Hercules_scars <- get_HRDs(Hercules_segs, chrominfo = chrominfo_grch38)
dfHercules_scars <- as.data.frame(Hercules_scars)

#Using previous criteria
Hercules_scars_prevdef <- get_HRDs(Hercules_segs, chrominfo = chrominfo_grch38, LOH_windos=c(15,1000), LST_segSize=10e6, LST_mindistance=3e6)
dfHercules_scars_prevdef <- as.data.frame(Hercules_scars_prevdef)

############## Merging scars levels from ovaHRDscar and previous definition ##############
colnames(dfHercules_scars_prevdef) <- c("pHRD_LOH","pLSTs","pnTAIs","pHRDsum")
pHRDsum <- dfHercules_scars_prevdef[,c(4)]

dfHercules_scars_bothdef <- cbind(dfHercules_scars,pHRDsum)

############## Merging with sample info and clinical information ##############
row.names(Hercules_sample_info) <- Hercules_sample_info$SampleID2

#Test to see if there is any merge
dfHercules_scars_status_info <- merge(dfHercules_scars_bothdef, Hercules_sample_info, by="row.names", all.x = TRUE)
notmerging <- dfHercules_scars_status_info[which(is.na(dfHercules_scars_status_info$Tissue)),1]
notmerging[!grepl("BDNA", notmerging)]
dim(dfHercules_scars_status_info)

colnames(Hercules_clinical_info)[1] <- "Patient"
dfHercules_scars_status_info_clin <- merge(dfHercules_scars_status_info, Hercules_clinical_info, by="Patient")
dim(dfHercules_scars_status_info_clin)

#This file will contain the amount of ovaHRDscars per sample
write.table(dfHercules_scars_status_info_clin, file=paste0(outputfolder, "Hercules_newscars_clinic.csv"), sep=",", row.names = FALSE)

##############################################################################################################################################
##############################################################################################################################################
########### A function used to select the average ovaHRDscar per patient taking for the average the selected tissues/samples #########################
##############################################################################################################################################
##############################################################################################################################################

#For each patient selects the samples for the corresponding selected_tissues.
#If the selected selected_tissues are not present, then take for the selectect tissues 2
#If not selected_tissues2 then mart it as NULL: selected_tissues2=NULL
#If none of the selected tissues are present, then take all the samples
#At the end calculate the average levels per patient of the selected tissues
#Return data frame with the average values per patient

average.levels.by.patient.tissues <- function(data, selected_tissues=c("OVA", "ADN"), selected_tissues2="OME"){
  selected_samples <- NULL
  for (patient in unique(data$Patient)){
    patient_samples <- data[data$Patient == patient,]
    patient_selectedsamples <- NULL
    if(any(patient_samples$Tissue %in% selected_tissues)){
      patient_selectedsamples <- patient_samples[patient_samples$Tissue %in% selected_tissues,]
    }
    if (length(selected_tissues2) > 0 & is.null(patient_selectedsamples) & any(patient_samples$Tissue %in% selected_tissues2)){
      print(paste0("Patient with samples from aditional tissues: ", patient))
      patient_selectedsamples <- patient_samples[patient_samples$Tissue %in% selected_tissues2,]
    }
    if (is.null(patient_selectedsamples)){
      print(paste0("Patient with no tissue of interest: ", patient))
      patient_selectedsamples <- patient_samples
    }
    selected_samples <- rbind(selected_samples, patient_selectedsamples)
  }
  #Getting the average value per patient
  dfHercules_scars_patient <- avg_by_sample(selected_samples)
  return(dfHercules_scars_patient)
}

##################################################################################################################################
########################################### Average info by patient ##############################################################
##################################################################################################################################

outputfolder = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Mainly for plots
setwd(outputfolder)
dfHercules_scars_status_info_clin <- read.table(file=paste0(outputfolder, "Hercules_newscars_clinic.csv"), sep=",", header = T)
BRCAmut <- scan(file="Z:/Documents/HRD_Hercules/input/BRCAmuted_samples.txt", what = "character")

#Select Illumina samples with purity >=0.2
dfHercules_scars_status_info_clin_no <- dfHercules_scars_status_info_clin[which(dfHercules_scars_status_info_clin$WGS.platform != "BGI" & dfHercules_scars_status_info_clin$Max.purity >= 0.20),]

#Selecting primary samples
dfHercules_scars_status_info_clin_pri  <- dfHercules_scars_status_info_clin_no[dfHercules_scars_status_info_clin_no$Sample.time == "primary",]


#Select OVA and ADN samples per patient, if not select then OME, if not OME then the rest of the samples.
#From the samples selected per patient, calculate the average value of ovaHRDscar
dfHercules_scars_patient <- average.levels.by.patient.tissues(dfHercules_scars_status_info_clin_pri)

#Changing tags according to the average
dfHercules_scars_patient$BRCAness <- c("1","0")[match(dfHercules_scars_patient$Patient %in% BRCAmut, c('TRUE', 'FALSE'))]
dfHercules_scars_patient$newHRDscars <- c("1","0")[match(dfHercules_scars_patient$HRDsum >= 54, c('TRUE', 'FALSE'))]
dfHercules_scars_patient$pHRDscars <- c("1","0")[match(dfHercules_scars_patient$pHRDsum >= 42, c('TRUE', 'FALSE'))]
dfHercules_scars_patient$TacaYaHRDscars <- c("1","0")[match(dfHercules_scars_patient$pHRDsum >= 63, c('TRUE', 'FALSE'))]

write.table(dfHercules_scars_patient, file=paste0(outputfolder, "Hercules_newscars_patient_OVA-ADN_HGSC.csv"), sep=",", row.names = FALSE)

###################################################################################################################################
###################################################################################################################################
################################################  Hercules tissue selection for survival comparison ###############################
###################################################################################################################################
###################################################################################################################################
library(reshape2)
library(dplyr)

dfHercules_scars_status_info_clin_pri2 <- dfHercules_scars_status_info_clin_pri
dfHercules_scars_status_info_clin_pri2$BRCAness <- c("1","0")[match(dfHercules_scars_status_info_clin_pri2$Patient %in% BRCAmut, c('TRUE', 'FALSE'))]
Hercules.patients.categories <- read.table(file="Z:/Documents/HRD_Hercules/Konstantinopoulos_piechart/Patient_categories.csv", sep = ",", header=T)

#The average ovaHRDscar value per patient, using all their samples
Hercules1 <- avg_by_sample(dfHercules_scars_status_info_clin_pri2)

#The average ovaHRDscar in OVA, if not OVA samples then the average of the remaining samples
Hercules2 <- average.levels.by.patient.tissues(dfHercules_scars_status_info_clin_pri2, selected_tissues = "OVA", selected_tissues2 = "NULL")

#The average ovaHRDscar in OVA and ADN, if not OVA/ADN samples then the average of the remaining samples
Hercules3 <- average.levels.by.patient.tissues(dfHercules_scars_status_info_clin_pri2, selected_tissues = c("OVA","ADN"), selected_tissues2 = "NULL")

#The average ovaHRDscar in OVA and ADN, if not OVA/ADN samples then the average of the remaining samples
Hercules4 <- average.levels.by.patient.tissues(dfHercules_scars_status_info_clin_pri2, selected_tissues = "OME", selected_tissues2 = c("OVA","ADN"))

#The average ovaHRDscar in OVA and ADN, if not OVA/ADN samples then the average of the remaining samples
Hercules5 <- average.levels.by.patient.tissues(dfHercules_scars_status_info_clin_pri2, selected_tissues = c("OVA", "ADN"), selected_tissues2 = "OME")

#This is to order the patients according to their categories.
#The categories were added according to the genomic alteration
#Order of the categories to appear in Figures:
categories.order <- c("BRCA1 somaticmut", "BRCA2 somaticmut", "CDK12 somaticmut", "FA somaticmut", "RADcore somaticmut",
                      "HR somaticmut", "EMSY amplification", "PTEN somaticmut", "MMR somaticmut", "NER somaticmut", "Other", "CCNE1 amplification" )
Hercules.patients.categories$categories <- factor(Hercules.patients.categories$categories,  levels=c(categories.order))
patient.order <- Hercules.patients.categories[order(Hercules.patients.categories$categories),1]
categories.order <- as.character(Hercules.patients.categories[order(Hercules.patients.categories$categories),2])


####Hercules info, function to change column names in input file
preproccessHercules <- function(datax){
  names(datax)[which(names(datax) == "Age.at.Diagnosis")] <- "age"
  return(datax)
}

####Hercules info, function to ignore patient that received PARPi ####
parp.ignore <- function(dataset){
  dataset.no.parp <- dataset[!dataset$Patient %in% c("H094"),]
  return(dataset.no.parp)
}

#For each of the test to make, Hercules1, Hercules2, etc...
#Run the pre-processing of data, the stratification of patients using the function getHRs and get the Hazard ratios
#Stratification will be stored in the data.frames Hercules1p.strat, Hercules2p.strat, etc...
#Hazard ratio values will be stored in the data.frames Hercules1.HR
for (x in 1:5){
  input <- paste0("Hercules",x)
  outputstrata <- paste0("Hercules",x,"p.strat")
  outputHR <- paste0("Hercules",x,".HR")
  preproccessed <- do.call(preproccessHercules, args = list(eval(parse(text=input))))
  assign(outputstrata,  do.call(getHRs, args=list(preproccessed, formulaHR ='Surv(PFI.time, Progression2)~ ', datalabel="Hercules", returndata = TRUE)))
  preproccessed.noparp <-   do.call(parp.ignore, args=list(preproccessed))
  assign(outputHR,  do.call(getHRs, args=list(preproccessed.noparp, formulaHR ='Surv(PFI.time, Progression2)~ Residual_tumor +', datalabel="Hercules", returndata = FALSE)))
}


###Merge stratification values in data frame
####The next lines are in order to generate Supplementary Figure4d
df <- data.frame(patient=factor(Hercules1p.strat$Patient), PatientCode=paste0("P",seq(1:nrow(Hercules1p.strat))), All=Hercules1p.strat$newHRDhigh, OVA=Hercules2p.strat$newHRDhigh,
                 OVA.ADN=Hercules3p.strat$newHRDhigh, OME_OVA.ADN=Hercules4p.strat$newHRDhigh, OVA.ADN_OME=Hercules5p.strat$newHRDhigh)
df <- melt(df)
df$value <- c("HRD","HRP")[match(df$value == 1, c('TRUE', 'FALSE'))]
#Order factor according to file Hercules.patients.categories
df$patient <- factor(df$patient, levels = patient.order)
df$PatientCode <- factor(df$PatientCode, levels=df$PatientCode[match(patient.order, df$patient)])

####Plotting Supplementary Figure4d
p <- ggplot(df, aes(y=PatientCode, x=variable, fill=value)) + geom_tile (color="black")
p <- p + scale_fill_manual(name = "Status", values=c("#FA5A41", "#34A0D3"))
p <- p + ylab("Patient") + xlab("")
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(size=rel(1), angle = 65, hjust = 1))
#p <- p + scale_y_discrete(labels=categories.order)
print(p)
ggsave(p, filename = paste0(outputfolder, "HR_status_Hercules-Tissue-comparision.svg"), height = 15, width =9, units = "cm")


####Merging Hazard ratio values in matrix for Supplementary Figure4e
datamatrixPFI <- cbind(Hercules1.HR[[2]],  Hercules1.HR[[1]][,c(1,4)], Hercules2.HR[[2]][,c(2,3)], Hercules2.HR[[1]][,c(1,4)],
                       Hercules3.HR[[2]][,c(2,3)], Hercules3.HR[[1]][,c(1,4)], Hercules4.HR[[2]][,c(2,3)], Hercules4.HR[[1]][,c(1,4)],
                       Hercules5.HR[[2]][,c(2,3)], Hercules5.HR[[1]][,c(1,4)])
datamatrixPFI[1,] <- c("", "N", "Prop" , "HR", "Pval", "N", "Prop", "HR", "Pval", "N", "Prop", "HR", "Pval", "N", "Prop", "HR", "Pval", "N", "Prop", "HR", "Pval")
datamatrixPFI[,1] <- c("","BRCAmut/del", "Telli2016", "Tacaya2020", "ovaHRDscar", "Telli2016", "Tacaya2020", "ovaHRDscar")
colnames(datamatrixPFI) <- NULL

#Plotting for Supplementary Figure4e
library("gridExtra")
tt2 <- ttheme_minimal(base_size = 7, padding = unit(c(1.6, 2.9), "mm"))
pdf(paste0(outputfolder, "Hercules_PFI_Tissue-comparision_HazardRatios.pdf"), width = 18, height = 12)       # Export PDF
grid.table(datamatrixPFI, theme=tt2)
dev.off()


###HR dot plots for Supplementary Figure4e
datamatrixPFI <- cbind(Hercules1.HR[[1]][c(5,8),c(1,4)], Hercules2.HR[[1]][c(5,8),c(1,4)], Hercules3.HR[[1]][c(5,8),c(1,4)])
df.ovaHRDscar <- data.frame(
  class = c(rep("Average",6),rep("OVA", 6), rep("OVA&ADN",6), rep("OME;OVA&ADN",6), rep("OVA&ADN;OME",6)),
  Method=c(rep(c("Telli2016", "Tacaya2020", "ovaHRDscar", "Telli2016.2", "Tacaya2020.2","ovaHRDscar.2"), 5)),
  HR = c(Hercules1.HR[[1]][c(3:8),c(1)], Hercules2.HR[[1]][c(3:8),c(1)], Hercules3.HR[[1]][c(3:8),c(1)], Hercules4.HR[[1]][c(3:8),c(1)], Hercules5.HR[[1]][c(3:8),c(1)]),
  pvalues=c(Hercules1.HR[[1]][c(3:8),c(4)], Hercules2.HR[[1]][c(3:8),c(4)], Hercules3.HR[[1]][c(3:8),c(4)], Hercules4.HR[[1]][c(3:8),c(4)], Hercules5.HR[[1]][c(3:8),c(4)])
)

df.ovaHRDscar$Method <- factor(df.ovaHRDscar$Method, levels=rev(c("Telli2016", "Tacaya2020", "ovaHRDscar", "Telli2016.2", "Tacaya2020.2","ovaHRDscar.2")))
df.ovaHRDscar$class <- factor(df.ovaHRDscar$class, levels=c("Average", "OVA", "OVA&ADN", "OME;OVA&ADN", "OVA&ADN;OME"))
df.ovaHRDscar$pvalues <- c(log10(df.ovaHRDscar$pvalues) * -1)

df.values <- data.frame(
  class = c("Average","Average","OVAthenOME", "OVAthenOME", "OMEthenOVA", "OMEthenOVA"),
  Samples=rep(c("Telli2016","ovaHRDscar"), 3),
  pvalues=c(Hercules1.HR[[1]][c(3,5),c(4)], Hercules3.HR[[1]][c(3,5),c(4)], Hercules4.HR[[1]][c(3,5),c(4)]),
  HR = c(Hercules1.HR[[1]][c(3,5),c(1)], Hercules3.HR[[1]][c(3,5),c(1)], Hercules4.HR[[1]][c(3,5),c(1)])
)

breaks = c(0.25, 0.35, 0.45, 0.55)

###Plotting dots for Supplementary Figure4e
p <- ggplot(df.ovaHRDscar, aes(x=class, y=Method))
p <- p + geom_point(aes(color=pvalues, size=HR)) + theme_classic()
p <- p + scale_color_gradientn(name="p.value", colours = c("blue4", "deepskyblue","red"),
                               values=c(0,0.60,1))
p <- p + scale_radius(range = c(7,2), limits = c(0.30, 0.7))
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1), angle = 45, hjust = 1), axis.text.y=element_text(size=rel(1)),
               legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)),
               axis.title=element_text(size=rel(1)))
p <- p + ylab("") + xlab("Selection criteria")
print(p)
ggsave(p, filename = paste0(outputfolder, "Hercules_HRdotspvalues.svg"), width = 11, height = 10, units = "cm")

###############################################################################################################################
###############################################################################################################################
########################################### Paired primary and relapse samples ################################################
###############################################################################################################################
###############################################################################################################################

library(tidyverse)
library(ggpubr)
library(rstatix)

dfHercules_scars_status_info_clin_pur_illumina <- read.table(file="C:/Users/fernpere/HRD/HRD_Hercules/Konstantinopoulos_piechart/Hercules_scars_mutation_classifications.csv", sep=",", header = T)
dfHercules_scars_status_info_clin_pur_illumina[dfHercules_scars_status_info_clin_pur_illumina$Tissue %in%  c("ADN","LN", "OTH","MES"),"Tissue"] <- "OTH"

#Add mutational categories per patient to each of their samples
Patient_category <- NULL
for (patient in unique(dfHercules_scars_status_info_clin_pur_illumina$Patient)){
      cat <- unique(dfHercules_scars_status_info_clin_pur_illumina[dfHercules_scars_status_info_clin_pur_illumina$Patient == patient, "categories"])
      if (length(cat)  == 1){
        Patient_category_aux <- data.frame(category_patient=cat, Patient=patient)
        Patient_category <- rbind(Patient_category, Patient_category_aux)
        next
      }
      cat <- cat[which(cat  != "Other")]
      if (length(cat) > 1){
          print (paste0("Multiple categories in patient:", patient))
      }
        Patient_category_aux <- data.frame(category_patient=cat, Patient=patient)
        Patient_category <- rbind(Patient_category, Patient_category_aux)
}


df <- merge(dfHercules_scars_status_info_clin_pur_illumina, Patient_category, by = "Patient")
df$mutation.present <- !(df$categories == "Other" & df$category_patient != "Other")

categories.order <- c("BRCA1 somaticmut", "BRCA2 somaticmut", "CDK12 somaticmut", "FA somaticmut", "RADcore somaticmut",
                      "HR somaticmut", "EMSY amplification", "PTEN somaticmut", "MMR somaticmut", "NER somaticmut", "CCNE1 amplification", "Other")
df$category_patient <- factor(df$category_patient, levels=c(categories.order))
df <- arrange(df, df$category_patient)
df$Patient <- factor(df$Patient, levels=unique(df$Patient))

anonymousID <- paste0("P",seq(1:length(unique(dfHercules_scars_status_info_clin_pur_illumina$Patient))))
patientID <- data.frame(Patient = unique(dfHercules_scars_status_info_clin_pur_illumina$Patient),
                        ID2 = factor(anonymousID, levels=anonymousID))


df3 <- merge(dfHercules_scars_status_info_clin_pur_illumina, patientID, by = "Patient")

#Making Figure4A
###Plot all those tissues with patients with samples of more than one time point
p <- ggplot(df3, aes(x=ID2, y=HRDsum, color=Sample.time))
p <- p + geom_hline(yintercept = 54, linetype="dashed", color = "darkred")
p <- p + theme(strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3), angle = 40, hjust = 1), axis.text.y=element_text(size=rel(1.6)), legend.text=element_text(size=rel(1)),
               strip.text.x = element_text(size=rel(1)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(2)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + geom_point(aes(y=HRDsum, color=Sample.time, shape=Tissue),  position= position_dodge(0.1),  size=3, alpha=0.9, show.legend = TRUE)
p <- p + ylab("ovaHRDscar") + scale_color_manual(name = "Time", values=c("orange3", "royalblue4", "darkgreen"))
p <- p + scale_shape_manual(values=c(19, 15, 17, 8, 7))
print(p)
ggsave(p, file=paste0(outputfolder, "Hercules_scars-patients-samples.svg"), width = 28, height = 11, units = "cm")


###### Pair test of matched tissues over different treatment time #####
#For Figure4c
patients <- unique(dfHercules_scars_status_info_clin_pur_illumina_paired$Patient)
paried_tissue_time <- NULL

for (x in patients){
  patient_samples <- dfHercules_scars_status_info_clin_pur_illumina_paired[dfHercules_scars_status_info_clin_pur_illumina_paired$Patient==x, ]
  patient_samples$Sample.time <- factor(patient_samples$Sample.time, levels=c("primary", "interval", "relapse"))
  patient_duplicated_tissues <- names(which(table(patient_samples$Tissue) > 1))
  for (t in patient_duplicated_tissues){
    patient_samples_tiss <- patient_samples[patient_samples$Tissue == t,]
    if (length(unique(patient_samples_tiss$Sample.time)) > 1){
      sample <- paste(x, t, sep="_")
      time_count <- 0
      for (j in sort(as.numeric(unique(patient_samples_tiss$Sample.time)))){
        time <- levels(patient_samples_tiss$Sample.time)[j]
        time_count <- time_count + 1
        HRDsum <- mean(patient_samples_tiss[(patient_samples_tiss$Sample.time == time), "HRDsum"])
        #sample_names <- patient_samples_tiss[(patient_samples_tiss$Sample.time == time), "Row.names"]
        data_sample <- data.frame(ovaHRDscar = HRDsum, sampleID=sample, Tx.time=time, tissue=t, Time.id=time_count)
        paried_tissue_time <- rbind(paried_tissue_time, data_sample)
      }
    }
  }
}

paried_tissue_time

stat.test <- paried_tissue_time  %>%
  wilcox_test(ovaHRDscar ~ Time.id, paired = TRUE, detailed = TRUE) %>%
  add_significance()
stat.test

#Plot the paired boxplots
#For Figure4c
b <- runif(nrow(paried_tissue_time), -0.01, 0.01)
p <- ggplot(paried_tissue_time)
p <- p + geom_boxplot(aes(x = as.numeric(Time.id), y = ovaHRDscar, group = Time.id), alpha = 0.8, colour="black")
p <- p +  geom_point(aes(x = as.numeric(Time.id) + b, y = ovaHRDscar, color=Tx.time, shape=tissue), size=4)
p <- p + geom_line(aes(x  = as.numeric(Time.id) + b, y = ovaHRDscar, group = sampleID))
p <- p + theme(strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1)), axis.text.y=element_text(size=rel(1.5)), legend.text=element_text(size=rel(1)),
               strip.text.x = element_text(size=rel(1)),
               legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(1.7)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + scale_x_continuous(breaks = c(1,2), labels = c("No Treatment", "Treatment")) + xlab("Time.id")
p <- p + scale_color_manual(name = "Time", values=c("orange3", "royalblue4", "darkgreen")) + scale_shape_manual(values=c(19, 17, 15, 7))
p <- p + ylim(30,103)
print(p)
ggsave(p, file=paste0(outputfolder, "Hercules_scars_pairedWtest.svg"), width = 10, height = 9, units = "cm")


####Paired samples
relapse_patients <- unique(dfHercules_scars_status_info_clin_pur_illumina[which(dfHercules_scars_status_info_clin_pur_illumina$Sample.time != "primary"),"Patient"])
primary_patients <- unique(dfHercules_scars_status_info_clin_pur_illumina[dfHercules_scars_status_info_clin_pur_illumina$Sample.time %in% "primary","Patient"])
primary_and_relapse_patients <- intersect(relapse_patients, primary_patients)
dfHercules_scars_status_info_clin_pur_illumina_paired <- dfHercules_scars_status_info_clin_pur_illumina[dfHercules_scars_status_info_clin_pur_illumina$Patient %in% primary_and_relapse_patients,]


###Plot all those tissues with patients with samples of more than one time point, color by purity
#For supplementary Figure4b
p <- ggplot(dfHercules_scars_status_info_clin_pur_illumina_paired, aes(x=Patient,y=HRDsum, color=Max.purity))
p <- p + geom_hline(yintercept = 54, linetype="dashed", color = "darkred")
p <- p + theme(strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1), angle = 40, hjust = 1), axis.text.y=element_text(size=rel(1)), legend.text=element_text(size=rel(1)),
               strip.text.x = element_text(size=rel(1)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(1.3)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + geom_point(aes(y=HRDsum, color=Max.purity, shape=Tissue), position= position_dodge(0.3), size= 3, alpha = 0.9, show.legend = TRUE)
p <- p + ylab("ovaHRDscar") + scale_shape_manual(values=c(19, 17, 15, 8, 7))
print(p)
ggsave(p, file=paste0(outputfolder, "Hercules_scars-purity.svg"), width = 15, height = 11, units = "cm")



#####Compare that the intra-patient variation is not higher than the inter-patient variation
#This if for Supplementary Figure 4a
patients <- unique(df3$Patient) #Patients with more than one sample
interpatient.comp <- NULL
for (p in 1:c(length(patients) - 1)){
  patient <- patients[p]
  patient.samples <- df3[df3$Patient == patient,]
  next.p <- p+1
  for (n in 1:nrow(patient.samples)){ #For each sample in the selected patient
    HRDvalue <- patient.samples[n,"HRDsum"]
    samplename1 <- strsplit(patient.samples[n,"Row.names"], "_")[[1]][2]
    ID2 <- as.character(unique(patient.samples$ID2)) ####This is to use the ID2 instead of patient name
    samplename1 <- paste(ID2, samplename1, sep="_")
    for (l in next.p:length(patients)){ #For each of the remaining patients
        n.p <- patients[l]
        n.p.samples <- df3[df3$Patient == n.p,]
        comparisions <- lapply(1:nrow(n.p.samples), function (x){
                        sample.x <- n.p.samples[x,]
                        HRdiff <- abs(HRDvalue - sample.x$HRDsum)
                        samplename2 <- strsplit(sample.x[,"Row.names"], "_")[[1]][2]
                        ID3 <- as.character(unique(sample.x$ID2)) ####This is to use the ID2 instead of patient name
                        samplename2 <- paste(ID3, samplename2, sep="_")
                        samplesnames <- paste(samplename1, samplename2, sep="-")
                        patient1 <- patient
                        patient2 <- n.p
                        res <- c(HRdiff, samplesnames, patient1, patient2)
                        names(res) <- c("HRdiff", "samplesnames", "Patient1", "Patient2")
                        return(res)
                  })
        res <- t(as.data.frame(comparisions, check.names = FALSE))
        res <- as.data.frame(res)
        rownames(res) <- res$samplenames
        interpatient.comp <- rbind(interpatient.comp, res)
    }
  }
}

interpatient.comp$HRdiff <- as.numeric(interpatient.comp$HRdiff)
samples_differences


HRDdifferences <- data.frame(condition=c(rep("Inter-patient", nrow(interpatient.comp)), rep("Intra-patient", nrow(samples_differences))),
                              ovaHRDscar=c(interpatient.comp$HRdiff, samples_differences$diffHRDscars),
                              ids=c(interpatient.comp$samplesnames, samples_differences$samples))
HRDdifferences$condition <- factor(HRDdifferences$condition, levels=c("Intra-patient","Inter-patient"))

u <- wilcox.test(samples_differences$diffHRDscars, interpatient.comp$HRdiff)
u$p.value

#Plotting Supplementary Figure 4a
p <- ggplot(HRDdifferences, aes(x = condition, y = ovaHRDscar))
p <- p + geom_violin(trim=FALSE, size=1) + geom_boxplot(alpha = 0.8, colour="black", outlier.fill = NULL, width=0.3)
p <- p + theme(strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(3), colour = "black"), axis.text.y=element_text(size=rel(2.5)), legend.text=element_text(size=rel(2)),
               strip.text.x = element_text(size=rel(2.5)),
               legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(3)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + ylim(c(0,110)) + ylab("Difference in ovaHRDscar") + xlab("") + scale_x_discrete(labels = c("Intra-patient", "Inter-patient"))
print(p)
ggsave(p, filename = paste0(outputfolder, "Intra_or_inter-patient_diff.svg"), width = 6, height = 7)


###Check if the  HRDscars difference intra patient is caused by purity
###For supplementary Figure 4c
patients <- names(which(table(df3$Patient) > 1)) #Patients with more than one sample
#For each of the patients, calculate the difference in ovaHRDscar (HRDsum) and Max.purity, for each pair of samples
samples_differences <- NULL
for (x in patients){
  comparision <- NULL
  patient_samples <- df3[df3$Patient == x,]
  patientdiffHRDsum <- NULL
  patientdiffpurity <- NULL
  patientcomparisions <- NULL
  for (sample in 1:c(nrow(patient_samples) - 1)){
    l <- sample + 1
    for (j in l:c(nrow(patient_samples))){
      diffHRDsum <- patient_samples[sample,"HRDsum"] - patient_samples[j, "HRDsum"]
      diffpurity <- patient_samples[sample,"Max.purity"] - patient_samples[j, "Max.purity"]
      samplename1 <- strsplit(patient_samples[sample,"Row.names"], "_")[[1]][2]
      samplename2 <- strsplit(patient_samples[j, "Row.names"], "_")[[1]][2]
      ID2 <- as.character(unique(patient_samples$ID2)) ####This is to use the ID2 instead of patient name
      samplesnames <- paste(ID2, samplename1, sep = "_") #This will generate an ID or tag for the difference, like: P4_pOVAL7-pPerL
      samplesnames <- paste(samplesnames, samplename2, sep="-")
      patientdiffHRDsum <- abs(c(patientdiffHRDsum, diffHRDsum))
      patientdiffpurity <- abs(c(patientdiffpurity, diffpurity))
      patientcomparisions <- c(patientcomparisions, samplesnames)
    }
  }
  #Add to data.frame the rest the difference of values for each sample
  comparision <- data.frame(patient=rep(x, length(patientdiffHRDsum)),
                            diffHRDscars = patientdiffHRDsum,
                            puritydiff = patientdiffpurity,
                            samples=patientcomparisions)

  samples_differences <- rbind(samples_differences, comparision)
}


library(ggrepel)
#P4 is an outlier, remove it from plotting
samples_differences <- samples_differences[!grepl("P4",samples_differences$samples),]

#Next lines are just to add the labels of the selected samples in the plot using ggrepel
#Otherwise all samples names are shown and that doesn't look nice
#This do not affect the number of data plotted or statistics
samples_differences2 <- samples_differences
samples_differences2$samples <- ""
selected <- which(samples_differences2$puritydiff > 0.6 | samples_differences2$diffHRDscars > 20)
samples_differences2[selected,"samples"] <- samples_differences[selected,"samples"]

##Plotting supplementary Figure 4c
p <- ggplot(samples_differences2, aes(diffHRDscars, puritydiff, label=samples)) +
  geom_text_repel(size=3) + geom_point()  +
  stat_cor(method = "pearson", label.x.npc="middle") +
  theme(panel.border = element_rect(linetype = "solid", fill = NA),
        axis.text.x=element_text(size=rel(1.5), colour = "black"), axis.text.y=element_text(size=rel(1.5)), legend.text=element_text(size=rel(2)),
        strip.text.x = element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(1.5)),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Difference in ovaHRDscar") + ylab("Difference in Purity")
print(p)
ggsave(p, file=paste0(outputfolder, "Hercules_scars_differenceHRDpurity.svg"), width = 14, height = 12, units = "cm")
