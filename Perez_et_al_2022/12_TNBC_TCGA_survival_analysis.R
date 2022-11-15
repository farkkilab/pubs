library(ggplot2)
library(gridExtra)

output.dir="Second_sample_set_noOutliers/"

#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")


setwd("C:/Users/fernpere/HRD/TCGA_analysis/TNBC/")
TNBCA.segs <- read.table(file="TNBCA_preprocessed_segs.tsv", sep="\t", header=TRUE)
TNBCA.segs <- TNBCA.segs[,c(1:5,6,8,9)]
TCGA.clinical.info <- read.csv(file="C:/Users/fernpere/HRD/TCGA_analysis/clinical/TCGA-CDR-SupplementalTableS1_page1_BRCA.tsv", header=TRUE, sep="\t")

BRCA.sample.sheet <- read.table(file="gdc_sample_sheet.2021-06-04_2.tsv", sep="\t", header = TRUE)
BRCA.sample.sheet <- BRCA.sample.sheet[,c(7,8,9)]
BRCA.sample.sheet$Sample.ID <- substr(BRCA.sample.sheet$Sample.ID, 1, 12)
#Only select primary samples
BRCA.sample.sheet <- BRCA.sample.sheet[grep("Primary", BRCA.sample.sheet$Sample.Type),c(1,3)]

#Scanning HRP and HRD samples
HRD.samples <- scan("HRD_samples.txt", what = "string")
HRP.samples <- scan("HRP_samples.txt", what = "string")

#Ignoring outliers
HRP.samples <- HRP.samples[which(!(HRP.samples %in% c("0d7cde44-ec86-415a-8667-0ea894d1c344","f1b4f790-083e-44f0-b924-06dc7c167a20")))]

#Reading drug info for samples
drug.info <- read.table(file="nationwidechildrens.org_clinical_drug_brca.txt", sep="\t", header=TRUE)
platinum.patients <- drug.info[grep("platin",drug.info$pharmaceutical_therapy_drug_name), "bcr_patient_barcode"]
chemo.patients <- drug.info[drug.info$pharmaceutical_therapy_type == "Chemotherapy","bcr_patient_barcode"]

#########################################################################################################################
####################################### Calculating scars for new values ################################################
#########################################################################################################################

TNBCA.scars <- get_HRDs(TNBCA.segs, chrominfo = chrominfo_grch38, LOH_windos=c(10,30), LST_segSize=5e6, LST_mindistance=2e6)
TNBCA.scars.id <- merge(TNBCA.scars, BRCA.sample.sheet, by.x="row.names", by.y="File.Name2")

#Using previous criteria
TNBCA.scars.prevdef <- get_HRDs(TNBCA.segs, chrominfo = chrominfo_grch38, LOH_windos=c(15,1000), LST_segSize=10e6, LST_mindistance=3e6)
TNBCA.scars.prevdef.id <- merge(TNBCA.scars.prevdef, BRCA.sample.sheet, by.x="row.names", by.y="File.Name2")

#Adding HR status
HR.status <- c("HRD","Undefined")[match(TNBCA.scars.id$Row.names %in% HRD.samples, c('TRUE', 'FALSE'))]
HR.status[which(TNBCA.scars.id$Row.names %in% HRP.samples)] <- "HRP"

#Merging previous and new scar values
TNBCA.scars.id$HR.status <- HR.status
TNBCA.scars.id$p.HRDsum <- TNBCA.scars.prevdef.id[,5]

#########################################################################################################################
####################################### Calculating cutoff for new values ###############################################
#########################################################################################################################

#Generating plots for Figure 5d

#Plotting distributions for HRDscar
p <- ggplot(TNBCA.scars.id, aes(x=HRDsum, fill=HR.status)) + geom_histogram(alpha=1, position="dodge", color="black", aes(y = ..density..), binwidth=7)
p <- p + geom_density(alpha=0.3)
p <- p + geom_vline(xintercept = 53, colour="red", linetype = "longdash", size=1)
p <- p + theme(axis.text=element_text(size=rel(1.2)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)),
               legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p +  guides(size = guide_legend(keyheight  = 3))
p <- p + ylab("Density \n") + xlab("\nNumber of HRD scars")
p <- p + scale_fill_manual(labels=c("HRD", "HRP", "Undefined"), values = c("#FA5A41", "#34A0D3", "#00b8384d") )
p <- p + xlim(-2,150) + ylim(0,0.035)
print(p)
ggsave(p, filename = paste0(output.dir, "Density_scars_HRD-HRP.svg"), width = 14, height = 12, units = "cm")


#Plotting distributions for Telli2016
p <- ggplot(TNBCA.scars.id, aes(x=p.HRDsum, fill=HR.status)) + geom_histogram(alpha=1, position="dodge", color="black", aes(y = ..density..), binwidth=7)
p <- p + geom_density(alpha=0.3)
p <- p + geom_vline(xintercept = 42, colour="red", linetype = "longdash", size=1)
p <- p + theme(axis.text=element_text(size=rel(1.2)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)),
               legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p +  guides(size = guide_legend(keyheight  = 3))
p <- p + ylab("Density \n") + xlab("\nNumber of HRD scars")
p <- p + scale_fill_manual(labels=c("HRD", "HRP", "Undefined"), values = c("#FA5A41", "#34A0D3", "#00b8384d") )
p <- p + xlim(-2,150) + ylim(0,0.035)
print(p)
ggsave(p, filename = paste0(output.dir, "Density_scars_HRD-HRP_Telli2016.svg"), width = 14, height = 12, units = "cm")



########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

#Function to get Balanced Accuracy values for Figure 5c and Figure 5d
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

#Performing bootstrapping of samples, to detect ideal cut-off value
cut.off.values <- c(25:95)
Accuracies <- sapply(cut.off.values, function(x){
  newHRDscores_HRDs <- TNBCA.scars.id[TNBCA.scars.id$HR.status %in% "HRD","HRDsum"]
  newHRDscores_HRPs <- TNBCA.scars.id[TNBCA.scars.id$HR.status %in% "HRP","HRDsum"]
  BA.values  <- NULL
  for (i in 1:10000){
    HRD.HRDsum.sampled <- sort(sample(newHRDscores_HRDs, 30, replace = TRUE))
    HRP.HRDsum.sampled <- sort(sample(newHRDscores_HRPs, 20, replace = TRUE))
    BA <- round(get_localdispertion_acc(HRD.HRDsum.sampled, HRP.HRDsum.sampled, status1="HRD", status2="HRP", cutoff=x),2)
    BA.values <- c(BA.values,BA)
  }
  mean.BA <- round(mean(BA.values),3)
  lower.CI.BA <- round(quantile(BA.values, probs = c(0.025)),4)
  upper.CI.BA <- round(quantile(BA.values, probs = c(0.975)),4)
  df.BA <- data.frame(mean= mean.BA, lower.CI=lower.CI.BA, upper.CI=upper.CI.BA)
  return(df.BA)
})

BA.values <- data.frame(Number.scars =cut.off.values,
                        Accuracy= unlist(Accuracies[1,]),
                        lower.CI=unlist(Accuracies[2,]), upper.CI=unlist(Accuracies[3,]))

BA.values[which.max(BA.values$Accuracy),]

#Plotting Figure 5c 
p <- ggplot(BA.values, aes(x=cut.off.values, y=Accuracy))
p <- p + geom_errorbar(aes(ymin=lower.CI, ymax=upper.CI), width=.1, colour="grey50")
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
p <- p + geom_vline(xintercept = 53, colour="red", linetype = "longdash", size=1)
p <- p + xlab("\nHR status cut-off point")
print(p)
ggsave(p, filename = paste0(output.dir, "Cut-off_HRDsatus.svg"), width = 10, height = 10, units = "cm")


#Getting the BA for Figure 5d
round(get_localdispertion_acc(TNBCA.scars.id[TNBCA.scars.id$HR.status %in% "HRD","HRDsum"],
                              TNBCA.scars.id[TNBCA.scars.id$HR.status %in% "HRP","HRDsum"],
                              status1="HRD", status2="HRP",
                              cutoff=53),2)

round(get_localdispertion_acc(TNBCA.scars.id[TNBCA.scars.id$HR.status %in% "HRD","p.HRDsum"],
                              TNBCA.scars.id[TNBCA.scars.id$HR.status %in% "HRP","p.HRDsum"],
                              status1="HRD", status2="HRP",
                              cutoff=42),2)

#########################################################################################################################
####################################### Survival analysis ###############################################################
#########################################################################################################################

TNBCA.scars.id <- TNBCA.scars.id[,c(1:8)]
TNBCA.scars.clin <- merge(TNBCA.scars.id, TCGA.clinical.info, by.x="Sample.ID", by.y="bcr_patient_barcode")

BRCAmutsamples <- scan("BRCAmut_samples.txt", what="character")
TNBCA.scars.clin$BRCAmut <- ifelse(TNBCA.scars.clin$Row.names %in% BRCAmutsamples,1,0)

#Ignoring samples in Stage I,X,Not available
samples.selected <- which(grepl("Stage II", TNBCA.scars.clin$ajcc_pathologic_tumor_stage) |
                          grepl("Stage III", TNBCA.scars.clin$ajcc_pathologic_tumor_stage) |
                          grepl("Stage IV", TNBCA.scars.clin$ajcc_pathologic_tumor_stage))

TNBCA.scars.clin <- TNBCA.scars.clin[samples.selected,]

TNBCA.scars.clin$tnbcHRDscar <- ifelse(TNBCA.scars.clin$HRDsum >= 53,1,0)
TNBCA.scars.clin$Telli2016 <- ifelse(TNBCA.scars.clin$p.HRDsum >= 42,1,0)

TNBCA.scars.clin$PFI.time <- as.numeric(TNBCA.scars.clin$PFI.time)
TNBCA.scars.clin$OS.time <- as.numeric(TNBCA.scars.clin$OS.time)

TNBCA.scars.clin$PFI.time.month <- TNBCA.scars.clin$PFI.time / 30.4
TNBCA.scars.clin$OS.time.month <- TNBCA.scars.clin$OS.time / 30.4

surv_objectPFI <- Surv(time = TNBCA.scars.clin$PFI.time.month, event = TNBCA.scars.clin$PFI)
surv_objectOS <- Surv(time = TNBCA.scars.clin$OS.time.month, event = TNBCA.scars.clin$OS)

#Plotting Figure 5e, 5f, 5g
plot <- make.merge.survplots(surv_objectPFI, TNBCA.scars.clin, variables = c("BRCAmut", "Telli2016" , "tnbcHRDscar"),
                             break_time=50, xlabplot="PFI time (months)", ylabplot="PFI probability", palette = c("#34A0D3", "#FA5A41"))
print(plot)
ggsave(plot, filename = paste0(output.dir, "Kaplan-Meir_PFI.svg"), width = 21, height = 12, units = "cm")

#Plotting Supplementary Figure 5c, 5d, 5e
plot <- make.merge.survplots(surv_objectOS, TNBCA.scars.clin, variables = c("BRCAmut", "Telli2016" , "tnbcHRDscar"),
                             break_time=50, xlabplot="OS time (months)", ylabplot="OS probability", palette = c("#34A0D3", "#FA5A41"))
print(plot)
ggsave(plot, filename = paste0(output.dir, "Kaplan-Meir_OS.svg"), width = 21, height = 12, units = "cm")

################################## Hazard ratios tables for TCGA dataset ########################################

#Next block will generate tables and dot plots for Supplementary Figures 5f, 5g
TNBCA.scars.clin$age <- as.numeric(TNBCA.scars.clin$age_at_initial_pathologic_diagnosis)

formulas <- c('Surv(PFI.time.month, PFI)~ ', 'Surv(OS.time.month, OS)~ age +')
outnames <- c('PFI','OS')
for (j in 1:length(formulas)){
  res.cox <- Cox_regresion_variables(TNBCA.scars.clin, c("BRCAmut", "Telli2016" , "tnbcHRDscar"), formula =formulas[j])
  res.cox$Hazzard_ratio <- paste0(res.cox[,c(2)], "(CI:", res.cox[,c(3)], "-", res.cox[,c(4)], ")")
  res.cox$N <- c(sum(TNBCA.scars.clin$BRCAmut == 1),
                 sum(TNBCA.scars.clin$Telli == 1),
                 sum(TNBCA.scars.clin$tnbcHRDscar == 1))
  res.cox$Prop <- c(sum(TNBCA.scars.clin$BRCAmut == 1)/nrow(TNBCA.scars.clin),
                 sum(TNBCA.scars.clin$Telli == 1)/nrow(TNBCA.scars.clin),
                 sum(TNBCA.scars.clin$tnbcHRDscar == 1)/nrow(TNBCA.scars.clin))
  res.cox$Prop <- paste0(round(res.cox$Prop,2) * 100, "%")
  res.cox$TestName <- rownames(res.cox)
  res.cox.t <- res.cox[,c(10,8,9,7,6)]
  colnames(res.cox.t)[c(4,5)] <- c("Hazard ratio","Pval")
  res.cox.t[1,1] <- "BRCAmut/del"
  
  #Table with Hazard ratios
  tt2 <- ttheme_minimal(base_size = 8, padding = unit(c(1.3, 2.9), "mm"))
  print(paste0(output.dir, "HazardR_scores_", outnames[j], ".pdf"))
  pdf(paste0(output.dir, "HazardR_scores_", outnames[j], ".pdf"), height=10, width=15)
  grid.table(res.cox.t, theme=tt2, rows=NULL)
  dev.off()

  #Generating color dots for the tables
  coxvalues <- res.cox
  coxvalues$p.value <- (-1 * log10(coxvalues$p.value))
  coxvalues$attribute <- factor(rownames(coxvalues), levels = rev(rownames(coxvalues)))
  coxvalues$xval <- rep(1,nrow(coxvalues))

  p <- ggplot(coxvalues, aes(x=xval, y=attribute))
  p <- p + geom_point(aes(color=p.value, size=HR)) + theme_classic()
  p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue4", "deepskyblue","red"),
                                 values=c(0,0.90,1), guide = guide_colourbar(direction = "horizontal"))
  p <- p + scale_radius(range = c(7,1), limits = c(0.25, 1.1), breaks = c(0.40, 0.50, 0.60, 0.7))
  p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
                 axis.text.x=element_blank(), axis.text.y=element_text(size=rel(1)),
                 legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)),
                 axis.title=element_text(size=rel(1)))
  p <- p + ylab("") + xlab("")
  ggsave(p, filename = paste0(output.dir, "Dotpvalue_", outnames[j] ,".svg"), width = 12, height = 5.5, units = "cm")
}