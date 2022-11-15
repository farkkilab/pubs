library(survival)
library(survminer)
library("extrafont")
library(gridExtra)
font_import()
loadfonts()
library(reshape2)

setwd("C:/Users/fernpere/ovaHRDscar_manuscript_scripts/")
#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")

#For reproducibility of results
set.seed(1979)  


################################################## Defining variables ######################################################

#Input common folder with files with scars info and clinical info
input.folder <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/"
#Output folder for plots and
output.folder <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/Survival_54/"

setwd(input.folder)

#OVA-TCGA scars and clinical info paths
TCGA.scars.clinical.new.file <-  "TCGA_HRD_newScars_clinical-info_HGSC.csv"
#OVA-TCGA previous version of scars
TCGA.scars.clinical.prev.file <- "TCGA_HRD_prevScars_clinical-info_HGSC.csv"
#PCAWG
PCAWG.scars.file <- "PCAWG_newScars_clinical-info_HGSC.csv"
#HERCULES
HERCULES.scars.file <- "Hercules_newscars_patient_OVA-ADN_HGSC.csv"

#Extra clinical info for PCAWG samples from the OV-AU cohort
PCAWG.clin.file <- "C:/Users/fernpere/HRD/TCGA_analysis/PCAW/donor.OV-AU.tsv"

#TCGA residual tumor after surgery, info for each sample
residual.file <- "C:/Users/fernpere/HRD/TCGA_analysis/clinical/TCGA-OV-Clinical-Table_S1.2.csv"

#List of patients that received platinum (cis-plating or carboplatin) in TCGA
platin.takers.list.file <- "C:/Users/fernpere/HRD/TCGA_analysis/clinical/Platin-receivers.txt"

#Driver mutations in PCAWG
mutations.driver.file <- "C:/Users/fernpere/HRD/TCGA_analysis/PCAW/TableS3_panorama_driver_mutations_ICGC_samples.public_OVA.tsv"

###################################### Loading input data and adjusting columns ###########################################

######Loading OVA-TCGA new scars info
TCGA.scars.clinical.new <- read.table(file=paste0(input.folder,TCGA.scars.clinical.new.file), sep=",", header=T, row.names=1)
names(TCGA.scars.clinical.new)[5] <- "age"
TCGA.scars.clinical.prev <- read.table(file=paste0(input.folder,TCGA.scars.clinical.prev.file), sep=",", header=T, row.names=1)
names(TCGA.scars.clinical.prev)[c(1,2,3,4)] <- c("pHRD_LOH", "pLSTs", "pnTAIs", "pHRDsum")
TCGA.scars.clinical.prev <- TCGA.scars.clinical.prev[,c(1:4)]
#Merging values
TCGA.scars.clinical <- cbind(TCGA.scars.clinical.new, TCGA.scars.clinical.prev)

#For this analysis we ignored samples without PFI info and in low_stages (Low stages are I and II)
TCGA.scars.clinical <- TCGA.scars.clinical[-c(which(is.na(TCGA.scars.clinical$PFI.time))),]
TCGA.scars.clinical <- TCGA.scars.clinical[which(TCGA.scars.clinical$clinical_stage != "Low stage"),]

#Make binarization for residual tumor in TCGA, 0= No residual tumor, 1 = residual tumor
residual_tumor_info <- read.csv(file=residual.file, header=T)
residual_tumor_info$TUMORRESIDUALDISEASE <- as.character(residual_tumor_info$TUMORRESIDUALDISEASE)
residual_tumor_info$TUMORRESIDUALDISEASE[which(residual_tumor_info$TUMORRESIDUALDISEASE == "No Macroscopic disease")] <- 0
residual_tumor_info$TUMORRESIDUALDISEASE[which(residual_tumor_info$TUMORRESIDUALDISEASE != "0" & residual_tumor_info$TUMORRESIDUALDISEASE !="")] <- 1
residual_tumor_info$TUMORRESIDUALDISEASE[which(residual_tumor_info$TUMORRESIDUALDISEASE == "")] <- NA
residual_tumor_info$TUMORRESIDUALDISEASE <- as.numeric(residual_tumor_info$TUMORRESIDUALDISEASE)

#Adding the BRCAness status column in OVA-TCGA patients:
TCGA.scars.clinical$BRCAness <- c(1,0)[match(grepl("BRCA",TCGA.scars.clinical$Germline_mutations)
                                             | grepl("BRCA",TCGA.scars.clinical$Somatic_mutations)
                                             | grepl("BRCA",TCGA.scars.clinical$CNV_status), c('TRUE', 'FALSE'))]

#Merging with residual tumor information for TCGA samples
row.names(residual_tumor_info) <- residual_tumor_info$BCRPATIENTBARCODE
TCGA_OVAscars_info_residual <- merge(TCGA.scars.clinical, residual_tumor_info, by="row.names")

#Select only patients that received platinum based chemotherapy in TCGA
platinum.patients <- scan(file=platin.takers.list.file, what="char")
TCGA_OVAscars_info_residual <- TCGA_OVAscars_info_residual[which(TCGA_OVAscars_info_residual$Row.names %in% platinum.patients),]
print(paste0("Number of patients selected from OVA-TCGA: ", nrow(TCGA_OVAscars_info_residual)))

#Validation set in TCGA correspond to those samples not used for the HRD-AIs optimization (Undefined)
TCGA_valset <- TCGA_OVAscars_info_residual[TCGA_OVAscars_info_residual$HRDstatus == "Undefined",]

#######Reading PCAWG info
PCAWGscar <- read.table(file=paste0(input.folder, PCAWG.scars.file), header=T, sep=",")
PCAW.PFS <- read.table(file = PCAWG.clin.file, header=T, sep="\t")
PCAWG.mutations <- read.table(mutations.driver.file, header=T, sep="\t")
PCAWG.TP53mutants <- PCAWG.mutations[which(PCAWG.mutations$gene == "TP53"), "sample_id"]

#Adjust data for survival status in PCAWG
PCAWGscar[which(PCAWGscar$donor_vital_status == "alive"),"donor_vital_status"] <- 0
PCAWGscar[which(PCAWGscar$donor_vital_status == "deceased"),"donor_vital_status"] <- 1
PCAWGscar$donor_vital_status <- as.numeric(PCAWGscar$donor_vital_status)

#Columns of interest from clinical data
PCAW.PFS <- PCAW.PFS[,c(1,7,8,11,12)]

#Adjusting time for progression and progression status in PCAWG
PCAWGscar_PFS  <- merge(PCAWGscar, PCAW.PFS, by="icgc_donor_id", all.x = TRUE)
PCAWGscar_PFS$donor_relapse_interval <- PCAWGscar_PFS$donor_relapse_interval/30.4 #For months
PCAWGscar_PFS$Progression <- rep(0, nrow(PCAWGscar_PFS))
PCAWGscar_PFS$Progression[PCAWGscar_PFS$disease_status_last_followup == "relapse" | PCAWGscar_PFS$disease_status_last_followup == "progression"] <- 1
#Selecting only high grade samples
PCAWGscar_HG_PFS <- PCAWGscar_PFS[which(PCAWGscar_PFS$tumour_grade >= 3),]
PCAWG.g2.samples <- PCAWGscar_PFS[which(PCAWGscar_PFS$tumour_grade == 2),2]
PCAWG.g2.samples.TP53mut <- PCAWG.g2.samples[PCAWG.g2.samples %in% PCAWG.TP53mutants]
PCAWG.g2.samples.TP53mut.rows <- PCAWGscar_PFS[PCAWGscar_PFS[,2] %in% PCAWG.g2.samples.TP53mut,]
PCAWGscar_HG_PFS <- rbind(PCAWGscar_HG_PFS, PCAWG.g2.samples.TP53mut.rows)
print(paste0("Number of patients selected from PCAWG: ", nrow(PCAWGscar_HG_PFS) ))

######Hercules info loading and editing
Hercules <- read.table(file=paste0(input.folder,HERCULES.scars.file), sep=",", header=T)
names(Hercules)[which(names(Hercules) == "Age.at.Diagnosis")] <- "age"
#Hercules$Residual_tumor[which(is.na(Hercules$Residual_tumor))] <- 1
print(paste0("Number of patients selected from HERCULES: ", nrow(Hercules)))

##################################################################################################################################
##################################################################################################################################
##################################### Get the Hazard ratios for PanCanAtlas-TCGA #################################################
##################################################################################################################################
##################################################################################################################################

#Just a test
Cox_regresion_variables(TCGA_OVAscars_info_residual, variables =c("TUMORRESIDUALDISEASE"), formula ='Surv(OS.time, OS)~ ')

########For TCGA, do cox regressions using the getHRDs function
#First doing Cox regressions analysis using co-variables, age for OS and TUMORRESIDUALDISEASE for PFI
#HR stands for Hazard ratios
#The function getHRs will do the COX and return HR, p.values, confidence intervals and other info by each of the predefined variables
TCGA.HR.OS <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(OS.time, OS)~ age +', datalabel="PanCanAtlas", BRCAness.comp=TRUE)
TCGA.HR.PFI <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(PFI.time, PFI)~ TUMORRESIDUALDISEASE +', datalabel="PanCanAtlas", BRCAness.comp=TRUE)
TCGA_valset.HR.OS <- getHRs(TCGA_valset, formulaHR ='Surv(OS.time, OS)~ age +', datalabel="PanCanAtlas", BRCAness.comp=TRUE)
TCGA_valset.HR.PFI <- getHRs(TCGA_valset, formulaHR ='Surv(PFI.time, PFI)~ TUMORRESIDUALDISEASE +', datalabel="PanCanAtlas", BRCAness.comp=TRUE)

#Same calculations but no co-variables used
TCGA.HR.OS2 <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(OS.time, OS)~ ', datalabel="PanCanAtlas", BRCAness.comp=TRUE)
TCGA.HR.PFI2 <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(PFI.time, PFI)~ ', datalabel="PanCanAtlas", BRCAness.comp=TRUE)
TCGA_valset.HR.OS2 <- getHRs(TCGA_valset, formulaHR ='Surv(OS.time, OS)~ ', datalabel="PanCanAtlas", BRCAness.comp=TRUE)
TCGA_valset.HR.PFI2 <- getHRs(TCGA_valset, formulaHR ='Surv(PFI.time, PFI)~ ', datalabel="PanCanAtlas", BRCAness.comp=TRUE)

########For PCAWG, do cox regressions using the getHRDs function
#First doing Cox regressions analysis using co-variables, age for OS and TUMORRESIDUALDISEASE for PFI
PCAWG.HR.OS <- getHRs(PCAWGscar_HG_PFS, formulaHR ='Surv(donor_survival_time, donor_vital_status)~ age +', datalabel="PCAWG", BRCAness.comp=TRUE)
PCAWG.HR.OS2 <- getHRs(PCAWGscar_HG_PFS, formulaHR ='Surv(donor_survival_time, donor_vital_status)~ ', datalabel="PCAWG", BRCAness.comp=TRUE)
#Not used PFI, we are not sure about treatment in PCAWG samples

#Ignoring patient that received PARP inhibitors in double blinded study
Hercules <- Hercules[Hercules$Patient != "H094",]

########For Hercules, do cox regressions using the getHRDs function
Hercules_HR.OS <- getHRs(Hercules, formulaHR ='Surv(OS.time, Survival2)~ age +', datalabel="Hercules", BRCAness.comp=TRUE)
Hercules_HR.PFI <- getHRs(Hercules, formulaHR ='Surv(PFI.time, Progression2)~ Residual_tumor +', datalabel="Hercules", BRCAness.comp=TRUE)
Hercules_HR.OS2 <- getHRs(Hercules, formulaHR ='Surv(OS.time, Survival2)~ ', datalabel="Hercules", BRCAness.comp=TRUE)
Hercules_HR.PFI2 <- getHRs(Hercules, formulaHR ='Surv(PFI.time, Progression2)~ ', datalabel="Hercules", BRCAness.comp=TRUE)

################################## Merging results for HR when using co-variables ########################################################
###########################################################################################################
################## A function for plotting dotplots for the Cox regression tables ########################

plot.dotplots <- function(matrix.values, datasets = c("TCGA","TCGA-training", "HERCULES"), analysis="OS_covariables"){
    matrix.values[6:9,1] <- c("Telli2016-2", "Telli2016-54-2", "Takaya2020-2","ovaHRDscar-2")
    pvalues.columns <- c(5, 9 , 13)
    HR.columns <- c(4, 8 , 12)

    for (h in 1:length(datasets)){
      df2 <- data.frame(
        xval = as.factor(rep(1, nrow(matrix.values))), #random value
        row = factor(matrix.values[,1], levels=rev(c("BRCAmut/del", "Telli2016", "Telli2016-54","Takaya2020", "ovaHRDscar", "Telli2016-2", "Telli2016-54-2", "Takaya2020-2", "ovaHRDscar-2"))),
        pval = (-1 * log10(as.numeric(matrix.values[,pvalues.columns[h]]))),
        HR = as.numeric(matrix.values[,HR.columns[h]])
      )
      
      print("Saving dot plot")
      p <- ggplot(df2, aes(x=xval, y=row))
      p <- p + geom_point(aes(color=pval, size=HR)) + theme_classic()
      p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue4", "deepskyblue","red"),
                                     values=c(0,0.90,1), guide = guide_colourbar(direction = "horizontal"))
      p <- p + scale_radius(range = c(5,1), limits = c(0.25, 1.1), breaks = c(0.40, 0.50, 0.60, 0.7, 0.8))
      p <- p + theme(axis.text=element_text(size=rel(0.7)), strip.placement = "outside",strip.background = element_blank(),
                     axis.text.x=element_blank(), axis.text.y=element_text(size=rel(1)),
                     legend.text=element_text(size=rel(0.8)), legend.title=element_text(size=rel(0.8)),
                     axis.title=element_text(size=rel(0.8)))
      p <- p + ylab("") + xlab("")
      print(p)
      ggsave(p, filename = paste0(output.folder, "/DotPlots/",analysis, "_", datasets[h],"_dotspvalues.svg"), width = 11, height = 5.25, units = "cm")
    }
}


################################## Merging results for HR when no using co-variables ########################################################
#########For OS

#Next table is used to generate dot plots (colored dots)
#First columns correspond to the whole OVA-TCGA columns
#Middle columns to the OVA-TCGA excluding the validation set
#Last columns to the PCAWG cox results
datamatrix.OS2 <- cbind(TCGA.HR.OS[[2]], TCGA.HR.OS[[1]][,c(1,4)],
                       TCGA_valset.HR.OS[[2]][,c(2,3)], TCGA_valset.HR.OS[[1]][,c(1,4)],
                       PCAWG.HR.OS[[2]][,c(2,3)], PCAWG.HR.OS[[1]][,c(1,4)])
datamatrix.OS2 <- datamatrix.OS2[-1,]
datamatrix.OS2[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")

#Next table is used to generate table with HR info
datamatrix.OS <- cbind(TCGA.HR.OS[[2]], TCGA.HR.OS[[1]][,c(5,4)],
                       TCGA_valset.HR.OS[[2]][,c(2,3)], TCGA_valset.HR.OS[[1]][,c(5,4)],
                       PCAWG.HR.OS[[2]][,c(2,3)], PCAWG.HR.OS[[1]][,c(5,4)])
datamatrix.OS <- datamatrix.OS[-1,]
datamatrix.OS[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")
colnames(datamatrix.OS) <- c("Names", "N", "Prop" , "Hazard ratio", "Pval",
                                "N", "Prop", "Hazard ratio", "Pval",
                                "N", "Prop", "Hazard ratio", "Pval")

#Storing results table in PDF
#For Figure3h
tt2 <- ttheme_minimal(base_size = 8, padding = unit(c(1.3, 2.9), "mm"))
pdf(paste0(output.folder, "HazardR_scores_OS.pdf"), height=10, width=15)
grid.table(datamatrix.OS, theme=tt2, rows=NULL)
dev.off()

#Saving dots for Figure3h
plot.dotplots(datamatrix.OS2, datasets = c("TCGA","TCGA-training", "PCAWG"))


#########For PFI, when covariables
#Next table is used to generate dot plots (colored dots)
datamatrix.PFI2 <- cbind(TCGA.HR.PFI[[2]], TCGA.HR.PFI[[1]][,c(1,4)],
                        TCGA_valset.HR.PFI[[2]][,c(2,3)], TCGA_valset.HR.PFI[[1]][,c(1,4)],
                        Hercules_HR.PFI[[2]][,c(2,3)], Hercules_HR.PFI[[1]][,c(1,4)])
datamatrix.PFI2 <- datamatrix.PFI2[-1,]
datamatrix.PFI2[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")

#Next table is used to generate table with HR info
datamatrix.PFI <- cbind(TCGA.HR.PFI[[2]],  TCGA.HR.PFI[[1]][,c(5,4)],
                        TCGA_valset.HR.PFI[[2]][,c(2,3)], TCGA_valset.HR.PFI[[1]][,c(5,4)],
                        Hercules_HR.PFI[[2]][,c(2,3)], Hercules_HR.PFI[[1]][,c(5,4)])
colnames(datamatrix.PFI) <- c("", "N", "Prop" , "Hazard ratio", "Pval",
                              "N", "Prop", "Hazard ratio", "Pval",
                              "N", "Prop", "Hazard ratio", "Pval")
datamatrix.PFI <- datamatrix.PFI[-1,]
datamatrix.PFI[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")

#Storing results table in PDF
#For Figure3d
pdf(paste0(output.folder, "HazardR_scores_PFI.pdf"), height=8, width=15)
grid.table(datamatrix.PFI, theme=tt2, rows=NULL)
dev.off()

##Dotplots for Figure3d
plot.dotplots(datamatrix.PFI2, datasets = c("TCGA","TCGA-training", "HERCULES"), analysis="PFI_covariables")


#########For OS when no covariables
datamatrix.OS2 <- cbind(TCGA.HR.OS2[[2]], TCGA.HR.OS2[[1]][,c(1,4)],
                        TCGA_valset.HR.OS2[[2]][,c(2,3)], TCGA_valset.HR.OS2[[1]][,c(1,4)],
                        PCAWG.HR.OS2[[2]][,c(2,3)], PCAWG.HR.OS2[[1]][,c(1,4)])
datamatrix.OS2 <- datamatrix.OS2[-1,]
datamatrix.OS2[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")


#Next table is used to generate table with HR info
datamatrix.OS <- cbind(TCGA.HR.OS2[[2]], TCGA.HR.OS2[[1]][,c(5,4)],
                        TCGA_valset.HR.OS2[[2]][,c(2,3)], TCGA_valset.HR.OS2[[1]][,c(5,4)],
                        PCAWG.HR.OS2[[2]][,c(2,3)], PCAWG.HR.OS2[[1]][,c(5,4)])
datamatrix.OS <- datamatrix.OS[-1,]
datamatrix.OS[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")
colnames(datamatrix.OS) <- c("Names", "N", "Prop" , "Hazard ratio", "Pval",
                              "N", "Prop", "Hazard ratio", "Pval",
                              "N", "Prop", "Hazard ratio", "Pval")


#Storing results table in PDF
#For Supplementary Figure3h
tt2 <- ttheme_minimal(base_size = 8, padding = unit(c(1.8, 2.9), "mm"))
pdf(paste0(output.folder, "HazardR_scores_OS2.pdf"), height=10, width=15)
grid.table(datamatrix.OS, theme=tt2, rows=NULL)
dev.off()

#Saving dots for figure Supplementary Figure3h
plot.dotplots(datamatrix.OS2, datasets = c("TCGA","TCGA-training", "PCAWG"), analysis="OS_Non-covariables")


#########For PFI, no covariables
#Next table is used to generate dot plots (colored dots)
datamatrix.PFI2 <- cbind(TCGA.HR.PFI2[[2]],  TCGA.HR.PFI2[[1]][,c(1,4)],
                         TCGA_valset.HR.PFI2[[2]][,c(2,3)], TCGA_valset.HR.PFI2[[1]][,c(1,4)],
                         Hercules_HR.PFI2[[2]][,c(2,3)], Hercules_HR.PFI2[[1]][,c(1,4)])
datamatrix.PFI2 <- datamatrix.PFI2[-1,]
datamatrix.PFI2[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")


#Next table is used to generate table with HR info
datamatrix.PFI <- cbind(TCGA.HR.PFI2[[2]],  TCGA.HR.PFI2[[1]][,c(5,4)],
                        TCGA_valset.HR.PFI2[[2]][,c(2,3)], TCGA_valset.HR.PFI2[[1]][,c(5,4)],
                        Hercules_HR.PFI2[[2]][,c(2,3)], Hercules_HR.PFI2[[1]][,c(5,4)])
colnames(datamatrix.PFI) <- c("", "N", "Prop" , "Hazard ratio", "Pval",
                              "N", "Prop", "Hazard ratio", "Pval",
                              "N", "Prop", "Hazard ratio", "Pval")
datamatrix.PFI <- datamatrix.PFI[-1,]
datamatrix.PFI[,1][c(1,3,6:9)] <- c("BRCAmut/del", "Telli2016-54", "Telli2016","Telli2016-54", "Takaya2020", "ovaHRDscar")


#Storing results table in PDF
#For Supplementary Figure3d
tt2 <- ttheme_minimal(base_size = 8, padding = unit(c(1.8, 2.9), "mm"))
pdf(paste0(output.folder, "HazardR_scores_PFI2.pdf"), height=10, width=14)
grid.table(datamatrix.PFI, theme=tt2, rows=NULL)
dev.off()

#Saving dots for figure Supplementary Figure3d
plot.dotplots(datamatrix.PFI2, datasets = c("TCGA","TCGA-training", "HERCULES"), analysis="PFI_Non-covariables")


############################################################################################################################
################################################ Kaplan meier Survival plots ###############################################
############################################################################################################################

#This will return the additional columns in to the input data.frame: BRCAness, ovaHRDscar, tacayaHRDhigh, etc...
TCGA_stratified <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(OS.time, OS)~ ', datalabel="PanCanAtlas", returndata = TRUE)
PCAWG_stratified <- getHRs(PCAWGscar_HG_PFS, formulaHR ='Surv(donor_survival_time, donor_vital_status)~ ', datalabel="PCAWG", returndata = TRUE)
Hercules_stratified <- getHRs(Hercules, formulaHR ='Surv(OS.time, Survival2)~ ', datalabel="Hercules", returndata = TRUE)


### TCGA PFI survival plots and tables###
surv_objectPFI <- Surv(time = TCGA_stratified$PFI.time, event = TCGA_stratified$PFI)
plot <- make.merge.survplots(surv_objectPFI, TCGA_stratified, variables = c("BRCAness","Telli2016", "ovaHRDscar"), max.time = 75,
                             break_time=15, xlabplot="PFI time (months)", ylabplot="PFI probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)
ggsave(plot, filename=paste0(output.folder, "KaplanMeier_PFI-TCGA.svg"), width = 25, height = 14, units = "cm")

### TCGA OS survival plots and tables###
surv_objectOS <- Surv(time = TCGA_stratified$OS.time, event = TCGA_stratified$OS)
plot <- make.merge.survplots(surv_objectOS, TCGA_stratified, variables = c("BRCAness","Telli2016", "ovaHRDscar"), max.time = 125,
                             break_time=25, xlabplot="OS time (months)", ylabplot="OS probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)
ggsave(plot, filename=paste0(output.folder, "KaplanMeier_OS-TCGA.svg"), width = 25, height = 14, units = "cm")


### PCAWG OS plots ####
surv_objectOS <- Surv(time = PCAWG_stratified$donor_survival_time, event = PCAWG_stratified$donor_vital_status)
surv_objectPFI <- Surv(time = PCAWG_stratified$donor_relapse_interval, event = PCAWG_stratified$Progression)
plot <- make.merge.survplots(surv_objectOS, PCAWG_stratified, variables = c("BRCAness","Telli2016", "ovaHRDscar"), max.time = 60,
                             break_time=15, xlabplot="OS time (months)", ylabplot="OS probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)
ggsave(plot, filename=paste0(output.folder, "KaplanMeier_OS-PCAWG.svg"), width = 25, height = 14, units = "cm")

### Hercules PFI plots ###
#Ignoring patient that received PARP inhibitors in double blinded clinical trial
Hercules_stratified <- Hercules_stratified[!Hercules_stratified$Patient %in% "H094",]

surv_objectPFI <- Surv(time = Hercules_stratified$PFI.time, event = Hercules_stratified$Progression2)
surv_objectOS <- Surv(time = Hercules_stratified$OS.time, event = Hercules_stratified$Survival2)
plot <- make.merge.survplots(surv_objectPFI, Hercules_stratified, variables = c("BRCAness","Telli2016", "ovaHRDscar"), max.time = 40,
                             break_time=10, xlabplot="PFI time (months)", ylabplot="PFI probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)
ggsave(plot, filename=paste0(output.folder, "KaplanMeier_PFI-HERCULES.svg"), width = 25, height = 14, units = "cm")

plot <- make.merge.survplots(surv_objectOS, Hercules_stratified, variables = c("BRCAness","Telli2016", "ovaHRDscar"),
                             break_time=10, xlabplot="OS time (months)", ylabplot="OS probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)

################################################ Survival plots for TCGA removing training set ############################################################


#This will return the additional columns: BRCAness, ovaHRDscar, tacayaHRDhigh, etc...
TCGA_stratified <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(OS.time, OS)~ ', datalabel="PanCanAtlas", returndata = TRUE)

#Selecting only samples with Not_clear HRD status (the validation set)
TCGA_valset <- TCGA_stratified[TCGA_stratified$HRDstatus == "Not_clear",]

###Do the PFI  analysis
surv_objectPFI  <- Surv(time = TCGA_valset$PFI.time, event = TCGA_valset$PFI)
plot <- make.merge.survplots(surv_objectPFI, TCGA_valset, variables = c("Telli2016", "ovaHRDscar"),
                             break_time=25, xlabplot="PFI time (months)", ylabplot="PFI probability", palette = c("#34A0D3" ,"#FA5A41"))
print(plot)

##################################################################################################
##################################################################################################
###################################Comparing median survival times ###############################
##################################################################################################
##################################################################################################


###Series of analysis to generate Supp Figure 3d and Supp Figure 3e###
###Will bootstrap dataset, do Mann-Whitney u-test to compare the difference of median survivals between HRD/HRP#####

#This will return the additional stratification columns in to the input data.frame: BRCAness, ovaHRDscar, tacayaHRDhigh, etc...
TCGA_stratified <- getHRs(TCGA_OVAscars_info_residual, formulaHR ='Surv(OS.time, OS)~ ', datalabel="PanCanAtlas", returndata = TRUE)
TCGA_valset_stratified <- getHRs(TCGA_valset, formulaHR ='Surv(OS.time, OS)~ ', datalabel="PanCanAtlas", returndata = TRUE)
PCAWG_stratified <- getHRs(PCAWGscar_HG_PFS, formulaHR ='Surv(donor_survival_time, donor_vital_status)~ ', datalabel="PCAWG", returndata = TRUE)
Hercules_stratified <- getHRs(Hercules, formulaHR ='Surv(OS.time, Survival2)~ ', datalabel="Hercules", returndata = TRUE)

#Concatenating stratified datasets
datasets.boots <- list(TCGA_stratified, TCGA_valset_stratified, PCAWG_stratified, Hercules_stratified)

#Changing colnames for survival info in PCAWG and HERCULEs dataset
colnames(datasets.boots[[3]])[colnames(datasets.boots[[3]]) == "donor_survival_time"] <- "OS.time"
colnames(datasets.boots[[3]])[colnames(datasets.boots[[3]]) == "donor_vital_status"] <- "OS"
colnames(datasets.boots[[3]])[colnames(datasets.boots[[3]]) == "Progression"] <- "PFI"
colnames(datasets.boots[[3]])[colnames(datasets.boots[[3]]) == "donor_relapse_interval"] <- "PFI.time"
colnames(datasets.boots[[4]])[colnames(datasets.boots[[4]]) == "Progression2"] <- "PFI"
colnames(datasets.boots[[4]])[colnames(datasets.boots[[4]]) == "Survival2"] <- "OS"

datasets.names <- c("OVA-TCGA","OVA-TCGA**","PCAWG","HERCULES")
names.surv.objs <- c("PFI", "OS")
res.values <- NULL
res.values2 <- NULL
res.values3 <- NULL

#This for loop will go by dataset by dataset
for (d in 1:length(datasets.boots)){
  data.full <- datasets.boots[[d]]
  boots <- 1000 #Number of bootstraps
  boots.median.pfi <- NULL
  boots.median.os <- NULL
  #Loop for boostraping, in each replicate calculate the median survival times
  for (i in 1:boots){
    vals <- sort(sample(nrow(data.full), replace = TRUE)) #selecting bootstrapped samples
    set <- data.full[vals,] #selecting bootstrapped samples
    #Calculating median survival times of independient variables of merged variables
    surv_objectPFI <- Surv(time = set$PFI.time, event = set$PFI) #Surv object with bootstrapped samples
    surv_objectOS <- Surv(time = set$OS.time, event = set$OS) #Surv object with bootstrapped samples
    l <- median.surv.times(surv_objectPFI, set, variables = c("BRCAness","Telli2016", "Telli2016.54", "ovaHRDscar")) #Median PFI times by variable
    m <- median.surv.times(surv_objectOS, set, variables = c("BRCAness","Telli2016", "Telli2016.54", "ovaHRDscar")) #Median OS times by variable
    boots.median.pfi <- rbind(boots.median.pfi, l) #Concatenating results in a dataframe
    boots.median.os <- rbind(boots.median.os, m)  #Concatenating results in a dataframe
  }

  surv_objectPFI <- Surv(time = data.full$PFI.time, event = data.full$PFI) #Surv object with all dataset
  surv_objectOS <- Surv(time = data.full$OS.time, event = data.full$OS) #Surv object with all dataset
  list.surv.obj <- list(surv_objectPFI, surv_objectOS)
  list.surv.boots <- list(boots.median.pfi, boots.median.os)

  for (j in 1:2){
    surv.boots <- list.surv.boots[[j]]
    medians.surv <- median.surv.times(list.surv.obj[[j]], data.full, variables = c("BRCAness","Telli2016", "ovaHRDscar","Telli2016.54")) #Median surv times for whole dataset
    medians.surv.prop <- medians.surv[medians.surv$marker.class == 1,] / medians.surv[medians.surv$marker.class == 0,]
    surv.boots.pos <- surv.boots[surv.boots$marker.class == 1,] #For HRD positives, for each bootstrap replicate
    surv.boots.neg <- surv.boots[surv.boots$marker.class == 0,] #For HRD negatives, for each bootstrap replicate
    surv.boots.fold <- surv.boots.pos/surv.boots.neg #Getting the ratio of HRD/HRP median survival times for plot
    surv.boots.fold <- surv.boots.fold[!is.na(surv.boots.fold[,2]),] #Ignoring NAs, column for Telli2016
    surv.boots.fold <- surv.boots.fold[!is.na(surv.boots.fold[,3]),] #Ignoring NAs, column for Telli2016-54
    surv.boots.fold <- surv.boots.fold[!is.na(surv.boots.fold[,4]),] #Ignoring NAs, column for ovaHRDscar
        CI.1.fold = paste0("CI:",round(quantile(surv.boots.fold[,2], 0.025, na.rm = TRUE),1),
                  ",", round(quantile(surv.boots.fold[,2], 0.975, na.rm = TRUE),1))
    CI.2.fold = paste0("CI:",round(quantile(surv.boots.fold[,3], 0.025, na.rm = TRUE),1),
                  ",", round(quantile(surv.boots.fold[,3], 0.975, na.rm = TRUE),1))
    
    
    #There can be some NAs, because for Telli2016 or ovaHRDscar, during bootstrap was not selected positive samples
    
     test.fold1 <- wilcox.test(surv.boots.fold[,3], surv.boots.fold[,2], alternative = "greater")
     p.value.fold1 <- signif(test.fold1$p.value, digits = 4)
     test.fold2 <- wilcox.test(surv.boots.fold[,4], surv.boots.fold[,2], alternative = "greater")
     p.value.fold2 <- signif(test.fold2$p.value, digits = 4)
     test.fold3 <- wilcox.test(surv.boots.fold[,4], surv.boots.fold[,3], alternative = "greater")
     p.value.fold3 <- signif(test.fold3$p.value, digits = 4)
    #Concatenating results, showing foldchanges using whole dataset, and confidence intervals from bootstraping and p-value from bootstraping
    df.iter <- data.frame(dataset = datasets.names[d],
                          comparison = names.surv.objs[j],
                          Diff.Telli2016=paste0(round(medians.surv[2,2] / medians.surv[1,2],1), "(", CI.1.fold, ")"),
                          Diff.Telli=paste0(round(medians.surv[2,3] / medians.surv[1,3],1), "(", CI.2.fold, ")"),
                          Diff.ovaHRDscar=paste0(round(medians.surv[2,4] / medians.surv[1,4],1), "(", CI.2.fold, ")"),
                          pval1= p.value.fold1, pval2= p.value.fold2, pval3= p.value.fold3)
    #Concatenating results, showing foldchanges using whole dataset, and confidence intervals from bootstraping and p-value from bootstraping
    df.iter2 <- data.frame(dataset = rep(datasets.names[d],2),
                           comparison = rep(names.surv.objs[j],2),
                           Diff=c(round(medians.surv[2,2] / medians.surv[1,2],1), round(medians.surv[2,3] / medians.surv[1,3],1)),
                           CI.up = c(round(quantile(surv.boots.fold[,2], 0.975),1), round(quantile(surv.boots.fold[,3], 0.975),1)),
                           CI.down = c(round(quantile(surv.boots.fold[,2], 0.025),1), round(quantile(surv.boots.fold[,3], 0.025),1)),
                           Algorithm = c("Telli2016","ovaHRDscar"))
    #Concatenating results, only foldchanges in bootstraping replicates
    df.iter3 <- data.frame(dataset = datasets.names[d],
                           comparison = names.surv.objs[j],
                           Differences=c(surv.boots.fold[,2],surv.boots.fold[,3],surv.boots.fold[,4]),
                           Algorithm = c(rep("Telli2016",length(surv.boots.fold[,2])), rep("Telli-54",length(surv.boots.fold[,3])), rep("ovaHRDscar",length(surv.boots.fold[,4]))))
    res.values <- rbind(res.values, df.iter)
    res.values2 <- rbind(res.values2, df.iter2)
    res.values3 <- rbind(res.values3, df.iter3)
  }
}

tt2 <- ttheme_minimal(base_size = 8, padding = unit(c(1.8, 2.9), "mm"))
pdf(paste0(output.folder, "Fold_medianSurv_table.pdf"), height=10, width=15)
grid.table(res.values, rows=NULL, theme=tt2)
dev.off()

res.values3$dataset <- factor(res.values3$dataset, levels=datasets.names)
res.values3$Algorithm <- factor(res.values3$Algorithm, levels=c("Telli2016", "Telli-54","ovaHRDscar"))

#Theme to be used in the Figures
mytheme <-  theme(strip.placement = "outside",strip.background = element_blank(),
                   axis.text.x=element_text(size=rel(1.2)), axis.text.y=element_text(size=rel(1.2)),
                   legend.text=element_text(size=rel(1.2)),
                   strip.text.x = element_text(size=rel(1.2)),
                   legend.title=element_text(size=rel(1.2)), axis.title=element_text(size=rel(1.2)),
                   panel.spacing = unit(1, "lines"),
                   legend.position = "right",
                   panel.border = element_rect(linetype = "solid", fill = NA),
                   panel.background = element_rect(fill = "white"),
                   panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                   panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank())


#Generating Supp Figure 3d
res.values3.os  <- res.values3[res.values3$dataset != "HERCULES",]
p <- ggplot(subset(res.values3.os, comparison == "OS"), aes(x=dataset, y=Differences, fill=Algorithm)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Cohort") + ylab("Fold change median OS: HRD/HRP") +
  scale_fill_manual(values=c('grey100','grey60','grey30')) +
  mytheme + ylim(0,2.6)
print(p)
ggsave(p, filename=paste0(output.folder, "Fold_medianOS_boxplot_Telli54.svg"), width = 12.5, height = 10, units = "cm")

#Generating Supp Figure 3e
res.values3.pfi  <- res.values3[res.values3$dataset != "PCAWG",]
p <- ggplot(subset(res.values3.pfi, comparison == "PFI"), aes(x=dataset, y=Differences, fill=Algorithm)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Cohort") + ylab("Fold change median PFI: HRD/HRP") +
  scale_fill_manual(values=c('grey100','grey60','grey30')) +
  mytheme + ylim(0,4)
print(p)
ggsave(p, filename=paste0(output.folder, "Fold_medianPFI_boxplot_Telli54.svg"), width = 12.5, height = 10, units = "cm")

###################################################################################################################################
###################################################################################################################################
################################################  CHORD comparison #####################################################
###################################################################################################################################
###################################################################################################################################

CHORD <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/PCAW/CHORD_PCAWG_ovary.csv", sep =",", header=TRUE)
CHORD <- CHORD[,c(2,8)]
CHORD$CHORD.status <- c(1,0)[match(CHORD[2] == "HR_deficient", c('TRUE', 'FALSE'))]
CHORD$CHORD.status[CHORD[2] == "cannot_be_determined"] <- NA


PCAWG_stratified <- getHRs(PCAWGscar_HG_PFS, formulaHR ='Surv(donor_survival_time, donor_vital_status)~ ', datalabel="PCAWG", returndata = TRUE)
PCAWG_stratified.chord <- merge(PCAWG_stratified, CHORD, by.x="sample_id", by.y="sample")
PCAWG_stratified.chord <- PCAWG_stratified.chord[!is.na(PCAWG_stratified.chord$CHORD.status),]

surv_objectOS <- Surv(time = PCAWG_stratified.chord$donor_survival_time, event = PCAWG_stratified.chord$donor_vital_status)

plot <- make.merge.survplots(surv_objectOS, PCAWG_stratified.chord, variables = c("CHORD.status", "ovaHRDscar"), max.time = 60,
                             break_time=15, xlabplot="OS time (months)", ylabplot="OS probability", palette = c("#34A0D3", "#FA5A41"))
print(plot)
ggsave(plot, filename = paste0(output.folder, "KaplanMeier_OS_CHORD-status_PCAWG.svg"), height = 12, width = 14, units = "cm")



t1 <- Cox.regresion.variables2(PCAWG_stratified.chord, variables =c("CHORD.status", "ovaHRDscar"), formula ='Surv(donor_survival_time, donor_vital_status)~ ')


coxvalues <- t1
coxvalues$p.value <- (-1 * log10(coxvalues$p.value))
coxvalues$attribute <- factor(rownames(coxvalues), levels = rev(rownames(coxvalues)))
coxvalues$xval <- rep(1,nrow(coxvalues))

p <- ggplot(coxvalues, aes(x=xval, y=attribute))
p <- p + geom_point(aes(color=p.value, size=HR)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue4", "deepskyblue","red"),
                               values=c(0,0.90,1), guide = guide_colourbar(direction = "horizontal"))
p <- p + scale_radius(range = c(7,1), limits = c(0.25, 1.1), breaks = c(0.40, 0.50, 0.60, 0.7, 0.8))
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_blank(), axis.text.y=element_text(size=rel(1)),
               legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)),
               axis.title=element_text(size=rel(1)))
p <- p + ylab("") + xlab("")
print(p)
ggsave(p, filename = paste0(output.folder, "Dotpvalue_CHROD.svg"), width = 11, height = 5.5, units = "cm")


##################################################################################################
##################################################################################################
####################Comparing median survival adjusting different cutoff values of ovaHRDscar ###########
##################################################################################################
##################################################################################################

###Series of analysis to generate Supp Figure 3d and Supp Figure 3e###
###Will bootstrap dataset, do t-test to compare the difference of median survivals between HRD/HRP#####

TCGA_OVAscars_info_residual.notrain <- TCGA_OVAscars_info_residual[which(TCGA_OVAscars_info_residual$HRDstatus == "Undefined"),]

#Adding datasets to list
datasets <- list(TCGA_OVAscars_info_residual,
                 TCGA_OVAscars_info_residual.notrain,
                 PCAWGscar_HG_PFS,
                 Hercules)

datasets.names <- c("OVA-TCGA","OVA-TCGA**", "PCAWG","HERCULES")
names.surv.objs <- c("PFI", "OS")

#Changing colnames for survival info in PCAWG and HERCULEs dataset
colnames(datasets[[3]])[colnames(datasets[[3]]) == "donor_survival_time"] <- "OS.time"
colnames(datasets[[3]])[colnames(datasets[[3]]) == "donor_vital_status"] <- "OS"
colnames(datasets[[3]])[colnames(datasets[[3]]) == "Progression"] <- "PFI"
colnames(datasets[[3]])[colnames(datasets[[3]]) == "donor_relapse_interval"] <- "PFI.time"
colnames(datasets[[4]])[colnames(datasets[[4]]) == "Progression2"] <- "PFI"
colnames(datasets[[4]])[colnames(datasets[[4]]) == "Survival2"] <- "OS"

boots <- 1000 #Number of bootstraps
res.values <- NULL #Dataframe storing main results
cutoff.values <- 45:60

for (cutoffval in cutoff.values){
  for (d in 1:length(datasets)){
    dat.stratified <- datasets[[d]]
    dat.stratified$ovaHRDscar <- ifelse(dat.stratified$HRDsum>=cutoffval, 1, 0)
    
    #This for loop will go by dataset by dataset
    boots.median.pfi <- NULL
    boots.median.os <- NULL
    #Loop for boostraping, in each replicate calculate the median survival times
    for (i in 1:boots){
      vals <- sort(sample(nrow(dat.stratified), replace = TRUE)) #selecting bootstrapped samples
      set <- dat.stratified[vals,] #selecting bootstrapped samples
      #Calculating median PFI and OS values
      surv_objectPFI <- Surv(time = set$PFI.time, event = set$PFI) #Surv object with bootstrapped samples
      surv_objectOS <- Surv(time = set$OS.time, event = set$OS) #Surv object with bootstrapped samples
      l <- median.surv.times(surv_objectPFI, set, variables = c("ovaHRDscar")) #Median PFI times by variable
      m <- median.surv.times(surv_objectOS, set, variables = c("ovaHRDscar")) #Median OS times by variable
      boots.median.pfi <- rbind(boots.median.pfi, l) #Concatenating results in a dataframe
      boots.median.os <- rbind(boots.median.os, m)  #Concatenating results in a dataframe
    }
      
    list.surv.boots <- list(boots.median.pfi, boots.median.os)
    
      for (j in 1:2){
        surv.boots <- list.surv.boots[[j]]
        surv.boots.pos <- surv.boots[surv.boots$marker.class == 1,1] #For HRD positives, for each bootstrap replicate
        surv.boots.neg <- surv.boots[surv.boots$marker.class == 0,1] #For HRD negatives, for each bootstrap replicate
        surv.boots.fold <- surv.boots.pos/surv.boots.neg #Getting the ratio of HRD/HRP median survival times for plot
        surv.boots.fold <- surv.boots.fold[!is.na(surv.boots.fold)] #Ignoring NAs
        #Concatenating results, only foldchanges in bootstraping replicates
        df.iter <- data.frame(dataset = datasets.names[d],
                               comparison = names.surv.objs[j],
                               Differences=c(surv.boots.fold),
                               Algorithm = rep("ovaHRDscar",length(surv.boots.fold)),
                               cutoff = cutoffval)
        res.values <- rbind(res.values, df.iter)
      }
  
  }
}

#Compare with Manwhitney U test each cutoff value vs the value of 54
df.result <- NULL
for (data.name in datasets.names){
    res.values.dat  <- res.values[res.values$dataset %in% data.name,]
    for (l in cutoff.values){
      res.values.dat.aux <- res.values.dat[res.values.dat$cutoff %in% c(l,54),]
      for(survclass in c("OS","PFI")){
        res.values.dat.class <- res.values.dat.aux[res.values.dat.aux$comparison %in% survclass,]
        surv.fold1 = res.values.dat.class[res.values.dat.class$cutoff == l,"Differences"]
        surv.fold2 = res.values.dat.class[res.values.dat.class$cutoff == 54,"Differences"]
        test.fold <- wilcox.test(surv.fold1, surv.fold2, alternative = "less")
        p.value.fold <- signif(test.fold$p.value, digits = 4)
        df.iter <- data.frame(dataset = data.name,
                              comparison = survclass,
                              cutoffsvals = paste0(l, "-54"),
                              pval= p.value.fold)
        df.result <- rbind(df.result, df.iter)
      }
    }
}

res.values$cutoff <- factor(res.values$cutoff, levels=unique(res.values$cutoff))
res.values$dataset <- factor(res.values$dataset, levels=c("OVA-TCGA","OVA-TCGA**", "PCAWG","HERCULES"))

res.values.pfi <- res.values[res.values$dataset != "PCAWG",]
p <- ggplot(subset(res.values.pfi, comparison == "PFI"), aes(x=dataset, y=Differences, fill=cutoff)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Cohort") + ylab("Fold change median PFI: HRD/HRP") + mytheme + ylim(0,3.85)
print(p)
ggsave(p, filename=paste0(output.folder, "Fold_medianPFI-ovaHRDscar_cutoffs_boxplot.svg"), width = 13, height = 11, units = "cm")

res.values.os <- res.values[res.values$dataset != "HERCULES",]
p <- ggplot(subset(res.values.os, comparison == "OS"), aes(x=dataset, y=Differences, fill=cutoff)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Cohort") + ylab("Fold change median OS: HRD/HRP") + mytheme + ylim(0,2.9)
print(p)
ggsave(p, filename=paste0(output.folder, "Fold_medianOS-ovaHRDscar_cutoffs_boxplot.svg"), width = 11, height = 11, units = "cm")
