##################################################################################################################################
########################################### ovaHRDscar using SNParrays #######################################
##################################################################################################################################
library("DescTools")

setwd("C:/Users/fernpere/ovaHRDscar_manuscript_scripts/")
#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")

#Reading input files
segments <- read.table(file="Z:/Documents/HRD_Hercules/input/TERVA_ASCAT_segments.txt", header = TRUE)
Hercules_clinical_info <- read.csv(file="Z:/Documents/HRD_Hercules/input/220129_Clinical_binary.csv", header = T, sep =",", row.names=1)
samplesheet <- read.table(file="Z:/Documents/HRD_Hercules/input/TERVA_ASCAT_samplesheet.csv", header = TRUE, sep =",", row.names = 1)
functionalHR <- read.table(file="Z:/Documents/HRD_Hercules/input/TERVA_functional_values2.csv", header = TRUE, sep =",")
genomic.scars <- read.table(file = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC/Hercules_newscars_clinic.csv", header=TRUE, sep=",")
genomic.alterations <- read.table(file="Z:/Documents/HRD_Hercules/input/selected_genes_summary_20201216.csv", sep=",", header=TRUE)
genomic.scars <- genomic.scars[genomic.scars$Max.purity >= 0.2 & genomic.scars$WGS.platform != "BGI",]
telliscores <- read.table(file="Z:/Documents/HRD_Hercules/input/TERVA_TelliHRDscore.csv", sep=",",header=TRUE)
telliscores <- telliscores[,c(1,5)]
names(telliscores) <- c("Sample","Telli2016")

output.folder <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/Survival_54/"

####Formatting Hercules information
#Remove NA values
segments <- segments[which(!is.na(segments[,5])),]

###Formatting input data
ploy = rep(2, nrow(segments))
terva.segs <- data.frame(SampleID = segments[,1],
                         Chromosome = paste("chr",segments[,2],sep=""),
                         Start_position = segments[,3],
                         End_position = segments[,4],
                         total_cn = (segments[,5] + segments[,6]),
                         A_cn  = segments[,5],
                         B_cn = segments[,6],
                         ploidy = ploy)

############## Calculating scars in input segments for Hercules ##############
terva.Prev.scars <- get_HRDs(terva.segs, chrominfo = chrominfo_grch38, LOH_windos = c(15,5000), LST_segSize = 10e6, LST_mindistance=3e6)
terva.Prev.scars <- as.data.frame(terva.Prev.scars)
Telli2016 <- terva.Prev.scars$HRDsum

terva.scars <- get_HRDs(terva.segs, chrominfo = chrominfo_grch38)
terva.scars <- as.data.frame(terva.scars)
terva.scars <- cbind(terva.scars, Telli2016)

terva.scars.names <- merge(terva.scars, samplesheet, by="row.names")
head(terva.scars.names)

#Ignore samples with noisy BAFplots
noisySamples <- c("OCP.12", "OCP.45.DUP.1", "OCP.45")
terva.scars.names <- terva.scars.names[!(terva.scars.names$Row.names %in% noisySamples),]

##########################Merge with functional scores #####################################
terva.scars.names[terva.scars.names$Row.names == "OCP.19", "Sample_Name"] <- "H289_pAdnR1_TR"
terva.scars.functional <- merge(terva.scars.names, functionalHR, by.x = "Sample_Name", by.y = "Sample.ID", all.x = TRUE)

#Ignore FFPE samples and Oncosys samples
terva.scars.functional <- terva.scars.functional[terva.scars.functional$Tissue.type == "FF",]
terva.scars.functional <- terva.scars.functional[!grepl("S", terva.scars.functional$Patient),]

############################################ Survival analysis #################################################
terva.scars.functional.clin <- merge(terva.scars.functional, Hercules_clinical_info, by.y = "Patient.card..Patient.cohort.code_Patient.Card", by.x = "Patient")
parp.samples <- c("H279","H243","H245") #These patient received PARP inhibitors as first line therapy in a double blinded clinical trial
terva.scars.functional.clin <- terva.scars.functional.clin[!(terva.scars.functional.clin$Patient %in% parp.samples),]

#This is to get the average values by patients
selected_samples <- NULL
tissue.of.interest <- c("Ova", "Adn")
for (patient in unique(terva.scars.functional.clin$Patient)){
  patient_samples <- terva.scars.functional.clin[terva.scars.functional.clin$Patient == patient,]
  patient_selectedsamples <- NULL
  if(any(patient_samples$Tissue %in% tissue.of.interest)){
    patient_selectedsamples <- patient_samples[patient_samples$Tissue %in% tissue.of.interest,]
  }else{
    if (any(patient_samples$Tissue == "Ome")){
      patient_selectedsamples <- patient_samples[patient_samples$Tissue == "Ome",]
    }
  }
  if (is.null(patient_selectedsamples)){
    print(paste0("Patient with no tissue of interest: ", patient))
    patient_selectedsamples <- patient_samples
  }
  selected_samples <- rbind(selected_samples, patient_selectedsamples)
}

write.table(selected_samples, file="D:/users/fperez/TERVA_selected_samples.csv", sep=",")

#Getting the average value per sample
terva.scars.functional.clin_patient <- avg_by_sample(selected_samples, columns_to_avg = c(7, 8))
terva.scars.functional.clin_patient$ovaHRDscar<- c("1","0")[match(terva.scars.functional.clin_patient$HRDsum >= 54, c(TRUE,FALSE))]
terva.scars.functional.clin_patient$Telli <- c("1","0")[match(terva.scars.functional.clin_patient$Telli2016 >= 42, c(TRUE,FALSE))]
terva.scars.functional.clin_patient$Tacaya <- c("1","0")[match(terva.scars.functional.clin_patient$Telli2016 >= 63, c(TRUE,FALSE))]

terva.scars.functional.clin_patient[terva.scars.functional.clin_patient$Patient %in% parp.samples,]

###Kaplan-Meier for Supp Figure3g
surv_objectPFI <- Surv(time = terva.scars.functional.clin_patient$PFI.time, event = terva.scars.functional.clin_patient$Progression2)
plot <- make.merge.survplots(surv_objectPFI, terva.scars.functional.clin_patient, variables = c("Telli" , "ovaHRDscar"), max.time=30,
                             break_time=5, xlabplot="PFS time (months)", ylabplot="PFS probability", palette = c("#34A0D3", "#FA5A41"))
print(plot)
ggsave(plot, filename = paste0(output.folder, "KaplanMeier_PFI_TERVA.svg"), height = 12, width = 14, units = "cm")


#For all samples
t1 <- Cox.regresion.variables2(terva.scars.functional.clin_patient, variables =c("Telli", "Tacaya", "ovaHRDscar"), formula ='Surv(PFI.time, Progression2)~ Residual_tumor +')

#Ignoring sample with BRCAmut
terva.scars.noBCRAmut <- terva.scars.functional.clin_patient[terva.scars.functional.clin_patient$Patient !=  "H197",]
t2 <- Cox.regresion.variables2(terva.scars.noBCRAmut, variables =c("Telli", "Tacaya", "ovaHRDscar"), formula ='Surv(PFI.time, Progression2)~ Residual_tumor + ')

#coxvalues <- rbind(t1, t2)
coxvalues <- t1
coxvalues$p.value <- (-1 * log10(coxvalues$p.value))
coxvalues$attribute <- factor(rownames(coxvalues), levels = rev(rownames(coxvalues)))
coxvalues$xval <- rep(1,nrow(coxvalues))

coxvalues

###Color dots for figure 
p <- ggplot(coxvalues, aes(x=xval, y=attribute))
p <- p + geom_point(aes(color=p.value, size=HR)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue4", "deepskyblue","red"),
                               values=c(0,0.90,1), guide = guide_colourbar(direction = "horizontal"))
p <- p + scale_radius(range = c(6,3), limits = c(0.08, 0.52))
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_blank(), axis.text.y=element_text(size=rel(1)),
               legend.text=element_text(size=rel(1)), legend.title=element_text(size=rel(1)),
               axis.title=element_text(size=rel(1)))
p <- p + ylab("") + xlab("")
print(p)
ggsave(p, filename = paste0(output.folder, "Dotpvalue_Terva.svg"), width = 11, height = 7, units = "cm")
