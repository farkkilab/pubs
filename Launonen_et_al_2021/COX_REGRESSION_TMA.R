##################################################################################################################
##################################### COX REGRESSION #############################################################
##################################################################################################################

#For Figure 2, 3 and Supplementary Figures 1, 2, 3 as well as Supplementary Tables 2 and 3

library(dplyr)
library(tidyr)
library(survival)
library(survminer)


ann <- read.csv("TMA_clinicaldata.csv")
ann <- ann[, c("Identifier","Category", "PFI_time")]

ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "short"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"
ann$Category <- as.character(ann$Category)
ann$Category[which(ann$`Category`== "HR")] <- "HRwt"
ann$Category[which(ann$`Category`== "BRCA1")] <- "BRCA1/2 mutated"
ann$Category[which(ann$`Category`== "BRCA2")] <- "BRCA1/2 mutated"


all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")

all_celltypes_24092020$GlobalCellType <- as.character(all_celltypes_24092020$GlobalCellType)

all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD11c+CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+CD11c+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD11c+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD11c+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD4"), "GlobalCellType"] <- "CD4+ Effector T-cells"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD8"), "GlobalCellType"] <- "CD8+ T-cells"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD163+Macrophages"), "GlobalCellType"] <- "CD163+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD11c+APC"), "GlobalCellType"] <- "CD11c+APC"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "FOXP3+CD4+Tregs"), "GlobalCellType"] <- "FOXP3+CD4+ T-regs"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "High-PDL1"), "GlobalCellType"] <- "Functional stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Non-proliferative_Stroma"), "GlobalCellType"] <- "Non-proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Low_eccentricity_medium_vimentin"), "GlobalCellType"] <- "Low eccentricity"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Proliferative_Stroma"), "GlobalCellType"] <- "Proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "High-proliferative_Stroma"), "GlobalCellType"] <- "High-proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Hyperfunctional epithelial"), "GlobalCellType"] <- "Functional epithelial"
all_celltypes_24092020 <- all_celltypes_24092020[-which(all_celltypes_24092020$GlobalCellType == "Negative"),]


all_clusters <- all_celltypes_24092020[, c(1, 46)]

metacluster_percentages <- all_clusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
metacluster_percentages <- as.data.frame(metacluster_percentages)

metacluster_percentages <- metacluster_percentages[, -c(3)]
metacluster_percentages$proportion <- metacluster_percentages$proportion*100


metacluster_percentages <- metacluster_percentages %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

metacluster_percentages <- as.data.frame(metacluster_percentages)
rownames(metacluster_percentages) <- metacluster_percentages[, 1]
metacluster_percentages <- metacluster_percentages[, -1]

metacluster_percentages[is.na(metacluster_percentages)] <- 0


ann <- ann[rownames(metacluster_percentages),]

metacluster_percentages$Patient <- rownames(metacluster_percentages)

median_expressions <- merge(metacluster_percentages, all_data_TMA, by.x="Patient", by.y="Identifier")


colnames(median_expressions) <- make.names(colnames(median_expressions))

median_expressions$Type <- as.character(median_expressions$Type)

median_expressions[which(median_expressions$Type == "BRCA1"), "Type"] <- "BRCA1/2 mutated"
median_expressions[which(median_expressions$Type == "BRCA2"), "Type"] <- "BRCA1/2 mutated"

median_expressions_abundance <- median_expressions


all_data_TMA_cat_abundance <- median_expressions_abundance
all_data_TMA_cat_BRCA_abundance <- median_expressions_abundance[which(median_expressions_abundance$Type == "BRCA1/2 mutated"),]

all_data_TMA_cat_BRCA_abundance <- median_expressions_abundance[which(median_expressions_abundance$Type == "BRCA1"| median_expressions_abundance$Type == "BRCA2"),]

all_data_TMA_cat_HR_abundance <-  median_expressions_abundance[which(median_expressions_abundance$Type == "HR"),]

for (i in seq_along(colnames(median_expressions_abundance))[c(6, 7, 11, 19)]){
  
  all_data_TMA_cat_abundance[, i] <- cut(all_data_TMA_cat_abundance[, i], breaks = c(-Inf, median(all_data_TMA_cat_abundance[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_BRCA_abundance[, i] <- cut(all_data_TMA_cat_BRCA_abundance[, i], breaks = c(-Inf, median(all_data_TMA_cat_BRCA_abundance[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_HR_abundance[, i] <- cut(all_data_TMA_cat_HR_abundance[, i], breaks = c(-Inf, median(all_data_TMA_cat_HR_abundance[, i]), Inf), labels=c("low", "high"))
  
}



profepi_surv <- all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Proliferating epithelial"),]
#find out median Ki67 expression in each patient
median_expressions <- profepi_surv %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                                      median_pSTAT1 = median(pSTAT1), median_P21 = median(P21), median_PD1=median(PD1), median_PDL1=median(PDL1), median_cCasp3=median(cCasp3))

median_expressions <- merge(median_expressions, all_data_TMA, by.x="Sample", by.y="Identifier")
median_expressions$Type <- as.character(median_expressions$Type)
median_expressions[which(median_expressions$Type == "BRCA1" | median_expressions$Type == "BRCA2"), "Type"] <- "BRCA1/2 mutated"

all_data_TMA_cat_BRCA <- median_expressions[which(median_expressions$Type== "BRCA1/2 mutated"),]

all_data_TMA_cat_BRCA <- median_expressions[which(median_expressions$Type== "BRCA1"| median_expressions$Type== "BRCA2"),]


all_data_TMA_cat_BRCA[which(all_data_TMA_cat_BRCA$median_Ki67 > 9.022178), "median_Ki67"] <- "high"
all_data_TMA_cat_BRCA[which(all_data_TMA_cat_BRCA$median_Ki67 <= 9.022178), "median_Ki67"] <- "low"

all_data_TMA_cat_HR <- median_expressions[which(median_expressions$Type == "HR"),]

all_data_TMA_cat_HR[which(all_data_TMA_cat_HR$median_Ki67 > 9.022178), "median_Ki67"] <- "high"
all_data_TMA_cat_HR[which(all_data_TMA_cat_HR$median_Ki67 <= 9.022178), "median_Ki67"] <- "low"





all_patients_cox <- all_data_TMA_cat_abundance


more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0



res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4..Effector.T.cells, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD8..T.cells, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ FOXP3.CD4..T.regs, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ IBA1.CD163..Macrophages, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)

#BRCAmut

all_patients_cox <- all_data_TMA_cat_BRCA_abundance

more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0



res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4..Effector.T.cells, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD8..T.cells, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ FOXP3.CD4..T.regs, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ IBA1.CD163..Macrophages, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



#HRwt

all_patients_cox <- all_data_TMA_cat_HR_abundance


more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4..Effector.T.cells, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD8..T.cells, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ FOXP3.CD4..T.regs, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ IBA1.CD163..Macrophages, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)

##Ki67


all_patients_cox <- all_data_TMA_cat_BRCA


more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0

res.cox <- coxph(Surv(PFI_time, status_PFI) ~ median_Ki67, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)

#HRwt

all_patients_cox <- all_data_TMA_cat_HR


more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ median_Ki67, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


###########################################################################################################
############################## MULTIVARIATE ###############################################################
###########################################################################################################


all_patients_cox <- all_data_TMA_cat_abundance


more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0



res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4..Effector.T.cells + CD8..T.cells + IBA1.CD163..Macrophages + FOXP3.CD4..T.regs + Type, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



all_patients_cox <- all_data_TMA_cat_BRCA_abundance
more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4..Effector.T.cells + IBA1.CD163..Macrophages + FOXP3.CD4..T.regs, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)




all_patients_cox_Ki67 <- all_data_TMA_cat_BRCA


all_patients_cox <- merge(all_patients_cox, all_patients_cox_Ki67[, c("Sample", "median_Ki67", "Type")], by.x="Patient",by.y="Sample")



more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0

#all_patients_cox <- all_patients_cox[which(all_patients_cox$Type.x == "BRCA1"),]

res.cox <- coxph(Surv(PFI_time, status_PFI) ~  CD4..Effector.T.cells +median_Ki67 + FOXP3.CD4..T.regs, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)

res.cox <- coxph(Surv(PFI_time, status_PFI) ~ CD4..Effector.T.cells + FOXP3.CD4..T.regs, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



############################################################################################################


#BRCAmut and HRwt

#OS and PFI and HR status

all_patients_cox <- all_data_TMA_cat_abundance

more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0

all_patients_cox$Type.x <- as.character(all_patients_cox$Type.x)
all_patients_cox[which(all_patients_cox$Type.x == "BRCA1/2 mutated"), "Type.x"] <- "mutBRCA"

res.cox <- coxph(Surv(PFI_time, status_PFI) ~ Type.x, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



all_patients_cox <- all_data_TMA_cat_abundance


all_patients_cox$Type.x <- as.character(all_patients_cox$Type.x)
all_patients_cox[which(all_patients_cox$Type.x == "BRCA1/2 mutated"), "Type.x"] <- "mutBRCA"

res.cox <- coxph(Surv(OS_time, status_OS) ~ Type.x, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)



################################################################################################################

#FOR IMMUNE AND TUMOR DIVERSITY
#get simpson diversity from "SIMPSON_CALC_TMA.R"

simpson <- data.frame(simpson_tumor, simpson.i)

median_SDI <- merge(simpson, all_data_TMA, by.x = "sample", by.y = "Identifier")


all_data_TMA_cat_abundance <- median_SDI
all_data_TMA_cat_BRCA_abundance <- median_SDI[which(median_SDI$Type.x == "BRCA1/2 mutated"),]
all_data_TMA_cat_HR_abundance<- median_SDI[which(median_SDI$Type.x == "HRwt"),]

for (i in seq_along(colnames(median_SDI))[c(2,3)]){
  #i=7  
  #all_data_TMA_cat_abundance <- median_expressions_abundance
  all_data_TMA_cat_abundance[, i] <- cut(all_data_TMA_cat_abundance[, i], breaks = c(-Inf, median(all_data_TMA_cat_abundance[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_BRCA_abundance[, i] <- cut(all_data_TMA_cat_BRCA_abundance[, i], breaks = c(-Inf, median(all_data_TMA_cat_BRCA_abundance[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_HR_abundance[, i] <- cut(all_data_TMA_cat_HR_abundance[, i], breaks = c(-Inf, median(all_data_TMA_cat_HR_abundance[, i]), Inf), labels=c("low", "high"))
  
}


all_patients_cox <- all_data_TMA_cat_abundance

more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0

res.cox <- coxph(Surv(PFI_time, status_PFI) ~ simpson_tumor, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ simpson_immune, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)

#BRCA

all_patients_cox <- all_data_TMA_cat_BRCA_abundance

more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0

res.cox <- coxph(Surv(PFI_time, status_PFI) ~ simpson_tumor, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ simpson_immune, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


#HRwt
all_patients_cox <- all_data_TMA_cat_HR_abundance

more.than.2.years <- which(all_patients_cox$PFI_time > 800)
all_patients_cox$PFI_time[more.than.2.years] <- 800
all_patients_cox$status_PFI[more.than.2.years] <- 0

res.cox <- coxph(Surv(PFI_time, status_PFI) ~ simpson_tumor, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)


res.cox <- coxph(Surv(PFI_time, status_PFI) ~ simpson_immune, data =  all_patients_cox)
summary(res.cox)
cox.zph(res.cox)

