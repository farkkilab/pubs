###########################################################################################################################################
############################# KAPLAN MEIER CURVES FOR IMMUNE CELL PROPORTIONS #############################################################
###########################################################################################################################################

#For Figures 2 and 3 and Supplementary Figures 2 and 3

library(survival)
library(survminer)
library(readxl)
library(data.table)
library(dplyr)

all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")

immune_clusters_sorted <- all_celltypes %>%
  group_by(Sample, GlobalCellType) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))


immune_clusters_sorted$freq <- immune_clusters_sorted$freq*100
immune_clusters_sorted <- immune_clusters_sorted[, -3]
immune_clusters_sorted_wide <- pivot_wider(immune_clusters_sorted, names_from="Sample", values_from="freq")
immune_clusters_sorted_wide <- as.data.frame(immune_clusters_sorted_wide)
rownames(immune_clusters_sorted_wide) <- immune_clusters_sorted_wide$GlobalCellType
immune_clusters_sorted_wide <- immune_clusters_sorted_wide[,-1]


immune_clusters_sorted_wide <- as.data.table(t(immune_clusters_sorted_wide), keep.colnames = T, keep.rownames = T)


rownames(immune_clusters_sorted_wide) <- immune_clusters_sorted_wide$rn
colnames(immune_clusters_sorted_wide)[1] <- "Patient"
immune_clusters_sorted_wide[is.na(immune_clusters_sorted_wide)] <- 0
immune_clusters_sorted_wide <- as.data.frame(immune_clusters_sorted_wide)




clinical_data_TMA <- read_excel("TMA_clinicaldata.xlsx")

#creating status_PFI and status_dead

clinical_data_TMA$status_PFI <- 1
clinical_data_TMA$status_PFI[which(clinical_data_TMA$`Recurrence`=="no recurrence")] <- 0

clinical_data_TMA$status_OS <- 1
clinical_data_TMA$status_OS[which(clinical_data_TMA$`Deceased`=="No")] <- 0
clinical_data_TMA$status_OS[which(clinical_data_TMA$`Deceased`=="Lost")] <- 0
#merging merged_data_TMA and clinical_data_TMA

all_data_TMA_surv <- clinical_data_TMA
#colnames(all_data_TMA_surv)[which(names(all_data_TMA_surv)=="PFI(days)")] <- "PFI_time"
#colnames(all_data_TMA_surv)[which(names(all_data_TMA_surv)=="OS(days)")] <- "OS_time"



all_data_TMA_surv[which(all_data_TMA_surv$Category == "BRCA1" | all_data_TMA_surv$Category == "BRCA2"), "Category"] <- "BRCA1/2 mutated"


median_expressions <- merge(immune_clusters_sorted_wide, all_data_TMA_surv[, c("status_PFI", "PFI_time", "status_OS", "OS_time", "Category", "Identifier")], by.x="Patient", by.y="Identifier")

theme <- theme(panel.border = element_rect(colour = "black", size=1.7, fill=NA), 
               panel.background = element_blank(), plot.title=element_text(hjust=0.5)) 


colnames(median_expressions)[c(3, 4:10, 14,16)] <- c("Bcells", "CD11cAPC", "CD11cCD163IBA1Macrophages", "CD11cIBA1Macrophages","CD163IBA1Macrophages", "CD163Macrophages","CD4", "CD8", "Tregs", "IBA1Macrophages")


names <-c("Patient","Apoptotic","Bcells", "CD11c+APC", "CD11c+CD163+IBA1+Macrophages", "CD11c+IBA1+Macrophages","CD163+IBA1+Macrophages", "CD163+Macrophages","CD4", "CD8", "EMT","Endothelia","Epithelial","Tregs", "Functional epithelial","IBA1+Macrophages")


for (i in seq_along(colnames(median_expressions))[c(3, 4:10, 14,16)]){
  
  
  median_all <- median_expressions
  
  median_expressions_BRCA <- median_expressions[which(median_expressions$Category== "BRCA1/2 mutated"),]
  
  
  median_expressions_HR <- median_expressions[which(median_expressions$Category == "HR"),]
  
  all_data_TMA_cat <- median_all
  all_data_TMA_cat[, i] <- cut(all_data_TMA_cat[, i], breaks = c(-Inf, median(all_data_TMA_cat[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_BRCA <- median_expressions_BRCA
  all_data_TMA_cat_BRCA[, i] <- cut(all_data_TMA_cat_BRCA[, i], breaks = c(-Inf, median(all_data_TMA_cat_BRCA[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_HR <- median_expressions_HR
  all_data_TMA_cat_HR[, i] <- cut(all_data_TMA_cat_HR[, i], breaks = c(-Inf, median(all_data_TMA_cat_HR[, i]), Inf), labels=c("low", "high"))
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat_BRCA)[i],", data = all_data_TMA_cat_BRCA)")
  eval(parse(text=formula))
  
  p1 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat_BRCA, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("BRCA1/2\n",names[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   ggtheme = theme, legend.title="", xlab="PFI(days)", ylab="Survival",
                   legend.labs=c("low %", "high %"))+ guides(colour=guide_legend(nrow=2))
  p1 <- ggpar(p1, 
              font.main = c(20),
              font.x = c(20),
              font.y = c(20),
              font.caption = c(16), 
              font.legend = c(20), 
              font.tickslab = c(16))
  
  print(p1)
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat_HR)[i],", data = all_data_TMA_cat_HR)")
  eval(parse(text=formula))
  
  p2 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat_HR, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("HRwt\n",names[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   ggtheme = theme, legend.title="", xlab="PFI(days)", ylab="Survival", legend.labs=c("low %", "high %"))+ guides(colour=guide_legend(nrow=2))
  p2 <- ggpar(p2, 
              font.main = c(20),
              font.x = c(20),
              font.y = c(20),
              font.caption = c(16), 
              font.legend = c(20), 
              font.tickslab = c(16))
  
  print(p2)
  
  
  
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat)[i],", data = all_data_TMA_cat)")
  eval(parse(text=formula))
  
  p3 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("All patients\n",names[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   ggtheme = theme, legend.title="", xlab="PFI(days)", ylab="Survival",legend.labs=c("low %", "high %"))+ guides(colour=guide_legend(nrow=2))
  p3 <- ggpar(p3, 
              font.main = c(20),
              font.x = c(20),
              font.y = c(20),
              font.caption = c(16), 
              font.legend = c(20), 
              font.tickslab = c(16))
  
  print(p3)
  
  
  
}



#######################################################################################################################################################################
################################ KAPLAN MEIER GRAPHS FOR IMMUNE AND TUMOR SDI #########################################################################################
#######################################################################################################################################################################
#simpson data from "SIMPSON_CALC_TMA.R"
simpson_data <- data.frame(simpson, simpson.i)
simpson_data <- merge(simpson_data, all_data_TMA_surv[, c("status_PFI", "PFI_time", "status_OS", "OS_time", "Category", "Identifier")], by.x="sample", by.y="Identifier")

theme <- theme(panel.border = element_rect(colour = "black", size=1.7, fill=NA), 
               panel.background = element_blank(), plot.title=element_text(hjust=0.5)) 


for (i in seq_along(colnames(simpson_data))[c(2:3)]){
  
  
  median_all <- simpson_data
  
  median_expressions_BRCA <- simpson_data[which(simpson_data$Category== "BRCA1/2 mutated"),]
  
  
  median_expressions_HR <- simpson_data[which(simpson_data$Category == "HRwt"),]
  
  all_data_TMA_cat <- median_all
  all_data_TMA_cat[, i] <- cut(all_data_TMA_cat[, i], breaks = c(-Inf, median(all_data_TMA_cat[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_BRCA <- median_expressions_BRCA
  all_data_TMA_cat_BRCA[, i] <- cut(all_data_TMA_cat_BRCA[, i], breaks = c(-Inf, median(all_data_TMA_cat_BRCA[, i]), Inf), labels=c("low", "high"))
  
  
  all_data_TMA_cat_HR <- median_expressions_HR
  all_data_TMA_cat_HR[, i] <- cut(all_data_TMA_cat_HR[, i], breaks = c(-Inf, median(all_data_TMA_cat_HR[, i]), Inf), labels=c("low", "high"))
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat_BRCA)[i],", data = all_data_TMA_cat_BRCA)")
  eval(parse(text=formula))
  
  p1 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat_BRCA, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("BRCA1/2\n",names[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   ggtheme = theme, legend.title="", xlab="PFI(days)", ylab="Survival",
                   legend.labs=c("low", "high"))+ guides(colour=guide_legend(nrow=2))
  p1 <- ggpar(p1, 
              font.main = c(20),
              font.x = c(20),
              font.y = c(20),
              font.caption = c(16), 
              font.legend = c(20), 
              font.tickslab = c(16))
  
  print(p1)
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat_HR)[i],", data = all_data_TMA_cat_HR)")
  eval(parse(text=formula))
  
  p2 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat_HR, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("HRwt\n",names[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   ggtheme = theme, legend.title="", xlab="PFI(days)", ylab="Survival", legend.labs=c("low", "high"))+ guides(colour=guide_legend(nrow=2))
  p2 <- ggpar(p2, 
              font.main = c(20),
              font.x = c(20),
              font.y = c(20),
              font.caption = c(16), 
              font.legend = c(20), 
              font.tickslab = c(16))
  
  print(p2)
  
  
  
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat)[i],", data = all_data_TMA_cat)")
  eval(parse(text=formula))
  
  p3 <- ggsurvplot(fit, size=1.5,data=all_data_TMA_cat, cencor=T, pval=T, pval.coord=c(1200, 0.80), 
                   title=paste("All patients\n",names[i],"and PFI"), risk.table=T, risk.table.height = 0.2, pval.size=6,
                   xlim = c(0, 2000), break.time.by=400, palette = c("deepskyblue", "sienna2"), risk.table.y.text=F,legend=c(0.15, 0.15),
                   ggtheme = theme, legend.title="", xlab="PFI(days)", ylab="Survival",legend.labs=c("low", "high"))+ guides(colour=guide_legend(nrow=2))
  p3 <- ggpar(p3, 
              font.main = c(20),
              font.x = c(20),
              font.y = c(20),
              font.caption = c(16), 
              font.legend = c(20), 
              font.tickslab = c(16))
  
  print(p3)
  
  
  
}






#######################################################################################################################################################################
################################ KAPLAN MEIER GRAPHS FOR KI67 EXPRESSION ##############################################################################################
#######################################################################################################################################################################


profepi_surv <- all_celltypes[which(all_celltypes$GlobalCellType == "Proliferating epithelial"),]

#find out median Ki67 expression in each patient
median_expressions <- profepi_surv %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                                      median_pSTAT1 = median(pSTAT1), median_P21 = median(P21), median_PD1=median(PD1), median_PDL1=median(PDL1), median_cCasp3=median(cCasp3))

median_expressions <- merge(median_expressions, all_data_TMA_surv, by.x="Sample", by.y="Identifier")

quantile(median_expressions$median_Ki67, 0.6667)
#9.022178

for (i in seq_along(colnames(median_expressions))[2]){
  
  
  all_data_TMA_cat_BRCA <- median_expressions[which(median_expressions$Category== "BRCA1/2 mutated"),]
  
  all_data_TMA_cat_BRCA[which(all_data_TMA_cat_BRCA$median_Ki67 > 9.022178), "median_Ki67"] <- "high"
  all_data_TMA_cat_BRCA[which(all_data_TMA_cat_BRCA$median_Ki67 <= 9.022178), "median_Ki67"] <- "low"
  
  all_data_TMA_cat_HR <- median_expressions[which(median_expressions$Category == "HR"),]
  
  all_data_TMA_cat_HR[which(all_data_TMA_cat_HR$median_Ki67 > 9.022178), "median_Ki67"] <- "high"
  all_data_TMA_cat_HR[which(all_data_TMA_cat_HR$median_Ki67 <= 9.022178), "median_Ki67"] <- "low"
  
  
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat_BRCA)[2],", data = all_data_TMA_cat_BRCA)")
  eval(parse(text=formula))
  
  p <- ggsurvplot(fit, data=all_data_TMA_cat_BRCA, cencor=T, pval=T, pval.coord=c(1000, 0.75), title=paste("BRCA1/2", colnames(all_data_TMA_cat_BRCA)[2],"PFI"), risk.table=T, risk.table.height = 0.3, xlim = c(0, 1500), break.time.by=500, conf.int = T)
  print(p)
  
  
 
  formula <- paste0("fit <- survfit(Surv(PFI_time,status_PFI) ~ ",colnames(all_data_TMA_cat_HR)[2],", data = all_data_TMA_cat_HR)")
  eval(parse(text=formula))
  
  p <- ggsurvplot(fit, data=all_data_TMA_cat_HR, cencor=T, pval=T, pval.coord=c(750, 0.75), title=paste("HRwt", colnames(all_data_TMA_cat_HR)[2],"PFI"), risk.table=T, risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, conf.int = T)
  print(p)
  
}




