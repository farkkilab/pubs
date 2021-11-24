################################################################################################################################################################
###################################### NEIGHBORHOOD PAIRS HEATMAP ##############################################################################################
################################################################################################################################################################

#For Supplementary Figure 5

library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(readxl)


setwd("interaction_frequencies")
temp= list.files(pattern=".csv")
myfiles = lapply(temp, read.csv)

for (i in c(1:112)){
  
  myfiles[[i]]$core <- paste(temp)[i]
}



library(dplyr)

all_fractions <- bind_rows(myfiles[[1]], myfiles[[2]], myfiles[[3]], myfiles[[4]],myfiles[[5]], myfiles[[6]], myfiles[[7]], myfiles[[8]],
                           myfiles[[9]], myfiles[[10]], myfiles[[11]], myfiles[[12]],myfiles[[13]], myfiles[[14]], myfiles[[15]], myfiles[[16]],
                           myfiles[[17]], myfiles[[18]], myfiles[[19]], myfiles[[20]],myfiles[[21]], myfiles[[22]], myfiles[[23]], myfiles[[24]],
                           myfiles[[25]], myfiles[[26]], myfiles[[27]], myfiles[[28]],myfiles[[29]], myfiles[[30]], myfiles[[31]], myfiles[[32]],
                           myfiles[[33]], myfiles[[34]], myfiles[[35]], myfiles[[36]],myfiles[[37]], myfiles[[38]], myfiles[[39]], myfiles[[40]],
                           myfiles[[41]], myfiles[[42]], myfiles[[43]], myfiles[[44]],myfiles[[45]], myfiles[[46]], myfiles[[47]], myfiles[[48]],
                           myfiles[[49]], myfiles[[50]], myfiles[[51]], myfiles[[52]],myfiles[[53]], myfiles[[54]], myfiles[[55]], myfiles[[56]],
                           myfiles[[57]], myfiles[[58]], myfiles[[59]],myfiles[[60]], myfiles[[61]], myfiles[[62]],myfiles[[63]], myfiles[[64]], myfiles[[65]], myfiles[[66]],
                           myfiles[[67]], myfiles[[68]], myfiles[[69]],myfiles[[70]], myfiles[[71]], myfiles[[72]],myfiles[[73]], myfiles[[74]], myfiles[[75]], myfiles[[76]],
                           myfiles[[77]], myfiles[[78]], myfiles[[79]],myfiles[[80]], myfiles[[81]], myfiles[[82]],myfiles[[83]], myfiles[[84]], myfiles[[85]], myfiles[[86]],
                           myfiles[[87]], myfiles[[88]], myfiles[[89]],myfiles[[90]], myfiles[[91]], myfiles[[92]],myfiles[[93]], myfiles[[94]], myfiles[[95]], myfiles[[96]],
                           myfiles[[97]], myfiles[[98]], myfiles[[99]],myfiles[[100]], myfiles[[101]], myfiles[[102]],myfiles[[103]], myfiles[[104]], myfiles[[105]], myfiles[[106]],
                           myfiles[[107]], myfiles[[108]], myfiles[[109]],myfiles[[110]],myfiles[[111]],myfiles[[112]])






tma_mapping <- read_excel("TMA_mapping.xlsx")
clinical_data_TMA <- read_excel("TMA_clinicaldata.xlsx")

tma_mapping <- merge(tma_mapping, clinical_data_TMA[, c(1:2)], by.x= "Patient", by.y="Identifier")
tma_mapping <- tma_mapping[, c(1, 3, 6)]

tma_mapping <- na.omit(tma_mapping)


all_fractions$core <- gsub("^.{0,14}", "", all_fractions$core)
all_fractions$core <- substr(all_fractions$core, 1, nchar(all_fractions$core) -4)


all_fractions <- merge(all_fractions, tma_mapping, by.x = "core", by.y="TMA1", all=T)

all_fractions[which(all_fractions$Category == "BRCA1"), "Category"] <- "BRCA1/2 mutated"
all_fractions[which(all_fractions$Category == "BRCA2"), "Category"] <- "BRCA1/2 mutated"


all_fractions_mean <- all_fractions %>% group_by(cluster, core) %>% summarise(Apoptotic = mean(Apoptotic),
                                                                                 B.cells = mean(B.cells),
                                                                                 CD11c.APC = mean(CD11c.APC),
                                                                                 CD11c.CD163.IBA1.Macrophages = mean(CD11c.CD163.IBA1.Macrophages),
                                                                                 CD163.IBA1.Macrophages = mean(CD163.IBA1.Macrophages),
                                                                                 CD163.Macrophages= mean(CD163.Macrophages),
                                                                                 CD4 = mean(CD4),
                                                                                 CD8 = mean(CD8),
                                                                                 EMT = mean(EMT),
                                                                                 Endothelia = mean(Endothelia),
                                                                                 Epithelial = mean(Epithelial),
                                                                                 FOXP3.CD4.Tregs = mean(FOXP3.CD4.Tregs),
                                                                                 High.PDL1 = mean(High.PDL1),
                                                                                 High.Vimentin = mean(High.Vimentin),
                                                                                 High.proliferative_Stroma = mean(High.proliferative_Stroma),
                                                                                 Low.vimentin = mean(Low.vimentin),
                                                                                 Low_eccentricity_medium_vimentin = mean(Low_eccentricity_medium_vimentin),
                                                                                 Mesenchymal = mean(Mesenchymal),
                                                                                 Non.proliferative_Stroma = mean(Non.proliferative_Stroma),
                                                                                 Proliferating.EMT = mean(Proliferating.EMT),
                                                                                 Proliferating.epithelial = mean(Proliferating.epithelial),
                                                                                 Proliferative_Stroma = mean(Proliferative_Stroma),
                                                                                 CD11c.IBA1.Macrophages = mean(CD11c.IBA1.Macrophages),
                                                                                 IBA1.Macrophages = mean(IBA1.Macrophages),
                                                                                 High.P21 = mean(High.P21),
                                                                                 Hyperfunctional.epithelial = mean(Hyperfunctional.epithelial))


all_fractions_mean <- merge(all_fractions_mean, tma_mapping[, c("Patient", "Category")], by= "Patient")


all_fractions_mean[which(all_fractions_mean$Category == "BRCA1"), "Category"] <- "BRCA1/2 mutated"
all_fractions_mean[which(all_fractions_mean$Category == "BRCA2"), "Category"] <- "BRCA1/2 mutated"

all_fractions_mean_unique <- unique(all_fractions_mean)


all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")
all_celltypes$GlobalCellType <- as.character(all_celltypes$GlobalCellType)

all_celltypes[which(all_celltypes$GlobalCellType == "CD11c+CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+CD11c+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD11c+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD11c+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD4"), "GlobalCellType"] <- "CD4+ Effector T-cells"
all_celltypes[which(all_celltypes$GlobalCellType == "CD8"), "GlobalCellType"] <- "CD8+ T-cells"
all_celltypes[which(all_celltypes$GlobalCellType == "CD163+Macrophages"), "GlobalCellType"] <- "CD163+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD11c+APC"), "GlobalCellType"] <- "CD11c+APC"
all_celltypes[which(all_celltypes$GlobalCellType == "FOXP3+CD4+Tregs"), "GlobalCellType"] <- "FOXP3+CD4+ T-regs"
all_celltypes[which(all_celltypes$GlobalCellType == "High-PDL1"), "GlobalCellType"] <- "Functional stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "Non-proliferative_Stroma"), "GlobalCellType"] <- "Non-proliferative Stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "Low_eccentricity_medium_vimentin"), "GlobalCellType"] <- "Low eccentricity"
all_celltypes[which(all_celltypes$GlobalCellType == "Proliferative_Stroma"), "GlobalCellType"] <- "Proliferative Stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "High-proliferative_Stroma"), "GlobalCellType"] <- "High-proliferative Stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "Hyperfunctional epithelial"), "GlobalCellType"] <- "Functional epithelial"

all_clusters <- all_celltypes[, c(2, 46)]

metacluster_percentages <- all_clusters %>% group_by(cores, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
metacluster_percentages <- as.data.frame(metacluster_percentages)

metacluster_percentages <- metacluster_percentages[, -c(3)]


metacluster_percentages <- metacluster_percentages %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

metacluster_percentages <- as.data.frame(metacluster_percentages)

metacluster_percentages[is.na(metacluster_percentages)] <- 0

#the rows of all_fractions are the centering cell types

#the columns are the neighboring cell types


colnames(metacluster_percentages)[3] <- "CD11c.APC"
colnames(metacluster_percentages)[6] <- "FOXP3.CD4.Tregs"
colnames(metacluster_percentages)[8] <- "High.PDL1"
colnames(metacluster_percentages)[9] <- "High.proliferative_Stroma"
colnames(metacluster_percentages)[7] <- "Hyperfunctional.epithelial"
colnames(metacluster_percentages)[11] <- "Proliferating.EMT"
colnames(metacluster_percentages)[12] <- "Proliferating.epithelial"
colnames(metacluster_percentages)[21] <- "CD11c.CD163.IBA1.Macrophages"
colnames(metacluster_percentages)[19] <- "CD11c.IBA1.Macrophages"
colnames(metacluster_percentages)[20] <- "CD163.IBA1.Macrophages"
colnames(metacluster_percentages)[17] <- "High.Vimentin"
colnames(metacluster_percentages)[18] <- "IBA1.Macrophages"
colnames(metacluster_percentages)[22] <- "Low.vimentin"
colnames(metacluster_percentages)[24] <- "Non.proliferative_Stroma"
colnames(metacluster_percentages)[25] <- "B.cells"
colnames(metacluster_percentages)[26] <- "CD163.Macrophages"
colnames(metacluster_percentages)[27] <- "High.P21"
colnames(metacluster_percentages)[14] <- "CD4"
colnames(metacluster_percentages)[15] <- "CD8"
colnames(metacluster_percentages)[13] <- "Proliferative_Stroma"
colnames(metacluster_percentages)[23] <- "Low_eccentricity_medium_vimentin"



all_fractions_merged_abundance <- merge(all_fractions, metacluster_percentages, by.x="core",by.y="cores")

#now for each row divide the proportion being neighbor with the proportion of abundance


all_fractions_merged_abundance$Apoptotic.x <- all_fractions_merged_abundance$Apoptotic.x/all_fractions_merged_abundance$Apoptotic.y
all_fractions_merged_abundance$B.cells.x <- all_fractions_merged_abundance$B.cells.x/all_fractions_merged_abundance$B.cells.y
all_fractions_merged_abundance$CD11c.APC.x <- all_fractions_merged_abundance$CD11c.APC.x/all_fractions_merged_abundance$CD11c.APC.y
all_fractions_merged_abundance$CD11c.CD163.IBA1.Macrophages.x <- all_fractions_merged_abundance$CD11c.CD163.IBA1.Macrophages.x/all_fractions_merged_abundance$CD11c.CD163.IBA1.Macrophages.y
all_fractions_merged_abundance$CD163.IBA1.Macrophages.x <- all_fractions_merged_abundance$CD163.IBA1.Macrophages.x/all_fractions_merged_abundance$CD163.IBA1.Macrophages.y
all_fractions_merged_abundance$CD163.Macrophages.x <- all_fractions_merged_abundance$CD163.Macrophages.x/all_fractions_merged_abundance$CD163.Macrophages.y
all_fractions_merged_abundance$CD4.x <- all_fractions_merged_abundance$CD4.x/all_fractions_merged_abundance$CD4.y
all_fractions_merged_abundance$CD8.x <- all_fractions_merged_abundance$CD8.x/all_fractions_merged_abundance$CD8.y
all_fractions_merged_abundance$EMT.x <- all_fractions_merged_abundance$EMT.x/all_fractions_merged_abundance$EMT.y
all_fractions_merged_abundance$Endothelia.x <- all_fractions_merged_abundance$Endothelia.x/all_fractions_merged_abundance$Endothelia.y
all_fractions_merged_abundance$Epithelial.x <- all_fractions_merged_abundance$Epithelial.x/all_fractions_merged_abundance$Epithelial.y
all_fractions_merged_abundance$FOXP3.CD4.Tregs.x <- all_fractions_merged_abundance$FOXP3.CD4.Tregs.x/all_fractions_merged_abundance$FOXP3.CD4.Tregs.y
all_fractions_merged_abundance$High.PDL1.x <- all_fractions_merged_abundance$High.PDL1.x/all_fractions_merged_abundance$High.PDL1.y
all_fractions_merged_abundance$High.Vimentin.x <- all_fractions_merged_abundance$High.Vimentin.x/all_fractions_merged_abundance$High.Vimentin.y
all_fractions_merged_abundance$High.proliferative_Stroma.x <- all_fractions_merged_abundance$High.proliferative_Stroma.x/all_fractions_merged_abundance$High.proliferative_Stroma.y
all_fractions_merged_abundance$Low.vimentin.x <- all_fractions_merged_abundance$Low.vimentin.x/all_fractions_merged_abundance$Low.vimentin.y
all_fractions_merged_abundance$Low_eccentricity_medium_vimentin.x <- all_fractions_merged_abundance$Low_eccentricity_medium_vimentin.x/all_fractions_merged_abundance$Low_eccentricity_medium_vimentin.y
all_fractions_merged_abundance$Mesenchymal.x <- all_fractions_merged_abundance$Mesenchymal.x/all_fractions_merged_abundance$Mesenchymal.y
all_fractions_merged_abundance$Non.proliferative_Stroma.x <- all_fractions_merged_abundance$Non.proliferative_Stroma.x/all_fractions_merged_abundance$Non.proliferative_Stroma.y
all_fractions_merged_abundance$Proliferating.EMT.x <- all_fractions_merged_abundance$Proliferating.EMT.x/all_fractions_merged_abundance$Proliferating.EMT.y
all_fractions_merged_abundance$Proliferating.epithelial.x <- all_fractions_merged_abundance$Proliferating.epithelial.x/all_fractions_merged_abundance$Proliferating.epithelial.y
all_fractions_merged_abundance$Proliferative_Stroma.x <- all_fractions_merged_abundance$Proliferative_Stroma.x/all_fractions_merged_abundance$Proliferative_Stroma.y
all_fractions_merged_abundance$CD11c.IBA1.Macrophages.x <- all_fractions_merged_abundance$CD11c.IBA1.Macrophages.x/all_fractions_merged_abundance$CD11c.IBA1.Macrophages.y
all_fractions_merged_abundance$IBA1.Macrophages.x <- all_fractions_merged_abundance$IBA1.Macrophages.x/all_fractions_merged_abundance$IBA1.Macrophages.y
all_fractions_merged_abundance$High.P21.x <- all_fractions_merged_abundance$High.P21.x/all_fractions_merged_abundance$High.P21.y
all_fractions_merged_abundance$Hyperfunctional.epithelial.x <- all_fractions_merged_abundance$Hyperfunctional.epithelial.x/all_fractions_merged_abundance$Hyperfunctional.epithelial.y

all_fractions_merged_abundance$Category <- as.character(all_fractions_merged_abundance$Category)
all_fractions_merged_abundance[which(all_fractions_merged_abundance$Category == "BRCA1"), "Category"] <- "BRCA1/2 mutated"
all_fractions_merged_abundance[which(all_fractions_merged_abundance$Category == "BRCA2"), "Category"] <- "BRCA1/2 mutated"



p_values <- all_fractions_merged_abundance %>% group_by(cluster)%>%
  summarise_each(funs(wilcox.test(.[Category == "BRCA1/2 mutated"], .[Category == "HR"])$p.value), vars = Apoptotic.x:Hyperfunctional.epithelial.x)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]

colnames(p_values) <- colnames(all_fractions_merged_abundance)[c(3:30)]


cores <- read_excel("C:/LocalData/ingamari/clinical_data/TMA_mapping.xlsx")
cores <- cores[, c(2, 5)]

data_for_hm <- merge(all_fractions_merged_abundance, cores, by.y="TMA1", by.x="core")




data_for_hm_mean <- data_for_hm %>% group_by(Patient.x, cluster) %>% summarise(CD163.IBA1.Macrophages = mean(CD163.IBA1.Macrophages.x),
                                                                               CD163.Macrophages = mean(CD163.Macrophages.x),
                                                                               CD4 = mean(CD4.x),
                                                                               CD8 = mean(CD8.x),
                                                                               Endothelia = mean(Endothelia.x),
                                                                               Epithelial = mean(Epithelial.x), FOXP3.CD4.Tregs = mean(FOXP3.CD4.Tregs.x),
                                                                               High.PDL1 = mean(High.PDL1.x), 
                                                                               Non.proliferative_Stroma = mean(Non.proliferative_Stroma.x),
                                                                               Proliferating.epithelial = mean(Proliferating.epithelial.x))



data_for_hm_mean <- pivot_longer(data_for_hm_mean, cols = 3:12, names_to = "Neighbor")

data_for_hm_mean[which(data_for_hm_mean$Neighbor == "Proliferating.epithelial"), "Neighbor"] <- "Proliferating epithelial"
data_for_hm_mean[which(data_for_hm_mean$Neighbor == "High.PDL1"), "Neighbor"] <- "High-PDL1"
data_for_hm_mean[which(data_for_hm_mean$Neighbor == "Non.proliferative_Stroma"), "Neighbor"] <- "Non-proliferative_Stroma"
data_for_hm_mean[which(data_for_hm_mean$Neighbor == "CD163.IBA1.Macrophages"), "Neighbor"] <- "CD163+IBA1+Macrophages"
data_for_hm_mean[which(data_for_hm_mean$Neighbor == "CD163.Macrophages"), "Neighbor"] <- "CD163+Macrophages"
data_for_hm_mean[which(data_for_hm_mean$Neighbor == "FOXP3.CD4.Tregs"), "Neighbor"] <- "FOXP3+CD4+Tregs"


data_for_hm_mean <- data_for_hm_mean[which(data_for_hm_mean$cluster == "CD4" |data_for_hm_mean$cluster == "CD8" |data_for_hm_mean$cluster == "FOXP3+CD4+Tregs" |data_for_hm_mean$cluster == "Proliferating epithelial" |data_for_hm_mean$cluster == "Epithelial" |data_for_hm_mean$cluster == "High-PDL1" |data_for_hm_mean$cluster == "Non-proliferative_Stroma" |data_for_hm_mean$cluster == "Endothelia" |data_for_hm_mean$cluster == "CD163+IBA1+Macrophages" |data_for_hm_mean$cluster == "CD163+Macrophages"),]

#WE ARE INTERESTED IN THE FOLLOWING PAIRS

#Proliferating epithelial - CD4
#Proliferating epithelial - CD8
#High-PDL1 - CD163+Macrophages
#Proliferating epithelial - CD163+IBA1+Macrophages
#Proliferating epithelial - Endothelia
#Proliferating epithelial - High-PDL1
#Proliferating epithelial - FOXP3+CD4+Tregs
#Non-proliferative_Stroma - FOXP3+CD4+Tregs
#CD4 - CD4
#Epithelial - CD8




data_for_hm_mean_agg_profepiCD4 <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_profepiCD4 = mean(value[cluster == "Proliferating epithelial" & Neighbor == "CD4"| cluster == "CD4" & Neighbor == "Proliferating epithelial"], na.rm = T))


data_for_hm_mean_agg_profepiCD8 <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_profepiCD8 = mean(value[cluster == "Proliferating epithelial" & Neighbor == "CD8"| cluster == "CD8" & Neighbor == "Proliferating epithelial"], na.rm = T))


data_for_hm_mean_agg_profepiCD163IBA1 <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_profepiCD163IBA1 = mean(value[cluster == "Proliferating epithelial" & Neighbor == "CD163+IBA1+Macrophages"| cluster == "CD163+IBA1+Macrophages" & Neighbor == "Proliferating epithelial"], na.rm = T))


data_for_hm_mean_agg_profepiTregs <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_profepiTregs = mean(value[cluster == "Proliferating epithelial" & Neighbor == "FOXP3+CD4+Tregs"| cluster == "FOXP3+CD4+Tregs" & Neighbor == "Proliferating epithelial"], na.rm = T))


data_for_hm_mean_agg_profepiEndothelia <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_profepiEndothelia = mean(value[cluster == "Proliferating epithelial" & Neighbor == "Endothelia"| cluster == "Endothelia" & Neighbor == "Proliferating epithelial"], na.rm = T))


data_for_hm_mean_agg_profepiFunctionalStroma <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_profepiFuncStroma = mean(value[cluster == "Proliferating epithelial" & Neighbor == "High-PDL1"| cluster == "High-PDL1" & Neighbor == "Proliferating epithelial"], na.rm = T))


data_for_hm_mean_agg_CD163FunctionalStroma <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_CD163FuncStroma = mean(value[cluster == "CD163+Macrophages" & Neighbor == "High-PDL1"| cluster == "High-PDL1" & Neighbor == "CD163+Macrophages"], na.rm = T))


data_for_hm_mean_agg_CD8Epithelial <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_CD8Epithelial = mean(value[cluster == "CD8" & Neighbor == "Epithelial"| cluster == "Epithelial" & Neighbor == "CD8"], na.rm = T))


data_for_hm_mean_agg_TregNonprofStroma <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_TregNonprofStroma = mean(value[cluster == "FOXP3+CD4+Tregs" & Neighbor == "Non-proliferative_Stroma"| cluster == "Non-proliferative_Stroma" & Neighbor == "FOXP3+CD4+Tregs"], na.rm = T))

data_for_hm_mean_agg_CD4CD4 <- data_for_hm_mean %>% group_by(Patient.x) %>% summarise(mean_CD4CD4 = mean(value[cluster == "CD4" & Neighbor == "CD4"| cluster == "CD4" & Neighbor == "CD4"], na.rm = T))

data_for_hm_mean_agg <- bind_cols(data_for_hm_mean_agg_CD163FunctionalStroma, data_for_hm_mean_agg_CD4CD4[, 2], data_for_hm_mean_agg_CD8Epithelial[, 2], data_for_hm_mean_agg_profepiCD163IBA1[, 2], data_for_hm_mean_agg_profepiCD4[, 2], data_for_hm_mean_agg_profepiCD8[, 2], data_for_hm_mean_agg_profepiEndothelia[, 2], data_for_hm_mean_agg_profepiFunctionalStroma[, 2], data_for_hm_mean_agg_profepiTregs[, 2], data_for_hm_mean_agg_TregNonprofStroma[, 2])




HRstatus <- read.csv("C:/LocalData/ingamari/clinical_data/all_data_TMA.csv")
HRstatus <- HRstatus[, c("Identifier", "Type", "PFI_time")]


HRstatus <- HRstatus[which(HRstatus$Identifier %in% data_for_hm_mean_agg$Patient.x),]

HRstatus[which(HRstatus$PFI_time > 365), "PFI_time"] <- "long"
HRstatus[which(HRstatus$PFI_time < 365), "PFI_time"] <- "short"
HRstatus[which(HRstatus$PFI_time == 56), "PFI_time"] <- "short"


#heatmap ja annotation!



colours<- list("Type"=c("BRCA1"="royalblue","BRCA2"="lightblue", "HR"="red3"), "PFI_time"=c("long" = "#3399CC", "short"="#003366"))



Rowann <- HeatmapAnnotation(df=HRstatus[, c(2:3)], which="row", col=colours,
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(1, "mm"))

#ht_opt(COLUMN_ANNO_PADDING=unit(2, "mm"), ROW_ANNO_PADDING=unit(1, "mm"))


data_for_hm_mean_agg_hm <- data_for_hm_mean_agg
data_for_hm_mean_agg_hm <- as.data.frame(data_for_hm_mean_agg_hm)
rownames(data_for_hm_mean_agg_hm) <- data_for_hm_mean_agg_hm$Patient.x
data_for_hm_mean_agg_hm <- data_for_hm_mean_agg_hm[, -1]

hmap <- Heatmap(scale(as.matrix(data_for_hm_mean_agg_hm)), clustering_method_rows ="ward.D2",name="Neighborhoods normalized", 
                row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                left_annotation = Rowann, row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "Interactions", row_title = "patients",
                height = unit(10.5, "cm") , width = unit(3, "cm"),border="white",
                rect_gp = gpar(col = "white", lwd = 2))


draw(hmap)



##################################################################################################################################################################
################################### HEATMAP OF CELL PERCENTAGES AND NEIGHBORHOODS COMBINED ######################################################################
##################################################################################################################################################################


combined_percentages_neighborhood <- cbind(data_for_hm_mean_agg, metacluster_percentages)
combined_percentages_neighborhood <- as.data.frame(combined_percentages_neighborhood)
combined_percentages_neighborhood <- combined_percentages_neighborhood[, -1]


hmap <- Heatmap(scale(as.matrix(combined_percentages_neighborhood)), clustering_method_rows ="ward.D2",name="Z-score", 
                row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                left_annotation = Rowann, row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "Interactions and cell percentages", row_title = "patients",
                height = unit(10.5, "cm") , width = unit(9, "cm"),border="white",
                rect_gp = gpar(col = "white", lwd = 2))


draw(hmap)



############################################################################################################################################################################
####################################### ENDOTHELIAL CELL NEIGHBORHOODS #####################################################################################################
############################################################################################################################################################################


data_for_hm_mean <- data_for_hm %>% group_by(Patient.x, cluster) %>% summarise(CD163.IBA1.Macrophages = mean(CD163.IBA1.Macrophages.x),
                                                                               CD163.Macrophages = mean(CD163.Macrophages.x),
                                                                               CD4 = mean(CD4.x),
                                                                               CD8 = mean(CD8.x),
                                                                               CD11c.APC = mean(CD11c.APC.x),
                                                                               IBA1.Macrophages = mean(IBA1.Macrophages.x), FOXP3.CD4.Tregs = mean(FOXP3.CD4.Tregs.x),
                                                                               CD11c.IBA1.Macrophages = mean(CD11c.IBA1.Macrophages.x), 
                                                                               CD11c.CD163.IBA1.Macrophages = mean(CD11c.CD163.IBA1.Macrophages.x),
                                                                               B.cells = mean(B.cells.x))



data_for_hm_mean <- pivot_longer(data_for_hm_mean, cols = 3:12, names_to = "Neighbor")


#select the cluster with endothelial cells

data_for_hm_endothalia <- data_for_hm_mean[which(data_for_hm_mean$cluster == "Endothelia"),]

data_for_hm_endothalia$Neighbor <- as.factor(data_for_hm_endothalia$Neighbor)
data_for_hm_endothalia$value <- as.numeric(data_for_hm_endothalia$value)
ggplot(data_for_hm_endothalia, aes(x=Neighbor, y=value)) + geom_boxplot() + 
  geom_jitter() + 
 coord_flip()


kruskal.test(value ~ Neighbor, data = data_for_hm_endothalia)
pairwise.wilcox.test(data_for_hm_endothalia$value, data_for_hm_endothalia$Neighbor,
                     p.adjust.method = "BH")

#IBA1+CD163+ macrophages were neighboring endothelial cells significantly more than B-cells, CD11c APCs, 
#IBA1+CD163+CD11c+ macrophages, IBA1+CD11c+ macrophages, CD163+ macrophages, CD4+ T-cells
#not however cds, Tregs or IBA1 macrophages










