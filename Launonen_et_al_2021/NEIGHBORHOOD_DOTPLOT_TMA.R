##########################################################################################################################
#################################### FOLD CHANGE DOTPLOT #################################################################
##########################################################################################################################

#For Supplementary Figure 5

library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(tidyr)

#all_fractions

setwd("interaction_frequencies")
temp= list.files(pattern=".csv")
myfiles = lapply(temp, read.csv)

for (i in c(1:112)){
  
  myfiles[[i]]$core <- paste(temp)[i]
}



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


all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")
all_celltypes$GlobalCellType <- as.character(all_celltypes$GlobalCellType)
all_celltypes$Subtype <- as.character(all_celltypes$Subtype)

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
library(tidyr)

metacluster_percentages <- metacluster_percentages %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

metacluster_percentages <- as.data.frame(metacluster_percentages)
metacluster_percentages[is.na(metacluster_percentages)] <- 0

#the rows of all_fractions are the centering cell types

#the columns are the neighboring cell types

#have to normalize by the neighboring cell type per core


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


#p-values

p_values <- all_fractions_merged_abundance %>% group_by(cluster)%>%
  summarise_each(funs(wilcox.test(.[Category == "BRCA1/2 mutated"], .[Category == "HR"])$p.value), vars = Apoptotic.x:Hyperfunctional.epithelial.x)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(all_fractions_merged_abundance)[c(3:30)]



#median normalized neighborhood fractions

mean_expression_new <- all_fractions_merged_abundance %>% 
  group_by(cluster, Category) %>%
  summarise(Apoptotic = median(Apoptotic.x, na.rm = TRUE),
            B.cells = median(B.cells.x, na.rm = TRUE),
            CD11c.APC = median(CD11c.APC.x, na.rm = TRUE),
            CD11c.CD163.IBA1.Macrophages = median(CD11c.CD163.IBA1.Macrophages.x, na.rm = TRUE),
            CD163.IBA1.Macrophages = median(CD163.IBA1.Macrophages.x, na.rm = TRUE),
            CD163.Macrophages= median(CD163.Macrophages.x, na.rm = TRUE),
            CD4 = median(CD4.x, na.rm = TRUE),
            CD8 = median(CD8.x, na.rm = TRUE),
            EMT = median(EMT.x, na.rm = TRUE),
            Endothelia = median(Endothelia.x, na.rm = TRUE),
            Epithelial = median(Epithelial.x, na.rm = TRUE),
            FOXP3.CD4.Tregs = median(FOXP3.CD4.Tregs.x, na.rm = TRUE),
            High.PDL1 = median(High.PDL1.x, na.rm = TRUE),
            High.Vimentin = median(High.Vimentin.x, na.rm = TRUE),
            High.proliferative_Stroma = median(High.proliferative_Stroma.x, na.rm = TRUE),
            Low.vimentin = median(Low.vimentin.x, na.rm = TRUE),
            Low_eccentricity_medium_vimentin = median(Low_eccentricity_medium_vimentin.x, na.rm = TRUE),
            Mesenchymal = median(Mesenchymal.x, na.rm = TRUE),
            Non.proliferative_Stroma = median(Non.proliferative_Stroma.x, na.rm = TRUE),
            Proliferating.EMT = median(Proliferating.EMT.x, na.rm = TRUE),
            Proliferating.epithelial = median(Proliferating.epithelial.x, na.rm = TRUE),
            Proliferative_Stroma = median(Proliferative_Stroma.x, na.rm = TRUE),
            CD11c.IBA1.Macrophages = median(CD11c.IBA1.Macrophages.x, na.rm = TRUE),
            IBA1.Macrophages = median(IBA1.Macrophages.x, na.rm = TRUE),
            High.P21 = median(High.P21.x, na.rm = TRUE),
            Hyperfunctional.epithelial = median(Hyperfunctional.epithelial.x, na.rm = TRUE))



median_expression_long_new <- pivot_longer(mean_expression_new, cols = 3:28, names_to = "Neighbour", values_to = "proportion")

median_expression_wide_new <- pivot_wider(median_expression_long_new, names_from = "Category", values_from = "proportion")

#fold change

median_expression_wide_new$fold_change <- median_expression_wide_new$`BRCA1/2 mutated`/median_expression_wide_new$HR
median_expression_wide_new$log_fold_change <- log2(median_expression_wide_new$fold_change)



fold_change_data_new <- median_expression_wide_new[, c("cluster", "Neighbour", "log_fold_change")]


log_fold_change_data_wide_new <- pivot_wider(fold_change_data_new, names_from = "cluster", values_from = "log_fold_change")

log_fold_change_data_wide_new <- as.data.frame(log_fold_change_data_wide_new)
rownames(log_fold_change_data_wide_new) <- log_fold_change_data_wide_new[, 1]

log_fold_change_data_wide_new <- log_fold_change_data_wide_new[, -c(1)]



log_fold_change_data_wide_new <- as.data.table(t(log_fold_change_data_wide_new), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide_new <- as.data.frame(log_fold_change_data_wide_new)
rownames(log_fold_change_data_wide_new) <- log_fold_change_data_wide_new[, 1]
log_fold_change_data_wide_new <- log_fold_change_data_wide_new[, -c(1)]

log_fd_heatmap_new <- log_fold_change_data_wide_new


fod_fd_hmap_sig <- log_fd_heatmap_new


fod_fd_hmap_sig <- as.data.frame(fod_fd_hmap_sig)

setDT(fod_fd_hmap_sig, keep.rownames = "cluster")[]
fod_fd_hmap_sig_long <- pivot_longer(fod_fd_hmap_sig, cols = 2:27, names_to = "Foldchange")
colnames(fod_fd_hmap_sig_long) <- c("cluster", "Neighbour", "log_foldchange")

setDT(p_values, keep.rownames = "cluster")[]
p_values_long <- pivot_longer(p_values, cols = 2:29, names_to = "pvalue")
colnames(p_values_long) <- c("cluster", "Neighbour", "pvalue")


p_values_long <- p_values_long[-which(p_values_long$cluster == "Other" | p_values_long$cluster == "Negative" | p_values_long$Neighbour == "Other" | p_values_long$Neighbour == "Negative"),]

p_values_long_fdr <- p.adjust(p_values_long$pvalue, method="BH")

p_values_long_fdr <- cbind(p_values_long[, c(1, 2)], p_values_long_fdr)
colnames(p_values_long_fdr)[3] <- "pvalue"

fod_fd_hmap_sig_long <- fod_fd_hmap_sig_long[-which(fod_fd_hmap_sig_long$cluster == "Other" | fod_fd_hmap_sig_long$cluster == "Negative"),]


poly.results <- cbind(fod_fd_hmap_sig_long, p_values_long_fdr)
poly.results <- poly.results[, c(1, 2, 3, 6)]


truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-5, range.lims=c(-4, 4), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 

library(ggplot2)
theme <- theme(panel.border = element_rect(colour = "black", size=1.7, fill=NA), 
               panel.background = element_blank(), plot.title=element_text(hjust=0.5)) 


y_order <- c("Proliferating epithelial","Epithelial","EMT","Proliferating EMT", "Mesenchymal","Apoptotic","Functional epithelial", "CD8+ T-cells", "CD4+ effector T-cells","FOXP3+CD4+ T-regs", "IBA1+CD163+ macrophages", "IBA1+CD163+CD11c+ macrophages","IBA1+ macrophages", "IBA1+CD11c+ macrophages","CD11c+ APC", "CD163+ macrophages","B-cells",  "Non-proliferative stroma","Low-vimentin","High-vimentin","Endothelia", "Proliferative stroma","High-proliferative stroma", "Functional stroma", "Low-eccentricity","High-P21")


x_order <- c("Proliferating epithelial","Epithelial","EMT","Proliferating EMT", "Mesenchymal","Apoptotic","Functional epithelial", "CD8+ T-cells", "CD4+ effector T-cells","FOXP3+CD4+ T-regs", "IBA1+CD163+ macrophages", "IBA1+CD163+CD11c+ macrophages","IBA1+ macrophages", "IBA1+CD11c+ macrophages","CD11c+ APC", "CD163+ macrophages","B-cells",  "Non-proliferative stroma","Low-vimentin","High-vimentin","Endothelia", "Proliferative stroma","High-proliferative stroma", "Functional stroma", "Low-eccentricity","High-P21")


to.plot[which(to.plot$cluster == "Hyperfunctional epithelial"), "cluster"] <- "Functional epithelial"
to.plot[which(to.plot$cluster == "CD4"), "cluster"] <- "CD4+ effector T-cells"
to.plot[which(to.plot$cluster == "CD8"), "cluster"] <- "CD8+ T-cells"
to.plot[which(to.plot$cluster == "FOXP3+CD4+Tregs"), "cluster"] <- "FOXP3+CD4+ T-regs"
to.plot[which(to.plot$cluster == "CD11c+IBA1+Macrophages"), "cluster"] <- "IBA1+CD11c+ macrophages"
to.plot[which(to.plot$cluster == "CD163+IBA1+Macrophages"), "cluster"] <- "IBA1+CD163+ macrophages"
to.plot[which(to.plot$cluster == "IBA1+Macrophages"), "cluster"] <- "IBA1+ macrophages"
to.plot[which(to.plot$cluster == "CD163+Macrophages"), "cluster"] <- "CD163+ macrophages"
to.plot[which(to.plot$cluster == "CD11c+CD163+IBA1+Macrophages"), "cluster"] <- "IBA1+CD163+CD11c+ macrophages"
to.plot[which(to.plot$cluster == "High-PDL1"), "cluster"] <- "Functional stroma"
to.plot[which(to.plot$cluster == "High-proliferative_Stroma"), "cluster"] <- "High-proliferative stroma"
to.plot[which(to.plot$cluster == "High-Vimentin"), "cluster"] <- "High-vimentin"
to.plot[which(to.plot$cluster == "Low_eccentricity_medium_vimentin"), "cluster"] <- "Low-eccentricity"
to.plot[which(to.plot$cluster == "Non-proliferative_Stroma"), "cluster"] <- "Non-proliferative stroma"
to.plot[which(to.plot$cluster == "Proliferative_Stroma"), "cluster"] <- "Proliferative stroma"
to.plot[which(to.plot$cluster == "CD11c+APC"), "cluster"] <- "CD11c+ APC"


to.plot[which(to.plot$Neighbour == "Hyperfunctional.epithelial"), "Neighbour"] <- "Functional epithelial"
to.plot[which(to.plot$Neighbour == "CD4"), "Neighbour"] <- "CD4+ effector T-cells"
to.plot[which(to.plot$Neighbour == "CD8"), "Neighbour"] <- "CD8+ T-cells"
to.plot[which(to.plot$Neighbour == "FOXP3.CD4.Tregs"), "Neighbour"] <- "FOXP3+CD4+ T-regs"
to.plot[which(to.plot$Neighbour == "CD11c.IBA1.Macrophages"), "Neighbour"] <- "IBA1+CD11c+ macrophages"
to.plot[which(to.plot$Neighbour == "CD163.IBA1.Macrophages"), "Neighbour"] <- "IBA1+CD163+ macrophages"
to.plot[which(to.plot$Neighbour == "IBA1.Macrophages"), "Neighbour"] <- "IBA1+ macrophages"
to.plot[which(to.plot$Neighbour == "CD163.Macrophages"), "Neighbour"] <- "CD163+ macrophages"
to.plot[which(to.plot$Neighbour == "CD11c.CD163.IBA1.Macrophages"), "Neighbour"] <- "IBA1+CD163+CD11c+ macrophages"
to.plot[which(to.plot$Neighbour == "High.PDL1"), "Neighbour"] <- "Functional stroma"
to.plot[which(to.plot$Neighbour == "High.proliferative_Stroma"), "Neighbour"] <- "High-proliferative stroma"
to.plot[which(to.plot$Neighbour == "High.Vimentin"), "Neighbour"] <- "High-vimentin"
to.plot[which(to.plot$Neighbour == "Low_eccentricity_medium_vimentin"), "Neighbour"] <- "Low-eccentricity"
to.plot[which(to.plot$Neighbour == "Non.proliferative_Stroma"), "Neighbour"] <- "Non-proliferative stroma"
to.plot[which(to.plot$Neighbour == "Proliferative_Stroma"), "Neighbour"] <- "Proliferative stroma"
to.plot[which(to.plot$Neighbour == "High.P21"), "Neighbour"] <- "High-P21"
to.plot[which(to.plot$Neighbour == "Low.vimentin"), "Neighbour"] <- "Low-vimentin"
to.plot[which(to.plot$Neighbour == "CD11c.APC"), "Neighbour"] <- "CD11c+ APC"
to.plot[which(to.plot$Neighbour == "B.cells"), "Neighbour"] <- "B-cells"
to.plot[which(to.plot$Neighbour == "Proliferating.EMT"), "Neighbour"] <- "Proliferating EMT"
to.plot[which(to.plot$Neighbour == "Proliferating.epithelial"), "Neighbour"] <- "Proliferating epithelial"



p_neighbours_HRstatus <- ggplot(to.plot,  aes(x=factor(Neighbour, levels=x_order), y=factor(cluster, levels=rev(y_order)))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="darkorange3", mid = "gray90", breaks=seq(-4,0.1,4), limits=c(-4, 4)) +
  scale_size_area("p-values", trans="log10",max_size = 3, breaks=c(1e-5, 1e-1, 0.05), limits=c(1e-5,0.05)) +
  theme(axis.text=element_text(size=rel(0.7)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"), panel.border = element_rect(colour = "black", size=1.5, fill=NA), aspect.ratio=1) + xlab("Neighbouring celltype") + ylab("Reference celltype") 




print(p_neighbours_HRstatus)
