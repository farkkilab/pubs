#for figure 3a
#T-cell adundancies per neighborhood pre and post chemo

#percentage of all cells per neighborhood
#percentage of all CD8+T-cells in here
library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(stats)

data_expr <- read.csv("E:/Sciset/data_030123.csv")

all_clusters_sorted <- data_expr %>%
  group_by(Sample_code, neighbordood_cluster2, GlobalCellType2) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

all_clusters_sorted$freq <- all_clusters_sorted$freq*100
all_clusters_sorted <- all_clusters_sorted[, -4]
all_clusters_sorted_wide <- pivot_wider(all_clusters_sorted, names_from="Sample_code", values_from="freq")
all_clusters_sorted_wide <- as.data.frame(all_clusters_sorted_wide)
all_clusters_sorted_wide <- all_clusters_sorted_wide[which(all_clusters_sorted_wide$GlobalCellType2 == "CD8.T.cells"),]
rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide$neighbordood_cluster2
all_clusters_sorted_wide <- all_clusters_sorted_wide[,-c(1, 2)]
all_clusters_sorted_wide <- as.data.table(t(all_clusters_sorted_wide), keep.colnames = T, keep.rownames = T)

rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide$rn
colnames(all_clusters_sorted_wide)[1] <- "Sample"
all_clusters_sorted_wide[is.na(all_clusters_sorted_wide)] <- 0
all_clusters_sorted_wide <- as.data.frame(all_clusters_sorted_wide)
rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide[, 1]
all_clusters_sorted_wide <- all_clusters_sorted_wide[, -1]


info <- read.csv("L:/ltdk_farkkila/Data/Sciset/clinical_data.csv")


all_clusters_sorted_wide <- all_clusters_sorted_wide[info$Sample_code,]


#heatmap

rownames(info) <- info$Sample_code
ann <- info[, c("HRD_status", "Stage", "CUD_Treatment.strategy")]
#rownames(ann) <- ann[, 1]
#ann <- ann[,-c(1)]

ann <- as.data.frame(ann)

colnames(ann)[3] <- "treatment"


#mean per pre and per post


all_clusters_sorted_wide <- cbind(all_clusters_sorted_wide, info)


mean_cells <- all_clusters_sorted_wide %>% group_by(Stage) %>% summarize(CD11c.myeloid = mean(CD11c.myeloid),
                                                                         CD8_CD4_T.cells = mean(CD8_CD4_T.cells),
                                                                         Desmin.positive = mean(Desmin.positive),
                                                                         EMT = mean(EMT),
                                                                         Epithelial = mean(Epithelial),
                                                                         epithelial_and_proliferating_epithelial = mean(epithelial_and_proliferating_epithelial),
                                                                         Epithelial_EMT = mean(Epithelial_EMT),
                                                                         Fibroblast = mean(Fibroblast),
                                                                         IBA1.CD163.Macrophages = mean(IBA1.CD163.Macrophages),
                                                                         Immune = mean(Immune),
                                                                         Macrophages = mean(Macrophages),
                                                                         Myofibroblast = mean(Myofibroblast),
                                                                         Proliferating.EMT = mean(Proliferating.EMT),
                                                                         Proliferating.epithelial = mean(Proliferating.epithelial),
                                                                         SMA.CD31.positive = mean(SMA.CD31.positive),
                                                                         SMA.Desmin.myofibroblast = mean(SMA.Desmin.myofibroblast),
                                                                         stroma = mean(stroma),
                                                                         `tumor-stroma-interface` = mean(`tumor-stroma-interface`))


mean_cells <- t(mean_cells)

colnames(mean_cells) <- c("interval", "primary")

mean_cells <- mean_cells[-1,]

mean_cells <- as.data.frame(mean_cells)
mean_cells$interval <- as.numeric(mean_cells$interval)
mean_cells$primary <- as.numeric(mean_cells$primary)
mean_cells <- mean_cells[, c(2, 1)]

mean_cells_log <- mean_cells
mean_cells$primary <- log(mean_cells$primary)
mean_cells$interval <- log(mean_cells$interval)

Heatmap(mean_cells, cluster_columns = FALSE, width = 1, height=5, border="white",
        rect_gp = gpar(col = "gray48", lwd = 2))



#also statistics

#paired wilcox test
#FDR correction

#all clusters sorted wide has the data

#keep paired 
all_clusters_sorted <- merge(all_clusters_sorted, info[, c("Sample_code", "Patient.x", "Stage", "Paired")])
freqs <- all_clusters_sorted[which(all_clusters_sorted$Paired == T),]

freqs <- freqs[which(freqs$GlobalCellType2 == "CD8.T.cells"),]
freqs <- pivot_wider(freqs[,-c(1,3)], values_from = freq, names_from = neighbordood_cluster2)
freqs[is.na(freqs)] <- 0
freqs_1 <- pivot_wider(freqs, names_from = Stage, values_from = unique(all_clusters_sorted$neighbordood_cluster2) )
freqs_1 <- as.data.frame(freqs_1)

for(i in unique(all_clusters_sorted$neighbordood_cluster2)){
  int <- freqs_1[,paste0(i, "_interval")]
  prim <- freqs_1[, paste0(i, "_primary")]
  print(i)
  print(wilcox.test(int, prim, p.adjust.method = 'none',paired = T, na.action(na.omit(int, prim))))
  
}


#as list
#then p.adjust

pvals <- c(0.03125, 0.03125, 0.4375, 0.03125, 0.03125, 0.4375, 0.03125, 0.09375, 0.09375, 0.0625, 0.03125, 0.5625, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125)
p_vals_adjusted <- p.adjust(pvals, method="BH")
neighborhoods <- unique(all_clusters_sorted$neighbordood_cluster2)

pvals_combined <- data.frame(p_vals_adjusted, neighborhoods)

#then add to heatmap 
#1) color codes for neighborhoods
#2) p-values

# colors = c("Proliferating.epithelial"= "#100EAF","epithelial_and_proliferating_epithelial"="#0A087F",
#            "Epithelial"="#29279F","Epithelial_EMT" ="#3F3596","EMT"="#261B83" ,"Proliferating.EMT" ="#110575", 
#            "tumor-stroma-interface" ="#3d5a80", "Desmin.positive"="#383f51","SMA.Desmin.myofibroblast"="#dddbf1",
#            "stroma"="#929487",  "SMA.CD31.positive" ="#D3D4D9", "Myofibroblast"="#d1beb0",  "Fibroblast"="#ab9f9d",
#            "Macrophage"="#247ba0","IBA1.CD163.Macrophages"= "#70c1b3","CD11c.myeloid"="#b2dbbf", 
#            "CD8_CD4_T.cells"="#ff1654","Immune"="#f3ffbd" )


colours<- list("Neighborhood"=c("Proliferating.epithelial"= "#100EAF","epithelial_and_proliferating_epithelial"="#0A087F",
                                "Epithelial"="#29279F","Epithelial_EMT" ="#3F3596","EMT"="#261B83" ,"Proliferating.EMT" ="#110575", 
                                "tumor-stroma-interface" ="#3d5a80", "Desmin.positive"="#383f51","SMA.Desmin.myofibroblast"="#dddbf1",
                                "stroma"="#929487",  "SMA.CD31.positive" ="#D3D4D9", "Myofibroblast"="#d1beb0",  "Fibroblast"="#ab9f9d",
                                "Macrophages"="#247ba0","IBA1.CD163.Macrophages"= "#70c1b3","CD11c.myeloid"="#b2dbbf", 
                                "CD8_CD4_T.cells"="#ff1654","Immune"="#f3ffbd" ))


ann <- data.frame(row.names(mean_cells))
colnames(ann) <- "Neighborhood"
Rowann <- HeatmapAnnotation(df=ann, which="row", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(1, "mm"))

rownames(pvals_combined) <- pvals_combined$neighborhoods
pvals_combined <- pvals_combined[rownames(mean_cells),]

ha = rowAnnotation(foo = anno_simple(pvals_combined$p_vals_adjusted, pch = 1, 
                                     pt_gp = gpar(col = "black"), pt_size = unit(-log(pvals_combined$p_vals_adjusted)*0.3, "snpc"))
)

#lgd3 = Legend(labels = -log(pvals_combined$p_vals_adjusted), legend_gp = gpar(fill = 7:9), title = "legend3")


Heatmap(mean_cells,name="log(CD8+T-cell %)", cluster_columns = FALSE, height = unit(8, "cm") , width = unit(1.3, "cm"), border="white",
        rect_gp = gpar(col = "gray48", lwd = 2), left_annotation = ha, right_annotation = Rowann)



#