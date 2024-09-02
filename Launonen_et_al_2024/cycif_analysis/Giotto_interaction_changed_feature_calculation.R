#Giotto Interaction changed feature calculation for S figure 3c-d

library(devtools)
library(remotes)
library(Giotto)
library(reticulate)

#install.packages('devtools')
#remotes::install_github("RubD/Giotto@cless") 
#write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
#Sys.which("make")
#install.packages('ClusterR')
#install.packages('ggraph')
#install.packages('reticulate')

installGiottoEnvironment()
Giotto:::checkGiottoEnvironment()
#remotes::install_github("drieslab/Giotto@suite")

use_python("C:\\PYTHON~1\\python.exe")


conda_list()
installGiottoEnvironment(force_miniconda = TRUE)

obs_all <- read.csv("P:/h345/afarkkilab/Projects/Sciset/data_20231125.csv", header=T)

my_instructions = createGiottoInstructions(python_path = 'C:\\PYTHON~1\\python.exe')


#make a expression matrix for aech sample and save

for (i in c("s01","s02", "s03","s04","s05","s06",'s07', "s08", "s09", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20", "s21", "s22", "s23")){
#i="s04"
s21 <- obs_all[which(obs_all$imageid ==i),]

#colnames(data)

raw_exprs = s21[, c(44, 3:34)]
rownames(raw_exprs) <- raw_exprs$ID
raw_exprs <- raw_exprs[, -c(1)]

#keep only markers of interest
raw_exprs <- raw_exprs[, c("Ki67", "MHCI", "MHCII", "TIM3", "pSTAT1", "TAZ", "HE4", "SNAT1", "Annexin", "FOXOA3")]

#rm(data)
#rm(raw_exprs)
#gc()
rownames(raw_exprs) <- paste0("cell_", rownames(raw_exprs))
write.table(raw_exprs, paste0("P:/h345/afarkkilab/Projects/Sciset/giotto/giottoICF/merged_myeloids_all2/raw_exprs_Giotto_",i,".csv"))
}

#make different categories: merge myeloids, merge stromal metaclusters, divide CD8+T-cells into subtypes

obs_all$Global_celltype_tumor_stroma_merged <- obs_all$GlobalCellType2
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "EMT" |obs_all$Global_celltype_tumor_stroma_merged == "Epithelial" |obs_all$Global_celltype_tumor_stroma_merged == "Proliferating.EMT" |obs_all$Global_celltype_tumor_stroma_merged == "Proliferating.epithelial" ), "Global_celltype_tumor_stroma_merged"] <- "Tumor"
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "Fibroblast" |obs_all$Global_celltype_tumor_stroma_merged == "SMA.CD31.positive.cell" |obs_all$Global_celltype_tumor_stroma_merged == "Myofibroblast" |obs_all$Global_celltype_tumor_stroma_merged == "Desmin.positive.cell" | obs_all$Global_celltype_tumor_stroma_merged == "SMA.Desmin.positive.cell" | obs_all$Global_celltype_tumor_stroma_merged == "Other" ), "Global_celltype_tumor_stroma_merged"] <- "Stroma"
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "IBA1.CD163.Macrophages" |obs_all$Global_celltype_tumor_stroma_merged == "CD163.Macrophages" |obs_all$Global_celltype_tumor_stroma_merged == "IBA1.CD11c.Macrophages" ), "Global_celltype_tumor_stroma_merged"] <- "Macrophages"

obs_all$Global_celltype_tumor_stroma_merged_CD8_cat <- obs_all$Global_celltype_tumor_stroma_merged
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "CD8.T.cells" & obs_all$pSTAT1 > 0.5), "Global_celltype_tumor_stroma_merged_CD8_cat"] <- "pSTAT1_CD8.T.cell"
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "CD8.T.cells" & obs_all$TIM3 > 0.5), "Global_celltype_tumor_stroma_merged_CD8_cat"] <- "TIM3_CD8.T.cell"
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "CD8.T.cells" & obs_all$pSTAT1 > 0.5 & obs_all$TIM3 > 0.5), "Global_celltype_tumor_stroma_merged_CD8_cat"] <- "TIM3_pSTAT1_CD8.T.cell"


obs_all$GlobalCellType2_CD8_cat <- obs_all$GlobalCellType2
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "CD8.T.cells" & obs_all$pSTAT1status== "pSTAT1pos"), "GlobalCellType2_CD8_cat"] <- "pSTAT1_CD8.T.cell"
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "CD8.T.cells" & obs_all$TIM3status== "TIM3pos"), "GlobalCellType2_CD8_cat"] <- "TIM3_CD8.T.cell"
obs_all[which(obs_all$Global_celltype_tumor_stroma_merged == "CD8.T.cells" & obs_all$pSTAT1status== "pSTAT1pos" & obs_all$TIM3status== "TIM3pos"), "GlobalCellType2_CD8_cat"] <- "TIM3_pSTAT1_CD8.T.cell"


#Run Giotto one sample at a time

unique(obs_all$Global_celltype_tumor_stroma_merged_CD8_cat)

#fit a loop
 
for (i in c("s01","s02", "s03","s04","s05","s06",'s07', "s08", "s09", "s11", "s12", "s13", "s14", "s15", "s16", "s17", "s18", "s19", "s20", "s21", "s22", "s23")){

s21 <- obs_all[which(obs_all$imageid == i),]

spatial_locs = s21[, c(44,35, 36)]
rownames(spatial_locs) = spatial_locs$ID
spatial_locs = spatial_locs[,-c(1)]
colnames(spatial_locs) <- c("X.X", "Y.Y")
spatial_locs$ID <- rownames(spatial_locs)
library(data.table)
spatial_locs <- as.data.table(spatial_locs)

expression <- readExprMatrix(paste0("P:/h345/afarkkilab/Projects/Sciset/giotto/giottoICF/merged_myeloids_all2/raw_exprs_Giotto_",i,".csv"), transpose = T)

my_giotto_object = createGiottoObject(expression = expression,
                                      spatial_locs = spatial_locs[,-c(3)], 
                                      instructions = my_instructions)



metadata <- s21[, c("GlobalCellType2", "Global_celltype_tumor_stroma_merged_CD8_cat", "neighbordood_cluster2", "ID", "X.Position", "Y.Position")]
colnames(metadata)[4] <- "ID"
metadata$ID <- paste0("cell_", metadata$ID)

codex_test<-addCellMetadata(my_giotto_object, new_metadata = metadata,
                           )

codex_test = createSpatialNetwork(gobject = codex_test, maximum_distance_delaunay = 45)

library(future)

ICFsForesHighGenes =  findInteractionChangedFeats(gobject = codex_test,
                                                  selected_feats = c("Ki67", "MHCI", "MHCII", "TIM3", "pSTAT1", "TAZ", "HE4", "SNAT1", "Annexin", "FOXOA3"),
                                                  spatial_network_name = 'Delaunay_network',
                                                  cluster_column = "Global_celltype_tumor_stroma_merged_CD8_cat",
                                                  diff_test = 'permutation',
                                                  expression_values = "raw",
                                                  adjust_method = 'fdr',
                                                  nr_permutations = 1000,
                                                  do_parallel = T)



write.csv(ICFsForesHighGenes[[1]], paste0("P:/h345/afarkkilab/Projects/Sciset/giotto/giottoICF/Giotto_merged_myeloids2_",i,".csv"))
}


