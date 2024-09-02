
#Giotto Interaction changed feature calculation for the validation IDS dataset

library(devtools)
library(remotes)
library(Giotto)
library(reticulate)

library(ggraph)
library(igraph)
installGiottoEnvironment()
Giotto:::checkGiottoEnvironment()
#remotes::install_github("drieslab/Giotto@suite")

#use_python("C:\\PYTHON~1\\python.exe")


#conda_list()
#installGiottoEnvironment(force_miniconda = TRUE)

obs_all <- read.csv('/media/oncosys/Expansion/cellcycle/cellcycle_phenotypes_TIM3.csv')

my_instructions = createGiottoInstructions(python_path = '/opt/conda/bin/python')


#make a expression matrix for aech sample and save

for (i in unique(obs_all$imageid)){
  #i="s04"
  s21 <- obs_all[which(obs_all$imageid ==i),]
  
  #colnames(data)
  
  raw_exprs = s21[, c(1, 3:5)]
  rownames(raw_exprs) <- raw_exprs$CellID
  raw_exprs <- raw_exprs[, -c(1)]
  
  #keep only markers of interest
  raw_exprs <- raw_exprs[, c('CD163', 'CD11c', 'TIM.3')]
  
  #rm(data)
  #rm(raw_exprs)
  #gc()
  rownames(raw_exprs) <- paste0("cell_", rownames(raw_exprs))
  write.table(raw_exprs, paste0("/media/oncosys/Expansion/cellcycle/Giotto/raw_exprs_Giotto_",i,".csv"))
}

library(dplyr)
#Run Giotto one sample at a time
#fit a loop

for (i in unique(obs_all$imageid)){
  #i="S050_iAdn"
  s21 <- obs_all[which(obs_all$imageid == i),]
  
  spatial_locs = s21[, c(1,6,7)]
  rownames(spatial_locs) = spatial_locs$CellID
  spatial_locs = spatial_locs[,-c(1)]
  colnames(spatial_locs) <- c("X.X", "Y.Y")
  spatial_locs$CellID <- rownames(spatial_locs)
  library(data.table)
  spatial_locs <- as.data.table(spatial_locs)
  
  #expression <- readExprMatrix(paste0("/media/oncosys/Expansion/cellcycle/Giotto/raw_exprs_Giotto_",i,".csv"), transpose = T)
  expression <- readExprMatrix(paste0("/home/oncosys/Desktop/cellcycle_testing/raw_exprs_Giotto_",i,".csv"), transpose = T)
  
  my_giotto_object = createGiottoObject(raw_exprs  = expression,
                                        custom_expr = expression,
                                        spatial_locs = spatial_locs[,-c(3)], 
                                        instructions = my_instructions)
  
  
  
  metadata <- s21[, c("phenotype", "CellID",'spatial_kmeans', "X_centroid", "Y_centroid")]
  #colnames(metadata)[4] <- "ID"
  metadata$CellID <- paste0("cell_", metadata$CellID)
  
  codex_test<-addCellMetadata(my_giotto_object, new_metadata = metadata,
  )
  
  #subset myeloid RCNs
  subset_cell_IDs = metadata[which(metadata$spatial_kmeans %in% c('13', '17', '9', '6', '16', '2', '12')),]$CellID
  #
  codex_test = subsetGiotto(codex_test, cell_ids = subset_cell_IDs)
  
  #spatial interaction changes functional markers
  
  #need delaunay triangulation
  
  codex_test = createSpatialNetwork(gobject = codex_test, maximum_distance_delaunay = 45)
  #plotStatDelaunayNetwork(gobject = codex_test, maximum_distance = 50, save_plot = F)
  
  #miten saada tasta network niinkuin ennen
  delaunay = codex_test@spatial_network$Delaunay_network$networkDT
  #nyt on from ja to
  graph <- graph_from_data_frame(data.frame(delaunay$from, delaunay$to))
  xpos <- delaunay$sdimx_begin
  ypos <- delaunay$sdimy_begin
  identity <- delaunay$from
  to_be_merged <- bind_cols(xpos, ypos, identity)
  colnames(to_be_merged)[c(1:3)] <- c('X_centroid', 'Y_centroid', 'CellID')
  #miks seuraavasta tulee 0 obs????
  #testaa CellID:lla
  s21$CellID <- paste0('cell_', s21$CellID)
  to_be_merged_merged <- merge(to_be_merged, s21, by=c('CellID'))
  #to_be_merged_merged <- merge(to_be_merged, s21, by=c('X_centroid', 'Y_centroid'))
  
  
  components <- components(graph)
  membership <- data.frame(components$membership)
  membership$CellID <- rownames(membership)
  to_be_merged_merged <- merge(to_be_merged_merged, membership, by='CellID')
  #write.csv(to_be_merged_merged, paste0( '/media/oncosys/Expansion/cellcycle/Giotto/',i,"_spatial_kmeans_13.csv"))
  write.csv(to_be_merged_merged, paste0( '/home/oncosys/Desktop/cellcycle_testing/',i,"_spatial_kmeans_all_mac.csv"))
  
  
  
  library(future)
  #findInteractionChangedFeats
  ICFsForesHighGenes =  findInteractionChangedGenes(gobject = codex_test,
                                                    spatial_network_name = 'Delaunay_network',
                                                    cluster_column = "phenotype",
                                                    diff_test = 'permutation',
                                                    expression_values = "custom",
                                                    adjust_method = 'fdr',
                                                    nr_permutations = 1000,
                                                    do_parallel = T)
  
  #save the interaction changed features
  write.csv(ICFsForesHighGenes[[1]], paste0("/home/oncosys/Desktop/cellcycle_testing/Giotto_phenotypes_",i,"_all_mac.csv"))
  #write.csv(ICFsForesHighGenes[[1]], paste0("/media/oncosys/Expansion/cellcycle/Giotto/Giotto_phenotypes_",i,"_all_mac.csv"))
}

