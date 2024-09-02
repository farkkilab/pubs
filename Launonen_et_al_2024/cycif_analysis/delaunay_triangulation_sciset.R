
#
obs_all <-read.csv('/home/oncosys/Desktop/data_20231125.csv')
install.packages('spatstat')
install.packages('igraph')
install.packages('dplyr')

library(spatstat)
library(igraph)
library(dplyr)

#for each RCN individually


for (i in c(unique(obs_all$imageid))){
  
  obs <- obs_all[which(obs_all$neighbordood_cluster2 == 'Immune'),]
  #obs <- obs_all[which(obs_all$neighbordood_cluster2 == 'IBA1.CD163.Macrophages'|obs_all$neighbordood_cluster2 == 'Macrophages'|obs_all$neighbordood_cluster2 == 'CD11c.myeloid'),]
  
  obs <- obs[which(obs$imageid == i), c("X.Position", 'Y.Position', 'ID')]
  #obs <- obs[, c("X.Position", 'Y.Position', 'CellID')]
  X <- obs

  gc()
  X <- ppp(X$X.Position, X$Y.Position,   window=owin(xrange=c(0, max(obs$X.Position)),yrange=c(0, max(obs$Y.Position))),marks=X$ID)
  
  Dnet <- delaunayNetwork(X)
  rm(X)
  gc()
  len <- lengths_psp(as.psp(Dnet))
  Net <- thinNetwork(Dnet, retainedges = (len <= 45))
  rm(len)
  rm(Dnet)
  gc()

  library(igraph)

  
  
  graph <- graph_from_data_frame(data.frame(Net$from, Net$to))
  

  
  xpos <- Net$lines$ends$x0
  ypos <- Net$lines$ends$y0
  identity <- Net$from
  
  library(dplyr)
  to_be_merged <- bind_cols(xpos, ypos, identity)
  colnames(to_be_merged)[c(1:3)] <- c('X.Position', 'Y.Position', 'identity')
  
  #merge with obs
  
  to_be_merged_merged <- merge(to_be_merged, obs, by=c('X.Position', 'Y.Position'))
  
  components <- components(graph)
  membership <- data.frame(components$membership)
  membership$identity <- rownames(membership)
  
  to_be_merged_merged <- merge(to_be_merged_merged, membership, by='identity')
  
  write.csv(to_be_merged_merged, paste0("/media/oncosys/T7/sciset/",i,'_Immune_0724.csv'))
}
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

