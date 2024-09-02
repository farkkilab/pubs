#delaunay triangulation for new set of IDS images

library(spatstat)
library(igraph)
library(dplyr)

obs_all <- read.csv('/media/oncosys/Expansion/cellcycle/cellcycle_phenotypes_RCNs.csv')

#merge RCNs that are of same category

obs_all[which(obs_all$spatial_kmeans == '5'),'spatial_kmeans'] <- '10'
obs_all[which(obs_all$spatial_kmeans == '12'),'spatial_kmeans'] <- '2'
obs_all[which(obs_all$spatial_kmeans == '13'),'spatial_kmeans'] <- '17'

#loop through RCNs
for (i in c('10', '2', '6', '9', '17')){
  for (j in unique(obs_all$imageid)){
obs <- obs_all[which(obs_all$spatial_kmeans == i),]
#sample by sample

obs <- obs[which(obs$imageid == j), c("X_centroid", 'Y_centroid', 'CellID')]

X <- obs
X <- ppp(X$X_centroid, X$Y_centroid,   window=owin(xrange=c(0, max(obs$X_centroid)),yrange=c(0, max(obs$Y_centroid))),marks=X$CellID)
Dnet <- delaunayNetwork(X)
len <- lengths_psp(as.psp(Dnet))
Net <- thinNetwork(Dnet, retainedges = (len <= 45))

#pdf(paste0(i,j,".pdf"))
#plot(Net)
#dev.off()

graph <- graph_from_data_frame(data.frame(Net$from, Net$to))
xpos <- Net$lines$ends$x0
ypos <- Net$lines$ends$y0
identity <- Net$from
to_be_merged <- bind_cols(xpos, ypos, identity)
colnames(to_be_merged)[c(1:3)] <- c('X_centroid', 'Y_centroid', 'identity')
to_be_merged_merged <- merge(to_be_merged, obs, by=c('X_centroid', 'Y_centroid'))
components <- components(graph)
membership <- data.frame(components$membership)
membership$identity <- rownames(membership)
to_be_merged_merged <- merge(to_be_merged_merged, membership, by='identity')
write.csv(to_be_merged_merged, paste0( i, j,".csv"))

}}


