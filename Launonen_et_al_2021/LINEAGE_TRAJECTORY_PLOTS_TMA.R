##############################################################################################################################
################################ LINEAGE TRAJECTORY PLOTS ####################################################################
##############################################################################################################################

#For Figure 3 and Supplementary Figure 1

#install.packages("devtools")
#devtools::install_git("https://git.embl.de/velten/STEMNET/", build_vignettes=TRUE)


require(STEMNET)
require(ggplot2)
require(cluster)
require(pheatmap)
require(gridExtra)
library(dplyr)
library(RColorBrewer)


set.seed(1234)


#For STEMNET, create a vector annotating each cell as part of one of the most mature populations, 
#or NA. Then, call stemnet.


#create a vector of annotations and a matrix of gene expressions 

all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")
all_cells_tumor <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]



stemnet_pop <- as.vector(all_cells_tumor$GlobalCellType)

#create few random NAs to train
ind <- floor(runif(30, min=0, max=80000))
stemnet_pop[ind] <- NA


geneExpression <- as.matrix(scale(all_cells_tumor[, c("pSTAT1", "Ki67", "PDL1", "P21", "Ecadherin", "vimentin", "cCasp3")]))

stemnet_result <- runSTEMNET(geneExpression, stemnet_pop)

#The output of STEMNET is a matrix of class probabilities:

CAPframe <- plot(stemnet_result)


mycolors <- rev(c("#386CB0", "#BEAED4", "#FDC086", "#BF5B17", "#7FC97F","#F0027F",  "#FFFF99"))

lineage_plot <- ggplot(CAPframe$PlotData, aes(x=CAPframe$PlotData$x, y=CAPframe$PlotData$y)) 
lineage_plot + geom_point(aes(color=all_cells_tumor$GlobalCellType), size=0.005) + 
  scale_color_manual(values=mycolors) + 
  theme_minimal() + theme(aspect.ratio = 1) + xlab("x") + ylab("y")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20))

lineage_plot + geom_point(aes(color=all_cells_tumor$Ki67), size=0.005) + 
  theme_minimal() + theme(aspect.ratio = 1)+ theme(aspect.ratio = 1) + xlab("x") + ylab("y")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20))+ 
  scale_color_continuous(name = "Ki67")

lineage_plot + geom_point(aes(color=all_cells_tumor$vimentin), size=0.005) + 
  theme_minimal() + theme(aspect.ratio = 1)+ theme(aspect.ratio = 1) + xlab("x") + ylab("y")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20))+ 
  scale_color_continuous(name = "Vimentin")

lineage_plot + geom_point(aes(color=all_cells_tumor$PDL1), size=0.005) + 
  theme_minimal() + theme(aspect.ratio = 1)+ theme(aspect.ratio = 1) + xlab("x") + ylab("y")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size=20))+ 
  scale_color_continuous(name = "PDL1")





#for immune cells

all_cells_immune <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Immune"),]



stemnet_pop_immune <- as.vector(all_cells_immune$GlobalCellType)

stemnet_pop_immune[which(stemnet_pop_immune == "IBA1+ Macrophages" | stemnet_pop_immune == "CD163+ Macrophages")]<-NA


geneExpression_immune <- as.matrix(scale(all_cells_immune[, c("FOXP3","pSTAT1", "Ki67", "PDL1", "P21","PD1","CD3d","CD4","CD8a","IBA1", "CD11c","CD163","CD20","cCasp3")]))


stemnet_result_immune <- runSTEMNET(geneExpression_immune, stemnet_pop_immune)

#The output of STEMNET is a matrix of class probabilities:

CAPframe_immune <- plot(stemnet_result_immune)

mycolors <- colorRampPalette(brewer.pal(10, "Paired"))(10)

lineage_plot_immune <- ggplot(CAPframe_immune$PlotData, aes(x=CAPframe_immune$PlotData$x, y=CAPframe_immune$PlotData$y)) 
lineage_plot_immune + geom_point(aes(color=all_cells_immune$GlobalCellType), size=0.005) + 
  scale_color_manual(values=mycolors) + 
  theme_minimal() + theme(aspect.ratio = 1, legend.position="none") + xlab("x") + ylab("y")+
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10))
#width=3



