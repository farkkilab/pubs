########################################################################################################################################
##################################### PCA PLOTS ########################################################################################
########################################################################################################################################

#For Figure 4 and Supplementary Figure 4

library(FactoMineR)
library(devtools)
library(factoextra)

#get metacluster percentages with HR status info for example from "HEATMAPS_TMA.R"
all_data <- metacluster_percentages

all_data <- all_data[, -c(1)]
colnames(all_data)[31] <- "Type"

pca_data <- all_data


#all celltypes
pca_data <- pca_data[, c(3:28)]

res.pca <- PCA(pca_data, scale.unit = T)

grp <- rownames(res.pca$var$coord)

#colors 12
colors <- c("mediumorchid3" ,"mediumorchid3","mediumorchid3","mediumorchid3" ,"darkgoldenrod1", "darkgoldenrod1", "mediumorchid3" ,"darkgoldenrod1", "darkgoldenrod1", "black", "darkgoldenrod1", "black")


theme <- theme(panel.border = element_rect(colour = "black", size=1.7, fill=NA), 
               panel.background = element_blank(), plot.title=element_text(hjust=0.5)) 


fviz_pca_var(res.pca, alpha.var = "contrib", col.var = grp, palette=colors, ggtheme=theme, select.var = list(contrib = 12), arrowsize=1.5, labelsize=5, col.circle = "White", col.quanti.sup = "white", repel=TRUE)+
  labs(title ="All patients")+ theme(legend.position="none", axis.title = element_text(size = 12),axis.text = element_text(size = 12))


#colors 12
colors <- c("black","mediumorchid3" ,"mediumorchid3","mediumorchid3","mediumorchid3" ,"black","darkgoldenrod1","black","mediumorchid3","black", "black","black")

fviz_pca_var(res.pca, alpha.var = "contrib", col.var = grp, palette=colors, ggtheme=theme, select.var = list(contrib = 12), arrowsize=1.5, labelsize=5, col.circle = "White", col.quanti.sup = "white", repel=TRUE, axes=c(2,3))+
  labs(title ="All patients")+ theme(legend.position="none", axis.title = element_text(size = 12),axis.text = element_text(size = 12))



#eigenvalues <- res.pca$eig
#head(eigenvalues[, 1:2])

grp <- rownames(res.pca$var$coord)
colors <- c("black","mediumorchid3" ,"mediumorchid3","mediumorchid3","mediumorchid3" ,"black","darkgoldenrod1","black","mediumorchid3","mediumorchid3","mediumorchid3","black", "darkgoldenrod1","black","black")

fviz_pca_ind(res.pca, label="none", habillage=all_data$Type, addEllipses = TRUE, axes=c(2,3),  select.var = list(contrib = 15))+
  labs(title ="All patients")+ scale_color_manual(values=c("royalblue", "red3"))+ theme(aspect.ratio=1) + xlim(-10, 10) + ylim(-10, 10)

grp <- rownames(res.pca$var$coord)
colors <- c("black","mediumorchid3" ,"mediumorchid3","mediumorchid3" ,"black","mediumorchid3" ,"darkgoldenrod1", "darkgoldenrod1", "mediumorchid3" ,"darkgoldenrod1", "darkgoldenrod1", "black", "darkgoldenrod1", "black", "darkgoldenrod1")

fviz_pca_ind(res.pca, label="none", habillage=all_data$Type, addEllipses = T,  select.var = list(contrib = 15), pointsize=3)+
  labs(title ="All patients")+ scale_color_manual(values=c("royalblue", "red3"))+ theme(aspect.ratio=1) + xlim(-10, 10) + ylim(-10, 10)


#do separately for BRCAmut and HRwt patients


pca_data <- all_data[which(all_data$Type == "BRCA1/2 mutated"),]
annot <- all_data[which(all_data$Type == "BRCA1/2 mutated"),]
rownames(pca_data) <- pca_data$Patient

pca_data <- pca_data[, c(3:28)]

res.pca <- PCA(pca_data, scale.unit = T)


grp <- rownames(res.pca$var$coord)

colors <- c("black","mediumorchid3","mediumorchid3","mediumorchid3","mediumorchid3","darkgoldenrod1", "darkgoldenrod1","darkgoldenrod1","mediumorchid3","mediumorchid3","darkgoldenrod1","darkgoldenrod1", "black","darkgoldenrod1","black")


colors <- c("mediumorchid3","mediumorchid3","mediumorchid3","mediumorchid3","darkgoldenrod1", "darkgoldenrod1","mediumorchid3","darkgoldenrod1","darkgoldenrod1", "black","darkgoldenrod1","black")




fviz_pca_var(res.pca, alpha.var = "contrib", col.var = grp, palette=colors, ggtheme=theme, arrowsize=1.5, col.circle = "White", col.quanti.sup = "white", select.var = list(contrib = 12)) +
  labs(title ="BRCAmut")


colors <- c("mediumorchid3","mediumorchid3","mediumorchid3","black","black","mediumorchid3","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","black", "black","black")


fviz_pca_var(res.pca, alpha.var = "contrib", labelsize=5,axes = c(2,3),col.var = grp, palette=colors, ggtheme=theme, arrowsize=1.5, col.circle = "White", col.quanti.sup = "white", select.var = list(contrib = 12), repel=TRUE)+ theme(legend.position="none",axis.title = element_text(size = 12),
                                                                                                                                                                                                                                        axis.text = element_text(size = 12)) +
  labs(title ="BRCAmut")

fviz_pca_ind(res.pca, label="none", addEllipses = TRUE,habillage=annot$PFI, axes = c(1,2), ggtheme=theme, pointsize = 3)+ labs(title="BRCA1/2mutated")+ scale_color_manual(values=c("#3399CC", "#003366"))+ theme(aspect.ratio=1) + ylim(c(-10, 10)) + xlim(c(-10, 10))

fviz_pca_ind(res.pca, label="none", axes=c(2,3),addEllipses = TRUE,habillage=annot$PFI, ggtheme=theme, pointsize = 3)+ labs(title="BRCA1/2 mutated")+ scale_color_manual(values=c("#3399CC", "#003366"))+ theme(aspect.ratio=1)

#HRwt

pca_data <- all_data[which(all_data$Type == "HRwt"),]
annot <- all_data[which(all_data$Type == "HRwt"),]
rownames(pca_data) <- pca_data$Patient

pca_data <- pca_data[, c(2:27)]

res.pca <- PCA(pca_data, scale.unit = T)

grp <- rownames(res.pca$var$coord)

colors <- c("black","mediumorchid3","mediumorchid3","mediumorchid3","mediumorchid3","darkgoldenrod1","black","darkgoldenrod1","darkgoldenrod1","mediumorchid3","mediumorchid3","mediumorchid3","darkgoldenrod1","black", "darkgoldenrod1")

colors <- c("black","mediumorchid3","mediumorchid3","mediumorchid3","darkgoldenrod1","black","darkgoldenrod1","darkgoldenrod1","mediumorchid3","mediumorchid3","black", "darkgoldenrod1")


fviz_pca_var(res.pca, alpha.var = "contrib", labelsize=5, col.var = grp, palette=colors, ggtheme=theme, arrowsize=1.5, col.circle = "White", col.quanti.sup = "white", repel=TRUE, select.var = list(contrib = 12)) + labs(title="HRwt") + theme(legend.position="none", axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                 axis.text = element_text(size = 12))


colors <- c("mediumorchid3","mediumorchid3","mediumorchid3","darkgoldenrod1","black","mediumorchid3","darkgoldenrod1","darkgoldenrod1","mediumorchid3","darkgoldenrod1","black", "darkgoldenrod1")


fviz_pca_var(res.pca, alpha.var = "contrib", labelsize=5,col.var = grp, palette=colors, ggtheme=theme, arrowsize=1.5, col.circle = "White", col.quanti.sup = "white", repel=TRUE, axes = c(2,3), select.var = list(contrib = 12)) + labs(title="HRwt")+ theme(legend.position="none", axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                              axis.text = element_text(size = 12))




fviz_pca_ind(res.pca, label="none", addEllipses = TRUE,habillage=annot$PFI, ggtheme=theme, pointsize = 3)+ labs(title="HRwt")+ scale_color_manual(values=c("#3399CC", "#003366"))+ theme(aspect.ratio=1)+ ylim(c(-13, 13)) + xlim(c(-13, 13))


fviz_pca_ind(res.pca, label="none", addEllipses = TRUE,habillage=annot$PFI, axes = c(2,3), ggtheme=theme, pointsize = 3)+ labs(title="HRwt")+ scale_color_manual(values=c("#3399CC", "#003366"))+ theme(aspect.ratio=1)+ ylim(c(-13, 13)) + xlim(c(-13, 13))




