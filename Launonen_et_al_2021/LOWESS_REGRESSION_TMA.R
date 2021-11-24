############################################################################################################################################
########################### CORRELATION SCATTER PLOTS WITH LOWESS REGRESSION ###############################################################
############################################################################################################################################

#For Figure 2 and Supplementary Figure 2

library(ggpubr)

immune_from_immune <- read.csv("immune_percentages_from_immune.csv")

ann <- read.csv("TMA_clinicaldata.csv")

#add simpson.i from "SIMPSON_CALC_TMA.R"

annotations_SDI <- cbind(ann, simpson.i)

immune_from_immune <- cbind(immune_from_immune, annotations_SDI)


colors <- c("blue", "red")
colors <- colors[as.numeric(immune_from_immune$Type)]

#lowess plots

#CD11c
plot(immune_from_immune$simpson_immune, immune_from_immune$CD11c..APC, col=colors, pch=16)
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "CD11c..APC"]), lty="dashed", lwd=2, col="blue")
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "HRwt"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "HRwt"), "CD11c..APC"]), lty="dashed", lwd=2, col="red")
lines(lowess(immune_from_immune$simpson_immune, immune_from_immune$CD11c..APC), lty="dashed", lwd=2)

#IBA1CD11c

plot(immune_from_immune$simpson_immune, immune_from_immune$IBA1.CD11c..Macrophages, col=colors, pch=16)
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "IBA1.CD11c..Macrophages"], f=10), lty="dashed", lwd=2, col="blue")
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "HRwt"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "HRwt"), "IBA1.CD11c..Macrophages"], f=10), lty="dashed", lwd=2, col="red")
lines(lowess(immune_from_immune$simpson_immune, immune_from_immune$IBA1.CD11c..Macrophages, f=10), lty="dashed", lwd=2)


#IBA1CD163

plot(immune_from_immune$simpson_immune, immune_from_immune$IBA1.CD163..Macrophages, col=colors, pch=16)
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "IBA1.CD163..Macrophages"], f=10), lty="dashed", lwd=2, col="blue")
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "HRwt"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "HRwt"), "IBA1.CD163..Macrophages"], f=10), lty="dashed", lwd=2, col="red")
lines(lowess(immune_from_immune$simpson_immune, immune_from_immune$IBA1.CD163..Macrophages, f=10), lty="dashed", lwd=2)


#IBA1CD163CD11c

plot(immune_from_immune$simpson_immune, immune_from_immune$IBA1.CD163.CD11c..Macrophages, col=colors, pch=16)
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "BRCA1/2 mutated"), "IBA1.CD163.CD11c..Macrophages"], f=10), lty="dashed", lwd=2, col="blue")
lines(lowess(immune_from_immune[which(immune_from_immune$Type == "HRwt"), "simpson_immune"], immune_from_immune[which(immune_from_immune$Type == "HRwt"), "IBA1.CD163.CD11c..Macrophages"], f=10), lty="dashed", lwd=2, col="red")
lines(lowess(immune_from_immune$simpson_immune, immune_from_immune$IBA1.CD163.CD11c..Macrophages, f=10), lty="dashed", lwd=2)





#width=3


#ggscatter original plot
ggscatter(immune_from_immune, x = "simpson_immune", y ="CD11c..APC" , 
          conf.int = F, size=1,
          cor.method = "spearman", fullrange=F,
          xlab = "Immune SDI", ylab = "CD11c+ APC", palette = c("blue", "red"), cor.coef.size = 2, color="Type"
)+
  stat_cor(aes(color = Type),label.y.npc = "top", method="spearman") + border("black") + grids(linetype = "dashed") +
  theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
        aspect.ratio=1) + stat_cor(method="spearman") + ylim(0, 40)
