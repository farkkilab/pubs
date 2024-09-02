#for Supplementary Figure 3e-i
library(ggpubr)
library(ggplot2)
library(dplyr)
#Use macrophages all from myelonet_heatmap

Macrophages <- read.csv("D:/Sciset/scimap/plots/dalaunay/Macrophage_blobs_size_mean_functional_marker_expression.csv")
Macrophages2 <- read.csv("D:/Sciset/scimap/plots/dalaunay/Immune_blobs_size_mean_functional_marker_expression.csv")
Macrophages3 <- read.csv("D:/Sciset/scimap/plots/dalaunay/IBA1CD163_Macrophage_blobs_size_mean_functional_marker_expression.csv")
Macrophages4 <- read.csv("D:/Sciset/scimap/plots/dalaunay/CD11c_blobs_size_mean_functional_marker_expression.csv")

Macrophages_all <- bind_rows(Macrophages, Macrophages1, Macrophages2, Macrophages3)

#size

my_comparisons <- list( c("CD11c", "IBA1CD163_Macrophages"), c("CD11c", "Immune"),c("CD11c", "Macrophages"), c("IBA1CD163_Macrophages", "Immune"),c("IBA1CD163_Macrophages", "Macrophages"),c("Immune", "Macrophages") )
ggplot(Macrophages_all, aes(x=factor(blob, levels = c("CD11c", "Immune",  "Macrophages", "IBA1CD163_Macrophages")), y=log(number), fill=blob)) + geom_violin() + geom_boxplot(width=0.3) + scale_fill_manual(values=c( "#d8315b", "#3e92cc","#fffaff",  "#6a2450")) + theme_bw() + ylab("log10 of myelonet size") + xlab("RCN")+stat_compare_means(comparisons=my_comparisons)

#summary statistics for the sizes of the myelonets

summarystat<- function(x) {
  z1 <- mean(x)
  z2 <- median(x)
  z3 <- sd(x)
  return(list(mean=z1, median=z2, sd=z3))
}

summarystat(Macrophages_all[which(Macrophages_all$blob == "CD11c"), "number"])

summarystat(Macrophages_all[which(Macrophages_all$blob == "Immune"), "number"])

summarystat(Macrophages_all[which(Macrophages_all$blob == "IBA1CD163_Macrophages"), "number"])

summarystat(Macrophages_all[which(Macrophages_all$blob == "Macrophages"), "number"])


ggplot(Macrophages_all, aes(x=number, fill=blob)) +
  geom_histogram(color="black")+
  geom_vline(aes(xintercept=median(number)), color="blue",
             linetype="dashed") + xlim(1, 300)

#confidence intervals

l.model <- lm(number ~ 1, Macrophages_all[which(Macrophages_all$blob == "CD11c"),])
confint(l.model, level=0.95)


l.model <- lm(number ~ 1, Macrophages_all[which(Macrophages_all$blob == "Macrophages"),])
confint(l.model, level=0.95)


l.model <- lm(number ~ 1, Macrophages_all[which(Macrophages_all$blob == "IBA1CD163_Macrophages"),])
confint(l.model, level=0.95)


l.model <- lm(number ~ 1, Macrophages_all[which(Macrophages_all$blob == "Immune"),])
confint(l.model, level=0.95)


#also, calculate the proportion of CD8+T-cells in different neighborhoods
#compare to other neighborhoods


data <- read.csv("E:/sciset/data_20231125.csv")

#percentage of cells per RCN

percentage_cells <- data %>% group_by(neighbordood_cluster2, GlobalCellType2) %>% summarize (n=n()) %>% mutate(freq=n/sum(n))
percentage_cells <- percentage_cells[which(percentage_cells$GlobalCellType2 == "CD8.T.cells"),]

ggplot(percentage_cells, aes(x=neighbordood_cluster2, y=freq)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

percentage_cells <- data %>% group_by(neighbordood_cluster2, GlobalCellType2) %>% summarize (n=n()) %>% mutate(freq=n/sum(n))

ggplot(percentage_cells[which(percentage_cells$GlobalCellType2 == "CD4.T.cells"),], aes(x=neighbordood_cluster2, y=freq)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


percentage_cells <- data %>% group_by(neighbordood_cluster2, GlobalCellType2) %>% summarize (n=n()) %>% mutate(freq=n/sum(n))

ggplot(percentage_cells[which(percentage_cells$GlobalCellType2 == "FOXP3.CD4.Tregs"),], aes(x=neighbordood_cluster2, y=freq)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


percentage_cells <- percentage_cells[which(percentage_cells$GlobalCellType2 == "CD8.T.cells" | percentage_cells$GlobalCellType2 == "CD4.T.cells" | percentage_cells$GlobalCellType2 == "FOXP3.CD4.Tregs"),]


ggplot(percentage_cells, aes(x=neighbordood_cluster2, y=freq, fill=GlobalCellType2)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


mean_cells <- percentage_cells %>% group_by(neighbordood_cluster2) %>% mutate(all=sum(freq))


unique(mean_cells[, c("neighbordood_cluster2", "all")])

#calculate per sample -> take the mean -> calculate CI

percentage_cells <- data %>% group_by(Sample_code, neighbordood_cluster2, GlobalCellType2) %>% summarize (n=n()) %>% mutate(freq=n/sum(n))

percentage_cells <- percentage_cells[which(percentage_cells$GlobalCellType2 == "CD8.T.cells" | percentage_cells$GlobalCellType2 == "CD4.T.cells" | percentage_cells$GlobalCellType2 == "FOXP3.CD4.Tregs"),]

mean_cells <- percentage_cells %>% group_by(Sample_code, neighbordood_cluster2) %>% mutate(all=sum(freq))

mean_cells <- unique(mean_cells[,c(1, 2, 6)])


l.model <- lm(all ~ 1, mean_cells[which(mean_cells$neighbordood_cluster2 == "Immune"),])
confint(l.model, level=0.95)
#0.1040863 0.1667547
#0.135

l.model <- lm(all ~ 1, mean_cells[which(mean_cells$neighbordood_cluster2 == "Macrophages"),])
confint(l.model, level=0.95)
#0.07481665 0.1371848
# 0.106

l.model <- lm(all ~ 1, mean_cells[which(mean_cells$neighbordood_cluster2 == "IBA1.CD163.Macrophages"),])
confint(l.model, level=0.95)
#0.05459872 0.09652183
#0.0756

l.model <- lm(all ~ 1, mean_cells[which(mean_cells$neighbordood_cluster2 == "CD11c.myeloid"),])
confint(l.model, level=0.95)
#0.1088247 0.1651053
#0.137

#then mean across samples

mean_c <- mean_cells %>% group_by(neighbordood_cluster2) %>% mutate(freq=mean(all))
mean_c <- unique(mean_c[,c(2, 4)])

#then plotting
mean_cells <- percentage_cells %>% group_by(neighbordood_cluster2, GlobalCellType2) %>% summarize(all=mean(freq))

ggplot(mean_cells, aes(x=factor(neighbordood_cluster2, levels=c("CD8_CD4_T.cells", "CD11c.myeloid", "Immune", "Macrophages", "Myofibroblast", "stroma", "tumor-stroma-interface", "SMA.CD31.positive", "SMA.Desmin.myofibroblast", "IBA1.CD163.Macrophages", "Desmin.positive", "Fibroblast", "Epithelial_EMT", "EMT", "Proliferating.EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Epithelial")), y=all, fill=factor(GlobalCellType2, levels=c("CD8.T.cells", "CD4.T.cells", "FOXP3.CD4.Tregs")))) + 
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_bw() + scale_fill_manual(values=c('#FF5A5F',"#57CC99","#BFD7EA" )) + xlab("RCN")


#boxplots of size comparisons myelonets

#calculate mean size per sample

#take paired samples
#plot


numbers <- Macrophages_all[which(Macrophages_all$number>9 & Macrophages_all$blob == "Macrophages"), ] %>% group_by(Sample_code) %>% summarize(mean_n = mean(number))

numbers <- numbers[which(numbers$Paired == T),]

numbers$Stage <- "chemo-naive"
numbers[grepl("i", numbers$Sample_code, fixed = TRUE), "Stage"] <- "chemo-exposed"

ggplot(numbers, aes(x=Stage, y=mean_n)) + geom_boxplot() + geom_point() + stat_compare_means(paired=T)




numbers <- Macrophages_all[which(Macrophages_all$number>9 & Macrophages_all$blob == "Immune"), ] %>% group_by(Sample_code) %>% summarize(mean_n = mean(number))

numbers <- numbers[which(numbers$Paired == T),]
numbers$Stage <- "chemo-naive"
numbers[grepl("i", numbers$Sample_code, fixed = TRUE), "Stage"] <- "chemo-exposed"

ggplot(numbers, aes(x=Stage, y=mean_n)) + geom_boxplot() + geom_point() + stat_compare_means(paired=T)

#


numbers <- Macrophages_all[which(Macrophages_all$number>9 & Macrophages_all$blob == "IBA1CD163_Macrophages"), ] %>% group_by(Sample_code) %>% summarize(mean_n = mean(number))

numbers <- numbers[which(numbers$Paired == T),]
numbers$Stage <- "chemo-naive"
numbers[grepl("i", numbers$Sample_code, fixed = TRUE), "Stage"] <- "chemo-exposed"

ggplot(numbers, aes(x=Stage, y=mean_n)) + geom_boxplot() + geom_point() + stat_compare_means(paired=T)


#

numbers <- Macrophages_all[which(Macrophages_all$number>9 & Macrophages_all$blob == "CD11c"), ] %>% group_by(Sample_code) %>% summarize(mean_n = mean(number))

numbers <- numbers[which(numbers$Paired == T),]

numbers$Stage <- "chemo-naive"
numbers[grepl("i", numbers$Sample_code, fixed = TRUE), "Stage"] <- "chemo-exposed"

ggplot(numbers, aes(x=Stage, y=mean_n)) + geom_boxplot() + geom_point() + stat_compare_means(paired=T)









