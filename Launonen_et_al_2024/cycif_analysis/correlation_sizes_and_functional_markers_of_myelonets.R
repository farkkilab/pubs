#sizes of myelonets and functional marker expressions for myelonets for further analysis

data <- read.csv("E:/Sciset/data_20230511.csv")

clinical <- read.csv("L:/ltdk_farkkila/Data/Sciset/clinical_data/clinical_data.csv")

data <- merge(clinical, data, by="imageid")
data[which(data$paired_patients == TRUE),]
obs_all <- data

#first RCN rich in Macrophages
#load in results from delaunay triangulation computed in cpu1

macrophage_blobs_all <- list()
for (i in c("s01","s02","s03","s04","s05","s06",'s07', "s08", "s09", "s11", "s12", "s13", "s14", "s15","s16","s17", "s18", "s19", "s20", "s21", "s22", "s23")){
  
  macrophage_blob <- read.csv(paste0("D:/sciset/delaunay/",i,"_macrophages.csv"))
  
  
  macrophage_blob <- unique(macrophage_blob[,-c(1)])
  colnames(macrophage_blob)[4] <- "ID"
  
  s01 <- obs_all[which(obs_all$imageid == i), c("ID", "pSTAT1", "TIM3", "MHCII","Ki67","MHCI","SNAT1","Annexin","Sample_code", "Patient", "GlobalCellType2")]
  
  macrophage_blob <- merge(macrophage_blob, s01, by="ID")
  
  macrophage_blobs_all[[i]] <- macrophage_blob
  
}


macrophage_blobs_all <- do.call(rbind, macrophage_blobs_all)

#then calculate sizes

blobs <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i),] %>% 
    group_by(components.membership) %>% 
    summarise(number = n())
  blob$Sample_code <- i
  
  blobs[[i]] <- blob
}


blobs <- do.call(rbind, blobs)

blobs <- blobs[-which(blobs$number<9),]

#
blobs$myelonet <- paste0(blobs$components.membership,".", blobs$Sample_code)
macrophage_blobs_all$myelonet <- paste0(macrophage_blobs_all$components.membership, ".", macrophage_blobs_all$Sample_code)

macrophage_blobs_CellIDs <- macrophage_blobs_all[which(macrophage_blobs_all$myelonet %in% blobs$myelonet), c("ID", "Sample_code","GlobalCellType2", "myelonet")]
write.csv(macrophage_blobs_CellIDs, "D:/sciset/delaunay/macrophage_myelonet_CellIDs.csv", row.names = F)

#sit boxplot mean size per sample pre ja post

#mean per sample

blobs_mean <- blobs[,-c(1)] %>%
  group_by(Sample_code)%>%
  summarize(mean_n = mean(number))


#all clinical

blobs_mean <- merge(blobs_mean, clinical[,c('Sample_code','Patient.x', 'Stage', 'HRD_status', 'paired_patient')], by="Sample_code")

blobs_mean <- blobs_mean[which(blobs_mean$åaired_patient == TRUE),]
blobs_mean$Stage <- factor(blobs_mean$Stage, levels=c("primary", "interval"))
pdf("E:/Sciset/scimap/plots/boxplots/boxplot_macrophage_blob_p0.3125.pdf")
ggplot(blobs_mean[which(blobs_mean$Paired == T),], aes(x=Stage, y=mean_n)) + geom_boxplot(aes(fill=Stage)) + geom_point() + scale_fill_manual(values=c("blue", "yellow")) + geom_line(aes(group=Patient.x))
dev.off()
print(wilcox.test(blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "primary"), "mean_n"], blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "interval"), "mean_n"], paired = T))

#sizes

mean(blobs_mean$mean_n)
median(blobs_mean$mean_n)
sd(blobs_mean$mean_n)

mean(blobs_mean[which(blobs_mean$Stage == "interval" & blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])

mean(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])

#functional markers of myeloid cells

blobs_m <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & (macrophage_blobs_all$GlobalCellType2 == "IBA1.CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "IBA1.CD11c.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD11c.myeloid")),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_MHCII = mean(MHCII),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_Annexin = mean(Annexin),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_m[[i]] <- blob
}

blobs_m <- do.call(rbind, blobs_m)

#functional markers of CD8T-cells inside the blobs
blobs_t <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & macrophage_blobs_all$GlobalCellType2 == "CD8.T.cells"),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_TIM3 = mean(TIM3),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_t[[i]] <- blob
}

blobs_t <- do.call(rbind, blobs_t)

colnames(blobs_t)[2:6] <- c("TIM3_cd8", "MHCI_cd8", "SNAT1_cd8", "Ki67_cd8", "pSTAT1_cd8")

merged_blobs <- merge(blobs, blobs_m, by=c("Sample_code", "components.membership"))
merged_blobs <- merge(merged_blobs, clinical, by="Sample_code")
merged_blobs <- merge(merged_blobs, blobs_t, by=c("Sample_code", "components.membership"), all.x = T)

#cell type proportions

blobs_p <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i),] %>% 
    group_by(components.membership, GlobalCellType2) %>% 
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  blob$Sample_code <- i
  
  blobs_p[[i]] <- blob
}

blobs_p <- do.call(rbind, blobs_p)
blobs_p <- blobs_p[,-3]
#tän pitäis mieluummin olla wide
blobs_p <- pivot_wider(blobs_p, names_from="GlobalCellType2", values_from="freq")
blobs_p <- as.data.frame(blobs_p)
blobs_p[is.na(blobs_p)] <- 0


#colnames(blobs_t)[2:6] <- c("TIM3_cd8", "MHCI_cd8", "SNAT1_cd8", "Ki67_cd8", "pSTAT1_cd8")

merged_blobs <- merge(merged_blobs, blobs_p, by=c("Sample_code", "components.membership"), all.x=T)
write.csv(merged_blobs, "E:/Sciset/scimap/plots/dalaunay/Macrophage_blobs_size_mean_functional_marker_expression_and_proportions.csv", row.names = F)


#hist(merged_blobs$number, breaks = 30)


ggplot(merged_blobs[which(merged_blobs$Stage == "interval"),], aes(x=mean_Annexin, y=log2(number), color=Stage)) + geom_point()+
  geom_smooth(method = "lm")

library(ggpubr)

ggscatter(merged_blobs, x = "pSTAT1_cd8", y ="number" , 
          conf.int = F, size=1,
          cor.method = "pearson", fullrange=F,
          ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
)+
  stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
  theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
        aspect.ratio=1) + stat_cor(method="spearman") + ylim(0, 1500)


merged_blobs <- merged_blobs[-which(merged_blobs$number < 10),]

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_size_Macrophages.pdf")
for (i in colnames(merged_blobs)[c(4:9, 22:26)]){
  p <- ggscatter(merged_blobs, x = i, y ="number" , 
                 conf.int = F, size=1,
                 cor.method = "pearson", fullrange=F,
                 ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
  )+
    stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
    theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
          aspect.ratio=1)
  
  print(p)
}

dev.off()



pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_func_markers_Macrophages.pdf")
for (i in colnames(merged_blobs)[c(4:9)]){
  for (j in colnames(merged_blobs)[c(22:26)]){
    p <- ggscatter(merged_blobs, x = i, y =j , 
                   conf.int = F, size=1,
                   cor.method = "pearson", fullrange=F,
                   palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
    )+
      stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
      theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
            aspect.ratio=1)
    
    print(p)
  }}

dev.off()


write.csv(merged_blobs, "E:/Sciset/scimap/plots/dalaunay/Macrophage_blobs_size_mean_functional_marker_expression.csv", row.names = F)


#then for individual myelonet RCNs


macrophage_blobs_all <- list()
for (i in c("s01","s02","s03","s04","s05","s06",'s07', "s08", "s09", "s11", "s12", "s13", "s14", "s15","s16","s17", "s18", "s19", "s20", "s21", "s22", "s23")){
  
  macrophage_blob <- read.csv(paste0("D:/sciset/delaunay/",i,"_IBA1CD163_macs.csv"))
  
  
  macrophage_blob <- unique(macrophage_blob[,-c(1)])
  colnames(macrophage_blob)[4] <- "ID"
  
  s01 <- obs_all[which(obs_all$imageid == i), c("ID", "pSTAT1", "TIM3", "MHCII","Ki67","MHCI","SNAT1","Annexin","Sample_code", "Patient", "GlobalCellType2")]
  
  #tää seuraava on vaan jos yli s18
  
  
  
  macrophage_blob <- merge(macrophage_blob, s01, by="ID")
  
  macrophage_blobs_all[[i]] <- macrophage_blob
  
}


macrophage_blobs_all <- do.call(rbind, macrophage_blobs_all)

blobs <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i),] %>% 
    group_by(components.membership) %>% 
    summarise(number = n())
  blob$Sample_code <- i
  
  blobs[[i]] <- blob
}


blobs <- do.call(rbind, blobs)

blobs <- blobs[-which(blobs$number<9),]
blobs$myelonet <- paste0(blobs$components.membership,".", blobs$Sample_code)
macrophage_blobs_all$myelonet <- paste0(macrophage_blobs_all$components.membership, ".", macrophage_blobs_all$Sample_code)

macrophage_blobs_CellIDs <- macrophage_blobs_all[which(macrophage_blobs_all$myelonet %in% blobs$myelonet), c("ID", "Sample_code","GlobalCellType2", "myelonet")]
write.csv(macrophage_blobs_CellIDs, "D:/sciset/delaunay/IBA1.CD163.macrophage_myelonet_CellIDs.csv", row.names = F)

#then boxplot mean size per sample pre ja post

#mean per sample

blobs_mean <- blobs[,-c(1)] %>%
  group_by(Sample_code)%>%
  summarize(mean_n = mean(number))


#add clinical

blobs_mean <- merge(blobs_mean, clinical[,c('Sample_code','Patient.x', 'Stage', 'HRD_status')], by="Sample_code")

blobs_mean <- blobs_mean[which(blobs_mean$paired_patient == TRUE),]
blobs_mean$Stage <- factor(blobs_mean$Stage, levels=c("primary", "interval"))
pdf("E:/Sciset/scimap/plots/boxplots/boxplot_IBA1CD163macrophage_blob_p_0.8438.pdf")
ggplot(blobs_mean[which(blobs_mean$Paired == T),], aes(x=Stage, y=mean_n)) + geom_boxplot(aes(fill=Stage)) + geom_point() + scale_fill_manual(values=c("blue", "yellow")) + geom_line(aes(group=Patient.x))
dev.off()
print(wilcox.test(blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "primary"), "mean_n"], blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "interval"), "mean_n"], paired = T))
#p=0.5625
#laske koko

mean(blobs_mean$mean_n)
median(blobs_mean$mean_n)
sd(blobs_mean$mean_n)

mean(blobs_mean[which(blobs_mean$Stage == "interval" & blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])

mean(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])



blobs_m <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & (macrophage_blobs_all$GlobalCellType2 == "IBA1.CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "IBA1.CD11c.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD11c.myeloid")),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_MHCII = mean(MHCII),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_Annexin = mean(Annexin),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_m[[i]] <- blob
}

blobs_m <- do.call(rbind, blobs_m)

blobs_t <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & macrophage_blobs_all$GlobalCellType2 == "CD8.T.cells"),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_TIM3 = mean(TIM3),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_t[[i]] <- blob
}

blobs_t <- do.call(rbind, blobs_t)

colnames(blobs_t)[2:6] <- c("TIM3_cd8", "MHCI_cd8", "SNAT1_cd8", "Ki67_cd8", "pSTAT1_cd8")

merged_blobs <- merge(blobs, blobs_m, by=c("Sample_code", "components.membership"))
merged_blobs <- merge(merged_blobs, clinical, by="Sample_code")
merged_blobs <- merge(merged_blobs, blobs_t, by=c("Sample_code", "components.membership"), all.x = T)

#


ggplot(merged_blobs[which(merged_blobs$Stage == "interval"),], aes(x=mean_Annexin, y=log2(number), color=Stage)) + geom_point()+
  geom_smooth(method = "lm")

library(ggpubr)

ggscatter(merged_blobs, x = "pSTAT1_cd8", y ="number" , 
          conf.int = F, size=1,
          cor.method = "pearson", fullrange=F,
          ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
)+
  stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
  theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
        aspect.ratio=1)


#loops

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_size_IBA1CD163.pdf")
for (i in colnames(merged_blobs)[c(4:9, 22:26)]){
  p <- ggscatter(merged_blobs, x = i, y ="number" , 
                 conf.int = F, size=1,
                 cor.method = "pearson", fullrange=F,
                 ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
  )+
    stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
    theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
          aspect.ratio=1)
  
  print(p)
}

dev.off()


merged_blobs <- merged_blobs[-which(merged_blobs$number < 10),]

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_func_markers_IBA1CD163.pdf")
for (i in colnames(merged_blobs)[c(4:9)]){
  for (j in colnames(merged_blobs)[c(22:26)]){
    p <- ggscatter(merged_blobs, x = i, y =j , 
                   conf.int = F, size=1,
                   cor.method = "pearson", fullrange=F,
                   palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
    )+
      stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
      theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
            aspect.ratio=1)
    
    print(p)
  }}

dev.off()

write.csv(merged_blobs, "E:/Sciset/scimap/plots/dalaunay/IBA1CD163_Macrophage_blobs_size_mean_functional_marker_expression.csv", row.names = F)


#CD11c+ myeloid

data <- read.csv("E:/Sciset/data_20230511.csv")

data <- merge(data, clinical, by="imageid")
data <- data[which(data$paired_patient == T),]
obs_all <- data

macrophage_blobs_all <- list()
for (i in c("s01","s02","s03","s04","s05","s06",'s07', "s08", "s09", "s11", "s12", "s13", "s14", "s15","s16","s17", "s18", "s19", "s20", "s21", "s22")){
  
  macrophage_blob <- read.csv(paste0("D:/sciset/delaunay/cd11c/",i,"_CD11c_0724.csv"))
  
  
  macrophage_blob <- unique(macrophage_blob[,-c(1)])
  #colnames(macrophage_blob)[3] <- "ID"
  
  s01 <- obs_all[which(obs_all$imageid == i), c("ID", "pSTAT1", "TIM3", "MHCII","Ki67","MHCI","SNAT1","Annexin","Sample_code", "Patient", "GlobalCellType2")]
  
  
  macrophage_blob <- merge(macrophage_blob, s01, by.x="ID", by.y="ID")
  
  macrophage_blobs_all[[i]] <- macrophage_blob
  
}


macrophage_blobs_all <- do.call(rbind, macrophage_blobs_all)

library(dplyr)
blobs <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i),] %>% 
    group_by(components.membership) %>% 
    summarise(number = n())
  blob$Sample_code <- i
  
  blobs[[i]] <- blob
}


blobs <- do.call(rbind, blobs)

blobs <- blobs[-which(blobs$number<9),]


blobs$myelonet <- paste0(blobs$components.membership,".", blobs$Sample_code)
macrophage_blobs_all$myelonet <- paste0(macrophage_blobs_all$components.membership, ".", macrophage_blobs_all$Sample_code)

macrophage_blobs_CellIDs <- macrophage_blobs_all[which(macrophage_blobs_all$myelonet %in% blobs$myelonet), c("ID", "Sample_code","GlobalCellType2", "myelonet")]

write.csv(macrophage_blobs_CellIDs, "D:/sciset/delaunay/CD11c_myelonet_CellIDs.csv", row.names = F)


#mean per sample

blobs_mean <- blobs[,-c(1)] %>%
  group_by(Sample_code)%>%
  summarize(mean_n = mean(number))


#lisää clinical

blobs_mean <- merge(blobs_mean, clinical[,c('Sample_code','Patient.x', 'Stage', 'HRD_status')], by="Sample_code")

blobs_mean$Paired <- FALSE
blobs_mean[which(blobs_mean$Patient.x == "H114" | blobs_mean$Patient.x == "H116" | blobs_mean$Patient.x == "H142" | blobs_mean$Patient.x == "H144" | blobs_mean$Patient.x == "H166" | blobs_mean$Patient.x == "H173"), "Paired"] <- TRUE

blobs_mean$Stage <- factor(blobs_mean$Stage, levels=c("primary", "interval"))
pdf("E:/Sciset/scimap/plots/boxplots/boxplot_CD11c_blob_0.563.pdf")
ggplot(blobs_mean[which(blobs_mean$Paired == T),], aes(x=Stage, y=mean_n)) + geom_boxplot(aes(fill=Stage)) + geom_point() + scale_fill_manual(values=c("blue", "yellow")) + geom_line(aes(group=Patient.x))
dev.off()
print(wilcox.test(blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "primary"), "mean_n"], blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "interval"), "mean_n"], paired = T))

#size

mean(blobs_mean$mean_n)
median(blobs_mean$mean_n)
sd(blobs_mean$mean_n)

mean(blobs_mean[which(blobs_mean$Stage == "interval" & blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])

mean(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])



blobs_m <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & (macrophage_blobs_all$GlobalCellType2 == "IBA1.CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "IBA1.CD11c.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD11c.myeloid")),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_MHCII = mean(MHCII),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_Annexin = mean(Annexin),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_m[[i]] <- blob
}

blobs_m <- do.call(rbind, blobs_m)

blobs_t <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & macrophage_blobs_all$GlobalCellType2 == "CD8.T.cells"),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_TIM3 = mean(TIM3),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_t[[i]] <- blob
}

blobs_t <- do.call(rbind, blobs_t)

colnames(blobs_t)[2:6] <- c("TIM3_cd8", "MHCI_cd8", "SNAT1_cd8", "Ki67_cd8", "pSTAT1_cd8")

merged_blobs <- merge(blobs, blobs_m, by=c("Sample_code", "components.membership"))
merged_blobs <- merge(merged_blobs, clinical, by="Sample_code")
merged_blobs <- merge(merged_blobs, blobs_t, by=c("Sample_code", "components.membership"), all.x = T)

#

ggplot(merged_blobs[which(merged_blobs$Stage == "interval"),], aes(x=mean_Annexin, y=log2(number), color=Stage)) + geom_point()+
  geom_smooth(method = "lm")

library(ggpubr)

ggscatter(merged_blobs, x = "pSTAT1_cd8", y ="number" , 
          conf.int = F, size=1,
          cor.method = "pearson", fullrange=F,
          ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
)+
  stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
  theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
        aspect.ratio=1)


#loops

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_size_CD11c_myeloid.pdf")
for (i in colnames(merged_blobs)[c(4:9, 22:26)]){
  p <- ggscatter(merged_blobs, x = i, y ="number" , 
                 conf.int = F, size=1,
                 cor.method = "pearson", fullrange=F,
                 ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
  )+
    stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
    theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
          aspect.ratio=1)
  
  print(p)
}

dev.off()


merged_blobs <- merged_blobs[-which(merged_blobs$number < 10),]

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_func_markers_CD11c_myeloid.pdf")
for (i in colnames(merged_blobs)[c(4:9)]){
  for (j in colnames(merged_blobs)[c(22:26)]){
    p <- ggscatter(merged_blobs, x = i, y =j , 
                   conf.int = F, size=1,
                   cor.method = "pearson", fullrange=F,
                   palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
    )+
      stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
      theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
            aspect.ratio=1)
    
    print(p)
  }}

dev.off()

write.csv(merged_blobs, "E:/Sciset/scimap/plots/dalaunay/CD11c_blobs_size_mean_functional_marker_expression.csv", row.names = F)


#mixed immune RCN

macrophage_blobs_all <- list()
for (i in c("s01","s02","s03","s04","s05","s06",'s07', "s08", "s09", "s11", "s12", "s13", "s14", "s15","s16","s17", "s18", "s19", "s20", "s21", "s22")){
  
  macrophage_blob <- read.csv(paste0("D:/sciset/delaunay/Immune/",i,"_Immune_0724.csv"))
  
  
  macrophage_blob <- unique(macrophage_blob[,-c(1)])
  #colnames(macrophage_blob)[4] <- "ID"
  
  s01 <- obs_all[which(obs_all$imageid == i), c("ID", "pSTAT1", "TIM3", "MHCII","Ki67","MHCI","SNAT1","Annexin","Sample_code", "Patient", "GlobalCellType2")]
  
  macrophage_blob <- merge(macrophage_blob, s01, by="ID")
  
  macrophage_blobs_all[[i]] <- macrophage_blob
  
}


macrophage_blobs_all <- do.call(rbind, macrophage_blobs_all)

blobs <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i),] %>% 
    group_by(components.membership) %>% 
    summarise(number = n())
  blob$Sample_code <- i
  
  blobs[[i]] <- blob
}


blobs <- do.call(rbind, blobs)

blobs <- blobs[-which(blobs$number<9),]

blobs$myelonet <- paste0(blobs$components.membership,".", blobs$Sample_code)
macrophage_blobs_all$myelonet <- paste0(macrophage_blobs_all$components.membership, ".", macrophage_blobs_all$Sample_code)

macrophage_blobs_CellIDs <- macrophage_blobs_all[which(macrophage_blobs_all$myelonet %in% blobs$myelonet), c("ID", "Sample_code","GlobalCellType2", "myelonet")]

write.csv(macrophage_blobs_CellIDs, "D:/sciset/delaunay/Immune_myelonet_CellIDs.csv", row.names = F)




#then boxplot mean size per sample pre ja post

#mean per sample

blobs_mean <- blobs[,-c(1)] %>%
  group_by(Sample_code)%>%
  summarize(mean_n = mean(number))


#add clinical

blobs_mean <- merge(blobs_mean, clinical[,c('Sample_code','Patient.x', 'Stage', 'HRD_status', 'paired_patient')], by="Sample_code")

blobs_mean <- blobs_mean[which(blobs_mean$Paired == T),]
blobs_mean$Stage <- factor(blobs_mean$Stage, levels=c("primary", "interval"))
pdf("E:/Sciset/scimap/plots/boxplots/boxplot_immune_blob_p_1.pdf")
ggplot(blobs_mean[which(blobs_mean$Paired == T),], aes(x=Stage, y=mean_n)) + geom_boxplot(aes(fill=Stage)) + geom_point() + scale_fill_manual(values=c("blue", "yellow")) + geom_line(aes(group=Patient.x))
dev.off()
print(wilcox.test(blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "primary"), "mean_n"], blobs_mean[which(blobs_mean$Paired == TRUE & blobs_mean$Stage == "interval"), "mean_n"], paired = T))
#p=1
#laske koko

mean(blobs_mean$mean_n)
median(blobs_mean$mean_n)
sd(blobs_mean$mean_n)

mean(blobs_mean[which(blobs_mean$Stage == "interval" & blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "interval"& blobs_mean$Paired == T), "mean_n"])

mean(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
median(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])
sd(blobs_mean[which(blobs_mean$Stage == "primary"& blobs_mean$Paired == T), "mean_n"])



blobs_m <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & (macrophage_blobs_all$GlobalCellType2 == "IBA1.CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "IBA1.CD11c.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD163.Macrophages" | macrophage_blobs_all$GlobalCellType2 == "CD11c.myeloid")),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_MHCII = mean(MHCII),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_Annexin = mean(Annexin),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_m[[i]] <- blob
}

blobs_m <- do.call(rbind, blobs_m)

blobs_t <- list()
for (i in unique(macrophage_blobs_all$Sample_code)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$Sample_code == i & macrophage_blobs_all$GlobalCellType2 == "CD8.T.cells"),] %>% 
    group_by(components.membership) %>% 
    summarise(mean_TIM3 = mean(TIM3),
              mean_MHCI = mean(MHCI),
              mean_SNAT1 = mean(SNAT1),
              mean_Ki67 = mean(Ki67),
              mean_pSTAT1 = mean(pSTAT1))
  
  blob$Sample_code <- i
  
  blobs_t[[i]] <- blob
}

blobs_t <- do.call(rbind, blobs_t)

colnames(blobs_t)[2:6] <- c("TIM3_cd8", "MHCI_cd8", "SNAT1_cd8", "Ki67_cd8", "pSTAT1_cd8")

merged_blobs <- merge(blobs, blobs_m, by=c("Sample_code", "components.membership"))
merged_blobs <- merge(merged_blobs, clinical, by="Sample_code")
merged_blobs <- merge(merged_blobs, blobs_t, by=c("Sample_code", "components.membership"), all.x = T)



ggplot(merged_blobs[which(merged_blobs$Stage == "interval"),], aes(x=mean_Annexin, y=log2(number), color=Stage)) + geom_point()+
  geom_smooth(method = "lm")

library(ggpubr)

ggscatter(merged_blobs, x = "pSTAT1_cd8", y ="number" , 
          conf.int = F, size=1,
          cor.method = "pearson", fullrange=F,
          ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
)+
  stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
  theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
        aspect.ratio=1)


#loops

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_size_immune.pdf")
for (i in colnames(merged_blobs)[c(4:9, 22:26)]){
  p <- ggscatter(merged_blobs, x = i, y ="number" , 
                 conf.int = F, size=1,
                 cor.method = "pearson", fullrange=F,
                 ylab = "size", palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
  )+
    stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
    theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
          aspect.ratio=1)
  
  print(p)
}

dev.off()

#merged_blobs$number <- as.numeric(merged_blobs$number)
#merged_blobs <- merged_blobs[-which(merged_blobs$number < 10),]

pdf("E:/Sciset/scimap/plots/dalaunay/cor_functional_markers_to_func_markers_immune.pdf")
for (i in colnames(merged_blobs)[c(4:9)]){
  for (j in colnames(merged_blobs)[c(22:26)]){
    p <- ggscatter(merged_blobs, x = i, y =j , 
                   conf.int = F, size=1,
                   cor.method = "pearson", fullrange=F,
                   palette = c("blue", "yellow"), cor.coef.size = 2, color="Stage"
    )+
      stat_cor(aes(color = Stage), method="pearson") + border("black") + grids(linetype = "dashed") +
      theme(legend.title=element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=10),
            aspect.ratio=1)
    
    print(p)
  }}

dev.off()


write.csv(merged_blobs, "E:/Sciset/scimap/plots/dalaunay/Immune_blobs_size_mean_functional_marker_expression.csv", row.names = F)
#


