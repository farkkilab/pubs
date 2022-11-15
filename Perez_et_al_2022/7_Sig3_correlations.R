#Script to calculate correlations between ovaHRDscars and Signature3 in PCAWG and TCGA
#Fernando Perez

library(ggplot2)
library(ggpubr)

##################### Defining variables ##################
currentdir <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Path to the repository script
outputdir <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Path to put the resultant images

#ovaHRDscar cutoff for HRD and HRP samples is 59 events
newcutoff <- 54
#Myriad cutoff from prevscars and Tacaya is 42 and 63
myriadcutoff <- 42
takayacutoff <- 63

setwd(currentdir)
#Defining inputs for OVA samples
sig3file <- "data/sig3_OVA-TCGA_Gulhan2019.csv"
TCGAprevHRDfile <- "TCGA_HRD_prevScars_clinical-info_HGSC.csv"
TCGAnewHRDfile <- "TCGA_HRD_newScars_clinical-info_HGSC.csv"
PCAWG_HRD <- "PCAWG_newScars_clinical-info_HGSC.csv"
pcawgSignatures <- "data/PCAWG_sigProfiler_SBS_signatures_OVA.csv"
pcaw_histology <- "data/pcawg_Supplementary_Table1_Ova.csv"


############################################################################################################################
################## Sig3 exposure (WES) and scars (SNP-arrays) correlation in TCGA  #########################################
############################################################################################################################

#Loading data
sig3 <- read.table(file=sig3file, header=T, sep=",", row.names = 1)
newHRDscores <- read.table(file=TCGAnewHRDfile, header=T, sep=",")
prevHRDscores <- read.table(file=TCGAprevHRDfile, header=T, sep=",")

#Merging scars data, previous definition and ovaHRDscar
newHRDsum <- newHRDscores[,c(2:5)]
prevHRDsum <- prevHRDscores[,c(5)]
HRDscores <- cbind(newHRDsum, prevHRDsum)
row.names(HRDscores) <- newHRDscores[,1]

#Merging input data
sig3HRD <- merge(sig3, HRDscores, by="row.names")

#Check pearson correlation between ovaHRDscar in SNP-arrays and Sig3 from WES
m <- cor.test(sig3HRD$exp_sig3, sig3HRD$HRDsum,  method = "pearson", use = "complete.obs")
m$p.value #p-value added to Figure

#Plotting correlations for Supplementary Figure 2i
p <- ggplot(sig3HRD,aes(HRDsum, exp_sig3)) + geom_point(size= 2)
p <- p + geom_smooth(method='glm') + stat_cor(method = "pearson", size=6) + scale_y_continuous(limit= c(0, max(sig3HRD$exp_sig3)))
p <- p + labs(x="ovaHRDscar (SNP-array)", y = "Sig3 exposure (WES)")
p <- p + theme(axis.text=element_text(size=rel(1.3)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)), legend.text=element_text(size=rel(1.3)),
               strip.text.x = element_text(size=rel(1.3)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(1.5)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, filename = paste0(outputdir, "Sig3-ovaHRDscar_correlation_TCGA.svg"), width = 10, height = 12, units = "cm")


#Check pearson correlation between Telli2016 in SNP-arrays and Sig3 from WES
t <- cor.test(sig3HRD$exp_sig3, sig3HRD$prevHRDsum,  method = "pearson", use = "complete.obs")
t$p.value #p-value added to Figure
#Plotting correlation of Sig3 and Telli2016 for Supplementary Figure2k
p <- ggplot(sig3HRD,aes(prevHRDsum, exp_sig3)) + geom_point(size= 2)
p <- p + geom_smooth(method='glm') + stat_cor(method = "pearson", size=6) + scale_y_continuous(limit= c(0, max(sig3HRD$exp_sig3)))
p <- p + labs(x="Telli2016 (SNP-array)", y = "Sig3 exposure (WES)")
p <- p + theme(axis.text=element_text(size=rel(1.3)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)), legend.text=element_text(size=rel(1.3)),
               strip.text.x = element_text(size=rel(1.3)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(1.5)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, filename = paste(outputdir, "Sig3-Telli2016_correlation_TCGA.svg", sep="/"), width = 10, height = 12, units = "cm")

###################################################################################################################################
########################################### sig3 status and HRD status (by scars) concordance in TCGA  ############################
###################################################################################################################################

#This function to calculate Cohen's Kappa, using as input list of samples classified as HRD  for two methods
kappafun <- function(pos1, pos2, neg1, neg2){
  a <- length(intersect(pos1, pos2)) #Both pos
  b <- length(intersect(pos1, neg2)) #First pos, second neg
  c <- length(intersect(pos2, neg1)) #Second pos, first neg
  d <- length(intersect(neg1, neg2)) #Both neg
  print(paste0("Both positive: ", a))
  print(paste0("First pos, second neg: ", b))
  print(paste0("Second pos, first neg: ", c))
  print(paste0("Both negative: ",d))

  p0 <- (a + d) / (a + b + c + d)
  print(p0)

  pPos <- ((a + b) / (a + b + c + d)) *  ((a + c) / (a + b + c + d))
  pNeg <- ((c + d) / (a + b + c + d)) *  ((b + d) / (a + b + c + d))
  pe = pPos + pNeg

  k = (p0 - pe)/ (1-pe)
  return(k)
}

#Loading data
sig3 <- read.table(file=sig3file, header=T, sep=",", row.names = 1)
newHRDscores <- read.table(file=TCGAnewHRDfile, header=T, sep=",")
prevHRDscores <- read.table(file=TCGAprevHRDfile, header=T, sep=",")


#Selecting only matched with ovaHRDscar info
sig3 <- sig3[row.names(sig3) %in% newHRDscores$Row.names,]

####Selecting the samples used in the Sig3 analysis
newHRDscores_subset <- newHRDscores[newHRDscores$Row.names %in% row.names(sig3),]
prevHRDscores_subset <- prevHRDscores[prevHRDscores$Row.names %in% row.names(sig3),]

#Selecting samples names with positive or negative status for Sig3, Telli2016 and ovaHRDscar
sig3positive_samples <- row.names(sig3[sig3$pass_mva_strict == 1,])
sig3negative_samples <- row.names(sig3[sig3$pass_mva_strict == 0,])

newdef_HRDpositiveSamples <- newHRDscores_subset[which(newHRDscores_subset$HRDsum  >= newcutoff),1]
newdef_HRDnegativeSamples <- newHRDscores_subset[which(newHRDscores_subset$HRDsum  < newcutoff),1]

myriad_HRDpositiveSamples <- prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  >= myriadcutoff),1]
myriad_HRDnegativeSamples <- prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  < myriadcutoff),1]

tacaya_HRDpositiveSamples <- prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  >= takayacutoff),1]
tacaya_HRDnegativeSamples <- prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  < takayacutoff),1]
#For Telli2016 but using the cutoff 54 as in ovaHRDscar
myriad_HRDpositiveSamples_54 <- prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  >= newcutoff),1]
myriad_HRDnegativeSamples_54 <- prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  < newcutoff),1]


#The next values correspond to  Figure2g
#Calculating Cohen's Kappa values for concordance of samples selected under each criteria
#For ovaHRDscar
kappafun(pos1=sig3positive_samples, pos2=newdef_HRDpositiveSamples, neg1=sig3negative_samples, neg2=newdef_HRDnegativeSamples)
#For Telli2016
kappafun(pos1=sig3positive_samples, pos2=myriad_HRDpositiveSamples, neg1=sig3negative_samples, neg2=myriad_HRDnegativeSamples)
#Extra (not in Figure) for Takaya2020
kappafun(pos1=sig3positive_samples, pos2=tacaya_HRDpositiveSamples, neg1=sig3negative_samples, neg2=tacaya_HRDnegativeSamples)
#Extra (not in Figure) for Telli2016 but using the cutoff 54 as in ovaHRDscar
kappafun(pos1=sig3positive_samples, pos2=myriad_HRDpositiveSamples_54, neg1=sig3negative_samples, neg2=myriad_HRDnegativeSamples_54)


length(sig3positive_samples) / nrow(newHRDscores_subset)

#Generating small pie-charts for description of abundance of HRD and HRP under each criteria
#For Figure2g
newHRDscores_subset$ovaHRDscar <- rep("HRP", nrow(newHRDscores_subset))
newHRDscores_subset[which(newHRDscores_subset$HRDsum  >= newcutoff),"ovaHRDscar"] <- "HRD"
table(newHRDscores_subset$ovaHRDscar) / length(newHRDscores_subset$ovaHRDscar)
#For ovaHRDscar
p <- ggplot(newHRDscores_subset, aes(x=factor(1), fill=ovaHRDscar))+ geom_bar(width = 1, colour="black")
p <- p + theme_void() + theme(legend.position = "none") + coord_polar("y")
p <- p + scale_fill_manual(values=c("#FA5A41", "#34A0D3"))
p
ggsave(p, filename = paste0(outputdir, "Piechart_HRD_HRP_ovaHRDcar.svg"))

#For Telli2016 for  Figure2g
prevHRDscores_subset$Telli2016status <- rep("HRP", nrow(prevHRDscores_subset))
prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  >= myriadcutoff),"Telli2016status"] <- "HRD"
table(prevHRDscores_subset$Telli2016status) / length(prevHRDscores_subset$Telli2016status)
p <- ggplot(prevHRDscores_subset, aes(x=factor(1), fill=Telli2016status))+ geom_bar(width = 1)
p <- p + theme_void() + theme(legend.position = "none") + coord_polar("y")
p <- p + scale_fill_manual(values=c("#FA5A41", "#34A0D3"))
p
ggsave(p, filename = paste0(outputdir, "Piechart_HRD_HRP_Telli2016.svg"))


#For Telli2016 for  Supp Figure2k
prevHRDscores_subset$Telli2016status <- rep("HRP", nrow(prevHRDscores_subset))
prevHRDscores_subset[which(prevHRDscores_subset$HRDsum  >= newcutoff),"Telli2016status"] <- "HRD"
table(prevHRDscores_subset$Telli2016status) / length(prevHRDscores_subset$Telli2016status)
p <- ggplot(prevHRDscores_subset, aes(x=factor(1), fill=Telli2016status))+ geom_bar(width = 1, colour="black")
p <- p + theme_void() + theme(legend.position = "none") + coord_polar("y")
p <- p + scale_fill_manual(values=c("#FA5A41", "#34A0D3"))
p
ggsave(p, filename = paste0(outputdir, "Piechart_HRD_HRP_Telli2016-54.svg"))


##############################################################################################################
######################### Sig3 and scars correlations in PCAWG samples  ######################################
##############################################################################################################

#Reading signature3 and scars
pcawg_sbs <- read.table(file=pcawgSignatures, sep=",", header=T)
HRDpcaw <- read.table(file=PCAWG_HRD, sep=",", header=T)

##Calculate SBS3 proportion in samples
i<-NULL
SBS3proportion <- NULL
for (i in 1:nrow(pcawg_sbs)){
  SBS3value <- pcawg_sbs[i,6]
  SBStotal <- sum(pcawg_sbs[i,c(4:68)])
  SBS3value <- SBS3value/SBStotal
  SBS3proportion <- c(SBS3proportion, SBS3value)
}
PCAWG_SBS3prop <- cbind(pcawg_sbs[,c(1:3)],SBS3proportion)

#Merging with HRD dataset
PCAWG_SBS3prop_HRD <- merge(HRDpcaw, PCAWG_SBS3prop, by.x = "icgc_specimen_id", by.y="Sample.Names")


#Finding correlation of Sig3 and ovaHRDscar
#For Figure2f:
p <- ggplot(PCAWG_SBS3prop_HRD,aes(HRDsum, SBS3proportion)) + geom_point(size= 2)
p <- p +  stat_cor(method = "pearson", size=5) + geom_smooth(method='glm')
p <- p + labs(x="ovaHRDscar (WGS)", y = "SBS3 proportion (WGS)")
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.6)), axis.text.y=element_text(size=rel(1.6)), legend.text=element_text(size=rel(1.6)),
               strip.text.x = element_text(size=rel(1.3)),
               legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(1.9)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, filename = paste(outputdir, "Sig3-ovaHRDscar_correlation_PCAWG.svg", sep="/"), width = 10, height = 12, units = "cm")


#Finding correlation of Sig3 and Telli2016
#For Supplementary Figure2j:
p <- ggplot(PCAWG_SBS3prop_HRD,aes(pHRDsum, SBS3proportion)) + geom_point(size= 2)
p <- p +  stat_cor(method = "pearson", size=5) + geom_smooth(method='glm')
p <- p + labs(x="Telli2016 (WGS)", y = "SBS3 proportion (WGS)")
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.6)), axis.text.y=element_text(size=rel(1.6)), legend.text=element_text(size=rel(1.3)),
               strip.text.x = element_text(size=rel(1.3)),
               legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(1.9)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, filename = paste(outputdir, "Sig3-Telli2016_correlation_PCAWG.svg", sep="/"), width = 10, height = 12, units = "cm")
