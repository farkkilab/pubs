library(ggplot2)
library(ggpubr)
library(rpart)
library(rpart.utils)
library(rpart.plot)


#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")

####### Defining variables ###########
outputfolder = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Mainly for plots

###### Reading input segments #######
segHRDTCGA <- read.table(file = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signatures/HRD_samples-TCGAsegmentsHGSC_2021.txt", header=T)
segHRPTCGA <- read.table(file = "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signatures/HRP_samples-TCGAsegmentsHGSC_2021.txt", header=T)

#Ploidy information
Ploidy_file <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signatures/samplesHR_ploidy-purity2.txt ", header = T, sep="\t")

################################################## Preparing input ############################################################

##Preprocesing the segments from TCGA portal
chrominfo = chrominfo_grch38
ploidy <- rep(2, nrow(segHRDTCGA))
segHRDTCGA <- cbind (segHRDTCGA, ploidy)
ploidy <- rep(2, nrow(segHRPTCGA))
segHRPTCGA <- cbind (segHRPTCGA, ploidy)

preprocessed_HRD_TCGA <- preparing.input(segHRDTCGA)
preprocessed_HRP_TCGA <- preparing.input(segHRPTCGA)

### Sample TCGA-13-1511 (HRP), is an HRP outlier, we removed it ###
preprocessed_HRP_TCGA <- preprocessed_HRP_TCGA[!preprocessed_HRP_TCGA[,1] %in%  "TCGA-13-1511",]

######################################## Calculating the nTAIs #############################################################
#Calculus of nTAIs in the GRCh38 datasets
res_ai <- calc.ai_new(seg = preprocessed_HRD_TCGA, chrominfo = chrominfo_grch38)
nTAIs_HRD_GRCh38_removedsmall <- res_ai[,1]

res_ai <- calc.ai_new(seg = preprocessed_HRP_TCGA, chrominfo = chrominfo_grch38)
nTAIs_HRP_GRCh38_removedsmall <- res_ai[,1]


res_ai<- calc.ai_new(seg = preprocessed_HRD_TCGA, chrominfo = chrominfo_grch38, shrink = FALSE, min.size=50)
nTAIs_HRD_GRCh38_notremoved <- res_ai[,1]

res_ai<- calc.ai_new(seg = preprocessed_HRP_TCGA, chrominfo = chrominfo_grch38, shrink = FALSE, min.size=50)
nTAIs_HRP_GRCh38_notremoved <- res_ai[,1]

#############################################################################################################################
######################################## Improvements for nTAIs #############################################################

counter_events <- function(dataf){
  samples <- unique(dataf$SampleID)
  count_selected <- sapply(samples, function(x){nrow(dataf[dataf$SampleID %in% x,])})
  return(count_selected)
}
#wilcox.test(counter_events(res_ai_HRD_TAIs), counter_events(res_ai_HRP_TAIs), alternative = "greater")

#Test for the next minimum AIs sizes the BA and mean difference between HRD and HRP samples
sizes <- c(0:20)
sizes_exp <- sizes * 1e6 #Convert to Mb

#Calculus of nTAIs in the GRCh38 datasets
ACC_bySegments <- NULL
pvalues <- NULL
for (s in sizes_exp){
  print(s)
  nTAIs_HRD <- calc.ai_new(seg = preprocessed_HRD_TCGA, chrominfo = chrominfo_grch38, min.size = s)[,1]
  nTAIs_HRP <- calc.ai_new(seg = preprocessed_HRP_TCGA, chrominfo = chrominfo_grch38, min.size = s)[,1]
  Utest <- wilcox.test(nTAIs_HRD, nTAIs_HRP, alternative = "greater")
  pval <- Utest$p.value
  pvalues <- c(pvalues, pval)
  acc_iteration <- get_dispertion_acc(nTAIs_HRD, nTAIs_HRP)
  ACC_bySegments <- c(ACC_bySegments, acc_iteration)
}

values_nTAIs <- data.frame(pvals = pvalues, BA=ACC_bySegments, minSize=sizes, row=rep("", length(sizes)))
values_nTAIs$minSize <- as.factor(values_nTAIs$minSize)
values_nTAIs$pvals <- (-1) * log10(values_nTAIs$pvals) #Pvals to log10

#For Sup. Figure 2f, bottom panel
##Plot circles to show the accuracy and pvalues of each iteration
p <- ggplot(values_nTAIs, aes(x=minSize, y=row))
p <- p + geom_point(aes(color=pvals, size=BA)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue3", "darkcyan", "red"), values=c(0,0.90,1), limits=c(3.8,12.5))
p <- p + scale_radius(breaks = c(0.75, 0.80, 0.85, 0.90), limits=c(0.58,0.91), range = c(2,7))
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)),
               legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.3)),
               axis.title=element_text(size=rel(1.5)))
p <- p +  guides(colour = guide_colourbar(order = 1), size = guide_legend(keyheight  = 1))
p <- p + ylab("") + xlab("\nMinimun size of allelic imbalances (Mb)")
print(p)
ggsave(p, filename = paste0(outputfolder,"nTAIs-dotplotsHGSC.svg"), width = 20, height = 9, units = "cm")


#For Sup. Figure 2f, upper panel
##Plot to show the pvalues of each iteration
p <- ggplot(values_nTAIs, aes(x = sizes, y = pvals)) +   geom_point(size=2) + geom_line()
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)),
               legend.text=element_text(size=rel(1.3)),
               strip.text.x = element_text(size=rel(1.1)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(1.5)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + ylab("-log10(p.val)") + xlab("\nMinimun size of allelic imbalances (Mb)")
print(p)
ggsave(p, filename=paste0(outputfolder, "Pvals_nTAIs_minsizes.svg"), width = 20, height = 7, units = "cm")


x = 0.0
n = 0.0001
for (j in 1:10000){
  x = x +  n
}
x = x - 10
print(x)
