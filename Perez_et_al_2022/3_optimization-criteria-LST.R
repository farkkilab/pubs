library(ggplot2)
library(rpart)
library(rpart.utils)
library(rpart.plot)
library(ggpubr)


#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")



####### Defining variables ###########
outputfolder = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines//" #Mainly for plots

###### Reading input #######
#This are GRCh38 segments for OVA-TCGA (PanCanAtlas)
segHRDTCGA <- read.table(file = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/HRD_samples-TCGAsegmentsHGSC_2021.txt", header=T)
segHRPTCGA <- read.table(file = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/HRP_samples-TCGAsegmentsHGSC_2021.txt", header=T)

#Plody information
Ploidy_file <- read.table(file="/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/samplesHR_ploidy-purity2.txt", header = T, sep="\t")

##Preprocesing the segments from TCGA portal
chrominfo = chrominfo_grch38
ploidy <- rep(2, nrow(segHRDTCGA))
segHRDTCGA <- cbind (segHRDTCGA, ploidy)
ploidy <- rep(2, nrow(segHRPTCGA))
segHRPTCGA <- cbind (segHRPTCGA, ploidy)

preprocessed_HRD_TCGA <- preparing.input(segHRDTCGA)
preprocessed_HRP_TCGA <- preparing.input(segHRPTCGA)
### Sample TCGA-13-1511 (HRP), is an HRP outlier, we excluded###
preprocessed_HRP_TCGA <- preprocessed_HRP_TCGA[!preprocessed_HRP_TCGA[,1] %in%  c("TCGA-13-1511"),]

#########################################################################################################################
#########################################################################################################################
#Plotting distribution of segments as in Popova et al. 2012 Figure 2A
#Proportion of segments equal or greater than a given segment size

##### Density calculation of events of HRD ####
Sizes <- c(seq(0.25,80, by=0.25))
all_samplesTCGAportal <- rbind(preprocessed_HRD_TCGA, preprocessed_HRP_TCGA)

#Proportion of segments of a given size
Proportions_segments_HRD <- segmentBySamples(preprocessed_HRD_TCGA, Sizes)
Proportions_segments_HRP <- segmentBySamples(preprocessed_HRP_TCGA, Sizes)
average_proportions_HRD <- log2(apply(Proportions_segments_HRD, FUN = mean, MARGIN = 1))
average_proportions_HRP <- log2(apply(Proportions_segments_HRP, FUN = mean, MARGIN = 1))

#Fitting a smoothing spline for the proportion of segments of given size
yHRD <- average_proportions_HRD[1:160]
xHRD <-  as.numeric(names(average_proportions_HRD)[1:160])
loHRD <- smooth.spline(xHRD, yHRD, spar=0.5)

yHRP <- average_proportions_HRP[1:160]
xHRP <-  as.numeric(names(average_proportions_HRP)[1:160])
loHRP <- smooth.spline(xHRP, yHRP, spar=0.5)

#Plotting differences in distribution of segments for Sup.Figure2b
svg(file=paste0(outputfolder, "HRP-HRD_segmentsProportions_average_HGSC.svg"), height = 6, width = 6, pointsize = 6)
plot(Sizes,average_proportions_HRD,xlim = c(0,28), ylim=c(-10.6,-4), pch=16, xaxt="n", main="",
     xlab="Segment size, Mb", ylab="Log2(Average proportions)", type="p", cex=0.8, col="darkred",
     cex.axis = 1.8, cex.lab=2.8)
points(Sizes,average_proportions_HRP, col="blue", pch=16, cex=0.8)
lines(loHRP,  col="blue")
lines(loHRD, col="darkred")
axis(side=1,at=seq(0,90,by=1), cex.axis=1.5)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
abline(v=2, lty=2)
dev.off()

#########################################################################################################################################################
################################################### Calculating best LSTs ###############################################################################
#########################################################################################################################################################

####### First generate vector with with size #########
chrominfo = chrominfo_grch38
LSTs_per_sampleTCGA_HRP <- LSTs(preprocessed_HRP_TCGA, chrominfo = chrominfo_grch38)
LSTs_per_sampleTCGA_HRD <- LSTs(preprocessed_HRD_TCGA, chrominfo = chrominfo_grch38)

####### Iteration of parameters that define an LST, get the number of LSTs under those parameters#########
mindistance_values <- c(1, 2, 3, 4) #In Mb, minimum AIs length for smoothing
segsizes_values <- c(1:20) #Mb minimum size of consecutive AIs after smoothing
tandemelements_Values = c(2,3) #Number of consecutive AIs after smoothing

chrominfo = chrominfo_grch38
#In the next pair of dataframes will be stored the amount of LSTs under the parameters s,m,t
#Just initialize them whit random numbers in the first two columns, letter will be removed the first first two columns
LSTs_all_parameters_HRP <- data.frame(x=rep(1,length(LSTs_per_sampleTCGA_HRP)), y=rep(1,length(LSTs_per_sampleTCGA_HRP)))
LSTs_all_parameters_HRD <- data.frame(x=rep(1,length(LSTs_per_sampleTCGA_HRD)), y=rep(1,length(LSTs_per_sampleTCGA_HRD)))
for (m in mindistance_values){
  m <- m * 1e6 #To Mb
  for (s in segsizes_values){
    s <- s * 1e6 #To Mb
    print(paste("distance=",m,"segzises=",s))
    for (t in tandemelements_Values){
        LSTs_per_sampleTCGA_HRP <- LSTs(preprocessed_HRP_TCGA, mindistance = m, segsizes = s, tandemelements = t)
        LSTs_all_parameters_HRP <- cbind(LSTs_all_parameters_HRP,LSTs_per_sampleTCGA_HRP)
        LSTs_per_sampleTCGA_HRD <- LSTs(preprocessed_HRD_TCGA, mindistance = m, segsizes = s, tandemelements = t)
        LSTs_all_parameters_HRD <- cbind(LSTs_all_parameters_HRD,LSTs_per_sampleTCGA_HRD)

    }
  }
}

#Generate column names
#Column names have the next structure: m1_t2_s5
#m means the minimum segment size for smoothing
#t number of tandem elements
#s minimum size of AI that is included in for a LSTs
column_names <- c("x","y")
for (m in mindistance_values){
  for (s in segsizes_values){
      for (t in tandemelements_Values){
      #Column names
      name <- paste(paste("m", m, sep=""), paste("t", t, sep=""), paste("s", s, sep=""), sep="_")
      print(name)
      column_names <- c(column_names, name)

    }
  }
}

######Add column and row names######
colnames(LSTs_all_parameters_HRP) <- column_names
LSTs_all_parameters_HRP <- LSTs_all_parameters_HRP[,-c(1,2)] #Removing random numbers
colnames(LSTs_all_parameters_HRD) <- column_names
LSTs_all_parameters_HRD <- LSTs_all_parameters_HRD[,-c(1,2)] #Removing random numbers

LSTs_all_parameters_HRD$samples <- row.names(LSTs_all_parameters_HRD)
LSTs_all_parameters_HRP$samples <- row.names(LSTs_all_parameters_HRP)

############ Plotting each distribution of segments in sample by HR status ###############
#Segments sizes per LSTs, considered with a minimum of 1Mb and two tandem segments
#No a plot in manuscript
m <- 1
t <- 2
average_proportions_HRD
columns_to_get <- NULL
for (s in segsizes_values){
    column <- paste(paste("m", m, sep=""), paste("t", t, sep=""), paste("s", s, sep=""), sep="_")
    columns_to_get <- c(columns_to_get, column)
}
for (sample in 1:length(row.names(LSTs_all_parameters_HRD))){
  if (sample == 1){
    values_sample <- LSTs_all_parameters_HRD[sample,columns_to_get]
    plot(segsizes_values, values_sample, type="l", lwd=1.2, ylim=c(4,86), xlim=c(5,20), pch=20, xlab="Segment size, Mb", ylab="State transitions")
  }else{
    values_sample <- LSTs_all_parameters_HRD[sample,columns_to_get]
    lines(segsizes_values, values_sample, type="l", lwd=1.2, pch=20)
  }

}
for (sample in 1:length(row.names(LSTs_all_parameters_HRP))){
    values_sample <- LSTs_all_parameters_HRP[sample,columns_to_get]
    lines(segsizes_values, values_sample, type="l", lwd=1.2, pch=20, col="blue")
}
legend(15,80, c("HRD","HRP"), col = c("black","blue"), lty = 1)

################# Calculus of p values and accuracy for the difference in abundance of LSTs between HRD and HRP########
######Function to use ############
calculate.pvalues.BA <- function(abundances.LSTs.HRD, abundances.LSTs.HRP, min.distances=mindistance_values, segments.sizes=segsizes_values, tandem.elements=2){
  values_by_mindist <- NULL
  for (m in min.distances){
    Pvalues_bySegments <- NULL
    ACC_bySegments <- NULL
    for (s in segments.sizes){
      column <- paste(paste("m", m, sep=""), paste("t", tandem.elements, sep=""), paste("s", s, sep=""), sep="_")
      print(column)
      all.pvals <- NULL
      all.BA <- NULL
      #for (i in 1:1000){ #This for loop was used for jacknifing
        #HRD.samples.extract <- nrow(abundances.LSTs.HRD)-1
        #HRP.samples.extract <- nrow(abundances.LSTs.HRP) - 1
        HRP.samples.extract <- nrow(abundances.LSTs.HRP)
        HRD.samples.extract <- nrow(abundances.LSTs.HRD)
        abundances.LSTs.HRD.sample <- sample(abundances.LSTs.HRD[,column], HRD.samples.extract, replace = FALSE) #This was used for jacknifing sampling
        abundances.LSTs.HRP.sample <- sample(abundances.LSTs.HRP[,column], HRP.samples.extract, replace = FALSE) #This was used for jacknifing sampling
        test_u <- wilcox.test(abundances.LSTs.HRD.sample, abundances.LSTs.HRP.sample, alternative = "greater")
        pvals <- test_u$p.value
        all.pvals <- c(all.pvals, pvals)
        acc_iteration <- get_dispertion_acc(abundances.LSTs.HRD.sample, abundances.LSTs.HRP.sample)
        all.BA <- c(all.BA, acc_iteration)
      #}
      Pvalues_bySegments <- c(Pvalues_bySegments, mean(all.pvals))
      ACC_bySegments <- c(ACC_bySegments, mean(all.BA))
    }

    j <- data.frame(pvalues=Pvalues_bySegments, BA=ACC_bySegments,
                    S_Mb = segments.sizes,
                    mdist_Mb = as.numeric(c(rep(m,length(Pvalues_bySegments)))),
                    values_by_tandem = rep(tandem.elements,length(Pvalues_bySegments)))

    values_by_mindist <- rbind(values_by_mindist, j)
  }
  values_by_mindist$logFDR <- (-1) * log10(values_by_mindist$pvalues)
  values_by_mindist$mdist_Mb <- as.factor(values_by_mindist$mdist_Mb)
  values_by_mindist$S_Mb <- as.factor(values_by_mindist$S_Mb)
  return(values_by_mindist)
}

#For each combination of parameters m l t#
BA.pvals.LSTs.t2 <- calculate.pvalues.BA(LSTs_all_parameters_HRD, LSTs_all_parameters_HRP, tandem.elements=2)
BA.pvals.LSTs.t3 <- calculate.pvalues.BA(LSTs_all_parameters_HRD, LSTs_all_parameters_HRP, tandem.elements=3, segments.sizes=c(1:12))

#Plotting accuracy and p.values of difference in abundance between HRD and HRP for two tandem AIs
#For Sup.Figure 2d
p <- ggplot(BA.pvals.LSTs.t2, aes(x=S_Mb, y=mdist_Mb))
p <- p + geom_point(aes(color=logFDR, size=BA)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue3", "darkcyan", "red"), values=c(0,0.90,1), limits=c(3.8,12.5))
p <- p + scale_radius(breaks = c(0.75, 0.80, 0.85, 0.90), limits=c(0.58,0.91), range = c(2,7))
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)),
               legend.text=element_text(size=rel(1.1)), legend.title=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.5)))
p <- p +  guides(colour = guide_colourbar(order = 1), size = guide_legend(keyheight  = 1))
p <- p + ylab("Space between AIs (Mb)\n") + xlab("\nMinimun size of AIs (Mb)")
print(p)
ggsave(p, filename = paste0(outputfolder, "LSTs-windows_pvalues_dotplot_segmentsHGSC_t2.svg"), width = 20, height = 9, units = "cm")
ggsave(p, filename = paste0(outputfolder, "LSTs-windows_pvalues_dotplot_segmentsHGSC_t2.png"), width = 20, height = 9, units = "cm")

#Plotting accuracy and p.values of difference in abundance between HRD and HRP  for three tandem AIs
#For Sup.Figure 2e
p <- ggplot(BA.pvals.LSTs.t3, aes(x=S_Mb, y=mdist_Mb))
p <- p + geom_point(aes(color=logFDR, size=BA)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue3", "darkcyan", "red"), values=c(0,0.90,1), limits=c(3.8,12.5))
p <- p + scale_radius(breaks = c(0.75, 0.80, 0.85, 0.90), limits=c(0.58,0.91), range = c(2,7))
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)),
               legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.3)),
               axis.title=element_text(size=rel(1.5)))
p <- p +  guides(colour = guide_colourbar(order = 1), size = guide_legend(keyheight  = 1))
p <- p + ylab("Space between AIs (Mb)\n") + xlab("\nMinimun size of AIs (Mb)")
print(p)
ggsave(p, filename = paste0(outputfolder, "LSTs-windows_pvalues_dotplot_segmentsHGSC_t3.svg"), width = 20, height = 9, units = "cm")
ggsave(p, filename = paste0(outputfolder, "LSTs-windows_pvalues_dotplot_segmentsHGSC_t3.png"), width = 20, height = 9, units = "cm")


#Best values for distinguishing LSTs are segments longer than 9Mb and a space (smoothing) of 1Mb
BA.pvals.LSTs.t2[which.max(BA.pvals.LSTs.t2$BA * BA.pvals.LSTs.t2$logFDR),]

#Printing the accuracy of the Telli2016 for separating HRD or HRPcolum
BA.pvals.LSTs.t2[BA.pvals.LSTs.t2$S_Mb == "10" & BA.pvals.LSTs.t2$mdist_Mb == "3",]

#Printing the accuracy of the of distinguished best cutoffs (values)
BA.pvals.LSTs.t2[BA.pvals.LSTs.t2$S_Mb == "12" & BA.pvals.LSTs.t2$mdist_Mb == "1",]

####### Boxplots to compare the separation between HRP and HRD samples using the new cutoff ##########
#The best value before reported was "m3_t2_s10"
status1 <- rep("HRD", length(LSTs_all_parameters_HRD[,"m3_t2_s10"]))
status2 <- rep("HRP", length(LSTs_all_parameters_HRP[,"m3_t2_s10"]))
Scars_definition <- rep("Previous", length(c(LSTs_all_parameters_HRD[,"m3_t2_s10"], LSTs_all_parameters_HRP[,"m3_t2_s10"])))
previosLSTS <- data.frame(HRDstatus=c(status1, status2), LSTs=c(LSTs_all_parameters_HRD[,"m3_t2_s10"], LSTs_all_parameters_HRP[,"m3_t2_s10"]), Scars_def=Scars_definition)

#The best value here identified was "m1_t2_s11"
Scars_definition2 <- rep("New proposal", length(c(LSTs_all_parameters_HRD[,"m1_t2_s9"], LSTs_all_parameters_HRP[,"m1_t2_s9"])))
newLSTS <- data.frame(HRDstatus=c(status1, status2), LSTs=c(LSTs_all_parameters_HRD[,"m1_t2_s9"], LSTs_all_parameters_HRP[,"m1_t2_s9"]), Scars_def=Scars_definition2)

LSTs_scars_comparision <- rbind(previosLSTS, newLSTS)


p <- ggplot(LSTs_scars_comparision, aes (x=HRDstatus, y=LSTs)) + geom_boxplot(alpha = 0.5, outlier.shape = NA)
p <- p + facet_wrap(~Scars_def)
p <- p + ylab("LSTs scars") + xlab("")
p <- p + theme(axis.text=element_text(size=rel(1.5)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)), legend.text=element_text(size=rel(1.4)),
               strip.text.x = element_text(size=rel(2.8)),
               legend.title=element_text(size=rel(2)), axis.title=element_text(size=rel(2.5)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + geom_point(aes(y=LSTs, color=HRDstatus), position= position_jitter(width= .3), size= 3, alpha = 0.7, show.legend = FALSE)
p <- p + ylim(0,50)
print(p)
ggsave(p, filename = "/home/fernpere/HRD/Figures/LSTs-previous_new_scars_boxplot.svg", width = 36, height = 20, units = "cm")
ggsave(p, filename = "/home/fernpere/HRD/Figures/LSTs-previous_new_scars_boxplot.png", width = 36, height = 20, units = "cm")


########################### Get cutoff for that separe LSTs for HRPs and HRDs ###############################
LSTs_per_sampleTCGA_HRP <- LSTs(preprocessed_HRP_TCGA, segsizes=9e6, mindistance=1e6, tandemelements=2)
LSTs_per_sampleTCGA_HRD <- LSTs(preprocessed_HRD_TCGA, segsizes=9e6, mindistance=1e6, tandemelements=2)

dA <- data.frame(LSTs=LSTs_per_sampleTCGA_HRD, anyvalue=rep(2,length(LSTs_per_sampleTCGA_HRD)), status=rep("HRD",length(LSTs_per_sampleTCGA_HRD)))
dB <- data.frame(LSTs=LSTs_per_sampleTCGA_HRP, anyvalue=rep(2,length(LSTs_per_sampleTCGA_HRP)), status=rep("HRP",length(LSTs_per_sampleTCGA_HRP)))
dz <- rbind(dA,dB)

fit <- rpart(status ~ ., data=dz, method='class', control=rpart.control(minsplit = 2, minbucket = 1, cp=0.001))

#Identify cutoff
value <- rpart.subrules.table(fit)[1,5]
cutoff <- as.numeric(value)
cutoff
#Cutoff is 13.5
