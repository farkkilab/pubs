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
outputfolder = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/" #Mainly for plots

###### Reading input segments #######
segHRDTCGA <- read.table(file = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/HRD_samples-TCGAsegmentsHGSC_2021.txt", header=T)
segHRPTCGA <- read.table(file = "/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/HRP_samples-TCGAsegmentsHGSC_2021.txt", header=T)

#Ploidy information
Ploidy_file <- read.table(file="/home/fernpere/HRD/TCGA_analysis/find_HRD_signatures/samplesHR_ploidy-purity2.txt ", header = T, sep="\t")
################################################## Preparing input ############################################################

##Preprocessing the segments from TCGA portal
chrominfo = chrominfo_grch38
ploidy <- rep(2, nrow(segHRDTCGA)) #Not really used this values, but necessary column in input
segHRDTCGA <- cbind (segHRDTCGA, ploidy)
ploidy <- rep(2, nrow(segHRPTCGA)) #Not really used this values, but necessary column in input
segHRPTCGA <- cbind (segHRPTCGA, ploidy)

preprocessed_HRD_TCGA <- preparing.input(segHRDTCGA)
preprocessed_HRP_TCGA <- preparing.input(segHRPTCGA)

####Make plot with the number of Allelic Imbalances per sample and pvalue, Sup.Figure2a
number_AI_HRD  <- table(preprocessed_HRD_TCGA$SampleID)
HRD_status <- rep("HRD", length(table(preprocessed_HRD_TCGA$SampleID)))
number_AI_HRP  <- table(preprocessed_HRP_TCGA$SampleID)
HRP_status <- rep("HRP", length(table(preprocessed_HRP_TCGA$SampleID)))

allelic_imbalances <- data.frame(status_HRD=c(HRD_status, HRP_status),
                 Alllelic_imbalances=c(number_AI_HRD, number_AI_HRP),
                 samples_id=c(names(table(preprocessed_HRD_TCGA$SampleID)), names(table(preprocessed_HRP_TCGA$SampleID))))

wilcox.test(number_AI_HRD, number_AI_HRP, alternative="greater", paired=FALSE)

#Sup.Figure2a
p <- ggplot(allelic_imbalances, aes (x=status_HRD , y=Alllelic_imbalances)) + geom_boxplot(alpha = 0.5, outlier.shape = NA)
p <- p + ylab("Number of allelic imbalances") + xlab("")
p <- p + ylim(0,670)
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.1)), axis.text.y=element_text(size=rel(1.1)), legend.text=element_text(size=rel(1.1)),
               strip.text.x = element_text(size=rel(1.1)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(1.3)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + geom_point(aes(y=Alllelic_imbalances, color=status_HRD), position= position_jitter(width= .3), size= 3, alpha = 0.7, show.legend = FALSE)
p <- p + scale_color_manual(values=c("#FA5A41", "#34A0D3"))
print(p)
ggsave(p, filename = paste0(outputfolder, "Raw-Allelic_imbalances-HRD_HRPs_HGSC.svg"), width = 9, height = 10, units = "cm")
ggsave(p, filename = paste0(outputfolder, "Raw-Allelic_imbalances-HRD_HRPs_HGSC.png"), width = 9, height = 10, units = "cm")

#### Remove chromosome LOH deletions ####
preprocessed_HRD_TCGA <- rm.chr.deletions(preprocessed_HRD_TCGA)
preprocessed_HRP_TCGA <- rm.chr.deletions(preprocessed_HRP_TCGA)

### Add the AIs size in the last column
sizes=(preprocessed_HRD_TCGA[,4] - preprocessed_HRD_TCGA[,3])
preprocessed_HRD_TCGA <- cbind(preprocessed_HRD_TCGA, sizes)

sizes=(preprocessed_HRP_TCGA[,4] - preprocessed_HRP_TCGA[,3])
preprocessed_HRP_TCGA <- cbind(preprocessed_HRP_TCGA, sizes)

### Sample TCGA-13-1511 (HRP), is an HRP outlier, we excluded it from subsequent analysis###
preprocessed_HRP_TCGA <- preprocessed_HRP_TCGA[!preprocessed_HRP_TCGA[,1] %in%  c("TCGA-13-1511"),]

##### Extract LOH events ######
preprocessed_HRD_TCGA <- preprocessed_HRD_TCGA[preprocessed_HRD_TCGA[,8] == 0 & preprocessed_HRD_TCGA[,7] != 0,,drop=F]
preprocessed_HRP_TCGA <- preprocessed_HRP_TCGA[preprocessed_HRP_TCGA[,8] == 0 & preprocessed_HRP_TCGA[,7] != 0,,drop=F]

sizes.HRP <- sapply(unique(preprocessed_HRP_TCGA$SampleID), function(x){
            sample.sizes <- preprocessed_HRP_TCGA[preprocessed_HRP_TCGA$SampleID %in% x,"sizes"]
            median(sample.sizes)/1e6
})

sizes.HRD <- sapply(unique(preprocessed_HRD_TCGA$SampleID), function(x){
                sample.sizes <- preprocessed_HRD_TCGA[preprocessed_HRD_TCGA$SampleID %in% x,"sizes"]
                median(sample.sizes)/1e6
})

HRP_status <- rep("HRP", length(unique(preprocessed_HRP_TCGA$SampleID)))

sizes.AIs <- data.frame(median.size=c(sizes.HRP,sizes.HRD), HR.status=c(HRP_status, HRD_status))

wilcox.test(sizes.HRP, sizes.HRD, alternative="greater", paired=FALSE)


p <- ggplot(sizes.AIs, aes (x=HR.status , y=median.size)) + geom_boxplot(alpha = 0.5, outlier.shape = NA)
p <- p + ylab("Median AIs length (Mb)") + xlab("")
#p <- p + ylim(0,670)
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.1)), axis.text.y=element_text(size=rel(1.1)), legend.text=element_text(size=rel(1.1)),
               strip.text.x = element_text(size=rel(1.1)),
               legend.title=element_text(size=rel(1.3)), axis.title=element_text(size=rel(1.3)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
p <- p + geom_point(aes(y=median.size, color=HR.status), position= position_jitter(width= .3), size= 3, alpha = 0.7, show.legend = FALSE)
p <- p + scale_color_manual(values=c("#FA5A41", "#34A0D3"))
print(p)

###################################################################################################################################
####################################### LOH size density inspection ###############################################################
###################################################################################################################################

### Plot LOH distributions per HR condition ###
segLOH1aux <-  preprocessed_HRD_TCGA
segLOH2aux <-  preprocessed_HRP_TCGA

segLOH1aux$Status <- rep("HRD", nrow(segLOH1aux))
segLOH2aux$Status <- rep("HRP", nrow(segLOH2aux))
segLOHconcat <- rbind(segLOH2aux, segLOH1aux)
segLOHconcat$sizes <- segLOHconcat$sizes/1e6


#Density comparision of LOH between HRD and HRP for Sup.Figure 2c
p <- ggplot(segLOHconcat, aes(x=sizes, fill=Status))
p <- p + geom_histogram(aes(y = ..density..), alpha=1, position="dodge", color="black", breaks=c(seq(0,80, by=2)))
p <- p + geom_density(alpha=0.4, outline.type="upper", kernel="biweight")
p <- p + scale_x_continuous(breaks=c(seq(0,80, by=10)), limits = c(-4,80), minor_breaks=c(seq(0,80, by=2)))
p <- p + theme_classic()
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.3)), axis.text.y=element_text(size=rel(1.3)),
               legend.text=element_text(size=rel(1.6)), legend.title=element_text(size=rel(1.6)),
               axis.title=element_text(size=rel(1.6)))
p <- p + guides(size = guide_legend(keyheight  = 1.2))
p <- p + ylab("Density \n") + xlab("\nSize of LOH (Mb)")
p <- p + scale_fill_discrete(labels=c("HRD", "HRP"))
p <- p + scale_fill_manual(values=c("#FA5A41", "#34A0D3"))
print(p)
ggsave(p, filename = paste0(outputfolder, "Density_LOH-events.svg"), width = 16, height = 14, units = "cm")
ggsave(p, filename = paste0(outputfolder, "Density_LOH-events.png"), width = 16, height = 14, units = "cm")

##################################################################################################################################################################
################################################## Inference of new cutoffs for HRD-LOH    #######################################################################
##################################################################################################################################################################

preprocessed_seg <- preprocessed_HRD_TCGA
preprocessed_seg2 <- preprocessed_HRP_TCGA

#Number of LoH events per window size
min_values <- c(0:20)
max_values <- c(25, 30, 40, 50, 400)
#400 as not limit, because there is not AIs bigger than 400Mbs (Max chr size is around 250Mb)

values.by.mindist <- NULL

#Iteration of possibles LOH sizes, calculus of BA and U.test p.values in each iteration
for (minval in min_values){
    BA_values <- NULL
    P_values <- NULL
    for (maxval in max_values){
        column <- paste(paste("Min", minval, sep=""), paste("Max", maxval, sep=""), sep="_")
        windows <-  c(minval,maxval)
        HRDs <- features.LOH(preprocessed_seg, windows)
        HRPs <- features.LOH(preprocessed_seg2, windows)
        BA <- get_dispertion_acc(HRDs[,1], HRPs[,1])
        BA_values <- c(BA_values, BA)
        test_u <- wilcox.test(HRDs[,1], HRPs[,1], alternative = "greater")
        pvals <- test_u$p.value
        P_values <- c(P_values, pvals)
    }

    Pvalues <- P_values

    j <- data.frame(pvals=Pvalues, BA=BA_values,
                    mdist_Mb = as.factor(as.numeric(c(rep(minval,length(P_values))))),
                    max_dist = as.factor(max_values))

     values.by.mindist <- rbind(values.by.mindist, j)
}


values.by.mindist$logFDR <- (-1) * log10(values.by.mindist$pvals)
y_labesl=max_values
y_labesl[length(y_labesl)] <- "No limit"

##Plot circles to show the accuracy and pvalues of each iteration
##For Figure2c
p <- ggplot(values.by.mindist, aes(x=mdist_Mb, y=max_dist))
p <- p + geom_point(aes(color=logFDR, size=BA)) + theme_classic()
p <- p + scale_color_gradientn(name="-log10(p.val)", colours = c("blue3", "darkcyan", "red"), values=c(0,0.90,1), limits=c(7,11.7))
p <- p + scale_radius(breaks = c(0.75, 0.80, 0.85, 0.90), limits=c(0.59,0.91), range = c(2,7))
p <- p + scale_y_discrete(labels=y_labesl)
p <- p + theme(axis.text=element_text(size=rel(1.1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.1)), axis.text.y=element_text(size=rel(1.1)),
               legend.text=element_text(size=rel(1.1)), legend.title=element_text(size=rel(1.3)),
               axis.title=element_text(size=rel(1.3)))
p <- p +  guides(colour = guide_colourbar(order = 1), size = guide_legend(keyheight  = 1))
p <- p + ylab("Maximun size of LOH (Mb)\n") + xlab("\nMinimun size of LOH (Mb)")
print(p)
ggsave(p, filename = paste0(outputfolder, "LOH-windows_pvalues_dotplot_HGSC.svg"), width = 20, height = 9, units = "cm")
ggsave(p, filename = paste0(outputfolder, "LOH-windows_pvalues_dotplot_HGSC.png"), width = 20, height = 9, units = "cm")


#Best value
max <- which.max(values.by.mindist[,c(2)] * values.by.mindist[,c(5)])
values.by.mindist[max,]

#Best cutoff found is LOH events from  10Mb to 50Mb

#Telli2016 BA:
values.by.mindist[which(values.by.mindist$mdist_Mb == 15 & values.by.mindist$max_dist == 400),]


###############Calculate improvement of accuracy between the current optimization and Telli206 ##############
#Previous score in Telli2016
windows <-  c(15,500)
HRDs <- features.LOH(preprocessed_seg, windows)
HRPs <- features.LOH(preprocessed_seg2, windows)
HRDs$HRDstatus <- rep("HRD", nrow(HRDs))
HRPs$HRDstatus <- rep("HRP", nrow(HRPs))
prevHRDscore <- rbind(HRDs, HRPs)
prevHRDscore$Scars_definition <- rep("Previous", nrow(prevHRDscore))
names(prevHRDscore)[1] <- "HRD.LOH_scars"


#New HRD score for LOH
windows <-  c(10,50)
HRDs <- features.LOH(preprocessed_seg, windows)
HRPs <- features.LOH(preprocessed_seg2, windows)
HRDs$HRDstatus <- rep("HRD", nrow(HRDs))
HRPs$HRDstatus <- rep("HRP", nrow(HRPs))
newHRDscore <- rbind(HRDs, HRPs)
newHRDscore$Scars_definition <- rep("New proposal", nrow(prevHRDscore))
names(newHRDscore)[1] <- "HRD.LOH_scars"

HRD_LOH_scars_comparision <- rbind(prevHRDscore[,c(-2)], newHRDscore[,c(-2)])

#Not in manuscript, plot for comparing the HRD and HRPs abudance of LOH using both comparisons
p <- ggplot(HRD_LOH_scars_comparision, aes (x=HRDstatus, y=HRD.LOH_scars)) + geom_boxplot(alpha = 0.5)
p <- p + facet_wrap(~Scars_definition)
p <- p + ylab("HRD-LOH scars") + xlab("")
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
p <- p + geom_point(aes(y=HRD.LOH_scars, color=HRDstatus),    position= position_jitter(width= .3), size= 3, alpha = 0.7, show.legend = FALSE)
p <- p + ylim(0,40)
print(p)
ggsave(p, filename = paste0(outputfolder, "HRD-LOH-previous_new_scars_boxplot.svg"), width = 36, height = 20, units = "cm")
ggsave(p, filename = paste0(outputfolder, "HRD-LOH-previous_new_scars_boxplot.png"), width = 36, height = 20, units = "cm")