library(ggplot2)
library(reshape)
library(ggridges)
library(parallel)
library(moments)
library(heatmaply)
require(gplots)
library(Matching)
library(kSamples)

#Importing other function inside the repository
extra.functions <- list.files("R/", full.names = TRUE)
sapply(extra.functions, source)
load("sysdata.rda")

##############################################################################################
######################## Declaring some functions ############################################
##############################################################################################

stats.by.tiss <- function(data, group.by="tissue", sampleid="SampleID"){
      group.vars <- sort(unique(data[,group.by]))
      skew.vals <- sapply(group.vars, function(x){
        sizes.x <- data[data[,group.by] == x,"sizes"]
        skew <- skewness(sizes.x)
      })
      median.counts <- sapply(group.vars, function(x){
        segs.x <- data[data[,group.by] == x,]
        counts.per.sample <- table(segs.x[,sampleid])
        #There are samples with size = 0, for those the count must be 0
        non.counts.samples <- segs.x[segs.x[,"sizes"] == 0,sampleid]
        counts.per.sample[names(counts.per.sample) %in% non.counts.samples] <- 0
        meadian.counts.group <- median(counts.per.sample)
        return(meadian.counts.group)
      })
      median.sizes <- sapply(group.vars, function(x){
        segs.x <- data[data[,group.by] == x,]
        samps <- unique(segs.x[,sampleid])
        median.size.vals <- sapply(samps, function(p){
                  median(segs.x[which(segs.x[,sampleid] == p),"sizes"])
        })
        median.sizes.tissue <- median(median.size.vals)
        return(median.sizes.tissue)
      })
      df.stats <- data.frame(skewness=skew.vals,
                             median.count=median.counts,
                             median.size=median.sizes,
                             cancer.type=group.vars)
      #Normalize by Z scores
      df.stats$skewness.z <- (df.stats$skewness - mean(df.stats$skewness))/sd(df.stats$skewness)
      df.stats$median.count.z <- (df.stats$median.count - mean(df.stats$median.count))/sd(df.stats$median.count)
      df.stats$median.size.z <- (df.stats$median.size - mean(df.stats$median.size))/sd(df.stats$median.size)
      return(df.stats)
}

#Count the number and median size of segments by sample, add final column with tissue
count.by.samp <- function(data, group.by="tissue", sampleid="SampleID"){
  group.vars <- sort(unique(data[,group.by]))
  number.segs <- lapply(group.vars, function(x){
    segs.x <- data[data[,group.by] == x,]
    counts.seg <- table(segs.x[,sampleid])
    median.size.vals <- sapply(names(counts.seg), function(p){
            median(segs.x[which(segs.x[,sampleid] == p),"sizes"])
    })
    cbind(counts.seg,median.size.vals)
  })
  df.counts <- NULL
  for (i in 1:length(number.segs)){
    df.aux <- data.frame(counts=number.segs[[i]][,1], median.size=number.segs[[i]][,2],
                         sample=rownames(number.segs[[i]]), tissue=rep(group.vars[i],length(number.segs[[i]][,1])))
    df.counts <- rbind(df.counts, df.aux)
  }
  return(df.counts)
}

############################################################################################################
#################################Identifiying BRCAmut samples and TNBCA ####################################
############################################################################################################

###First block of code is used to identify the IDs of TNBCA BRCA1/2 mutants and BRCA1/2wt
###There are two main type of IDs "0298819c-c447-44f3-9067-433801c5838a" and "TCGA-EW-A1P1"
##It was necessary to link some tables from TCGA to identify the samples with BRCA1/2 mutations reported by Knijnenburg et al 2018
##List of samples is further used

setwd("D:/users/fperez/HRD/TCGA_analysis/Segments_allCancers/")
output.folder <- "C:/Users/fernpere/HRD/Figures/Individual_figures_section_AIs/"

#Long sample ID `0298819c-c447-44f3-9067-433801c5838a` and tissue of origin
tcga_histology <- read.table(file="gdc_sample-tissue.txt", sep="\t", header=TRUE) #File made by hand using the filenames

#Sample sheet files downloaded from GDC, those correspond to all AI files in BRCA and OVA. Was added a column at the left
#With the file name as File.Name2
#The triple negative BRCA (TNBA) ids were taken from the GDC file: nationwidechildrens.org_clinical_patient_brca.txt
BRCA.sample.sheet <- read.table(file="gdc_sample_sheet.2021-06-04_2.tsv", sep="\t", header = TRUE)
OVA.sample.sheet <- read.table(file="gdc_sample_sheet.2021-06-02_2.tsv", sep="\t", header=TRUE)

#List of TNBCA samples, taken from Lehmann et al 2011
TNBCA <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/TNBC/12885_2020_6600_MOESM1_ESM_TNBC_ids.csv", sep=",", header=TRUE)


#########OVA HGSC samples, getting ID of BRCA1/2 mutants and BRCAwt###########
#List of HGSC samples used in the rest of the analysis. Next table has the manual annotation if sample is HRD or HRP and BRCA1/2 mutations
ovasamples <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/TCGA_HRD_newScars_clinical-info_HGSC.csv", sep=",", header=TRUE)
BRCA.germline <- grep("BRCA",ovasamples$Germline_mutations)
BRCA.alterations <- grep("BRCA",ovasamples$Biological_status)
ova.brcamut.samples <- ovasamples[unique(c(BRCA.germline, BRCA.alterations)),"Row.names"] #HGSC brcamut

#Only selected as BRCAwt those samples with full info
not.full.info.samples <- which(grepl("NA",ovasamples$Somatic_mutations) | grepl("NA",ovasamples$CNV_status) | grepl("NA",ovasamples$Hypermethylation_status))
BRCAwt.rows <- unique(c(BRCA.germline, BRCA.alterations, not.full.info.samples))
ova.brcawt.samples <- ovasamples[-BRCAwt.rows,"Row.names"]

########Get the samples in all TCGA with BRCA1/2 mutations/deletions/hypermethylation status according to the DDR-PanCan paper from Cell, 2018
#Next table from Knijnenburg et al 2018
DDR.status <- read.table(file="C:/Users/fernpere/HRD/TCGA_analysis/DDR/TCGA_DDR_Gene_summary.csv", sep =",", header = TRUE)

#Now identifying the ID's BRCA1/2 mutants
BRCA.status <- DDR.status[which(DDR.status$Gene.Symbol == "BRCA1" | DDR.status$Gene.Symbol == "BRCA2"),]
gene.names <- BRCA.status[,2]
BRCA.status <- BRCA.status[,-c(1,2)]
BRCA.status <- data.frame(t(BRCA.status))
colnames(BRCA.status) <- gene.names
BRCAmutants <- row.names(BRCA.status[which(BRCA.status$BRCA1 == 1 | BRCA.status$BRCA2 == 1),])
BRCAmutants <- gsub('[.]', '-', substr(BRCAmutants, 1, 12)) #This vector will be used in next block

#Samples with full info in the DDR, those are will be used to identify nonBRCAmut
BRCA1.info <- which(!is.na(BRCA.status$BRCA1))
BRCA2.info <- which(!is.na(BRCA.status$BRCA2))
ddr.full.info <- row.names(BRCA.status)[intersect(BRCA1.info, BRCA2.info)]
ddr.full.info <- substr(ddr.full.info, 1, 12)
ddr.full.info <- gsub('[.]', '-', ddr.full.info)

#####BRCAgermline mutations for all TCGA, extracted from PanCanAtlas info
BRCAgermlinemut <- scan(file="BRCAmut_BRCAsamples", what = "character")

######################## Get TNBCA samples long-IDs ###########################
#Merge with tissue name and sample names from TCGA for BRCA samples
hist.merged <- merge(tcga_histology, BRCA.sample.sheet, by.x = "sampleID", by.y = "File.Name2")
hist.merged$sampleName <- substr(hist.merged$Sample.ID, 1, 12)
#Select only Primary Tumors and TNBC
hist.merged <- hist.merged[which(grepl("Primary Tumor", hist.merged$Sample.Type)),]
TNBCA.ids <- merge(TNBCA, hist.merged, by.x="BARCODE", by.y="sampleName")

#Distinguish those TNBCA with BRCAmut or BRCAwt
TNBCA.ids.brcasomaticmut <- TNBCA.ids[which(TNBCA.ids$BARCODE %in% BRCAmutants), "sampleID"]
TNBCA.ids.brcagermut <- TNBCA.ids[which(TNBCA.ids$BARCODE %in% BRCAgermlinemut), "sampleID"]
TNBCA.ids.brcamut <- unique(c(TNBCA.ids.brcasomaticmut, TNBCA.ids.brcagermut))

nonBRCAmut <- TNBCA.ids$BARCODE[which(!(TNBCA.ids$sampleID %in% TNBCA.ids.brcamut))]
TNBCA.brcawt.longids <- nonBRCAmut[which(nonBRCAmut %in% ddr.full.info)]
TNBCA.brcawt <- TNBCA.ids[TNBCA.ids$BARCODE %in% TNBCA.brcawt.longids,"sampleID"]

###########Get BRCAmut OVA samples files ids#############
hist.merged2 <- merge(tcga_histology, OVA.sample.sheet, by.x = "sampleID", by.y = "File.Name2")
hist.merged2 <- hist.merged2[which(grepl("Primary Tumor", hist.merged2$Sample.Type)),]
hist.merged2$sampleName <- substr(hist.merged2$Sample.ID, 1, 12)

non.hgsc.ova <- hist.merged2[!(hist.merged2$sampleName %in% ovasamples$Row.names),"sampleID"] #Used to ignore other OV in TCGA non HGSC in next block
OVA.hgsc.ids <- hist.merged2[hist.merged2$sampleName %in% ovasamples$Row.names,"sampleID"]

OVA.brca.mut <- hist.merged2[which(hist.merged2$sampleName %in% ova.brcamut.samples),"sampleID"]
OVA.brca.wt <- hist.merged2[which(hist.merged2$sampleName %in% ova.brcawt.samples),"sampleID"]


##################################################################################################
##################################################################################################
############################# Allelic imbalance identification ###################################
##################################################################################################
##################################################################################################
#Next block perform the Pan-cancer analysis of AIs included in Figure1 and Supp Figure1

#Reading input
tcga_segments_raw <- read.table(file="segments_TCGA.txt", sep="\t", header=TRUE)

tcga_segments <- tcga_segments_raw

tcga.segs <- data.frame(SampleID = tcga_segments[,1],
                        Chromosome = tcga_segments[,2],
                        Start_position = tcga_segments[,3],
                        End_position = tcga_segments[,4],
                        total_cn = tcga_segments[,5],
                        A_cn  = tcga_segments[,6],
                        B_cn = tcga_segments[,7],
                        ploidy=rep(2,nrow(tcga_segments)))

tcga.segs <- preparing.input(tcga.segs)
tcga.segs <- rm.chr.deletions(tcga.segs)
sizes=(tcga.segs[,4] - tcga.segs[,3])
sizes <- sizes/1e6
tcga.segs$sizes <- sizes

###################################################################################################
################################## Plots for total AIs ############################################
###################################################################################################

#Selecting AIs longer than 3Mb and shorter than 50Mbs
tcga.segs.set  <- tcga.segs[tcga.segs$sizes > 3 & tcga.segs$sizes < 50, ]

#Adding the cancer type in the last column
tcga.segs.info <- merge(tcga.segs.set, tcga_histology, by.x="SampleID", by.y="sampleID")

#Selecting cancer types (tissues) with more than 255 samples
tissues.counts <- sort(table(tcga_histology$tissue))
tissues.selected <- names(tissues.counts[which(tissues.counts > 255)])
tcga.segs.info.set <- tcga.segs.info[tcga.segs.info$tissue %in% tissues.selected,]

#Excluding non-HGSC from OV
tcga.segs.info.set<- tcga.segs.info.set[which(!tcga.segs.info.set$SampleID %in% non.hgsc.ova),]

#Selecting only LOH events for some analysis
tcga.loh.info.set <- tcga.segs.info.set[tcga.segs.info.set$B_cn == 0,]

###################################################################################################
################################## Heatmaps for total AIs and LOH #################################
###################################################################################################

#Counting the number of AIs by sample
l.loh <- count.by.samp(tcga.loh.info.set)
j.AIs <- count.by.samp(tcga.segs.info.set)

#Order of tissues for the Figures
order.tissues <- rev(c("BLCA", "STAD", "LUSC", "BRCA", "LUAD",
                        "OV", "THCA", "KIRP", "KIRC", "LGG" , "UCEC", "LIHC", "CESC" ,
                       "COAD", "PRAD", "HNSC", "SKCM", "GBM"))

order.tissues2 <- rev(c("OV", "LUSC", "BRCA", "BLCA", "LUAD", "STAD",
                        "HNSC", "CESC", "LIHC", "SKCM", "COAD",  "LGG" , "GBM",
                        "PRAD", "UCEC", "KIRC", "THCA", "KIRP"))
l.loh$tissue<- factor(l.loh$tissue, levels=order.tissues2)
j.AIs$tissue<- factor(j.AIs$tissue, levels=order.tissues)

###There are some samples without LOH, for those will be added 0 counts and median size of 0
non.loh.samples <- j.AIs[!(j.AIs$sample %in% l.loh$sample),]
non.loh.samples[,c(1,2)] <- 0
l.loh <- rbind(l.loh, non.loh.samples)

#Barplot with total number of samples per cancer type
#For Figure 1a
p <- ggplot(j.AIs, aes(x=tissue)) + geom_bar(aes(y = (..count..)), fill="cadetblue4") +
      coord_flip() + xlab("") + ylab("Number of patients")
p <- p + theme(text=element_text(size=9,  family="Arial"),
               panel.border = element_blank(),
               axis.text.x=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.5)),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, file=paste0(output.folder, "Barplots_AIs_PanCan.svg"), width = 3, height = 10, units = "cm")
dev.off()

#Horizontal boxplot, number of AIs per type of cancer
#For Figure 1b
p <- ggplot(j.AIs, aes(x=tissue, y=counts)) + geom_boxplot(fill="antiquewhite3", outlier.shape = NA) +
  coord_flip() + xlab("") + ylab("Number of AIs per sample") + ylim(0,250)
p <- p + theme(text=element_text(size=9,  family="Arial"),
               panel.border = element_blank(),
               axis.text.x=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, file=paste0(output.folder, "Boxplots_AIs-counts_PanCan.svg"), width = 6, height = 8, units = "cm")


#Horizontal boxplots, median size per type of cancer
#For Figure 1c
p <- ggplot(j.AIs, aes(x=tissue, y=median.size)) + geom_boxplot(fill="antiquewhite3", outlier.shape = NA) +
  coord_flip() + xlab("Type of cancer") + ylab("Median AIs length per sample (Mb)") + ylim(0,40)
p <- p + theme(text=element_text(size=9,  family="Arial"),
               panel.border = element_blank(),
               axis.text.x=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, file=paste0(output.folder, "Boxplots_AIs-lenght_PanCan.svg"), width = 6, height = 8, units = "cm")

#Horizontal boxplot, number of LOH per type of cancer
#For Supplementary Figure 1a
p <- ggplot(l.loh, aes(x=tissue, y=counts)) + geom_boxplot(fill="antiquewhite3", outlier.shape = NA) +
  coord_flip() + xlab("") + ylab("Number of LOH per sample") + ylim(0,75)
p <- p + theme(text=element_text(size=9,  family="Arial"),
               panel.border = element_blank(),
               axis.text.x=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, file=paste0(output.folder, "Boxplots_LOH-counts_PanCan.svg"), width = 6, height = 8, units = "cm")


#Horizontal boxplots, median size per type of cancer
#For Supplementary Figure 1b
p <- ggplot(l.loh, aes(x=tissue, y=median.size)) + geom_boxplot(fill="antiquewhite3", outlier.shape = NA) +
  coord_flip() + xlab("Type of cancer") + ylab("Median LOH length per sample (Mb)") + ylim(0,40)
p <- p + theme(text=element_text(size=9,  family="Arial"),
               panel.border = element_blank(),
               axis.text.x=element_text(size=rel(1.2)),
               axis.title=element_text(size=rel(1.2)),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, file=paste0(output.folder, "Boxplots_LOH-lenght_PanCan.svg"), width = 6, height = 8, units = "cm")


####Heatmaps with AIs median sizes and skewness

#Including for the calculus of skewness that the samples with 0 LOH, had a length of 0 for the LOH
tcga.non.loh.info <- data.frame(SampleID=non.loh.samples$sample, Chromosome=rep("chr0", nrow(non.loh.samples)),
                         Start_position=rep(0, nrow(non.loh.samples)), End_position=rep(0, nrow(non.loh.samples)),
                         total_cn=rep(0, nrow(non.loh.samples)), A_cn=rep(0, nrow(non.loh.samples)),
                         A_cn.1=rep(0, nrow(non.loh.samples)), B_cn=rep(0, nrow(non.loh.samples)),
                         ploidy=rep(0, nrow(non.loh.samples)), V10=rep(0, nrow(non.loh.samples)),
                         sizes=rep(0, nrow(non.loh.samples)), tissue=non.loh.samples$tissue)

tcga.loh.info.set <- rbind(tcga.loh.info.set, tcga.non.loh.info)

#Calculating of stats for LOH
stats.loh <- stats.by.tiss(tcga.loh.info.set)
colnames(stats.loh) <- c("LOH.swkewness.raw","LOH.med.count.raw",
                         "LOH.med.size.raw", "cancer.type", "LOH.swkewness", "LOH.med.count", "LOH.med.size")

#Calculating of stats for all AIs
stats.segments <- stats.by.tiss(tcga.segs.info.set)
colnames(stats.segments) <- c("Seg.swkewness.raw","Seg.med.count.raw",
                              "Seg.med.size.raw", "cancer.type", "Seg.swkewness", "Seg.med.count", "Seg.med.size")


#Heatmap for Figure 1d
#Then reduce lwid to 2, to reduce dendogram size
par(cex.main=1, cex.lab=0.7, cex.axis=0.7)
stats.segs.m <- as.matrix(stats.segments[,c(5,6,7)])
svg(file=paste0(output.folder, "Matrix_AIs_stats_PanCan.svg"), height = 7, width = 3.5, family = "Arial")
heatmap.2(stats.segs.m, trace="none", key=TRUE, dendrogram = "row", Colv=FALSE, cexCol = 0.6, cexRow = 1, density.info="none",
          keysize = 1.5, lwid=c(2,3.5), key.par = list(cex=0.5), notecol="black",
          col =  colorRampPalette(c("darkblue","white","darkred"))(100))
dev.off()
dev.off()

#Heatmap for Supplementary Figure 1c
#Then reduce lwid to 2, to reduce dendogram size
par(cex.main=1, cex.lab=0.7, cex.axis=0.7)
stats.segs.m <- as.matrix(stats.loh[,c(5,6,7)])
svg(file=paste0(output.folder, "Matrix_LOH_stats_PanCan.svg"), height = 5.4, width = 4, family = "Arial")
heatmap.2(stats.segs.m, trace="none", key=TRUE, dendrogram = "row", Colv=FALSE, cexCol = 0.6, cexRow = 1, density.info="none",
          keysize = 0.8, lwid=c(2,2), key.par = list(cex=0.5), notecol="black",
          col =  colorRampPalette(c("darkblue","white","darkred"))(100))
dev.off()
dev.off()

###########################################################################################################################
###########################################################################################################################
################# OVA and TNBC, dot plots. calculating wilcoxon and foldchange difference in medians ######################
###########################################################################################################################
###########################################################################################################################

#Next block of code correspond to the calculation of AIs differences between HGSC and TNBC

#List of IDs in groups according to labels.list
#First element are all samples, second only the BRCAmutatns and third the BRCAwt
list.samples <- list(c(TNBCA.ids$sampleID, OVA.hgsc.ids),
                     c(TNBCA.ids.brcamut, OVA.brca.mut),
                     c(TNBCA.brcawt, OVA.brca.wt))
labels.list <- c("OVA-TNBCall", "OVA-TNBC-BRCAmut","OVA-TNBC-nonBRCA")

#Selecting TNBCA and OVA AIs and LOH
j.AIs.OVA.TNBCA.all <- j.AIs[j.AIs$sample %in% TNBCA.ids$sampleID | j.AIs$sample %in% OVA.hgsc.ids,]
LOH.OVA.TNBCA.all.pre <- l.loh[l.loh$sample %in% TNBCA.ids$sampleID | l.loh$sample %in% OVA.hgsc.ids,]

#There are some samples without LOH events, adding those with a 0 count and 0 median size
LOH.absent.TNBCA <- TNBCA.ids$sampleID[which(!(TNBCA.ids$sampleID %in% LOH.OVA.TNBCA.all.pre$sample))]
LOH.absent.OV <- OVA.hgsc.ids[which(!(OVA.hgsc.ids %in% LOH.OVA.TNBCA.all.pre$sample))]
LOH.absent.df.TNBCA <- data.frame(sample=LOH.absent.TNBCA, tissue=rep("BRCA", length(LOH.absent.TNBCA)),
                                  counts=rep(0, length(LOH.absent.TNBCA)), median.size=rep(0, length(LOH.absent.TNBCA)))
LOH.absent.df.OV <- data.frame(sample=LOH.absent.OV, tissue=rep("OV", length(LOH.absent.OV)),
                               counts=rep(0, length(LOH.absent.OV)), median.size=rep(0, length(LOH.absent.OV)))
LOH.OVA.TNBCA.all <- rbind(LOH.OVA.TNBCA.all.pre, LOH.absent.df.TNBCA, LOH.absent.df.OV)
LOH.OVA.TNBCA.all$tissue <- factor(LOH.OVA.TNBCA.all$tissue, levels=c("OV", "BRCA"))

#Comparisons results in U test
results.comparisons <- NULL
#Next for loop will produce the violin plots for Figure 1 between TNBC and OVA, using the samples selected in the list

for (set in 1:length(list.samples)){
  set.samples <- list.samples[[set]]
  label.set <- labels.list[set]
  j.AIs.OVA.TNBCA.set <- j.AIs[which(j.AIs$sample %in% set.samples), ]
  LOH.OVA.TNBCA.set <- LOH.OVA.TNBCA.all[which(LOH.OVA.TNBCA.all$sample %in% set.samples), ]

  #OV and BRCA samples, count difference AIs
  p <- ggplot(j.AIs.OVA.TNBCA.set , aes(x=tissue, y=counts)) +
    xlab("") + ylab("Number of AIs per sample") + geom_violin() +
    geom_boxplot(outlier.shape = NA, width=0.28) + ylim(0,280)
  p <- p + theme(text=element_text(size=9,  family="Arial"),
                 axis.text.x=element_text(size=rel(1.2)),
                 axis.title=element_text(size=rel(1.2)),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, size=1))
  #print(p)
  ggsave(p, file=paste0(output.folder, "Boxplots_", label.set, "_counts.svg"), width = 6, height = 8, units = "cm")

  #OV and BRCA samples, count difference LOH
  p <- ggplot(LOH.OVA.TNBCA.set, aes(x=tissue, y=counts)) +
    xlab("") + ylab("Number of LOH per sample") + geom_violin() +
    geom_boxplot(outlier.shape = NA, width=0.28) + ylim(0,130)
  p <- p + theme(text=element_text(size=9,  family="Arial"),
                 axis.text.x=element_text(size=rel(1.2)),
                 axis.title=element_text(size=rel(1.2)),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, size=1))
  #print(p)
  ggsave(p, file=paste0(output.folder, "Boxplots_", label.set, "_counts_LOH.svg"), width = 6, height = 8, units = "cm")

  #OV and BRCA samples, length difference AIs
  p <- ggplot(j.AIs.OVA.TNBCA.set , aes(x=tissue, y=median.size)) +
    xlab("") + ylab("Median length of AIs per sample (Mb)") + geom_violin() +
    geom_boxplot(outlier.shape = NA, width=0.28) + ylim(5, 22)
  p <- p + theme(text=element_text(size=9,  family="Arial"),
                 axis.text.x=element_text(size=rel(1.2)),
                 axis.title=element_text(size=rel(1.2)),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, size=1))
  #print(p)
  #ggsave(p, file=paste0(output.folder, "Boxplots_", label.set, "_lengths.svg"), width = 6, height = 8, units = "cm")

  #OV and BRCA samples, length difference LOH
  p <- ggplot(LOH.OVA.TNBCA.set , aes(x=tissue, y=median.size)) +
    xlab("") + ylab("Median length of AIs per sample (Mb)") + geom_violin() +
    geom_boxplot(outlier.shape = NA, width=0.28) + ylim(5, 22)
  p <- p + theme(text=element_text(size=9,  family="Arial"),
                 axis.text.x=element_text(size=rel(1.2)),
                 axis.title=element_text(size=rel(1.2)),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, size=1))
  #print(p)
  #ggsave(p, file=paste0(output.folder, "Boxplots_", label.set, "_lengths_LOH.svg"), width = 6, height = 8, units = "cm")

  w.test1 <- wilcox.test(j.AIs.OVA.TNBCA.set[j.AIs.OVA.TNBCA.set$tissue == "OV","counts"],
                         j.AIs.OVA.TNBCA.set[j.AIs.OVA.TNBCA.set$tissue == "BRCA","counts"])

  w.test2 <- wilcox.test(LOH.OVA.TNBCA.set[LOH.OVA.TNBCA.set$tissue == "OV","counts"],
                         LOH.OVA.TNBCA.set[LOH.OVA.TNBCA.set$tissue == "BRCA","counts"])

  df.values <- data.frame(pval.diff=c(w.test1$p.value, w.test2$p.value),
                          measure=c("TotalAIs", "LOH"), class=c(label.set,label.set))
  results.comparisons <- rbind(results.comparisons, df.values)
}
print("P value for the comparisions:")
print(results.comparisons)

###################################################################################################
################################## By window difference between OV and TNBCA ######################
###################################################################################################
#Next block of code for producing the Figure 1h

#####Calculus of difference of abundance in AIs of given sizes. Between 3Mb to 10Mb for example
sizes <- c(3, 10, 16, 22, 30, 50)

#List of IDs in groups according to labels.list
list.samples2 <- list(TNBCA.ids$sampleID, OVA.hgsc.ids,
                     TNBCA.ids.brcamut, OVA.brca.mut,
                     TNBCA.brcawt, OVA.brca.wt)
labels.list2 <- c("All.samples", "All.samples", "BRCAmut","BRCAmut", "BRCAwt", "BRCAwt")

diff.window.all  <- NULL

#For each of the elements in the list.data will evaluate AIs of the corresponding sizes
#After foor loop will be generated data frame with fold changes and p-values of differences between HGSC and TNBC
for (d in seq(1,length(list.samples2),by=2)){
  typedata <- labels.list2[d]
  #From the raw segments, get the corresponding samples
  set.samples.tnbc <- list.samples2[[d]]
  set.samples.ova <- list.samples2[[d + 1]]
  data.AIs <- tcga.segs.info.set[which(tcga.segs.info.set$SampleID %in% c(set.samples.tnbc, set.samples.ova)), ]
  data.LOH <- tcga.loh.info.set[which(tcga.loh.info.set$SampleID %in% c(set.samples.tnbc, set.samples.ova)), ]

  for (j in 1:(length(sizes) - 1)){
    print(sizes[j])
    #Selecting AIs of the corresponding size
    selected.by.size.AIs <- data.AIs[data.AIs$sizes >= sizes[j] & data.AIs$sizes < sizes[j+1],]
    selected.by.size.loh <- data.LOH[data.LOH$sizes >= sizes[j] & data.LOH$sizes < sizes[j+1],]

    #Finding absent LOH samples, those with no LOH of the selected size
    absent.samples.tnbc <- set.samples.tnbc[which(!(set.samples.tnbc %in% selected.by.size.loh[which(selected.by.size.loh$tissue == "BRCA"),"SampleID"]))]
    absent.samples.ova <- set.samples.ova[which(!(set.samples.ova %in% selected.by.size.loh[which(selected.by.size.loh$tissue == "OV"),"SampleID"]))]

    #Finding absent AIs samples, those with no AIs of the selected size
    absent.samples.tnbc.ais <- set.samples.tnbc[which(!(set.samples.tnbc %in% selected.by.size.AIs[which(selected.by.size.AIs$tissue == "BRCA"),"SampleID"]))]
    absent.samples.ova.ais <- set.samples.ova[which(!(set.samples.ova %in% selected.by.size.AIs[which(selected.by.size.AIs$tissue == "OV"),"SampleID"]))]

    #Count the number of AIs per sample
    l.loh.set <- count.by.samp(selected.by.size.loh)
    j.AIs.set <- count.by.samp(selected.by.size.AIs)

    #Adding 0 to the samples without LOH of the given size
    if (length(absent.samples.ova) > 0){
        LOH.absent.df.OV <- data.frame(counts=rep(0, length(absent.samples.ova)), median.size=rep(0, length(absent.samples.ova)),
                                   sample=absent.samples.ova, tissue=rep("OV", length(absent.samples.ova)))
    }else{
        LOH.absent.df.OV <- NULL
    }
    if (length(absent.samples.tnbc) > 0){
        LOH.absent.df.TNBCA <- data.frame(sample=absent.samples.tnbc, tissue=rep("BRCA", length(absent.samples.tnbc)),
                                        counts=rep(0, length(absent.samples.tnbc)), median.size=rep(0, length(absent.samples.tnbc)))
    }else{
      LOH.absent.df.TNBCA <- NULL
    }
    l.loh.set <- rbind(l.loh.set, LOH.absent.df.TNBCA, LOH.absent.df.OV)

    #Adding 0 to the samples without AIs of the given size
    if (length(absent.samples.ova.ais) > 0){
      AIs.absent.df.OV <- data.frame(counts=rep(0, length(absent.samples.ova.ais)), median.size=rep(0, length(absent.samples.ova.ais)),
                                     sample=absent.samples.ova.ais, tissue=rep("OV", length(absent.samples.ova.ais)))
    }else{
      AIs.absent.df.OV <- NULL
    }
    if (length(absent.samples.tnbc.ais) > 0){
      AIs.absent.df.TNBCA <- data.frame(sample=absent.samples.tnbc.ais, tissue=rep("BRCA", length(absent.samples.tnbc.ais)),
                                        counts=rep(0, length(absent.samples.tnbc.ais)), median.size=rep(0, length(absent.samples.tnbc.ais)))
    }else{
      AIs.absent.df.TNBCA <- NULL
    }
    j.AIs.set <- rbind(j.AIs.set, AIs.absent.df.TNBCA, AIs.absent.df.OV)


    #Get FoldChanges for the abundance of LOH in HGSC compared to TNBC
    med.diff1 <- median(l.loh.set[l.loh.set$tissue == "OV","counts"]) /
      median(l.loh.set[l.loh.set$tissue == "BRCA","counts"])
    med.diff1 <- med.diff1 - 1

    #Get FoldChanges for the abundance of AIs in HGSC compared to TNBC
    med.diff2 <- median(j.AIs.set[j.AIs.set$tissue == "OV","counts"]) /
      median(j.AIs.set[j.AIs.set$tissue == "BRCA","counts"])
    med.diff2 <- med.diff2 - 1

    #Get p.values for wilcoxon U test
    w.test1 <- wilcox.test(l.loh.set[l.loh.set$tissue == "OV","counts"],
                           l.loh.set[l.loh.set$tissue == "BRCA","counts"])
    pval1 <- w.test1$p.value

    w.test2 <- wilcox.test(j.AIs.set[j.AIs.set$tissue == "OV","counts"],
                           j.AIs.set[j.AIs.set$tissue == "BRCA","counts"])
    pval2 <- w.test2$p.value

    diff.window <- data.frame(foldchange = c(med.diff1,med.diff2), pval=c(pval1, pval2),
                              window.size=sizes[j], type=c("LOH","Total.AIs"),
                              sample.type=c(typedata, typedata))

    diff.window.all <- rbind(diff.window.all, diff.window)
  }
}

dim(diff.window.all)
diff.window.all$pval <- -1 * log10(diff.window.all$pval)
diff.window.all$window.size <- factor(diff.window.all$window.size, levels=sizes)
diff.window.all$type  <- factor(diff.window.all$type , levels=c("Total.AIs", "LOH"))

#Plotting results for Figure 1h
p <- ggplot(diff.window.all, aes(x=type, y=window.size))
p <- p + geom_point(aes(color=foldchange, size=pval)) + theme_classic()
p <- p + scale_color_gradientn(name="Foldchange", colours = c("darkorchid4", "cornflowerblue", "gold3"), values=c(0,0.4,1),limits = c(-0.35, 0.35))
p <- p + scale_radius(range = c(4,7), name="-log10(p.val)", limits = c(1.301, 4), breaks = c(1.5,2.5,3.5)) + facet_wrap(~ sample.type)
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1), angle = 45, vjust = 1, hjust = 1), axis.text.y=element_text(size=rel(1)),
               legend.text=element_text(size=rel(1.1)), legend.title=element_text(size=rel(1)),
               axis.title=element_text(size=rel(1)))
p <- p + ylab("AIs range size (Mb)\n") + xlab("")
print(p)
ggsave(p, filename = paste0(output.folder, "Dotplots_difference-TNBC-OVA_byWindow.svg"), width = 10, height = 8.5, units = "cm")
