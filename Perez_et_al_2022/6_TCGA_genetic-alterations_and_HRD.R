#Script to merge different pathways mutations for samples in TCGA, according to the classification made by  Konstantinopoulos
library(ggplot2)
library(scales)
library(dplyr)
library(colorspace)

###################### Defining variables ############################
#File with scars and RAD50s, BRCAs description
scars_clin <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/TCGA_HRD_newScars_clinical-info_HGSC.csv"

#Somatic mutations in other pathways
somatic <- "Somatic_mutations/mc3.v0.2.8.PUBLIC_OVA-KonstantinopoulosGenes2.csv"

#Gene deletions or amplifications in other pathways
CNVs <- "CNVs2/all_thresholded_ova-KonstantinopoulosGenes.tsv"

#Pathways of interest
pathways <- "utilities/Genes_of_interest_annotations.txt"

#Output folder were to store intermediate table and final plots
outputfolder <- "C:/Users/fernpere/HRD/TCGA_analysis/find_HRD_signaturesHGSC_germlines/"

#Cutoff point for HRD and HRP in ovaHRDscar
ovaHRDscar.cutoff <- 54

#Path with mutation files
setwd("C:/Users/fernpere/HRD/TCGA_analysis/")

######################## Reading input files ########################
scars_BRCARADinfo <- read.table(file=scars_clin, header=T, sep=",")
mutations <- read.csv(file=somatic, header=T, sep=",")
delitionsdup <- read.table(file=CNVs, header=T, sep="\t")
genePaths <- read.table(pathways)
colnames(genePaths) <- c("Hugo_Symbol", "Pathway")

##### Function for shortening the sample IDs: example: TCGA.04.1331.01A to TCGA-04-1331 #####

#When the samples names are the rownames
shortnames <- function(dataframe){
  renames <- NULL
  for (samplename in rownames(dataframe)){
    name_strings <- strsplit(samplename, split = '[.]')[[1]]
    renames <- c(renames, paste(name_strings[1], name_strings[2], name_strings[3], sep ="-"))
  }
  rownames(dataframe) <- renames
  return(dataframe)
}

#When the samples names are the first column
shortnames2 <- function(dataframe){
  renames <- NULL
  for (samplename in dataframe[,1]){
    name_strings <- strsplit(samplename, split = '[-]')[[1]]
    renames <- c(renames, paste(name_strings[1], name_strings[2], name_strings[3], sep ="-"))
  }
  dataframe[,1] <- renames
  return(dataframe)
}

#######################################################################################################################
#################################################### Process CNVs #####################################################
#######################################################################################################################

###Function to get the samples with at least a gene lost in the gene list indicated
sampleswithdeletions <- function(deletionmatrix, genelist){
  sample_list <- NULL
  for (i in rownames(deletionmatrix)){
      genedeleted <- which(deletionmatrix[i,genelist] == "-2")
      if (length(genedeleted) > 0){
        sample_list <- c(i, sample_list)
    }
  }
  return(sample_list)
}

####Editing input matrix
delitionsdupT <- t(delitionsdup)
delitionsdupTsample <- delitionsdupT[c(-1,-2,-3),]
colnames(delitionsdupTsample) <- c(delitionsdup$Gene.Symbol)
delitionsdupTsample <- shortnames(delitionsdupTsample)

### Getting geneloss categories
CCNE1amplificationSamples <- names(which(delitionsdupTsample[,"CCNE1"] == " 2"))
EMSYamplificationSamples <- names(which(delitionsdupTsample[,"C11orf30"] == " 2"))
PTENgenelossSamples <- names(which(delitionsdupTsample[,"PTEN"] == "-2"))
CDK12genelossSamples <- names(which(delitionsdupTsample[,"CDK12"] == "-2"))
BRCA1genelossSamples <- names(which(delitionsdupTsample[,"BRCA1"] == "-2"))
BRCA2genelossSamples <- names(which(delitionsdupTsample[,"BRCA2"] == "-2"))

HRgenestoignore <- c("BRCA1", "BRCA2", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D","RAD52", "RAD54B", "RAD54L")  #This genes are counted in a different category
HRgenes <- c(genePaths[genePaths$Pathway == "HR",1])
HRgenes <- setdiff(HRgenes, HRgenestoignore)
FAgenes <- c(genePaths[genePaths$Pathway == "FA",1])
MMRgenes <- c(genePaths[genePaths$Pathway == "MMR",1])
NERgenes <- c(genePaths[genePaths$Pathway == "NER",1])

HRgenelossSamples <- sampleswithdeletions(delitionsdupTsample, HRgenes)
FAgenelossSamples <- sampleswithdeletions(delitionsdupTsample, FAgenes)
MMRgenelossSamples <- sampleswithdeletions(delitionsdupTsample, MMRgenes)
NERgenelossSamples <- sampleswithdeletions(delitionsdupTsample, NERgenes)

delitionsdupTsample_RAD50s <- delitionsdupTsample[,grep("RAD5", colnames(delitionsdupTsample))]
RAD50genelossSamples <- sampleswithdeletions(delitionsdupTsample_RAD50s, colnames(delitionsdupTsample_RAD50s))

#### Convert in to FALSE or NEGATIVE geneloss/gain, in the order of the samples in the matrix scars_BRCARADinfo
CCNE1amplification <- scars_BRCARADinfo$Row.names %in% CCNE1amplificationSamples
EMSYamplification <- scars_BRCARADinfo$Row.names %in% EMSYamplificationSamples
PTENgeneloss <- scars_BRCARADinfo$Row.names %in% PTENgenelossSamples
CDK12geneloss <- scars_BRCARADinfo$Row.names %in% CDK12genelossSamples
BRCA1geneloss <- scars_BRCARADinfo$Row.names %in% BRCA1genelossSamples
BRCA2geneloss <- scars_BRCARADinfo$Row.names %in% BRCA2genelossSamples


HRgeneloss <- scars_BRCARADinfo$Row.names %in% HRgenelossSamples
FAgeneloss <- scars_BRCARADinfo$Row.names %in% FAgenelossSamples
MMRgeneloss <- scars_BRCARADinfo$Row.names %in% MMRgenelossSamples
NERgeneloss <- scars_BRCARADinfo$Row.names %in% NERgenelossSamples
RAD50sgeneloss <- scars_BRCARADinfo$Row.names %in% RAD50genelossSamples

#######################################################################################################################
############################################## Process somatic mutations ##############################################
#######################################################################################################################

#Select only pathogenic variants
mutations_pathogenic <- mutations[mutations$ACMG_classification %in% "Pathogenic",]

#Shortening samples names
mutations_pathogenic <- shortnames2(mutations_pathogenic)

#Add pathway names
mutations_pathogenic_pathway <- merge(mutations_pathogenic, genePaths, by="Hugo_Symbol", all.x=TRUE)

#Extract samples names with relevant mutations
BRCA1somaticmutSamples <- mutations_pathogenic[mutations_pathogenic$Hugo_Symbol %in% "BRCA1", 1]
BRCA2somaticmutSamples <- mutations_pathogenic[mutations_pathogenic$Hugo_Symbol %in% "BRCA2", 1]
CDK12somaticmutSamples <- mutations_pathogenic[mutations_pathogenic$Hugo_Symbol %in% "CDK12", 1]
RAD50somaticmutSamples <- mutations_pathogenic[which(grepl("RAD5", mutations_pathogenic$Hugo_Symbol)), 1]
PTENsomaticmutSamples <- mutations_pathogenic[mutations_pathogenic$Hugo_Symbol %in% "PTEN", 1]

mutations_pathogenic_pathway_noBRCAsRADs <- mutations_pathogenic_pathway[!(mutations_pathogenic_pathway$Hugo_Symbol %in% HRgenestoignore),]

FAsomaticmutSamples <- mutations_pathogenic_pathway_noBRCAsRADs[which(mutations_pathogenic_pathway_noBRCAsRADs$Pathway == "FA"),2]
HRsomaticmutSamples <- mutations_pathogenic_pathway_noBRCAsRADs[which(mutations_pathogenic_pathway_noBRCAsRADs$Pathway == "HR"),2]
MMRsomaticmutSamples <- mutations_pathogenic_pathway_noBRCAsRADs[which(mutations_pathogenic_pathway_noBRCAsRADs$Pathway == "MMR"),2]
NERsomaticmutSamples <- mutations_pathogenic_pathway_noBRCAsRADs[which(mutations_pathogenic_pathway_noBRCAsRADs$Pathway == "NER"),2]

#### Convert in to FALSE or NEGATIVE mutations, in the order of the samples in the matrix scars_BRCARADinfo
BRCA1somaticmut <- scars_BRCARADinfo$Row.names %in% BRCA1somaticmutSamples
BRCA2somaticmut <- scars_BRCARADinfo$Row.names %in% BRCA2somaticmutSamples
CDK12somaticmut <- scars_BRCARADinfo$Row.names %in% CDK12somaticmutSamples
RAD50somaticmut <- scars_BRCARADinfo$Row.names %in% RAD50somaticmutSamples
PTENsomaticmut <- scars_BRCARADinfo$Row.names %in% PTENsomaticmutSamples

FAsomaticmut <- scars_BRCARADinfo$Row.names %in% FAsomaticmutSamples
HRsomaticmut <- scars_BRCARADinfo$Row.names %in% HRsomaticmutSamples
MMRsomaticmut <- scars_BRCARADinfo$Row.names %in% MMRsomaticmutSamples
NERsomaticmut <- scars_BRCARADinfo$Row.names %in% NERsomaticmutSamples

#######################################################################################################################
########################## Find additional categories according in the scars_BRCARADinfo file #########################
#######################################################################################################################

BRCA1germlinemut <- grepl("BRCA1", scars_BRCARADinfo$Germline_mutations)
BRCA2germlinemut <- grepl("BRCA2", scars_BRCARADinfo$Germline_mutations)
BRCA1promotermethylated <- grepl("BRCA1", scars_BRCARADinfo$Hypermethylation_status)
RAD51Cpromotermethylated <- grepl("RAD51C", scars_BRCARADinfo$Hypermethylation_status)

#######################################################################################################################
########################## Merging all information in the scars_BRCARADinfo file #######################################
#######################################################################################################################


scars_BRCARADinfo_KonstantinopoulosInfo <- cbind(scars_BRCARADinfo,
                                                  BRCA1germlinemut, BRCA1somaticmut, BRCA1geneloss,
                                                  BRCA2germlinemut, BRCA2somaticmut, BRCA2geneloss,
                                                  BRCA1promotermethylated, CDK12somaticmut, CDK12geneloss,
                                                  RAD51Cpromotermethylated, FAgeneloss, FAsomaticmut,
                                                  RAD50somaticmut, RAD50sgeneloss, HRgeneloss, HRsomaticmut,
                                                  EMSYamplification, PTENgeneloss, PTENsomaticmut,
                                                  CCNE1amplification, MMRsomaticmut, MMRgeneloss,
                                                  NERsomaticmut, NERgeneloss)

###Generate category for each sample according to the order in the Pie Chart (Figure 2)
categories <- NULL
for (i in rownames(scars_BRCARADinfo_KonstantinopoulosInfo)){
      category <- 0
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$CCNE1amplification){
        category <- "CCNE1 amplification"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$NERsomaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$NERgeneloss){
        category <- "NER somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$MMRsomaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$MMRgeneloss){
          category <- "MMR somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$PTENgeneloss | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$PTENsomaticmut){
        category <- "PTEN somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$EMSYamplification){
        category <- "EMSY amplification"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$HRsomaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$HRgeneloss){
        category <- "HR somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$RAD50somaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$RAD50sgeneloss){
        category <- "RADcore somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$FAsomaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$FAgeneloss){
        category <- "FA somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$RAD51Cpromotermethylated){
        category <- "RAD51C prommethyl"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$CDK12somaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$CDK12geneloss){
        category <- "CDK12 somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA1promotermethylated){
        category <- "BRCA1 prommethyl"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA2somaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA2geneloss){
        category <- "BRCA2 somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA2germlinemut){
        category <- "BRCA2 germlinemut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA1somaticmut | scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA1geneloss){
        category <- "BRCA1 somaticmut"
      }
      if (scars_BRCARADinfo_KonstantinopoulosInfo[i,]$BRCA1germlinemut){
        category <- "BRCA1 germlinemut"
      }
      if (category == 0){
        category <- "Other"
      }
      categories <- c(categories, category)
}

(table(categories) * 100) / length(categories)

scars_BRCARADinfo_KonstantinopoulosCategories <- cbind(scars_BRCARADinfo_KonstantinopoulosInfo, categories)

columnstoremove <- c("CNV_status", "Hypermethylation_status", "Somatic_mutations", "Germline_mutations", "CNVstatus_HRgenes")
columnstoremovePos <- which(colnames(scars_BRCARADinfo_KonstantinopoulosCategories) %in% columnstoremove)

scars_BRCARADinfo_KonstantinopoulosCategories <- scars_BRCARADinfo_KonstantinopoulosCategories[,-columnstoremovePos]

#Write file in table
colnames(scars_BRCARADinfo_KonstantinopoulosCategories)[1] <- "Sample"

#Select only high grade samples
scars_BRCARADinfo_KonstantinopoulosCategories_HG <- scars_BRCARADinfo_KonstantinopoulosCategories[!is.na(scars_BRCARADinfo_KonstantinopoulosCategories$histological_grade),]

write.table(scars_BRCARADinfo_KonstantinopoulosCategories_HG, file="Categories_Konstantinopoulos_HR_HGSC_germline.csv", sep=",", row.names = FALSE)

#########################################################################################################################################
###################################################### Generate plots ###################################################################
#########################################################################################################################################

# Import data
df <- read.table(file="Categories_Konstantinopoulos_HR_HGSC_germline.csv", sep=",", header = TRUE)
#Set manual order to later label/order pie chart plot
categories.order <- c("BRCA1 germlinemut", "BRCA1 somaticmut", "BRCA2 germlinemut", "BRCA2 somaticmut", "BRCA1 prommethyl",
                      "CDK12 somaticmut", "RAD51C prommethyl", "FA somaticmut", "RADcore somaticmut",
                      "HR somaticmut", "EMSY amplification", "PTEN somaticmut", "MMR somaticmut", "NER somaticmut", "CCNE1 amplification", "Other")

# Force order of categories
df$categories <- factor(df$categories, levels=c(categories.order))

# Sort values by HRDsum within categories levels
df <- arrange(df, df$categories, -HRDsum)

# Create steps for the slices of the pie chart, one for each sample.
df$ymax <- 1:nrow(df)
df$ymin <- 0:(nrow(df)-1)

#For each category calculate the median and mean of ovaHRDscar values
#Also adding ymax and ymin values for the pie-chart sizes
categories.info <- df %>%
  group_by(categories) %>%
  summarise(MeanHRD = mean(HRDsum), MedianHRD= median(HRDsum), ymax = max(ymax), ymin = min(ymin))

# Create a special color palette using scale_fill_continuous_diverging
diverging_hcl(4, palette = "Blue-Red", h1=233, h2=10, register = "Palette.plot")

# Generating Figure 2f
p <- ggplot(df) +
  geom_rect(aes(fill=HRDsum, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  geom_rect(data= categories.info, aes(fill=MeanHRD, ymax=ymax, ymin=ymin, xmax=3, xmin=0), color='white', size=0.1) +
  scale_fill_continuous_diverging(palette="Palette.plot", p1 = 0.3, mid = ovaHRDscar.cutoff) +
  coord_polar(theta="y") +
  theme_void()
p
ggsave(p, filename = paste0(outputfolder, "Konstantinopoulos_PieChart_HR.svg"), height = 10, width = 10, units = "cm")


#####Calculating the differences in abundance of categories vs CCNE
#Generating boxplot for Supplementary Figure 2g
p <- ggplot(df, aes (x=reorder(categories,-HRDsum,FUN=median), y=HRDsum))
p <- p + geom_point(aes(y=HRDsum, color=categories),    position= position_jitter(width= .2), size= 3, alpha = 0.5, show.legend = FALSE)
p <- p + labs(x="", y = "ovaHRDscar") + stat_summary(fun.y=median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width=.6, size = 1)
p <- p + geom_hline(yintercept=ovaHRDscar.cutoff, linetype="dashed", color = "blue")
p <- p + theme(axis.text=element_text(size=rel(1)), strip.placement = "outside",strip.background = element_blank(),
               axis.text.x=element_text(size=rel(1.2), angle = 40, hjust = 1), axis.text.y=element_text(size=rel(1.2)),
               legend.text=element_text(size=rel(1)), strip.text.x = element_text(size=rel(1.2)),
               legend.title=element_text(size=rel(1)), axis.title=element_text(size=rel(1.4)),
               panel.spacing = unit(1, "lines"),
               panel.border = element_rect(linetype = "solid", fill = NA),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "lightgrey"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank())
print(p)
ggsave(p, filename = paste0(outputfolder, "Categories_Konstantinopoulos_HR_newScars.svg"), width = 15, height = 10, units = "cm")

round(table(df$categories) * 100 / nrow(df))


#Calculating U test p values of difference of abundance in categories vs CCNE1
CCNE1values <- df[df$categories == "CCNE1 amplification","HRDsum"]
y = "HR somaticmut"
Uresult <- wilcox.test(df[df$categories == y, "HRDsum"], CCNE1values, alternative="greater")

pvaluesUtest <- sapply(unique(df$categories), function(x){ Uresult <- wilcox.test(df[df$categories == x, "HRDsum"], CCNE1values)
                                            pval <- Uresult$p.value
                                            return(pval)})
names(pvaluesUtest) <- unique(df$categories)
sort(signif(pvaluesUtest, 3)) #p-values for Supplementary Figure 2g
round(pvaluesUtest, 4)  #p-values for Supplementary Figure 2g
