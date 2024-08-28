library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(GeoDiff)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(data.table)
library(cowplot)
library(preprocessCore)
library(Biobase)
library(reshape2)

library(umap)
library(Rtsne)

library(clusterProfiler)
library(msigdbr)
library(progeny)
library(reshape2)
library(biomaRt)
library(GSVA)
library(ggpubr)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'

dcc_path <- dir(file.path(data_dir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
pkc_path <- file.path(data_dir, 'metadata/Hs_R_NGS_WTA_v1.0.pkc')
anno_path <- file.path(data_dir, 'metadata/dcc_metadata_all.xlsx')

output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'
output_rds_path <- file.path(output_dir, 'geomx_qc_neggeo_ntc.RDS')

imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
main_var <- "Annotation_cell" # legend in sankey, 


# create dirs and source functions ----------------------------------------

dir.create(output_dir, showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'qc'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx dataset ------------------------------------------------------

geomx_obj <- readNanoStringGeoMxSet(dccFiles = dcc_path, 
                                    pkcFiles = pkc_path, # this goes into fData() - features (probes) annotation
                                    phenoDataFile = anno_path, # this goes into pData() - protocol (samples) annotation
                                    phenoDataSheet = "Sheet1",
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("Aoi", "Roi"),
                                    experimentDataColNames = c("Panel")) # TODO dunno if this is needed

pkcs <- annotation(geomx_obj)
modules <- gsub(".pkc", "", pkcs)
# explore
# View(assayData(geomx_obj)$exprs)
# dim(assayData(geomx_obj)$exprs)
# 
# View(pData(geomx_obj))
# View(pData(protocolData(geomx_obj)))
# View(fData(geomx_obj))
# featureType(geomx_obj)
# annotation(geomx_obj)
# 
# View(summary(geomx_obj, MARGIN = 1)) # for probes
# View(summary(geomx_obj, MARGIN = 2)) # for segments (rois)

print(paste('dim of raw dataset is: '))
print(dim(geomx_obj))

# manualfix of NTC --------------------------------------------------------

# manually add messed up info for NTC to sData()

sdt <- sData(geomx_obj)
sdt$NTC_ID <- apply(sdt, 1, function(x){
  # one NTC/batch
  if(x[['batch_nr']] %in% c(1,2,3,4,7,8)){
    ntc <- sdt$dcc_filename[sdt$`Slide Name` == 'No Template Control' & sdt$batch_nr == x[['batch_nr']]]
  } else{
    ntc <- NA
  }
  
  # the same NTC for batch 4 and 6
  if(x[['batch_nr']] == 6){
    ntc <- sdt$dcc_filename[sdt$`Slide Name` == 'No Template Control' & sdt$batch_nr == 4]
  }
  
  # 2 different NTC for batch 5
  if(x[['batch_nr']] == 5 & grepl('-E-', x[['dcc_filename']])){
    ntc <- sdt$dcc_filename[sdt$`Slide Name` == 'No Template Control' & sdt$batch_nr == 5 & grepl('-E-', sdt$dcc_filename)]
  } else if(x[['batch_nr']] == 5 & grepl('-B-', x[['dcc_filename']])){
    ntc <- sdt$dcc_filename[sdt$`Slide Name` == 'No Template Control' & sdt$batch_nr == 5 & grepl('-B-', sdt$dcc_filename)]
  }
  
  return(ntc)
})

#TODO check in manual if it really is Deduplicatedreads for NTC count
sdt$NTC <- apply(sdt, 1, function(x){
  ntc_cnt <- sdt$DeduplicatedReads[sdt$dcc_filename == x[['NTC_ID']]]
})

#add to protocolData
identical(rownames(protocolData(geomx_obj)@data), sdt$dcc_filename)
protocolData(geomx_obj)@data[, c("NTC_ID", "NTC")] <- sdt[, c("NTC_ID", "NTC")]

#change 'Area' and 'Nuclei' colnames for correct qc flags
pData(geomx_obj) <- dplyr::rename(pData(geomx_obj), 'area' = 'Area', 'nuclei' = 'Nuclei')

# make overall sankey plot ------------------------------------------------
count_segments <- geomx_obj@phenoData@data[main_var != 'NA' & !is.na(main_var), ]

plot_sankey(count_segments, imp_vars, main_var, 
            file.path(output_dir, 'qc/sankey_slides.png'))


# set and plot basic qc parameters ----------------------------------------
# Shift 0 counts to one -needed for  Q3 norm (but not 100% sure why)
geomx_obj <- shiftCountsOne(geomx_obj, useDALogic = TRUE)

qc_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10, 1 in log scale)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100) 
       minArea = 1000)         # Minimum segment area (5000)

# set up qc flags for segments
geomx_obj <- setSegmentQCFlags(geomx_obj, qcCutoffs = qc_params)

# rmv NTC segments
#TODO check if this is not messing up with latter functions
geomx_obj <- geomx_obj[, !(geomx_obj$`Slide Name` == 'No Template Control')]

qc_results_segment <- protocolData(geomx_obj)[["QCFlags"]]
qc_summary <- qc_summarize(qc_results_segment)

print('qc summary:')
print(qc_summary)

# plot qc histograms
# duplicated cause you have to iterate trough 2 lists of names
QC_histogram(sData(geomx_obj), "Trimmed (%)", "Segment", qc_params[["percentTrimmed"]], 
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_trim.png'))

QC_histogram(sData(geomx_obj), "Stitched (%)", "Segment", qc_params[["percentStitched"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_stich.png'))

QC_histogram(sData(geomx_obj), "Aligned (%)", "Segment", qc_params[["percentAligned"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_align.png'))

QC_histogram(sData(geomx_obj), "Saturated (%)", "Segment", qc_params[["percentSaturation"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_satur.png'))

QC_histogram(sData(geomx_obj), "area", "Segment", qc_params[["minArea"]], 
             scale_trans = "log10", file.path(output_dir, 'qc/qc_hist_area.png'))

QC_histogram(sData(geomx_obj), "nuclei", "Segment", qc_params[["minNuclei"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_nuclei.png'))


# negative geometric means ------------------------------------------------

# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(geomx_obj), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 

protocolData(geomx_obj)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(geomx_obj)[, negCols] <- sData(geomx_obj)[["NegGeoMean"]]

for(ann in paste0("NegGeoMean_", modules)) {
  QC_histogram(sData(geomx_obj), ann, "Segment", 2, scale_trans = "log10",
               file.path(output_dir, 'qc/qc_hist_neggeomean.png')) #TODO? why exactly thr = 2?
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(geomx_obj) <- pData(geomx_obj)[, !colnames(pData(geomx_obj)) %in% negCols]

# background modelling based on negative probes ---------------------------

sum(fData(geomx_obj)$Negative) # same as negativeControlSubset(geomx_obj)

# fit poisson distribution
geomx_obj <- fitPoisBG(geomx_obj)

summary(pData(geomx_obj)$sizefact)
summary(fData(geomx_obj)$featfact[fData(geomx_obj)$Negative])

# diagnose Poisson model
set.seed(123)
geomx_diag <- diagPoisBG(geomx_obj)

notes(geomx_diag)$disper

# If the dispersion is >2, one of these factors might be present in the data. 
# We can check for outlier ROIs. People can choose to set outliers to be missing 
# values and rerun the Poisson Background model.

length(which(assayDataElement(geomx_diag, "low_outlier") == 1, arr.ind = TRUE))
length(which(assayDataElement(geomx_diag, "up_outlier") == 1, arr.ind = TRUE))

# Or if a batch effect is assumed, the poisson model can be adjusted to take 
# different groups into account. Here we are grouping the ROIs by slide.

geomx_obj <- fitPoisBG(geomx_obj, groupvar = "Slide Name")

set.seed(123)
geomx_diag <- diagPoisBG(geomx_obj, split = TRUE)

notes(geomx_diag)$disper_sp

# remove flagged segments -------------------------------------------------

table(sData(geomx_obj)$NTC)

qc_results_segment$qc_status <- apply(qc_results_segment, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

geomx_obj <- geomx_obj[, qc_results_segment$qc_status == "PASS", ]

# remove segments with neggeomean < 1.5
geomx_obj <- geomx_obj[, sData(geomx_obj)[["NegGeoMean"]] >= 1.5]

print(paste('dim after removing bad quality segments: '))
print(dim(geomx_obj))

# rmv of probes based on geometric mean and grubbs test -------------------

# the geometric mean of that probe’s counts from all segments divided by the geometric mean 
# of all probe counts representing the target from all segments is less than 0.1
# the probe is an outlier according to the Grubb’s test in at least 20% of the segments
# typically don't change this parameters:

geomx_obj <- setBioProbeQCFlags(geomx_obj, 
                                qcCutoffs = list(minProbeRatio = 0.1,
                                                 percentFailGrubbs = 20), 
                                removeLocalOutliers = TRUE)

qc_results_probe <- fData(geomx_obj)[["QCFlags"]]

# summarise probe qc results
qc_probe_df <- data.frame(Passed = sum(rowSums(qc_results_probe[, -1]) == 0),
                          Global = sum(qc_results_probe$GlobalGrubbsOutlier),
                          Local = sum(rowSums(qc_results_probe[, -2:-1]) > 0
                                      & !qc_results_probe$GlobalGrubbsOutlier))

# retain only probes that passed qc (globally)
geomx_obj <- 
  subset(geomx_obj, 
         fData(geomx_obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geomx_obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

print(paste('dim after removing bad quality probes (globally): '))
print(dim(geomx_obj))

# aggregate probes to features --------------------------------------------

# collapse features to targets
geomx_obj <- aggregateCounts(geomx_obj)

print(paste('dim after collapsing features to targets: '))
print(dim(geomx_obj))

# filter based on LOQ per segment and per gene ----------------------------

# The LOQ is calculated based on the distribution of negative control probes 
# and is intended to approximate the quantifiable limit of gene expression per segment. 
# More stable in larger segments. May not be as accurate in segments with low negative probe counts (ex: <2).
# typically use 2 geometric SD above the geometric mean as the LOQ, which is reasonable for most studies. 
# recommend that a minimum LOQ of 2 be used if the LOQ calculated in a segment is below SD threshold.
# typically don't change this parameters

loq_cutoff <- 2
loq_min <- 2

gene_detect_thr <- 0.1 # segment is removed if <10% of genes > LOQ
segment_detect_rate_thr <- 0.01 # genes are removed if its expr > LOQ in less than 1% of segments

# Calculate LOQ for each segment
LOQ <- data.frame(row.names = colnames(geomx_obj))

for(module in modules){
  LOQ[, module] <-
    pmax(loq_min,
         pData(geomx_obj)[, paste0("NegGeoMean_", module)] * # coming from aggregate_counts
           pData(geomx_obj)[, paste0("NegGeoSD_", module)] ^ loq_cutoff)
}

pData(geomx_obj)$LOQ <- LOQ

# calculate if expr > LOQ per each gene per segment
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(geomx_obj)$Module == module
  Mat_i <- t(esApply(geomx_obj[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering 
LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]

# Save detection rate information to pheno data
# how many genes have been detected in each segment  
pData(geomx_obj)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(geomx_obj)$GeneDetectionRate <- pData(geomx_obj)$GenesDetected / nrow(geomx_obj)

print(paste("median gene nr is: ", as.character(median(pData(geomx_obj)$GenesDetected))))
print(paste("mean gene nr is: ", as.character(mean(pData(geomx_obj)$GenesDetected))))
print(paste("median gene detection rate is: ", as.character(median(pData(geomx_obj)$GeneDetectionRate))))

sapply(imp_vars, function(vname){
  plot_detection_rate(pData(geomx_obj), vname, 
                      file.path(output_dir, 'qc', paste0('gene_detect_rate_', vname, '.png')))
})

# filter out segments with too low gene detection rate
geomx_obj <- geomx_obj[, pData(geomx_obj)$GeneDetectionRate >= gene_detect_thr]

print(paste('dim after removing segments based on LOQ: '))
print(dim(geomx_obj))

# save to probe data in how many segments the given gene was detected
LOQ_Mat <- LOQ_Mat[, colnames(geomx_obj)]
fData(geomx_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(geomx_obj)$DetectionRate <- fData(geomx_obj)$DetectedSegments / nrow(pData(geomx_obj))
LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]

# plot detection rate per gene
plot_gene_detection_rate(fData(geomx_obj), file.path(output_dir, 'qc/gene_detection_rate.png'))

# manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative") # 1 bcs already collapsed to targets
neg_probes <- unique(negativeProbefData$TargetName)

# filter out genes detected > LOQ in less then thr nr of segments (1% for now)
geomx_obj <- 
  geomx_obj[fData(geomx_obj)$DetectionRate >= segment_detect_rate_thr |
              fData(geomx_obj)$TargetName %in% neg_probes, ]

print(paste('dim after removing genes based on LOQ: '))
print(dim(geomx_obj))


# additional background modelling for genes -------------------------------

geomx_obj <- BGScoreTest(geomx_obj)

sum(fData(geomx_obj)[["pvalues"]] < 1e-3, na.rm = TRUE)

# save geomx object after QC ----------------------------------------------

print(paste("median gene nr is: ", as.character(median(pData(geomx_obj)$GenesDetected))))
print(paste("mean gene nr is: ", as.character(mean(pData(geomx_obj)$GenesDetected))))
print(paste("median gene detection rate is: ", as.character(median(pData(geomx_obj)$GeneDetectionRate))))

saveRDS(geomx_obj, file = output_rds_path)

