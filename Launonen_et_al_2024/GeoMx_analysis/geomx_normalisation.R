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
library(DESeq2)

library(umap)
library(Rtsne)


# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'
input_rds_path <- file.path(output_dir, 'geomx_qc.RDS')
output_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
main_var <- "Annotation_cell" # legend in sankey, 

umap_vars <- c(imp_vars, "Patient", "Sample")
# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'umap_tsne', 'all'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'umap_tsne', 'tumor'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'umap_tsne', 'stroma'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load qc geomx data ------------------------------------------------------

geomx_obj <- readRDS(input_rds_path)

# Q3 normalisation --------------------------------------------------------

negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative") # 1 bcs already collapsed to targets
neg_probes <- unique(negativeProbefData$TargetName)

plot_q3_stats(geomx_obj, main_var, file.path(output_dir, 'qc/q3_stats.png'))


geomx_obj <- normalize(geomx_obj ,
                       norm_method = "quant", 
                       desiredQuantile = .75,
                       toElt = "q3_norm")

# quantile normalisation --------------------------------------------------
# from
# https://github.com/LevivanHijfte/NanoString_normalization_methods/blob/main/Data_preprocessing.R

norm.quantile = normalize.quantiles(as.matrix(geomx_obj@assayData$exprs))
dimnames(norm.quantile) = dimnames(geomx_obj@assayData$exprs)

# DESeq2 normalisation ----------------------------------------------------

#change to integers
expr_int <- apply(geomx_obj@assayData$exprs, c (1, 2), function (x) {(as.integer(x))})

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = expr_int,
                              colData = sData(geomx_obj),
                              design= ~ Segment ) # TODO examine eg if add Annotation_cell or NACT status?

# normalise
dds <- estimateSizeFactors(dds)

# sizeFactors(dds)[1:10] # have a look at size factors

#make df
deseq2_norm_counts <- counts(dds, normalized=TRUE)
dimnames(deseq2_norm_counts) = dimnames(geomx_obj@assayData$exprs)

# add quantile and dseq2 to geomx_obj -------------------------------------

# hacking GeoMx class object 
# TODO this is experimental - newassay is not identical and it may cause problems
# if so, store this in another mtx and use when needed
newassay <- new.env(parent=geomx_obj@assayData)
newassay$exprs <- geomx_obj@assayData$exprs
newassay$q3_norm <- geomx_obj@assayData$q3_norm
newassay$quant_norm <- norm.quantile
newassay$deseq2_norm <- deseq2_norm_counts

geomx_obj@assayData <- newassay

# plot effects of normalisation -------------------------------------------

plot_norm_effect(exprs(geomx_obj)[,1:10], 'Raw Counts', file.path(output_dir, 'qc/norm_raw.png'))


plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "q3_norm"),
                 'Q3 normalised', file.path(output_dir, 'qc/norm_q3.png'))

# TODO I don't like sth with this plot, why all outliers are the same in each segment?
plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "quant_norm"),
                 'Quantile normalised', file.path(output_dir, 'qc/norm_quant.png'))

# super similar to Q3 :0
plot_norm_effect(assayDataElement(geomx_obj[,1:10], elt = "deseq2_norm"),
                 'DESeq2 normalised', file.path(output_dir, 'qc/norm_deseq2.png'))


# make UMAP and t-SNE -----------------------------------------------------

# TODO change for any segment type
# divide for tumor and stroma and do dimentionality reduction for all
geomx_obj_tumor <- geomx_obj[, geomx_obj@phenoData@data$Segment == "tumor"]
geomx_obj_stroma <- geomx_obj[, geomx_obj@phenoData@data$Segment == "stroma"]

geomx_list <- list(all = geomx_obj, tumor = geomx_obj_tumor, stroma = geomx_obj_stroma)

geomx_list_dim_red <- lapply(1:length(geomx_list), function(n){
  
  geomx <- geomx_list[[n]]
  
  # run UMAP and tSNE on Q3 and quantile norm
  for(norm in c('q3_norm', 'quant_norm', 'deseq2_norm')){
    # update defaults for umap to contain a stable random_state (seed)
    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    
    umap_out <-
      umap(t(log2(assayDataElement(geomx , elt = norm))),  
           config = custom_umap)
    
    # save UMAP1 and 2 results to pData
    pData(geomx)[, c(paste0("UMAP1_", norm), paste0("UMAP2_", norm))] <- umap_out$layout[, c(1,2)]
    
    # set the seed for tSNE as well
    set.seed(42) 
    tsne_out <-
      Rtsne(t(log2(assayDataElement(geomx , elt = norm))),
            perplexity = ncol(geomx)*.15)
    
    # save tSNE1 and 2 results to pData
    pData(geomx)[, c(paste0("tSNE1_", norm), paste0("tSNE2_", norm))] <- tsne_out$Y[, c(1,2)]
  }
  
  # generate umap and tsne plots and color by variables
  for(method in c('UMAP', 'tSNE')){
    for(norm in c('q3', 'quant', 'deseq2')){
      for(color_var in umap_vars){
        plot_umap_tsne(pData(geomx), method_type = method, 
                       norm_type = norm, color_var = color_var,
                       output_name = file.path(output_dir, 'umap_tsne', names(geomx_list)[n], 
                                               paste0(method, '_', norm, '_', color_var, '.pdf')))
      }
    }
  }
  
  return(geomx)
})

# update objects
geomx_obj <- geomx_list_dim_red[[1]]
# geomx_qc_tumor <- geomx_list_dim_red[[2]]
# geomx_qc_stroma <- geomx_list_dim_red[[3]]

rm(geomx_list)
rm(geomx_list_dim_red)

# save geomx as RDS

saveRDS(geomx_obj, file = output_rds_path)
