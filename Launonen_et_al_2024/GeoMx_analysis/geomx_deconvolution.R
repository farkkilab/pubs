library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(SpatialDecon)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(reshape2)
library(Seurat)
library(tibble)
library(BayesPrism)
library(biomaRt)

# define variables --------------------------------------------------------

data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')
output_rds_path <- file.path(output_dir, 'geomx_qc_norm_deconv.RDS')

scrna_ref_path <- '/media/iganiemi/T7-iga/st/data/scrna/vaharautio_scrnaseq_dataset_downsampled_for_iga_processed.RDS'
scrna_ref_cleaned_path <- file.path(output_dir,'deconvolution', 'scrna_ref_cleaned.RDS')

norm_type <- 'q3_norm'
scrna_anno <- 'mid_lvl_ct' # either 'cell_type' or 'mid_lvl_ct'
ct_nr_thr <- 20 # best 20 or 45 to rmv cell states not abundant enough in scrnaseq

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'prism'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# prepare scrnaseq reference dataset --------------------------------------
geomx_obj <- readRDS(input_rds_path)
scrna_ref_obj <- readRDS(scrna_ref_path)


# repair synonymuous gene names
length(rownames(geomx_obj@assayData$exprs))
length(rownames(scrna_ref_obj@assays$RNA@data))
length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))

adjusted_genes_rna <- adjust_synonym_genes(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data))

#  make a new assay with renamed genes
RNA_common_genes <- scrna_ref_obj@assays$RNA
RNA_common_genes@counts@Dimnames[[1]] <- adjusted_genes_rna
RNA_common_genes@data@Dimnames[[1]] <- adjusted_genes_rna
scrna_ref_obj@assays$RNA_common_genes <- RNA_common_genes

length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))
length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA_common_genes@data)))

# select tumor ct
scrna_ref_obj@meta.data$cell_type <- ifelse(scrna_ref_obj@meta.data$cell_type == 'Epithelial cells', 
                                            'tumor', scrna_ref_obj@meta.data$cell_type)

# cell states - clustering tumor cells by patient
scrna_ref_obj@meta.data$cell_state <- ifelse(scrna_ref_obj@meta.data$cell_type == 'tumor', 
                                             paste0('tumor_', scrna_ref_obj@meta.data$patient), 
                                             scrna_ref_obj@meta.data$cell_type)

# fix cell type labels
scrna_ref_obj@meta.data$mid_lvl_ct <- ifelse(scrna_ref_obj@meta.data$mid_lvl_ct == 'Plasma cells', 
                                             'Bcells', scrna_ref_obj@meta.data$mid_lvl_ct)
scrna_ref_obj@meta.data$mid_lvl_ct <- ifelse(scrna_ref_obj@meta.data$mid_lvl_ct == 'Classical monocytes', 
                                             'Macrophages', scrna_ref_obj@meta.data$mid_lvl_ct)
scrna_ref_obj@meta.data$mid_lvl_ct <- ifelse(scrna_ref_obj@meta.data$mid_lvl_ct == 'Epithelial cells', 
                                             'tumor', scrna_ref_obj@meta.data$mid_lvl_ct)


# QC of cell states

plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
              input.labels=scrna_ref_obj@meta.data$cell_state,
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs",
              cexRow=0.6, cexCol=0.6,
              margins=c(6,6))

dev.off()

plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
              input.labels=scrna_ref_obj@meta.data$mid_lvl_ct,
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)

dev.off()

# check genes outliers
scrna_stat <- plot.scRNA.outlier(
  input=t(scrna_ref_obj@assays$RNA_common_genes@data), #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting.
)

View(scrna_stat)

# filter out outlier genes
scrna_filt <- cleanup.genes (input=t(scrna_ref_obj@assays$RNA_common_genes@data),
                             input.type="count.matrix",
                             species="hs", 
                             gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                             exp.cells=5)

dim(t(scrna_ref_obj@assays$RNA_common_genes@data))
dim(scrna_filt)

# geomx doesn't have to be filtered since later on they took only intersection of genes

# check expr concordance for different gene types
plot.bulk.vs.sc (sc.input = scrna_filt, bulk.input = geomx_raw)

# subset to protein coding genes
scrna_filt_pc <-  select.gene.type(scrna_filt, gene.type = "protein_coding")

#  make a new assay with filtered genes
RNA_common_genes_filt_pc <- scrna_ref_obj@assays$RNA_common_genes
RNA_common_genes_filt_pc@counts <- RNA_common_genes_filt_pc@counts[rownames(RNA_common_genes_filt_pc@counts) %in% colnames(scrna_filt_pc),  ]
RNA_common_genes_filt_pc@data <- RNA_common_genes_filt_pc@data[rownames(RNA_common_genes_filt_pc@data) %in% colnames(scrna_filt_pc),  ]
scrna_ref_obj@assays$RNA_common_genes_filt_pc <- RNA_common_genes_filt_pc

# save adjusted scRNAseq file
saveRDS(scrna_ref_obj, file = scrna_ref_cleaned_path)


# deconvolution by bayesprism ---------------------------------------------
# BayesPrism
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(input_rds_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)

# remove cells from cell states with < thr cells in ref scrnaseq
ct_freq <- as.data.frame(table(scrna_ref_obj@meta.data$cell_state))
low_ct_cells <- scrna_ref_obj@meta.data$cell_name[scrna_ref_obj@meta.data$cell_state %in% 
                                                    as.character(ct_freq$Var1[ct_freq$Freq < ct_nr_thr])]

scrna_ref_obj <- scrna_ref_obj[, !colnames(scrna_ref_obj) %in% low_ct_cells]

# make a prism object
prism_obj <- new.prism(
  reference=t(scrna_ref_obj@assays$RNA_common_genes_filt_pc@data), 
  mixture=t(geomx_obj@assayData$exprs),
  input.type="count.matrix", 
  cell.type.labels = scrna_ref_obj@meta.data[[scrna_anno]], 
  cell.state.labels = scrna_ref_obj@meta.data$cell_state,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run bayesprism
bprism_res <- run.prism(prism = prism_obj, n.cores=18)

# save res
saveRDS(bprism_res, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '_', ct_nr_thr, '.RDS')))


# deconvolution by SpatialDecon -------------------------------------------
# from
# https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(input_rds_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)

# filter geomx object from low complexity genes
geomx_stat <- plot.bulk.outlier(
  bulk.input=t(geomx_obj@assayData$exprs),#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=t(scrna_ref_obj@assays$RNA_common_genes@data), #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=scrna_ref_obj@meta.data$cell_type,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
)

View(geomx_stat)

geomx_stat_to_rm <- geomx_stat[ rowSums(geomx_stat[, -c(1,2)]) >= 1, ]
geomx_filtered <- geomx_obj[!(rownames(geomx_obj) %in% geomx_stat_to_rm),  ]


featureType(geomx_obj) <- "Target"
sampleNames(geomx_obj) <- sData(geomx_obj)[['dcc_filename']]

featureType(geomx_filtered) <- "Target"
sampleNames(geomx_filtered) <- sData(geomx_filtered)[['dcc_filename']]

# prepare cell profile matrix from reference scRNAseq

# format annotations
scrna_anno_dt <- scrna_ref_obj@meta.data[, c('cell_name', scrna_anno)]
rownames(scrna_anno_dt) <- NULL
colnames(scrna_anno_dt) <- c('cell_name', 'cell_type')

custom_oc_mtx <- create_profile_matrix(mtx = scrna_ref_obj@assays$SCT@data,            # cell x gene count matrix
                                       cellAnnots = scrna_anno_dt,  # cell annotations with cell type and cell name as columns
                                       cellTypeCol = "cell_type",  # column containing cell type
                                       cellNameCol = "cell_name",           # column containing cell ID/name
                                       matrixName = "oc_scrnaseq_ref_cell_type_sct", # name of final profile matrix
                                       outDir = output_dir,                    # path to desired output directory, set to NULL if matrix should not be written
                                       normalize = FALSE,                # Should data be normalized?
                                       minCellNum = 50,                   # minimum number of cells of one type needed to create profile, exclusive
                                       minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                       scalingFactor = 1,                # what should all values be multiplied by for final matrix
                                       discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

# run extended SpatialDecon with custom oc mtx ----------------------------

sd_res_custom <- runspatialdecon(object = geomx_obj,
                                    norm_elt = norm_type,                # normalized data
                                    raw_elt = "exprs",                      # expected background counts for every data point in norm
                                    X = custom_oc_mtx,                            # safeTME matrix, used by default
                                    #cell_counts = geomx_obj$Nuclei,      # nuclei counts, used to estimate total cells
                                    #is_pure_tumor = geomx_obj$istumor,   # identities of the Tumor segments/observations
                                    n_tumor_clusters = 5)               # how many distinct tumor profiles to append to safeTME

saveRDS(sd_res_custom, file = file.path(output_dir, 'deconvolution', 'spatial_decon', paste0('sd_res_', scrna_anno, 
                                  '_filt_geomx.rds')))


