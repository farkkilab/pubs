library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(GSVA)
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(DESeq2)
library(msigdbr)
library(tibble)

# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

input_bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
deconv_type <- 'mid_lvl_ct' # either mid_lvl_ct or cell_type

input_sd_deconv_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')

sig_additional_path <- '/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/additional_signatures_macro_tcells_msigdb_filt.csv'

# if not doing all hall_cp
selected_sig_path <- file.path('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/immune_signatures_selected_names.csv')

######
imp_vars <- c("Segment", "Annotation_cell", "NACT status", "PFS") # vals used for sankey, detection rate plots, 
gsva_vars <- c(imp_vars, 'dcc_filename', 'Patient') 

norm_type <- 'q3_norm' # either 'q3_norm' or 'quant_norm'

do_gsva_hal_cp_all <- FALSE # whether do gsva on all hallmark and cp paths from msigdb

# make dirs and source functions ------------------------------------------

dir.create(file.path(output_dir, 'gsva'), showWarnings = T, recursive = T)
#dir.create(file.path(output_dir, 'progeny'), showWarnings = T, recursive = T)

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

# load geomx obj from rds -------------------------------------------------

geomx_obj <- readRDS(input_rds_path)

# read expression mtx
expr_mtx <- assayDataElement(geomx_obj, elt = norm_type)

# load deconvoluted signal for macrophages and tcells ---------------------

if(deconv_type == 'mid_lvl_ct'){
  ct_names <- ct_names <- c("Tcells", "Bcells", "Fibroblasts", "NKcells", "Macrophages", "DCs", "tumor", "Endothelial cells")
} else if(deconv_type == 'cell_type'){

  ct_names <- c("Fibroblasts","Macrophages", "tumor", "Endothelial cells", "Classical monocytes",
                "Tem/Trm cytotoxic Tcells", "CD16+ NK cells", "pDC", "Plasma cells", "DC1", "Regulatory Tcells",
                "Naive B cells", "Migratory DCs", "NK cells", "Type 17 helper T cells", "Memory B cells",
                "Tcm/Naive helper Tcells", "CD16- NK cells")
}

deconv_res <- readRDS(input_bp_deconv_path)

# extract coeff of variation per cell type
# mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
# histogram suggests 0.2 as thr
deconv_ct_list <- lapply(ct_names, function(ct_name){
  cell_frac_cv <- as.data.frame(deconv_res@posterior.theta_f@theta.cv)
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- BayesPrism::get.exp(bp=deconv_res,
                                      state.or.type="type",
                                      cell.name=ct_name)
  
  deconv_ct <- deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ]
  deconv_ct <- varianceStabilizingTransformation(round(t(deconv_ct))) # normalisation
  
  return(deconv_ct)
})

# GSVA on all Hallmark + CP + additional -------------------------------------
# do GSVA on all Hallmark + CP from msigDB or on selected Hal + CP + additional pathways

# prepare msigdb signatures list
msigdb_df <- msigdbr(species = "Homo sapiens")
msigdb_df <- filter(msigdb_df, gs_cat == 'H' | 
                      gs_subcat %in% c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME', 'CP:PID', 'CP:WIKIPATHWAYS', 'GO:BP'))

msigdb_df$gene_symbol_adj <- adjust_synonym_genes(rownames(geomx_obj), msigdb_df$gene_symbol)

hal_cp_list <- lapply(unique(msigdb_df$gs_name), function(x){
  gs <- filter(msigdb_df, gs_name == x)
  gs_genes <- unique(gs$gene_symbol_adj)
})

names(hal_cp_list) <- unique(msigdb_df$gs_name)

# filter list to selected pathways
if(!do_gsva_hal_cp_all){
  selected_sig <- fread(selected_sig_path)
  hal_cp_list <- hal_cp_list[names(hal_cp_list) %in% selected_sig$pathway]
  out_name <- 'selected_and_additional'
} else{
  out_name <- 'hal_cp_full_and_additional'
}

# prepare additional signatures list
sig_list_additional <- as.list(fread(sig_additional_path))
sig_list_additional <- lapply(sig_list_additional, function(l){l[l !=""]})
sig_list_additional <- lapply(sig_list_additional, function(x){
  adjust_synonym_genes(rownames(geomx_obj), x)})

sig_list_all <- c(hal_cp_list, sig_list_additional)

expr_list <- deconv_ct_list
expr_list[[length(expr_list) + 1]] <- expr_mtx
names(expr_list) <- c(paste0('deconv_', ct_names, '_', deconv_type), 'all')

gsva_list_long <- lapply(1:length(expr_list), function(x){
  # do gsva
  gsva <- gsva(gsvaParam(expr_list[[x]], sig_list_all, kcdf="Gaussian", minSize = 5))
  # do ssgsea
  #gsva <- gsva(ssgseaParam(expr_list[[x]], sig_list_all, minSize = 5, normalize = F))
  
  # adjust df and save
  gsva_long <- melt(gsva)
  colnames(gsva_long) <- c('pathway','dcc_filename', 'gsva_score')
  gsva_long$expr_signal <- names(expr_list)[x]
  gsva_long <- left_join(gsva_long, pData(geomx_obj)[gsva_vars])
  
  fwrite(gsva_long, file.path(output_dir, 'gsva', paste0('ssgsea_notnorm_', names(expr_list)[x], '_', out_name,  '.csv')))
  
  return(gsva_long)
})