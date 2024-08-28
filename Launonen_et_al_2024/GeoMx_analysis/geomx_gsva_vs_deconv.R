library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ComplexHeatmap)
library(GSVA)
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(DESeq2)
library(msigdbr)
library(tibble)
library(ggcorrplot)
library(circlize)

# get variables -----------------------------------------------------------
data_dir <- '/media/iganiemi/T7-iga/st/data/geomx/nact_experiment/'
output_dir <- '/media/iganiemi/T7-iga/st/geomx-processing/results/nact2'

input_rds_path <- file.path(output_dir, 'geomx_qc_norm.RDS')

input_gsva_all_path <- file.path(output_dir, 'gsva', 'gsva_all_additional_msigdb_filt.csv')
input_gsva_deconv_macro_path <- file.path(output_dir, 'gsva', 'gsva_deconv_Macrophages_mid_lvl_ct_selected_and_additional.csv')
input_gsva_deconv_tcell_path <- file.path(output_dir, 'gsva', 'gsva_deconv_Tcells_mid_lvl_ct_selected_and_additional.csv')

input_bp_deconv_path <- file.path(output_dir, 'deconvolution', 'bp', 'bp_res_mid_lvl_ct_45.RDS')
input_sd_deconv_path <- file.path(output_dir, 'deconvolution', 'sd', 'sd_res_mid_lvl_ct_nofilt.rds')

# if not doing all hall_cp
selected_sig_path <- file.path('/media/iganiemi/T7-iga/st/geomx-processing/data/signatures/immune_signatures_selected_forpaper_less_names2.csv')

source('/media/iganiemi/T7-iga/st/geomx-processing/src/geomx_utils.R')

lm_dirname <- 'lm_gsva_deconv_fin'

dir.create(file.path(output_dir, 'gsva', lm_dirname), showWarnings = T, recursive = T)


# read files --------------------------------------------------------------

selected_sig <- fread(selected_sig_path)

ct_names <- c("Tcells", "Macrophages", "tumor")

sd_deconv <- readRDS(input_sd_deconv_path)
sd_deconv <- data.frame(pData(sd_deconv)[, 'prop_of_all'])
sd_deconv <- sd_deconv[, ct_names]

bp_deconv <- readRDS(input_bp_deconv_path)
bp_deconv <- BayesPrism::get.fraction(bp=bp_deconv,
                         which.theta="final",
                         state.or.type="type")
bp_deconv <- bp_deconv[, ct_names]

# iterate through 3 files
# iterate trough bp/sd deconv
gsva_type_list <- list('all' = input_gsva_all_path, 'deconv_macro' = input_gsva_deconv_macro_path,
                       'deconv_tcell' = input_gsva_deconv_tcell_path)

deconv_list <- list('bp' = bp_deconv, 'sd' = sd_deconv)

sapply(1:length(gsva_type_list), function(x){
  gsva <- fread(gsva_type_list[[x]])
  
  gsva <- gsva[gsva$pathway %in% selected_sig$pathway, ]
  
  #rename the pathways
  gsva <- left_join(gsva, selected_sig[, c('pathway', 'path_shortname')])
  gsva$pathway <- gsva$path_shortname
  selected_sig$pathway <- selected_sig$path_shortname
  
  gsva_type_name <- names(gsva_type_list)[x]
  print('XXXXXX')
  print(gsva_type_name)
  
  sapply(1:length(deconv_list), function(x){
    deconv <- deconv_list[[x]]
    deconv_name <- names(deconv_list)[x]
    print(deconv_name)
    
    dir.create(file.path(output_dir, 'gsva', lm_dirname, gsva_type_name, deconv_name), showWarnings = T, recursive = T)
    
    deconv <- as.data.frame(deconv)
    deconv <- tibble::rownames_to_column(deconv, 'dcc_filename')
    
    gsva_df <- left_join(gsva, deconv)
    # iterate trough all, pre, post, short, long, 
    
    gsva_list <- list('all' = gsva_df,
                      'pre' = gsva_df[gsva_df$`NACT status` == 'pre', ],
                      'post' = gsva_df[gsva_df$`NACT status` == 'post', ],
                      'long' = gsva_df[gsva_df$`NACT status` == 'post' & gsva_df$PFS == 'Long', ],
                      'short' = gsva_df[gsva_df$`NACT status` == 'post' & gsva_df$PFS == 'Short', ])
    
    
    sapply(1:length(gsva_list), function(x){
      gsva_inp <- gsva_list[[x]]
      gsva_name <- names(gsva_list)[x]
      print(gsva_name)
      
      # make scatterplot with lm
      gsva_lm <- lapply(unique(as.vector(gsva_inp$pathway)), function(path_name){
        
        lapply(unique(as.vector(gsva_inp$Segment)), function(seg){
          gsva_path <- gsva_inp[gsva_inp$pathway == path_name & gsva_inp$Segment == seg, ]
          
          lapply(colnames(deconv)[-1], function(ct){
            # fit lm with ct fraction as explanatory var
            lm_res <- lm(gsva_score~get(ct),data=gsva_path)
            lm_coef <- summary(lm_res)$coefficients[2]
            lm_rsq <- summary(lm_res)$adj.r.squared
            
            cor_pe <- cor.test(gsva_path$gsva_score, gsva_path[[ct]], method = 'pearson')
            
            gsva_lm_res <- list(pathway = path_name, deconv_ct = ct, segment = seg,
                                lm_coef = lm_coef, lm_rsq = lm_rsq,
                                cor_pe_coeff = cor_pe$estimate, cor_pe_pval = cor_pe$p.value)
            
            png(file = file.path(output_dir, 'gsva', lm_dirname, gsva_type_name, deconv_name,
                                 paste0('scatter_', path_name, '_', ct, '_', seg, '.png')))
            plot(gsva_path[[ct]], gsva_path$gsva_score, xlab = path_name, ylab = ct, sub = paste0('R2 = ', lm_rsq))
            abline(lm(gsva_score~get(ct),data=gsva_path),col='red')
            dev.off()
            
            return(gsva_lm_res)
          })
        })
      })

      gsva_lm <- unlist(gsva_lm, recursive = F)
      gsva_lm <- unlist(gsva_lm, recursive = F)
      gsva_lm_df <- rbindlist(gsva_lm, fill=TRUE)


      fwrite(gsva_lm_df, file.path(output_dir, 'gsva', lm_dirname, gsva_type_name, deconv_name,
                                   paste0('lm_gsva_', gsva_name, '.csv')))

      # make a corrplot
      lapply(unique(gsva_lm_df$segment), function(seg){
        gsva_lm_df_seg <- gsva_lm_df[gsva_lm_df$segment == seg, ]
        
        gsva_corr_coeff <- as.data.frame(dcast(gsva_lm_df_seg, pathway ~ deconv_ct,
                                               value.var = 'cor_pe_coeff'))
        
        gsva_corr_pval <- as.data.frame(dcast(gsva_lm_df_seg, pathway ~ deconv_ct,
                                              value.var = 'cor_pe_pval'))
        
        sapply(unique(selected_sig$path_type), function(path_type){
          print(path_type)
          gsva_corr_coeff_path_type <- dplyr::filter(gsva_corr_coeff, pathway %in%
                                                       selected_sig$pathway[selected_sig$path_type == path_type]) %>%
            tibble::column_to_rownames('pathway')
          
          gsva_corr_pval_path_type <- dplyr::filter(gsva_corr_pval, pathway %in%
                                                      selected_sig$pathway[selected_sig$path_type == path_type]) %>%
            tibble::column_to_rownames('pathway')
          
          gsva_corr_coeff_path_type[gsva_corr_pval_path_type >= 0.05] <- 0 #rmv insignificant correlations
          
          col_fun <- colorRamp2(c(-1, 0, 1), hcl_palette = 'RdBu', reverse = TRUE)
          
          path_type_heat <- Heatmap(as.matrix(gsva_corr_coeff_path_type),
                                    height = unit(6, "cm") , width = unit(6, "cm"),border="white",
                                    rect_gp = gpar(col = "white", lwd = 2), name=path_type,
                                    cluster_columns = F, cluster_rows= T,
                                    column_title_gp = gpar(fontsize = 14),
                                    show_heatmap_legend = T, col = col_fun,
                                    column_names_gp = gpar(fontsize = 12),
                                    row_names_gp = gpar(fontsize = 10),
                                    column_names_rot = 45,
                                    column_title = paste(gsva_name, seg),
                                    heatmap_legend_param = list(title = 'correlation'))
          
          pdf(file=file.path(output_dir, 'gsva', lm_dirname, gsva_type_name, deconv_name,
                             paste0('corrplot_gsva_', gsva_name,'_', path_type, '_', seg, '.pdf')),
              width = 8, height = 6)
          draw(path_type_heat, heatmap_legend_side = "left")
          dev.off()
          
        })
      })
    })
  })
})
