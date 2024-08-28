plot_sankey <- function(data, variables_to_plot, fill_var, output_name){
  count_mat <-   data %>%
    group_by_at(variables_to_plot) %>% 
    summarise(n = n())
  
  test_gr <- gather_set_data(count_mat, 1:length(variables_to_plot))
  
  
  test_gr$x <- mapvalues(test_gr$x, 
                         from=seq(1:length(variables_to_plot)), 
                         to=variables_to_plot)
  test_gr$x <- factor(test_gr$x,
                      levels = variables_to_plot)
  
  # plot Sankey
  ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = get(fill_var)), alpha = 0.5, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.2) +
    geom_parallel_sets_labels(color = "white", size = 3) +
    theme_classic(base_size = 17) + 
    theme(legend.position = "bottom",
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank()) +
    scale_y_continuous(expand = expansion(0)) + 
    scale_x_discrete(expand = expansion(0)) +
    labs(x = "", y = "") +
    annotate(geom = "segment", x = 4.25, xend = 4.25,
             y = 20, yend = 120, lwd = 2) +
    annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
             hjust = 0.5, label = "100 segments")
  
  ggsave(output_name, width = 2000, height = 2000, unit = 'px')
}

###############################################
qc_summarize <- function(QCResults){
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, colnames(QCResults)]),
                           Warning = colSums(QCResults[, colnames(QCResults)]))
  
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  
  QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))
  
  return(QC_Summary)
}

################################################
# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL,
                         output_name = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
  
  print(output_name)
  ggsave(output_name, width = 2000, height = 1000, unit='px')
}

#######################################
plot_detection_rate <- function(segment_data, fill_var, output_name){
  segment_data$DetectionThreshold <- 
    cut(segment_data$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
  
  # stacked bar plot of different cut points (1%, 5%, 10%, 15%)
  ggplot(segment_data,
         aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = get(fill_var))) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")
  
  ggsave(output_name, width = 2000, height = 2000, unit='px')
}

############################################
plot_gene_detection_rate <- function(gene_data, output_name){
  plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
  plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(gene_data$DetectionRate >= x)}))
  plot_detect$Rate <- plot_detect$Number / nrow(gene_data)
  rownames(plot_detect) <- plot_detect$Freq
  
  ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")
  
  ggsave(output_name, width = 2000, height = 1500, unit='px')
}

#######################################################

plot_q3_stats <- function(geomx_obj, ann_of_interest, output_name){
  Stat_data <- 
    data.frame(row.names = colnames(exprs(geomx_obj)),
               Segment = colnames(exprs(geomx_obj)),
               Annotation = pData(geomx_obj)[, ann_of_interest],
               Q3 = unlist(apply(exprs(geomx_obj), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = exprs(geomx_obj)[neg_probes, ])
  
  Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  plt1 <- ggplot(Stat_data_m,
                 aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")
  
  plt2 <- ggplot(Stat_data,
                 aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")
  
  plt3 <- ggplot(Stat_data,
                 aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
  
  btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                       rel_widths = c(0.43,0.57))
  plt_all <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
  
  ggsave(output_name, width=2000, height=1500, unit='px')
}

###########################################################
# plot normalisation effects

plot_norm_effect <- function(expr_data, norm_name, output_name){
  png(filename=output_name, width=1000, height=750, units="px")
  
  boxplot(expr_data,
          col = "#9EDAE5", main = norm_name,
          log='y', names = seq(1:ncol(expr_data)), xlab = "Segment",
          ylab = norm_name)
  
  dev.off()
}

############################################################
plot_umap_tsne <- function(pheno_data, method_type = c('UMAP', 'tSNE'), 
                           norm_type = c('q3', 'quant'), color_var, shape_var = 'Segment',
                           output_name){
  ggplot(pheno_data,
         aes(x = get(paste0(method_type, '1_', norm_type, '_norm')), 
             y = get(paste0(method_type, '2_', norm_type, '_norm')), 
             color = get(color_var), shape = get(shape_var))) +
    geom_point(size = 3) +
    xlab(paste0(method_type, '1_', norm_type, '_norm')) +
    ylab(paste0(method_type, '2_', norm_type, '_norm')) +
    scale_color_discrete(name = color_var) + 
    scale_shape_discrete(name = shape_var) + 
    theme_bw()
  
  ggsave(output_name, width = 2000, height = 1500, unit='px', device='pdf')
}

############################################################
calculate_ora <- function(gene_vect, bcg_gene_vect, msigdb_df, padj = 0.1){
  ora <- enricher(
    gene = gene_vect,
    pvalueCutoff = padj, # Can choose a FDR cutoff
    pAdjustMethod = "BH", 
    universe = bcg_gene_vect, 
    TERM2GENE = dplyr::select(msigdb_df, gs_name, human_gene_symbol)
  )
  
  if(!is.null(ora)){
    ora_df <- data.frame(ora@result) %>%
      filter(p.adjust <= padj)
  } else{
    ora_df <- data.frame()
  }
  return(ora_df)
}

##########################################################
plot_volcano_deg <- function(results, plot_name, top_n_lab, group_pos, group_neg){
  # Categorize Results based on P-value & FDR for plotting
  results$Color <- "NS or FC < 0.5"
  results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results$Color[results$FDR < 0.05] <- "FDR < 0.05"
  results$Color[results$FDR < 0.001] <- "FDR < 0.001"
  results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
  results$Color <- factor(results$Color,
                          levels = c("NS or FC < 0.5", "P < 0.05",
                                     "FDR < 0.05", "FDR < 0.001"))
  
  # pick top genes for either side of volcano to label
  # order genes for convenience:
  results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
  top_g <- c()
  for(cond in c("tumor", "stroma")) {
    ind <- results$Segment == cond
    top_g <- c(top_g,
               results[ind, 'Gene'][
                 order(results[ind, 'invert_P'], decreasing = TRUE)[1:top_n_lab]],
               results[ind, 'Gene'][
                 order(results[ind, 'invert_P'], decreasing = FALSE)[1:top_n_lab]])
  }
  top_g <- unique(unlist(top_g))
  results <- results[, -'invert_P'] # remove invert_P from matrix
  
  # Graph results
  volc <- ggplot(results,
                 aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                     color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = paste("Enriched in ",  group_neg, " <- log2(FC) -> Enriched in ", group_pos),
         y = "Significance, -log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results, (Gene %in% top_g) & (FDR < 0.05) & (Estimate > 0.5 | Estimate < -0.5)),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom") +
    facet_wrap(~Segment, scales = "fixed") +
    ggtitle(paste(plot_name, group_pos, group_neg))
  
  ggsave(file.path(output_dir, 'dge', paste0('volc_', plot_name, '_', group_pos, '_', group_neg,  '.png')), width = 4000, height = 2000, unit='px')
}

################################################################
# change ensembl/entrez into gene names for nested list of genes
# def
# inp
# args
# outp
gene_2names <- function(gene_inp_list, conv = c('ens', 'entrez'), type = 'list'){
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  #TODO make it for df if needed and also the other way around
  # gene_df_names <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id'),
  #                           filters = 'ensembl_gene_id',
  #                           values = as.character(unlist(gene_list)),
  #                           mart = ensembl)
  # 
  # 
  # 
  # marker_ind <- left_join(marker_ind, marker_ind_names)
  # rm(marker_ind_names)
  # 
  if(conv == 'ens'){
    conv_name <- 'ensembl_gene_id'
  } else if(conv == 'entrez'){
    conv_name <- 'entrezgene_id'
  } else{stop()}
  
  if(type == 'list'){
    
    gene_list <- lapply(gene_inp_list, function(x){
      gene_names <- getBM(attributes=c('external_gene_name', conv_name),
                          filters = conv_name,
                          values = x,
                          mart = ensembl)
      
      gene_names <- gene_names$external_gene_name
    })
    return(gene_list)
  }
}

############################################################
# make boxplot for pathway
pathway_boxplot <- function(df, pathway_colname, score_colname, color_colname, facet_var,
                            plot_title, output_path, statistic_test="t_test", ymin=-1, ymax=1.4,
                            manual_colours = c("#F8766D", "#00BA38", "#619CFF", "#C77CFF")){
  # per Anno cell type
  gsva_boxpl <- ggplot(data = df, aes(x = get(pathway_colname), y = get(score_colname), fill = get(color_colname))) +
    #geom_boxplot() +
    geom_violin() +
    geom_point(position= position_jitterdodge(dodge.width = 1, jitter.width= .3, jitter.height = 0),
               size= 0.2, alpha = 0.6) +
    stat_summary(fun = "mean", geom = "point", colour = "red", position = position_dodge(0.9), size=0.3) +
    geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = TRUE, size = 0.2, label.size = 2.8) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size = 5)) +
    ggtitle(plot_title)+
    xlab(pathway_colname) +
    ylab(paste0(score_colname)) +
    guides(fill=guide_legend(title=color_colname)) +
    scale_fill_manual(values=manual_colours) #+
    #scale_y_continuous(trans='log10')
    #ylim(ymin, ymax)
  
  if(length(facet_var) == 1){
    gsva_boxpl <- gsva_boxpl +
      facet_wrap(~get(facet_var), scales = "fixed", dir="v")
  } else if(length(facet_var) == 2){
    gsva_boxpl <- gsva_boxpl +
      facet_wrap(get(facet_var[1])~get(facet_var[2]), scales = "fixed", dir="v", nrow=2)
  } else if(length(facet_var) > 2){
    stop('only 1 or 2 variables for facet')
  }
  
  pdf(file= output_path, width=8, height=5)
  plot(gsva_boxpl)
  dev.off()
  #ggsave(output_path, height = 2000, width = 3000, unit = 'px', device='pdf')
}

############################################
#############################################
# adjust gene names in external vector to geomx names using ensembl synonyms
# input: 2 vectors with HGSC gene names. names in 2nd vector will be adjusted to the 1st one
# returns 2nd vector with common genes names if possible
# it doesnt have to be geomx, any vector is fine, but keep it to avoid confusion
adjust_synonym_genes <- function(geomx_gene_names, gene_vector){
  gene_vector <- unlist(gene_vector)
  geo_non_ex <- setdiff(gene_vector, geomx_gene_names)
  
  if(length(geo_non_ex) > 0){
    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    geo_non_ex_syn <- getBM(attributes = c('external_gene_name', 'external_synonym'),
                            filters = 'external_gene_name',
                            values = geo_non_ex,
                            mart = ensembl)
    
    
    geo_syn_in_gene_vector <- filter(geo_non_ex_syn, external_synonym %in% geomx_gene_names & 
                                       !(external_synonym %in% gene_vector)) %>%
      distinct(external_gene_name, .keep_all = T) %>% # it'll remove a handful of weird genes with multiple synonyms simultaneously present in scrna, may be ignored
      distinct(external_synonym, .keep_all = T)
    
    if(nrow(geo_syn_in_gene_vector) > 0){
      common_genes <- sapply(gene_vector, function(x){
        if(x %in% geo_syn_in_gene_vector$external_gene_name){
          gname <- geo_syn_in_gene_vector$external_synonym[geo_syn_in_gene_vector$external_gene_name == x]
        } else{
          gname <- x
        }
        return(gname)
      })
      
      stopifnot(length(common_genes) == length(gene_vector))
      return(common_genes)
    } else{
      return(gene_vector)
    }
  } else{
    return(gene_vector)
  }
}
