make_sample_lr_prod_activity_plots_v2 <- function(prioritization_tables, prioritized_tbl_oi, widths = NULL) {
    
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    
    sample_data = prioritization_tables$sample_prioritization_tbl %>% 
        dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
        dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                      lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
        dplyr::arrange(receiver) %>% 
        dplyr::group_by(receiver) %>% 
        dplyr::arrange(sender, .by_group = TRUE)
    
    group_data = prioritization_tables$group_prioritization_tbl %>% 
        dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
        dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, ligand_receptor_lfc_avg, 
                        activity, activity_scaled, direction_regulation, fraction_ligand_group, 
                        prioritization_score, scaled_avg_exprs_ligand) %>% 
        dplyr::filter(id %in% sample_data$id) %>% 
        dplyr::arrange(receiver) %>% 
        dplyr::group_by(receiver) %>% 
        dplyr::arrange(sender, .by_group = TRUE)
    
    group_data = group_data %>% 
        dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
    
    keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
    
    names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
    
    
    group_data_sorted <- group_data %>% 
        arrange(desc(activity_scaled)) %>% 
        mutate(lr_interaction = factor(lr_interaction, 
                                       levels = rev(unique(lr_interaction[order(activity_scaled, decreasing = T)]))))
    
    sample_data$lr_interaction <- factor(x = sample_data$lr_interaction, levels = levels(group_data_sorted$lr_interaction))
    
    # Plot
    p2 = group_data_sorted %>%
        ggplot(aes(direction_regulation, lr_interaction, fill = activity_scaled)) + 
        geom_tile(color = "whitesmoke") + 
        facet_grid(sender_receiver ~ group, scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.title = element_blank(), 
              axis.text.y = element_text(face = "bold.italic", size = 9), 
              axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
              strip.text.x.top = element_text(angle = 0), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.x = unit(0.2, "lines"), 
              panel.spacing.y = unit(0.25, "lines"), 
              strip.text.x = element_text(size = 10, color = "black", face = "bold"), 
              strip.text.y = element_blank(), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(fill = "Scaled Ligand\nActivity in Receiver")
    max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
    custom_scale_fill = scale_fill_gradientn(colours = c("white", 
                                                         RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]), 
                                             values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1), 
                                             limits = c(-1 * max_activity, max_activity))
    p2 = p2 + custom_scale_fill
    
    p3 = group_data_sorted %>% 
        ggplot(aes(direction_regulation, lr_interaction, fill = activity)) + 
        geom_tile(color = "whitesmoke") + 
        facet_grid(sender_receiver ~ group, scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + theme(axis.ticks = element_blank(), 
                              axis.title = element_blank(), axis.title.y = element_blank(), 
                              axis.text.y = element_blank(), axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
                              strip.text.x.top = element_text(angle = 0), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                              panel.spacing.x = unit(0.2, "lines"), panel.spacing.y = unit(0.25, "lines"), 
                              strip.text.x = element_text(size = 10, color = "black", face = "bold"), 
                              strip.text.y = element_blank(), 
                              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(fill = "Ligand\nActivity in Receiver")
    max_activity = (group_data$activity) %>% max()
    min_activity = (group_data$activity) %>% min()
    custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white", high = "darkorange", midpoint = 0)
    p3 = p3 + custom_scale_fill
    
    p1 = sample_data %>% ggplot(aes(sample, lr_interaction, 
                                    color = scaled_LR_pb_prod, size = keep_sender_receiver)) + 
        geom_point() + facet_grid(sender_receiver ~ group, scales = "free", 
                                  space = "free", switch = "y") + scale_x_discrete(position = "top") + 
        theme_light() + theme(axis.ticks = element_blank(), 
                              axis.title = element_blank(), axis.text.y = element_text(face = "bold.italic", size = 9), 
                              axis.text.x = element_text(size = 9,angle = 90, hjust = 0), 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), panel.spacing.x = unit(0.4,"lines"), 
                              panel.spacing.y = unit(0.25, "lines"), 
                              strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0), 
                              strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0), 
                              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled L-R\npseudobulk exprs product", 
             size = "Sufficient presence\nof sender & receiver") + 
        scale_size_manual(values = keep_sender_receiver_values)
    max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
    custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                              values = c(0, 0.35, 0.485,0.5, 0.515, 0.65, 1), 
                                              limits = c(-1 * max_lfc, max_lfc))
    p1 = p1 + custom_scale_fill
    
    if (!is.null(widths)) {
        p = patchwork::wrap_plots(p1, p2, p3, nrow = 1, guides = "collect", widths = widths)
    }
    else {
        p = patchwork::wrap_plots(p1, p2, p3, nrow = 1, guides = "collect", 
                                  widths = c(sample_data$sample %>% unique() %>% length(), 
                                             3 * (sample_data$group %>% unique() %>% length()), 
                                             3 * (sample_data$group %>% unique() %>% length())))
    }
    
    return(p)
}