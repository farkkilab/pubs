get_top_n_lr_pairs_v2 <- function(prioritization_tables, top_n, groups_oi = NULL, senders_oi = NULL, 
                                  receivers_oi = NULL, rank_per_group = TRUE) 
{
    prioritization_tbl_oi = prioritization_tables$group_prioritization_tbl %>% 
        dplyr::filter(group == top_group & fraction_expressing_ligand_receptor > 0) %>% 
        dplyr::filter((activity_scaled > 0 &direction_regulation == "up")|(activity_scaled < 0 & direction_regulation == "down"))%>% 
        dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score)
    
    if (!is.null(groups_oi)) {
        prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(group %in% 
                                                                            groups_oi)
    }
    if (!is.null(senders_oi)) {
        prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(sender %in% 
                                                                            senders_oi)
    }
    if (!is.null(receivers_oi)) {
        prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(receiver %in% 
                                                                            receivers_oi)
    }
    if (rank_per_group == TRUE) {
        prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::group_by(group) %>% 
            dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>% 
            dplyr::filter(prioritization_rank <= top_n)
    }
    else {
        prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>% 
            dplyr::filter(prioritization_rank <= top_n)
    }
    return(prioritization_tbl_oi)
}