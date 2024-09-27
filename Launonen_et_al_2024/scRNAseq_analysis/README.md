# hgsc_tme

This repository is complementary to the publication:

Launonen IM, Erkan EP, Niemiec I, et al. Chemotherapy induces myeloid-driven spatial T-cell exhaustion in ovarian cancer. Preprint. bioRxiv. 2024;2024.03.19.585657. Published 2024 Mar 20. doi:10.1101/2024.03.19.585657

The repository contains the scripts to reproduce the scRNA-seq analyses and generate the relevant figures in the manuscript.

Description of files
1. infer_cd8_tcell_phenotypes.Rmd: Used to infer CD8+ T cell phenotypes, using the reference atlas published by Andreatta et al. (2021). The reference atlas is available at: https://figshare.com/articles/dataset/ProjecTILs_human_reference_atlas_of_CD8_tumor-infiltrating_T_cells_version_1/21931875/2 (Figure 4f).
2. 03_analyze_lr_interactions.qmd: Used to analyze ligand-receptor interactions in the scRNA-seq data (Figure 4d).
3. compare_CD163_ITGAX_module_scores.Rmd: Used to compare CD163 and ITGAX module scores in macrophages between clinical groups (Figure 4g-h).
4. compare_cell_abundance_with_propeller.Rmd: Used to compare cell type proportions (Supplementary Figure 1f, Figure 1f).
