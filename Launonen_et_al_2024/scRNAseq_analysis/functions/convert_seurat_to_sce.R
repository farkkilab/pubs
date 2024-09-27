# Define function to convert Seurat to SCE
ConvertSeuratToSCE <- function(input, cell_idents) {
    
    keep_cells <- cell_idents
    input <- subset(x = input, idents = keep_cells)
    
    # Convert to SCE
    sce <- Seurat::as.SingleCellExperiment(input, assay = "RNA")
    
    # Make names for SCE
    sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
    SummarizedExperiment::colData(sce)$sample = SummarizedExperiment::colData(sce)$sample %>% make.names()
    SummarizedExperiment::colData(sce)$class = SummarizedExperiment::colData(sce)$class %>% make.names()
    SummarizedExperiment::colData(sce)$cell_type = SummarizedExperiment::colData(sce)$cell_type %>% make.names()
    
    return(sce)
}