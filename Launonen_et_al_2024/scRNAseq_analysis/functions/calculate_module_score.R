CalculateModuleScore <- function(seurat_object, cell_type, gene_sets) {
    
    # Change default idents
    Idents(seurat_object) <- "cell_type"
    
    # Subset cell type(s) of interest
    seurat_subset <- subset(x = seurat_object, idents = cell_type)
    
    # Calculate module scores
    seurat_subset <- UCell::AddModuleScore_UCell(obj = seurat_subset, features = gene_sets)
    
    return(seurat_subset)
}
