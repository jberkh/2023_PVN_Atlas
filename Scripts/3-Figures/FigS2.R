# Imports
set.seed(1, sample.kind = "Rejection")
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(Matrix)
  library(here)
  PROJDIR <- here()
  
  library(Seurat)
  library(reticulate)
  use_condaenv("2023_PVN_Atlas")
  
  #Package w/ own functions
  library(RExtra)
  library(SeuratExtra)
})

# Get data

# Get dataset
adata <- anndata::read_h5ad(here("Output/Data/2-Integration/3-Excluded.h5ad"))

df <- CreateSeuratObject(as.matrix(t(adata$X)), meta.data = adata$obs) %>%
  NormalizeData()

df[["mde"]] <- adata$obsm$X_mde %>% 
  (\(x) {row.names(x) <- colnames(df); x})(.) %>%
  CreateDimReducObject(., key = "MDE", global = T)

df[["vae"]] <- adata$obsm$X_vae %>% 
  (\(x) {row.names(x) <- colnames(df); x})(.) %>%
  CreateDimReducObject(., key = "VAE", global = F)

df@reductions$umap <- df@reductions$mde 

Idents(df) <- df$annotated
levels(df) <- sort(levels(df))

# Process data
df <- df %>%
  NormalizeData() %>%
  ScaleData() 

# Panel A
SFig2A <- function(df) {
  p <- DimPlot(df) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom"
    )
  return(p)
}
SFig2A(df)

# Save PNG
png(here("Output/Figures/FigS2A.png"), 
    res = 450, height = 18, width = 16, units = "cm")
SFig2A(df)
dev.off()

