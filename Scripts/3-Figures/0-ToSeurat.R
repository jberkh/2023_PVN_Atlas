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

# Get dataset
adata <- anndata::read_h5ad(here("Output/Data/2-Integration/4-Final.h5ad"))

df <- CreateSeuratObject(as.matrix(t(adata$X)), meta.data = adata$obs) %>%
  NormalizeData()

df[["mde"]] <- adata$obsm$X_mde %>% 
  (\(x) {row.names(x) <- colnames(df); x})(.) %>%
  CreateDimReducObject(., key = "MDE", global = T)

df[["vae"]] <- adata$obsm$X_vae %>% 
  (\(x) {row.names(x) <- colnames(df); x})(.) %>%
  CreateDimReducObject(., key = "VAE", global = F)

idents <- as.character(df$annotated)
idents[is.na(idents)] <- "Beine2022"
Idents(df) <- idents
levels(df) <- sort(levels(df))

VariableFeatures(df) <- adata$uns$highly_variable

saveRDS(df, here("Output/Data/2-Integration/4-Final.rds"))
