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

VariableFeatures(df) <- adata$uns$highly_variable

# Revision change: Small error in integration pipeline was corrected
# This lead to Oxt-Adarb2 clearly displaying two subclusters
# The subclusters are separated here 
df <- df %>%
  FindNeighbors(reduction = "vae", dims = 1:20, k.param = 50L) %>%
  FindSubCluster("Oxt-Adarb2") %>%
  RenameIdents("Oxt-Adarb2_0" = "Oxt-Adarb2") %>%
  RenameIdents("Oxt-Adarb2_1" = "Ghrh-Adarb2")
levels(df) <- sort(levels(df))

saveRDS(df, here("Output/Data/2-Integration/4-Final.rds"))
