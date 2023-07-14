library(tidyverse)
library(Matrix)
library(vroom)
library(here)

library(Seurat)
library(anndata)
library(reticulate)
use_condaenv("2023_PVN_Atlas")

library(RExtra)
library(SeuratExtra)

# Get data
df <- readRDS(here("Data/Steuernagel2022/CELLxGENE_HypoMap.rds"))


# Extract elements for new Seurat Object
data <- df@assays$RNA@counts
var <- df@assays$RNA@meta.features
obs <- df@meta.data
umap <- df@reductions$umap


# EnsemblID -> Symbol
var$feature_name <- gsub("_ENS.*$", "", var$feature_name) %>%
  StandardiseGeneNames()
data@Dimnames[[1]] <- var$feature_name


# Remove duplicate rows
# Inspecting duplicates yields >99% similarity between all duplicate genes
data <- data[!duplicated(var$feature_name),]
var <- var[!duplicated(var$feature_name),]
row.names(var) <- var$feature_name


# Create Seurat
df <- CreateSeuratObject(data, "HypoMap", meta.data = obs)
df@assays$RNA@meta.features <- var[1]
df@reductions$umap <- umap


# Metadata
names(df@meta.data)[18] <- "pctMito_RNA"


# Save
saveRDS(df, here("Output/Data/1-FromRaw/Steuernagel2022.rds"))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Steuernagel2022.h5ad"))