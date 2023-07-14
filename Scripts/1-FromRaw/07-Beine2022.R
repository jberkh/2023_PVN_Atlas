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
df <- CreateSeuratObject(Read10X(here("Data/Beine2022/Sample1"))) %>%
  merge(CreateSeuratObject(Read10X(here("Data/Beine2022/Sample2")))) %>%
  merge(CreateSeuratObject(Read10X(here("Data/Beine2022/Sample3"))))
Idents(df) <- (df$orig.ident <- "Beine2022")

# QC/Filtering
df$pctRibo_RNA <- PercentageFeatureSet(df, pattern = "^Rp[sl]")
df$pctMito_RNA <- PercentageFeatureSet(df, pattern = "^mt.")
df <- df[,df$nCount_RNA > 2000 & df$nFeature_RNA > 1000 & df$pctMito_RNA < 15]
VlnPlot(df, features = c("nCount_RNA", "nFeature_RNA", "pctRibo_RNA", "pctMito_RNA"), ncol = 4) & scale_y_log10()


# Basic clustering to annotate classes
df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('nCount_RNA','pctRibo_RNA','pctMito_RNA')) %>%
  RunPCA() %>%
  RunUMAP(dims=1:25)  %>%
  FindNeighbors(dims=1:25) %>%
  FindClusters(resolution = 2)

df <- df[,Idents(df) == 23]

Idents(df) <- "Beine2022"
df$assay <- "10x 3' v3"
df$Dataset <- "Beine2022"
df$Sample_ID <- "Beine2022"
df$assay <- "10x 3' v3"
df$disease <- "normal"
df$organism <- "Mus musculus"
df$sex <- "male"
df$suspension_type <- "nucleus"

#Save as rds
saveRDS(df,here('Output/Data/1-FromRaw/Beine2022.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Beine2022.h5ad"))