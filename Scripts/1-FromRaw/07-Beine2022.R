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
df <- (\(df) {
  df1 <- Read10X(here("Data/Beine2022/Sample1"))
  row.names(df1) <- StandardiseGeneNames(row.names(df1))
  
  df2 <- Read10X(here("Data/Beine2022/Sample2"))
  row.names(df2) <- row.names(df1)
  
  df3 <- Read10X(here("Data/Beine2022/Sample3"))
  row.names(df3) <- row.names(df1)
  
  df1 <- MergeRowsByGeneNames(df1)
  df2 <- MergeRowsByGeneNames(df2)
  df3 <- MergeRowsByGeneNames(df3)
  
  df <- CreateSeuratObject(df1) %>%
    merge(CreateSeuratObject(df2)) %>%
    merge(CreateSeuratObject(df3))
  
  Idents(df) <- (df$orig.ident <- "Beine2022")

  return(df)
})()


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

df <- df[,Idents(df) == 22]

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
saveRDS(df,here('Output/Data/1-FromRaw/Beine2022.new.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Beine2022.new.h5ad"))
