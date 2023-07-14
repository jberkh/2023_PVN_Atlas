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
system('sed -i "s/^Ctrl_GCTATGTAATAA/\\tCtrl_GCTATGTAATAA/g" Data/Son2021/GSM5076376_Sim1td_Cont_DGE.txt')
data <- vroom(here('Data/Son2021/GSM5076376_Sim1td_Cont_DGE.txt'), delim = "\t", row_names = T) %>%
  as.matrix() %>% Matrix()


# Gene names to most current MGI symbols
row.names(data) <- StandardiseGeneNames(row.names(data))


# Remove duplicate rows
# Inspecting duplicates yields ~95% similarity between all duplicate genes w/ 94% sparisity
data <- MergeRowsByGeneNames(data)

# Create Seurat
df <- CreateSeuratObject(data, 'Son2021')
Idents(df) <- (df$orig.ident <- "Son2021")


# QC/Filtering
df$pctRibo_RNA <- PercentageFeatureSet(df, pattern = "^Rp[sl]")
df$pctMito_RNA <- PercentageFeatureSet(df, pattern = "^mt.")
df <- df[,df$nCount_RNA > 1000 & df$nFeature_RNA > 500 & df$pctMito_RNA < 10]
VlnPlot(df, features = c("nCount_RNA", "nFeature_RNA", "pctRibo_RNA", "pctMito_RNA"), ncol = 4) & scale_y_log10()


# Basic clustering to annotate classes
df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('nCount_RNA','pctRibo_RNA','pctMito_RNA')) %>%
  RunPCA() %>%
  RunUMAP(dims=1:25)  %>%
  FindNeighbors(dims=1:15) %>%
  FindClusters(resolution = 1)

df$C7_named <- mapvalues(Idents(df), 0:12, c("C7-1: GLU" * 10, "NA", "C7-1: GLU" * 2))


# Standardise metadata
Idents(df) <- (df$Dataset <- df$orig.ident)
df$SRA_ID <- "SRR13687659-62"
df$Sample_ID <- df$Dataset + "_" + df$SRA_ID
df$assay <- "Drop-seq"
df$disease <- "normal"
df$organism <- "Mus musculus"
df$sex <- "male"
df$suspension_type <- "cell"
df@meta.data <- df@meta.data[-6:-7]


#Save as rds
saveRDS(df,here('Output/Data/1-FromRaw/Son2021.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Son2021.h5ad"))