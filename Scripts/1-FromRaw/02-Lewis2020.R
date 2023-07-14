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
data <- vroom(here('Data/Lewis2020/expression_matrix.tsv'), delim = "\t", row_names = T) %>%
  as.matrix() %>% Matrix()
cells <- vroom(here('Data/Lewis2020/cellAnnotation.tsv'), delim = '\t', row_names = T)
genes <- vroom(here('Data/Lewis2020/featureAnnotation.tsv'), delim = '\t', row_names = T) %>%
  (\(x) data.frame(
    "EnsemblID" = x[,1],
    "Gene" = StandardiseGeneNames(x[,5]))
  )()

# FKPM -> Raw counts
data <- t(data) %>% 
  (\(x) x / rowSums(x))() %>% 
  (\(x) x * cells$total_mass)() %>% t()
data@x <- as.numeric(as.integer(data@x))
data@Dimnames <- list(genes$Gene, cells$cell_id)


# Remove duplicate rows
# Inspecting duplicates yields ~75% similarity between all duplicate genes
data <- MergeRowsByGeneNames(data)
genes <- genes[!duplicated(genes$Gene),]


# Create Seurat
df <- CreateSeuratObject(data, 'Lewis2020')
df@assays$RNA@meta.features <- cbind(df@assays$RNA@meta.features, genes)
Idents(df) <- (df$orig.ident <- "Lewis2020")
df$SRA_ID <- "NA"


# QC/Filtering
df$pctRibo_RNA <- PercentageFeatureSet(df, pattern = "^Rp[sl]")
df$pctMito_RNA <- PercentageFeatureSet(df, pattern = "^mt.")
df <- df[,df$nCount_RNA > 1000 & df$nFeature_RNA > 1000 & df$pctMito_RNA < 10]
VlnPlot(df, features = c("nCount_RNA", "nFeature_RNA", "pctRibo_RNA", "pctMito_RNA"), ncol = 4) & scale_y_log10()


# Basic clustering to annotate classes
df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('nCount_RNA','pctRibo_RNA','percent_mito')) %>%
  RunPCA() %>%
  RunUMAP(dims=1:25)  %>%
  FindNeighbors(dims=1:25) 

df$C7_named <- "C7-1: GLU"


# Standardise metadata
Idents(df) <- (df$Dataset <- df$orig.ident)
df$Sample_ID <- "Lewis2020_NO_SRA"
df$assay <- "SMARTseq2"
df$disease <- "normal"
df$organism <- "Mus musculus"
df$sex <- "male"
df$suspension_type <- "cell"
df$Fluorogold <- cells[colnames(df),"Fluorogold"]
df$Lewis <- cells[colnames(df),"CellType"]

#Save as rds
saveRDS(df,here('Output/Data/1-FromRaw/Lewis2020.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Lewis2020.h5ad"))