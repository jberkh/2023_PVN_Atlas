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
data <- vroom(here("Data/Xu2020/GSE148568_compiled_data.csv.gz"), row_names = T)
meta <- vroom(here("Data/Xu2020/GSE148568_cell_metadata_after_qc.csv.gz"), row_names = T)
row.names(meta) <- meta$sample_id


# Gene names to most current MGI symbols
data <- Matrix(as.matrix(data[,colnames(data) %in% row.names(meta)]))
row.names(data) <- StandardiseGeneNames(row.names(data))


# Remove duplicate rows
# Inspecting duplicates yields ~86% similarity between all duplicate genes w/ 89% sparisity
data <- MergeRowsByGeneNames(data)


# Create Seurat
df <- CreateSeuratObject(data, 'Xu2020')
Idents(df) <- (df$orig.ident <- "Xu2020")


# QC/Filtering
df$pctRibo_RNA <- PercentageFeatureSet(df, pattern = "^Rp[sl]")
df$pctMito_RNA <- PercentageFeatureSet(df, pattern = "^mt.")
VlnPlot(df, features = c("nCount_RNA", "nFeature_RNA", "pctRibo_RNA", "pctMito_RNA"), ncol = 4) & scale_y_log10()


# Basic clustering to annotate classes
df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('nCount_RNA','pctRibo_RNA','pctMito_RNA')) %>%
  RunPCA() %>%
  RunUMAP(dims=1:25)  %>%
  FindNeighbors(dims=1:15) %>%
  FindClusters(resolution = 2) %>%
  FindSubCluster("4", resolution = 0.75) %>%
  FindSubCluster("5", resolution = 0.75) %>%
  RenameIdents("4_0" = "4") %>%
  RenameIdents("4_1" = "12") %>%
  RenameIdents("5_0" = "5") %>%
  RenameIdents("5_1" = "13")


df$C7_named <- c("C7-1: GLU" * 9, "C7-2: GABA", "C7-1: GLU" * 4) %>% 
  mapvalues(Idents(df), 0:13, .)
df$Xu <- c("Oxt (Magno)", "Avp (Magno)", "Oxt (Magno)", 
           "Sst", "Crh-Scgn", "Avp (Parvo)", "Trh-Reln", 
           "Crh-Reln", "Trh-Bdnf", "Dlx1", "Avp (Magno)", 
           "Tcf7l2-Shox2", "Trh-Nfix", "Oxt (Parvo)"
           ) %>% mapvalues(Idents(df), 0:13, .)


# Standardise metadata
Idents(df) <- (df$Dataset <- df$orig.ident)
df$SRA_ID <- mapvalues(meta$plate_label, unique(meta$plate_label), "SRR" + 11539717:11539726)
df$Sample_ID <- df$Dataset + "_" + df$SRA_ID
df$assay <- "SMART-SCRB"
df$disease <- "normal"
df$organism <- "Mus musculus"
df$sex <- "male"
df$suspension_type <- "cell"
df@meta.data <- df@meta.data[-6:-7]


#Save as rds
saveRDS(df,here('Output/Data/1-FromRaw/Xu2020.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Xu2020.h5ad"))