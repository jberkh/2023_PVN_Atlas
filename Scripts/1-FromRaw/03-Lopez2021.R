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
cells <- vroom(here('Data/Lopez2021/GSM4914022_PVN_CTRL_barcodes.tsv.gz'), delim='\t',col_names = F)[[1]] %>% 
  as.character() %>% (\(x) x + "_SRR13081815")()
genes<- vroom(here('Data/Lopez2021/GSM4914022_PVN_CTRL_genes.tsv.gz'), delim='\t',col_names = c('EnsemblID','Gene'))
data <- readMM(here('Data/Lopez2021/GSM4914022_PVN_CTRL_matrix.mtx.gz'))
data@Dimnames <- list(genes$EnsemblID, cells)
df <- data

cells <- vroom(here('Data/Lopez2021/GSM4914023_PVN_STRESS_barcodes.tsv.gz'), delim='\t',col_names = F)[[1]] %>% 
  as.character() %>% (\(x) x + "_SRR13081816")()
genes<- vroom(here('Data/Lopez2021/GSM4914023_PVN_STRESS_genes.tsv.gz'), delim='\t',col_names = c('EnsemblID','Gene'))
data <- readMM(here('Data/Lopez2021/GSM4914023_PVN_STRESS_matrix.mtx.gz'))
data@Dimnames <- list(genes$EnsemblID, cells)
df <- cbind(df, data)


# EnsemblID -> Symbol
genes$Gene <- StandardiseGeneNames(genes$Gene)
df@Dimnames[[1]] <- genes$Gene


# Remove duplicate rows
# Inspecting duplicates yields ~93% similarity between all duplicate genes w/ ~97% sparsity
df <- MergeRowsByGeneNames(df)
genes <- genes[!duplicated(genes$Gene),]


# Create Seurat
df <- CreateSeuratObject(df, 'Lopez2021')
df@assays$RNA@meta.features <- df@assays$RNA@meta.features %>% 
  cbind(.,"EnsemblID" = genes$EnsemblID)


# QC/Filtering
df$pctRibo_RNA <- PercentageFeatureSet(df, pattern = "^Rp[sl]")
df$pctMito_RNA <- PercentageFeatureSet(df, pattern = "^mt.")
df <- df[,df$nCount_RNA > 1500 & df$nFeature_RNA > 750 & df$pctMito_RNA < 10]
VlnPlot(df, features = c("nCount_RNA", "nFeature_RNA", "pctRibo_RNA", "pctMito_RNA"), ncol = 4) & scale_y_log10()


# Basic clustering to annotate classes
df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('nCount_RNA','pctRibo_RNA','pctMito_RNA')) %>%
  RunPCA() %>%
  RunUMAP(dims=1:25)  %>%
  FindNeighbors(dims=1:25) %>%
  FindClusters(resolution = 0.25) %>%
  FindSubCluster(5, resolution = 0.25) %>%
  RenameIdents("5_0" = "5") %>%
  RenameIdents("5_1" = "10")
  
df$C7_named <- c("C7-2: GABA", "C7-4: Oligo+Precursor", "C7-1: GLU", 
                 "C7-3: Astro-Ependymal", "C7-7: Vascular", "C7-2: GABA", 
                 "C7-3: Astro-Ependymal", "C7-5: Immune", "C7-4: Oligo+Precursor", 
                 "C7-4: Oligo+Precursor", "C7-1: GLU", "C7-2: GABA", 
                 "C7-3: Astro-Ependymal", "C7-7: Vascular", "C7-7: Vascular", 
                 "C7-4: Oligo+Precursor", "C7-5: Immune"
                 ) %>% mapvalues(Idents(df), 0:16, .)


# Standardise metadata
Idents(df) <- (df$Dataset <- df$orig.ident)
df$SRA_ID <- str_remove(colnames(df),"^[ACGT]*-[0-9]*_")
df$Sample_ID <- df$Dataset + "_" + df$SRA_ID
df$assay <- "10x 3' v2"
df$disease <- "normal"
df$organism <- "Mus musculus"
df$sex <- "male"
df$suspension_type <- "cell"
df$Stress <- ifelse(df$SRA_ID == "SRR13081815", 0, 1)
df@meta.data <- df@meta.data[-6:-7]


# Save as rds
saveRDS(df,here('Output/Data/1-FromRaw/Lopez2021.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Lopez2021.h5ad"))