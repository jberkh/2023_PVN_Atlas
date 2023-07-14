library(tidyverse)
library(Matrix)
library(vroom)
library(here)

library(Seurat)

library(RExtra)
library(SeuratExtra)

# Tao = DMV
df <- vroom(here("Data/Tao2021/GSE172411_dge.tsv.gz"), row_names = T) %>%
  CreateSeuratObject(., "Tao2021") %>%
  (\(df) {df[,df$nFeature_RNA > 5000]})(.)


df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures(n.features = 1500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:19) %>%
  FindNeighbors(reduction = "pca", dims = 1:19) %>%
  FindClusters(resolution = 1.2) %>%
  RenameIdents("0" = "Nppb") %>%
  RenameIdents("1" = "Trpv1") %>%
  RenameIdents("2" = "Pdyn") %>%
  RenameIdents("3" = "Grp") %>%
  RenameIdents("4" = "Atf3") %>%
  RenameIdents("5" = "Gm4881") %>%
  RenameIdents("6" = "Hypoglossal") %>%
  RenameIdents("7" = "Cck")

# Remove hypoglossal nucleus -> Not DMV
df <- df[,Idents(df) != "Hypoglossal"]

saveRDS(df, here("Output/Data/1-FromRaw/Tao2021.rds"))
