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
data <- vroom(here('Data/Short2021/GSE174085_new_allcells111620.matrix.txt.gz'), delim='\t', col_names = T, row_names = T) %>%
  as.matrix() %>% Matrix()
cells <- colnames(data)
genes <- row.names(data) %>%
  (\(x) data.frame(
    "EnsemblID" = str_remove(x, "\\.[0-9]*-.*$"),
    "Gene" = StandardiseGeneNames(str_remove(x, "^ENSMUSG[0-9]*\\.[0-9]*-")))
   )()
data@Dimnames[[1]] <- genes$Gene


if (F){
  library(biomaRt)
  # Smartseq2 gene length normalization
  # Get mean gene lengths from biomaRt
  Len <- useMart("ensembl",dataset="mmusculus_gene_ensembl",
                 host = "https://apr2022.archive.ensembl.org")
  Len <- getBM(attributes=c("ensembl_gene_id","transcript_length"), mart = Len)
  Len <- filter(Len, ensembl_gene_id %in% row.names(data)) %>%
    group_by(ensembl_gene_id) %>%
    summarise("ensembl_gene_id" = unique(ensembl_gene_id),
              "transcript_length" = mean(transcript_length)) %>%
    filter(ensembl_gene_id != "")
  
  #Filter genes without gene length found
  data <- data[row.names(data) %in% Len$ensembl_gene_id,]
  Len <- as.numeric(mapvalues(row.names(data), Len$ensembl_gene_id, Len$transcript_length))
  
  data <- data / Len
  data@x <- as.numeric(as.integer(median(Len) * data@x))
  }


# Remove duplicate rows
# Inspecting duplicates yields ~73% similarity between all duplicate genes
data <- MergeRowsByGeneNames(data)
genes <- genes[!duplicated(genes$Gene),]


# Create Seurat
df <- CreateSeuratObject(data, 'Short2021')
df$SRA_ID <- "SRR14467974-8575"


# QC/Filtering
df$pctRibo_RNA <- PercentageFeatureSet(df, pattern = "^Rp[sl]")
df$pctMito_RNA <- PercentageFeatureSet(df, pattern = "^mt.")
df <- df[,df$nCount_RNA > 1000 & df$nFeature_RNA > 500]
VlnPlot(df, features = c("nCount_RNA", "nFeature_RNA", "pctRibo_RNA", "pctMito_RNA"), ncol = 4) & scale_y_log10()


# Basic clustering to annotate classes
df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c("nFeature_RNA", 'pctMito_RNA','pctRibo_RNA')) %>%
  RunPCA() %>%
  RunUMAP(dims=1:25)


# Dataset too small to annotate accurately
df$C7_named <- "NA"


# Standardise metadata
Idents(df) <- (df$Dataset <- df$orig.ident)
df$Sample_ID <- df$Dataset + "_" + df$SRA_ID
df$assay <- "SMARTseq2"
df$disease <- "normal"
df$organism <- "Mus musculus"
df$sex <- "male"
df$suspension_type <- "cell"


#Save as rds
saveRDS(df,here('Output/Data/1-FromRaw/Short2021.rds'))
AnnData(
  X = t(df@assays$RNA@counts), 
  obs = df@meta.data, 
  var = df@assays$RNA@meta.features) %>%
  write_h5ad(., here("Output/Data/1-FromRaw/Short2021.h5ad"))