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
  
  #Package w/ own functions
  library(RExtra)
  library(SeuratExtra)
})


# Get data
df <- readRDS(here("Output/Data/2-Integration/4-Final.rds"))
df <- df[,Idents(df) != "Beine2022"]

theme_violin <- function() {
  return(theme(
    axis.title   = element_blank(),
    axis.text.x  = element_blank(),
    plot.title   = element_blank(),
    plot.margin = unit(c(0,1,0,1), "mm"),
    panel.grid   = element_blank(),
    legend.position = "none"
  ))
}
as_dataframe <- function(df, genes, subtypes) {
  df <- t(df[genes,]@assays$RNA@data) %>%
    as.data.frame() %>%
    cbind(., "ident" = as.character(Idents(df))) %>%
    cbind(., "dataset" = df$Dataset) %>%
    filter(ident %in% subtypes)
  return(df)
}

pallette  <- c("None" = "#AAAAAA",
               "Avp"  = "#E69F00", "Crh" = "#F0E442", 
               "Ghrh" = "#D55E00", "Oxt" = "#0072B2", 
               "Penk" = "#56B4E9", "Sst" = "#CC79A7",
               "Trh"  = "#009E73")

# Panel A
Fig2A <- function(df) {
  df$Neuropeptide <- str_remove(Idents(df), "-.*$")
  df$Neuropeptide[df$Neuropeptide == "Brs3"] <- "Trh"
  df$Neuropeptide[df$Neuropeptide %in% c("Asb4","Slc17a6")] <- "None"
  df$Neuropeptide <- factor(df$Neuropeptide, c("None", sort(unique(df$Neuropeptide))[-4]))
  p <- DimPlot(df, group.by = "Neuropeptide", reduction = "mde") +
    scale_color_manual(values = pallette) +
    guides(color = guide_legend(byrow = T, override.aes = list(size = 3))) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
      )
  return(p)
}
Fig2A(df)

# Panel B
Fig2B <- function(df) {
  Neuropeptides <- c("None","Avp" * 3,"Trh","Crh" * 2, "Ghrh","Oxt" * 2,"Penk","None","Sst" * 2,"Trh" * 5)
  Pallette  <- pallette
  VarFeat <- VariableFeatures(df)
  Expr <- log1p(AverageExpression(df, "RNA", VarFeat, slot = "data")$RNA)
  Cor <- cor(Expr, Expr)
  Dend <- dendsort::dendsort(hclust(dist(Cor)))
  
  p <- Heatmap(
    Cor,
    
    left_annotation = rowAnnotation(
      df = data.frame(Neuropeptide = Neuropeptides), 
      col = list("Neuropeptide" = Pallette)),
    row_title = "Clusters",
    cluster_rows = Dend,
    
    top_annotation = columnAnnotation(
      df = data.frame(Neuropeptide = Neuropeptides),
      col = list("Neuropeptide" = Pallette),
      show_legend = F),
    column_title = "Clusters",
    cluster_columns = Dend,
    
    heatmap_legend_param = list(title = "Correlation")
  )
  return(p)
}
Fig2B(df)

# Violins
# Panel C
Fig2C <- function(df) {
  Avp <- c("Avp-Adarb2","Avp-Tac1","Avp-Th")
  Oxt <- c("Oxt-Adarb2","Oxt-Foxp1")
  isAvpOxt <- Idents(df) %in% c(Avp,Oxt)
  Markers <- c("Calb1","Kcnmb4","Foxp1")
  ColorAvp <- as.character(pallette["Avp"])
  ColorOxt <- as.character(pallette["Oxt"])
  p <- df[,isAvpOxt] %>%
    VlnPlot(., Markers, ncol = 1) &
    scale_fill_manual(values = c(ColorAvp * 3, ColorOxt * 2)) &
    theme_bw() & 
    theme_violin()
  
  return(p)
}
Fig2C(df)


# Panel D
Fig2D <- function(df) {
  Crh <- c("Crh-Adarb2","Crh-Nr3c1")
  isCrh <- Idents(df) %in% Crh
  Markers <- c("Agtr1a","Avp","Nr3c1","Scgn")
  Color <- as.character(pallette["Crh"])
  p1 <- df[,isCrh] %>%
    VlnPlot(., Markers[1], ncol = 1) &
    scale_fill_manual(values = Color * 4) &
    scale_y_continuous(name = "Expression", limits = c(0,3)) &
    theme_bw() & 
    theme_violin()
  p2 <- df[,isCrh] %>%
    VlnPlot(., Markers[2], ncol = 1) &
    scale_fill_manual(values = Color * 4) &
    scale_y_continuous(name = "Expression", limits = c(0,6)) &
    theme_bw() & 
    theme_violin()
  p3 <- df[,isCrh] %>%
    VlnPlot(., Markers[3:4], ncol = 1) &
    scale_fill_manual(values = Color * 4) &
    scale_y_continuous(name = "Expression", limits = c(0,3)) &
    theme_bw() & 
    theme_violin()
  p <- wrap_plots(p1,p2,p3, heights = c(1,1,2))
  
  return(p)
}
Fig2D(df)

# Panel E
Fig2E <- function(df) {
  Trh <- c("Asb4-Adarb2","Trh-Ghrh","Trh-Nfib","Trh-Tac1","Trh-Omp","Trh-Ucn3")
  isTrh <- Idents(df) %in% Trh
  Markers <- c("Agtr1a","Thrb")
  Color <- as.character(pallette["Trh"])
  p <- df[,isTrh] %>%
    VlnPlot(., Markers, ncol = 1) &
    scale_fill_manual(values = Color * 6) &
    theme_bw() &
    theme_violin()
  
  return(p)
}
Fig2E(df)

# Save .pngs
png(here("Output/Figures/Fig2A.png"), 
    res = 450, height = 17, width = 15, units = "cm")
Fig2A(df)
dev.off()

png(here("Output/Figures/Fig2B.png"), 
    res = 300, height = 20, width = 22.5, units = "cm")
Fig2B(df)
dev.off()

png(here("Output/Figures/Fig2C.png"), 
    res = 300, height = 10, width = 25, units = "cm")
Fig2C(df)
dev.off()

png(here("Output/Figures/Fig2D.png"), 
    res = 450, height = 40/3, width = 10.75, units = "cm")
Fig2D(df)
dev.off()

png(here("Output/Figures/Fig2E.png"), 
    res = 450, height = 20/3, width = 29.5, units = "cm")
Fig2E(df)
dev.off()
