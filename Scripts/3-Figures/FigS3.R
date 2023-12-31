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
df <- readRDS(here("Output/Data/1-FromRaw/Lewis2020.rds")) %>%
  (\(df) {Idents(df) <- df$Lewis; df})(.) %>%
  FindSubCluster("Magnocellular") %>%
  RenameIdents("Magnocellular_0" = "Magnocellular 1") %>%
  RenameIdents("Magnocellular_1" = "Magnocellular 2") %>%
  (\(df) {levels(df) <- sort(levels(df)); df})(.)

# Panel A (UMAP)
FigS3A <- function(df) {
  p <- DimPlot(df, pt.size = 1.75, shape.by = "Fluorogold") + 
    labs(title = "Lewis et al. (2020)", shape = "Fluorogold", color = "Types") +
    scale_color_discrete(type = list(c("#E69F00", "#0072B2", "#009E73"))) +
    scale_shape_discrete(labels = c("FG-","FG+")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 12)
    )
  return(p)
}
FigS3A(df)

# Panel A (DotPlot)
FigS3B <- function(df) {
  p <- DotPlot(df, "RNA", c("Avp","Oxt","Calb1","Kcnmb4","Foxp1")) + 
    labs(title = "Lewis et al. (2020)") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size = 11),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      legend.title = element_text(size = 11)
    )
  return(p)
}
FigS3B(df)

# Save PNG
png(here("Output/Figures/FigS3A.png"), 
    res = 450, height = 16, width = 18, units = "cm")
FigS3A(df)
dev.off()

png(here("Output/Figures/FigS3B.png"), 
    res = 300, height = 12, width = 15, units = "cm")
FigS3B(df)
dev.off()

