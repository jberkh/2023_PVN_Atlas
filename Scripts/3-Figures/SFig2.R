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
  library(reticulate)
  use_condaenv("2023_PVN_Atlas")
  
  #Package w/ own functions
  library(RExtra)
  library(SeuratExtra)
})

# Get data
df <- readRDS(here("Output/Data/2-Integration/4-Final.rds"))

# Panel A-D
SFig2 <- function(df) {
  Markers = c("Slc17a6","Slc32a1","Gad"+1:2)
  p <- FeaturePlot(df, Markers, reduction = "mde", ncol = 2) &
    scale_color_gradient(low = "#D3D3D3", high = "#100DFE", limits = c(0,3)) &
    theme_bw() &
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "italic"),
      )
  return(p)
}
SFig2(df)


# Save PNG
png(here("Output/Figures/SFig2.png"), 
    res = 450, height = 25, width = 27.5, units = "cm")
SFig2(df)
dev.off()

