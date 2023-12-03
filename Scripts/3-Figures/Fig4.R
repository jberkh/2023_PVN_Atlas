# Imports
set.seed(1, sample.kind = "Rejection")
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
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
df <- df[,df$Dataset != "Beine2022"]

# Panel A
Fig4A <- function(df) {
  Markers = c("Mc4r")
  p <- VlnPlot(df, Markers, ncol = 1) &
    scale_y_continuous(limits = c(0,2)) &
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, color = "white")
    )
  return(p)
}
Fig4A(df)

# Panel B
Fig4B <- function(df) {
  levels(df) <- sort(levels(df), T)
  Markers = c("Mc4r", "Brs3", "Npy1r", "Npy5r", "Vgf", "Asb4")
  
  p <- DotPlot(df, "RNA", Markers) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  return(p)
}
Fig4B(df)


# Save PNG
png(here("Output/Figures/Fig4A.png"), 
    res = 450, height = 15, width = 15, units = "cm")
Fig4A(df)
dev.off()

png(here("Output/Figures/Fig4B.png"), 
    res = 450, height = 12, width = 15, units = "cm")
Fig4B(df)
dev.off()

