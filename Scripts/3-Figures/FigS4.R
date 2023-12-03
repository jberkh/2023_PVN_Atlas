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
# Panel B
FigS4A <- function(df) {
  Markers = c("Avp","Tac1","Th","Ddc")
  p <- FeaturePlot(df, Markers, reduction = "mde", combine = F)
  p <- patchwork::wrap_plots(p, ncol = 2) + 
    plot_layout(guides = "collect") &
    scale_color_gradient(low = "#D3D3D3", high = "#100DFE", na.value = "#100DFE", 
                         limits = c(0,4), labels = c(0:3, "4+")) &
    theme_bw() &
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "italic"),
      legend.position = "right",
    )
  return(p)
}
FigS4A(df)

# Panel B
FigS4B <- function(df) {
  df <- df[,Idents(df) %in% c("Sst-Sfrp2","Sst-Calb2")]
  Markers = c("Sst","Calb1","Calb2","Kcnip4","Gpr101","Sfrp2","Ar","Npy1r","Mc4r")
  p1 <- VlnPlot(df, Markers[1], ncol = 1) &
    theme_bw() &
    theme_violin()
  p2 <- VlnPlot(df, Markers[2:8], ncol = 1) &
    theme_bw() &
    theme_violin()
  p3 <- VlnPlot(df, Markers[9], ncol = 1) &
    scale_y_continuous(limits = c(0,3)) &
    theme_bw() &
    theme_violin()
  p <- wrap_plots(p1,p2,p3, heights = c(0.25,1.75,0.25))
  return(p)
}
FigS4B(df)

# Save PNG
png(here("Output/Figures/FigS4A.png"), 
    res = 450, height = 20, width = 22, units = "cm")
FigS4A(df)
dev.off()

png(here("Output/Figures/FigS4B.png"), 
    res = 300, height = 25, width = 12, units = "cm")
FigS4B(df)
dev.off()

