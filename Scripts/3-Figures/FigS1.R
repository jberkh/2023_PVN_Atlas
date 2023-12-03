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

# Panel B
FigS1A <- function(df) {
  df$Dataset <- factor(df$Dataset, levels = sort(unique(df$Dataset)))
  p <- DimPlot(df, group.by = "Dataset", reduction = "mde") +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.text = element_text(size = 14),
    )
  return(p)
}
FigS1A(df)

# Panel B
FigS1B <- function(df) {
  Markers = c("Slc17a6","Slc32a1","Gad"+1:2)
  p <- FeaturePlot(df, Markers, reduction = "mde", combine = F)
  p <- patchwork::wrap_plots(p, ncol = 1) + 
    plot_layout(guides = "collect") &
    scale_color_gradient(low = "#D3D3D3", high = "#100DFE", limits = c(0,3)) &
    theme_bw() &
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "italic"),
      legend.position = "right",
      )
  return(p)
}
FigS1B(df)


# Panel C
FigS1C <- function(df) {
  theme_custom <- function() {
    return(theme(
      plot.title = element_text(face = "bold.italic", size = 18, hjust = 0.5),
      panel.grid = element_blank(),
    ))
  }
  theme_violin <- function() {
    return(theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 12)
    ))
  }
  Trh <- c("Brs3-Adarb2","Trh-Ghrh","Trh-Nfib","Trh-Omp","Trh-Tac1","Trh-Ucn3")
  
  p1 <- FeaturePlot(df, "Trh", reduction = "mde")$data %>% 
    ggplot() +
    ggtitle("") +
    geom_point(aes(MDE_1, MDE_2, color = Trh, alpha = Trh), size = 1) +
    scale_color_gradient(low = "#AAAAAA", high = "#0000FF", name = "Normalized\nExpression") +
    scale_alpha_continuous(range = c(0.1, 1), guide =  "none")+
    theme_bw() + 
    theme_custom()
  
  p2 <- VlnPlot(df[,Idents(df) %in% Trh], "Trh", split.by = "assay") &
    theme_bw() & 
    theme_custom() & 
    theme_violin()
  
  p <- wrap_plots(p1,p2, heights = c(20,6)) &
    theme(
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12),
    )
  
  return(p)
}
FigS1C(df)

# Save PNG
png(here("Output/Figures/FigS1A.png"), 
    res = 450, height = 20, width = 24, units = "cm")
FigS1A(df)
dev.off()

png(here("Output/Figures/FigS1B.png"), 
    res = 450, height = 30, width = 10, units = "cm")
FigS1B(df)
dev.off()

png(here("Output/Figures/FigS1C.png"), 
    res = 450, height = 23, width = 21, units = "cm")
FigS1C(df)
dev.off()
