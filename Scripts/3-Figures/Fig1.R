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
df <- df[, Idents(df) != "Beine2022"]


# Panel B 
Fig1B <- function(df) {
  # Maximize color difference of adjecent clusters
  levels(df) <- c(
    "Crh-Adarb2", 
    "Avp-Tac1", 
    "Oxt-Adarb2", 
    "Trh-Nfib", 
    "Trh-Tac1",
    "Avp-Adarb2",
    "Sst-Sfrp2", 
    "Trh-Ghrh",
    "Slc17a6-Adarb2",
    "Oxt-Foxp1", 
    "Asb4-Adarb2", 
    "Crh-Nr3c1",
    "Trh-Omp", 
    "Avp-Th", 
    "Brs3-Adarb2",
    "Trh-Ucn3", 
    "Penk-Adarb2", 
    "Sst-Calb2")
  p <- DimPlot(df, reduction = "mde") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.text = element_text(size = 14),
      legend.position = "none",
    )
  return(p)
}
Fig1B(df)


# Panel C
Fig1C <- function(df) {
  theme_custom <- function() {
    return(theme(
      axis.title   = element_blank(),
      axis.text.x  = element_blank(),
      plot.title   = element_blank(),
      plot.margin = unit(c(0,1,0,1), "mm"),
      panel.grid   = element_blank(),
      legend.position = "none"
    ))
  }
  Markers <- c("Sim1","Otp",
               "Avp","Crh","Oxt","Penk","Sst","Trh",
               "Asb4","Adarb2","Tac1","Th","Brs3","Nr3c1","Foxp1",
               "Calb2","Sfrp2","Ghrh","Nfib","Omp","Ucn3")
  
  p <- VlnPlot(df, Markers, pt.size = 0, ncol = 1) & 
    theme_bw() & 
    theme_custom()
  return(p)
}
Fig1C(df)


# Panel D 
Fig1D <- function(df) {
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
Fig1D(df)


# Save .pngs ----

# Fig 1A -> Inkscape

png(here("Output/Figures/Fig1B.png"), 
    res = 450, height = 20, width = 20, units = "cm")
Fig1B(df)
dev.off()

png(here("Output/Figures/Fig1C.png"), 
    res = 300, height = 30, width = 18.5, units = "cm")
Fig1C(df)
dev.off()

png(here("Output/Figures/Fig1D.png"), 
    res = 450, height = 23, width = 21, units = "cm")
Fig1D(df)
dev.off()
