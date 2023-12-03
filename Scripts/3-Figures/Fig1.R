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
df <- df[, Idents(df) != "Beine2022"]


# Panel B 
Fig1C <- function(df) {
  # Maximize color difference of adjacent clusters
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
    "Ghrh-Adarb2",
    "Trh-Ucn3", 
    "Penk-Adarb2", 
    "Sst-Calb2",
    "Brs3-Adarb2")
  p <- DimPlot(df, reduction = "mde") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.text = element_text(size = 14),
      legend.position = "none",
    )
  return(p)
}
Fig1C(df)


# Panel C
Fig1D <- function(df) {
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
               "Avp","Crh","Ghrh","Oxt","Penk","Sst","Trh",
               "Asb4","Adarb2","Tac1","Th","Brs3","Nr3c1","Foxp1",
               "Calb2","Sfrp2","Nfib","Omp","Ucn3")
  
  p <- VlnPlot(df, Markers, pt.size = 0, ncol = 1) & 
    theme_bw() & 
    theme_custom()
  return(p)
}
Fig1D(df)


# Panel D 
Fig1E <- function(df) {
  levels(df) <- sort(levels(df), T)
  Receptors <- c("Avpr1a", "Avpr1b", "Avpr2", "Crhr"+1:2, "Ghrhr", 
                 "Opr" + c("m","k","d") + 1,"Oxtr", "Sstr" + 1:5, "Trhr")
  p <- DotPlot(df, "RNA", Receptors) +
    scale_x_discrete(name = "Neuropeptide receptors") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    )
  return(p)
}
Fig1E(df)


# Save .pngs ----
# Fig 1A -> Inkscape
# Fig 1B -> Inkscape
png(here("Output/Figures/Fig1C.png"), 
    res = 450, height = 20, width = 20, units = "cm")
Fig1C(df)
dev.off()

png(here("Output/Figures/Fig1D.png"), 
    res = 300, height = 30, width = 18.5, units = "cm")
Fig1D(df)
dev.off()

png(here("Output/Figures/Fig1E.png"), 
    res = 450, height = 20, width = 24, units = "cm")
Fig1E(df)
dev.off()
