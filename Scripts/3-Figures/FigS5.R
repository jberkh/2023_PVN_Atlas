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

# Panel A
FigS5A <- function(df) {
  theme_violin <- function() {
    return(theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank()
    ))
  }
  
  df <- df[,df$Dataset != "Beine2022"]
  
  Markers = c("Adarb2","Reln","Snca","Ntng1")
  p <- VlnPlot(df, Markers, ncol = 1) + 
    plot_layout(guides = "collect") &
    scale_y_continuous(limits = c(0,4)) &
    theme_bw() &
    theme_violin()
  return(p)
}
FigS5A(df)

# Panel B
FigS5B <- function(df) {
  id <- as.character(Idents(df))
  Idents(df) <- str_remove(id, "^.*-")
  df <- df %>%
    RenameIdents("Tac1" = "Adarb2-",
                 "Th" = "Adarb2-", 
                 "Nr3c1" = "Adarb2-", 
                 "Foxp1" = "Adarb2-", 
                 "Calb2" = "Adarb2-",
                 "Sfrp2" = "Adarb2-", 
                 "Ghrh" = "Adarb2-", 
                 "Nfib" = "Adarb2-", 
                 "Omp" = "Adarb2-", 
                 "Ucn3" = "Adarb2-") %>%
    RenameIdents("Adarb2" = "Adarb2+") %>%
    RenameIdents("Beine2022" = "Beine et al.\n(2022)")
  levels(df) <- sort(levels(df))
  
  p1 <- FeaturePlot(df, "Cacna1g", reduction = "mde") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, color = "white"),
      panel.grid = element_blank()
    )
  p2 <- VlnPlot(df, "Cacna1g") +
    theme_bw() +
    theme(
      plot.title = element_blank()
    )
  
  p <- patchwork::wrap_plots(
    p1, p2, ncol = 1, heights = c(8,2)
  )
  return(p)
}
FigS5B(df)


# Panel C
FigS5C <- function(df) {
  DMV <- readRDS(here("Output/Data/1-FromRaw/Tao2021.rds"))
  Receptors <- c("Avpr1a", "Avpr1b", "Avpr2", "Crhr"+1:2, "Ghrhr", 
                 "Opr" + c("m","k","d") + 1,"Oxtr", "Sstr" + 1:5, "Trhr")
  p <- c("Avpr1a", "Avpr1b", "Avpr2", "Crhr"+1:2, "Trhr", "Sstr" + 1:5, "Oxtr") %>%
    VlnPlot(readRDS(here("Output/Data/1-FromRaw/Tao2021.rds")), ., stack = T) & 
    theme_bw() & theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 14, face = "italic", angle = 90)
    )
  p <- DotPlot(DMV, "RNA", Receptors) +
    labs(title = "DMV Neuropeptide Receptors") +
    scale_x_discrete(name = "Neuropeptide Receptors") +
    scale_y_discrete(name = "Tao et al. (2022)\nClusters") +
    theme_bw()  +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      
    )
  return(p)  
}
FigS5C(df)


# Save PNG
png(here("Output/Figures/FigS5A.png"), 
    res = 450, height = 15, width = 30, units = "cm")
FigS5A(df)
dev.off()

png(here("Output/Figures/FigS5B.png"), 
    res = 450, height = 16, width = 16, units = "cm")
FigS5B(df)
dev.off()

png(here("Output/Figures/FigS5C.png"), 
    res = 450, height = 15, width = 20, units = "cm")
FigS5C(df)
dev.off()
