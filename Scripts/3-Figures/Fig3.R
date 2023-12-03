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

# Panel A
Fig3A <- function(df) {
  theme_custom <- function() {return(
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ))
  }
  palette <- c(
    "Adarb2-" = "grey70",
    "Adarb2+" = "#00BFC4",
    "Beine et al.\n(2022)" = "#FF2222"
  )
  df <- df %>% 
    RenameIdents("Avp-Tac1" = "Adarb2-") %>% 
    RenameIdents("Avp-Th" = "Adarb2-") %>% 
    RenameIdents("Crh-Nr3c1" = "Adarb2-") %>% 
    RenameIdents("Sst-Sfrp2" = "Adarb2-") %>% 
    RenameIdents("Sst-Calb2" = "Adarb2-") %>% 
    RenameIdents("Trh-Ghrh" = "Adarb2-") %>% 
    RenameIdents("Trh-Nfib" = "Adarb2-") %>% 
    RenameIdents("Trh-Omp" = "Adarb2-") %>% 
    RenameIdents("Trh-Tac1" = "Adarb2-") %>% 
    RenameIdents("Trh-Ucn3" = "Adarb2-") %>% 
    RenameIdents("Oxt-Foxp1" = "Adarb2-") %>% 
    
    RenameIdents("Asb4-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Avp-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Brs3-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Crh-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Ghrh-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Oxt-Adarb2" = "Adarb2+") %>%
    RenameIdents("Slc17a6-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Penk-Adarb2" = "Adarb2+") %>% 
    
    RenameIdents("Beine2022" = "Beine et al.\n(2022)") %>%
    (\(df) {levels(df) <- c("Adarb2-","Adarb2+","Beine et al.\n(2022)"); df})(.)
  
  p <- DimPlot(df, reduction = "mde") +
    scale_color_manual(name = "Cluster", values = palette) +
    theme_bw() + 
    theme_custom() 
  
  return(p)
  
}
Fig3A(df)

# Panel B
Fig3B <- function(df) {
  palette <- c(
    "Adarb2-" = "grey70",
    "Adarb2+" = "#00BFC4",
    "Beine et al.\n(2022)" = "#FF2222"
  )
  df <- df %>% 
    RenameIdents("Avp-Tac1" = "Adarb2-") %>% 
    RenameIdents("Avp-Th" = "Adarb2-") %>% 
    RenameIdents("Crh-Nr3c1" = "Adarb2-") %>% 
    RenameIdents("Sst-Sfrp2" = "Adarb2-") %>% 
    RenameIdents("Sst-Calb2" = "Adarb2-") %>% 
    RenameIdents("Trh-Ghrh" = "Adarb2-") %>% 
    RenameIdents("Trh-Nfib" = "Adarb2-") %>% 
    RenameIdents("Trh-Omp" = "Adarb2-") %>% 
    RenameIdents("Trh-Tac1" = "Adarb2-") %>% 
    RenameIdents("Trh-Ucn3" = "Adarb2-") %>% 
    RenameIdents("Oxt-Foxp1" = "Adarb2-") %>% 
    
    RenameIdents("Asb4-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Avp-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Brs3-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Crh-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Ghrh-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Oxt-Adarb2" = "Adarb2+") %>%
    RenameIdents("Slc17a6-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Penk-Adarb2" = "Adarb2+") %>% 
    
    RenameIdents("Beine2022" = "Beine et al.\n(2022)") %>%
    (\(df) {levels(df) <- c("Adarb2-","Adarb2+","Beine et al.\n(2022)"); df})(.)
    
  p <- VlnPlot(df, c("Adarb2","Reln","Snca","Ntng1"), ncol = 2) &
    scale_fill_manual(values = palette) &
    theme_bw() &
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 24),
      axis.text    = element_text(size = 22),
      plot.title   = element_text(size = 30, face = "italic", hjust = 0.5),
      panel.grid   = element_blank(),
      panel.border = element_rect(linewidth = 1.5),
      legend.position = "none"
    )
    
  return(p)
}
Fig3B(df)

# Panel E
Fig3E <- function(df) {
  theme_custom <- function() {return(
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5),
      panel.grid = element_blank(),
      strip.text = element_text(size = 14, face = "italic")
    ))
  }
  palette <- c(
    "Adarb2-" = "grey70",
    "Adarb2+" = "#00BFC4",
    "Beine et al.\n(2022)" = "#FF2222"
  )
  df <- df %>% 
    RenameIdents("Avp-Tac1" = "Adarb2-") %>% 
    RenameIdents("Avp-Th" = "Adarb2-") %>% 
    RenameIdents("Crh-Nr3c1" = "Adarb2-") %>% 
    RenameIdents("Sst-Sfrp2" = "Adarb2-") %>% 
    RenameIdents("Sst-Calb2" = "Adarb2-") %>% 
    RenameIdents("Trh-Ghrh" = "Adarb2-") %>% 
    RenameIdents("Trh-Nfib" = "Adarb2-") %>% 
    RenameIdents("Trh-Omp" = "Adarb2-") %>% 
    RenameIdents("Trh-Tac1" = "Adarb2-") %>% 
    RenameIdents("Trh-Ucn3" = "Adarb2-") %>% 
    RenameIdents("Oxt-Foxp1" = "Adarb2-") %>% 
    
    RenameIdents("Asb4-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Avp-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Brs3-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Crh-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Oxt-Adarb2" = "Adarb2+") %>%
    RenameIdents("Slc17a6-Adarb2" = "Adarb2+") %>% 
    RenameIdents("Penk-Adarb2" = "Adarb2+") %>% 
    
    RenameIdents("Beine2022" = "Beine et al.\n(2022)") %>%
    (\(df) {levels(df) <- c("Adarb2-","Adarb2+","Beine et al.\n(2022)"); df})(.)
  
  PlotCluster <- function(df) {
    GetData <- function(gene) {
      df = FeaturePlot(df, gene, reduction = "mde")$data
      names(df)[4] = "expr"
      df$gene = gene
      df$expr = df$expr / max(df$expr)
      return(df)
    }
    df <- GetData("Avp") %>%
      rbind(., GetData("Crh")) %>%
      rbind(., GetData("Oxt")) %>%
      rbind(., GetData("Penk")) %>%
      mutate(., gene = gene + "\n")
    levels(df$ident) <- "Adarb2+\n"
    df$gene <- factor(df$gene, levels = c(
      "Avp\n","Crh\n","Oxt\n","Penk\n"
    ))
    
    p <- ggplot(df) +
      geom_point(aes(MDE_1,MDE_2, color = expr, alpha = expr), size = 0.75) +
      scale_x_continuous(limits = c(1.75, 4.5)) +
      scale_y_continuous(limits = c(-1, 4)) +
      scale_color_gradient(low = "#CCCCCC", high = "#2222FF", 
                           limits = c(0, 1), na.value = "#2222FF",
                           name = "Normalized\nExpression") + 
      scale_alpha(range = c(0.5, 1), guide = "none") +
      facet_grid(gene ~ ident) +
      theme_bw() + 
      theme_custom() 
    
    return(p)
  }
  p <- wrap_plots(
    PlotCluster(df[,Idents(df) == "Adarb2+"]) +
      theme(
        legend.position = "none",
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,5), "mm")
      ),
    PlotCluster(df[,Idents(df) == "Beine et al.\n(2022)"]) +
      theme(
        axis.text.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
      ),
    ncol = 2
  )
  return(p)
}
Fig3E(df)

# Save .pngs ----
png(here("Output/Figures/Fig3A.png"), 
    res = 300, height = 17.5, width = 19, units = "cm")
Fig3A(df)
dev.off()

png(here("Output/Figures/Fig3B.png"), 
    res = 300, height = 30, width = 30, units = "cm")
Fig3B(df)
dev.off()

png(here("Output/Figures/Fig3E.png"), 
    res = 300, height = 24, width = 20, units = "cm")
Fig3E(df)
dev.off()

png(here("Output/Figures/Fig3D.png"), 
    res = 300, height = 12, width = 31, units = "cm")
Fig3D(df)
dev.off()
