# ———————————
# FIGURE_06.R — R SCRIPT FOR GENERATING FIGURE 6 PANELS
# ———————————

library(tidyverse) # v.2.0.0
library(ggrepel) # v.0.9.8

# DGE = Differential gene expression 
# DTE = Differential transcript expression 
# DTU = Differential transcript usage
# GvG = GABA vs GLU
# NvN = Neurons (GABA+GLU) vs Nonneurons (OLIG)

deg_results <- list(
  DLPFC_GvG = readRDS("DLPFC_DGE_GvG_Summary.rds"),
  DLPFC_NvN = readRDS("DLPFC_DGE_NvN_Summary.rds"),
  OFC_GvG = readRDS("OFC_DGE_GvG_Summary.rds"),
  OFC_NvN = readRDS("OFC_DGE_NvN_Summary.rds")  
)

det_results <- list(
  DLPFC_GvG = readRDS("DLPFC_DTE_GvG_Summary.rds"),
  DLPFC_NvN = readRDS("DLPFC_DTE_NvN_Summary.rds"),
  OFC_GvG = readRDS("OFC_DTE_GvG_Summary.rds"),
  OFC_NvN = readRDS("OFC_DTE_NvN_Summary.rds")  
)

dtu_results <- list(
  DLPFC_GvG = readRDS("DLPFC_DTU_GvG_Summary.rds"),
  DLPFC_NvN = readRDS("DLPFC_DTU_NvN_Summary.rds"),
  OFC_GvG = readRDS("OFC_DTU_GvG_Summary.rds"),
  OFC_NvN = readRDS("OFC_DTU_NvN_Summary.rds")  
)

# ———————————
# FIGS. 6A-6B
# ———————————
plot_fig_6A_6B <- function(comparison) {
  
  comparison_df <- deg_results[[comparison]] %>% filter(significance=="Significant")
  
  deg_prop <- comparison_df %>%
    group_by(novel_iso_per_gene) %>%
    summarise(n_genes=n_distinct(gene.id), .groups="drop") %>%
    mutate(pct=n_genes / sum(n_genes) * 100) %>%
    filter(!is.na(novel_iso_per_gene))

  deg_prop$novel_iso_per_gene <- factor(deg_prop$novel_iso_per_gene, levels=c("0","1","2","3","4",">4"))
  
  ggplot(deg_prop, aes(x=factor(novel_iso_per_gene), y=pct)) +
    geom_col(width=0.7) +
    scale_y_continuous(labels=function(x) paste0(x, "%")) +
    labs(x="# NOVEL ISO PER GENE", y="PROPORTION OF DEGs (%)", title=comparison) +    
    theme_bw() 

}

plot_fig_6A_6B("DLPFC_GvG")
plot_fig_6A_6B("DLPFC_NvN")
plot_fig_6A_6B("OFC_GvG")
plot_fig_6A_6B("OFC_NvN")


# ——————
# FIG.6C
# ——————
plot_fig_6C <- function(comparison) {
  
  region_df <- det_results[[comparison]]
  region_df <- region_df %>% filter(significance=="Significant")
  
  isoform_prop <- region_df %>%
    filter(structural_group %in% c("Known", "Novel")) %>%
    group_by(structural_group) %>%
    summarise(n_isoforms=n_distinct(isoform.id), .groups="drop") %>%
    mutate(pct=n_isoforms / sum(n_isoforms) * 100, group="All isoforms")
  
  ggplot(isoform_prop, aes(x=group, y=pct, fill=structural_group)) +
    geom_col(width=0.5) +
    geom_text(aes(label=paste0(n_isoforms, "\n(", round(pct, 1), "%)")), 
        position=position_stack(vjust=0.5), size=4) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x=NULL, y="ISOFORMS (%)", fill="STRUCTURAL GROUP", title=comparison) +
    theme_classic(base_size=14) 
}

plot_fig_6C("DLPFC_GvG")
plot_fig_6C("DLPFC_NvN")
plot_fig_6C("OFC_GvG")
plot_fig_6C("OFC_NvN")


# ——————
# FIG.6D
# ——————
plot_fig_6D <- function(condition) {
  
  DEG <- deg_results[[condition]]
  DET <- det_results[[condition]]
  deg_genes <- DEG %>% filter(significance=="Significant") %>% distinct(gene.id)
  det_genes <- DET %>% filter(significance=="Significant") %>% distinct(gene.id)
  
  overlap <- full_join(mutate(deg_genes, DEG=TRUE), mutate(det_genes, DET=TRUE), by="gene.id") %>%
    mutate(
      DEG = ifelse(is.na(DEG), FALSE, DEG),
      DET = ifelse(is.na(DET), FALSE, DET),
      category = case_when(
        DEG & DET ~ "DEG + DET",
        DEG & !DET ~ "DEG only",
        !DEG & DET ~ "DET only"
      )
    )
  
  sig_gene_counts <- overlap %>%
    group_by(category) %>%
    summarise(n_genes=n_distinct(gene.id), .groups="drop")

  sig_gene_counts$category <- factor(sig_gene_counts$category, levels=c("DEG only","DET only","DEG + DET"))
  
  ggplot(sig_gene_counts, aes(x=condition, y=n_genes, fill=category)) +
    geom_col() +
    geom_text(aes(label=n_genes), position=position_stack(vjust=0.5), size=4) +
    labs(x=NULL, y="# SIGNIFICANT GENES", fill="SIGNIFICANCE", title=condition) +
    theme_bw()

}

plot_fig_6D("DLPFC_GvG")
plot_fig_6D("DLPFC_NvN")
plot_fig_6D("OFC_GvG")
plot_fig_6D("OFC_NvN")


# ——————
# FIG.6E
# ——————
plot_fig_6E <- function(dataset_name, data_list = dtu_results) {
  
  region_df <- data_list[[dataset_name]]
  region_df$plot_group <- ifelse(
    region_df$significance == "Significant",
    as.character(region_df$structural_group),
    "NonSignificant"
  )
  
  region_df$plot_group <- factor(
    region_df$plot_group, levels=c("Known", "Novel", "Other", "NonSignificant")
  )
  
  region_df$plot_group_order <- factor(
    region_df$plot_group, levels=c("NonSignificant", "Known", "Other", "Novel")
  )
  
  DF_plot <- region_df[order(DF$plot_group_order), ]
  
  top_genes <- bind_rows(
    region_df %>%
      filter(significance=="Significant", stageR.padj == 0) %>%
      mutate(direction=ifelse(dIF < 0, "neg", "pos")) %>%
      group_by(direction) %>% arrange(desc(abs(dIF))) %>% slice_head(n=5) %>% ungroup(),
    
    region_df %>%
      filter(significance=="Significant", stageR.padj > 0, stageR.padj < 1e-10) %>%
      mutate(direction=ifelse(dIF < 0, "neg", "pos")) %>%
      group_by(direction) %>% arrange(stageR.padj) %>% slice_head(n=15) %>% ungroup(),
    
    region_df %>%
      filter(significance=="Significant", stageR.padj >= 1e-10) %>%
      mutate(direction=ifelse(dIF < 0, "neg", "pos")) %>%
      group_by(direction) %>% arrange(stageR.padj) %>% slice_head(n=30) %>% ungroup()
  ) %>%
    distinct(isoform.id, .keep_all=TRUE)
  
  DF_plot <- na.omit(DF_plot)

  p <- ggplot(DF_plot, aes(x = dIF, y = -log10(stageR.padj), color = plot_group)) +
    geom_point(alpha = 0.75, size = 3) +
    
    geom_text_repel(
      data = top_genes,
      aes(label = gene_symbol),
      size = 4,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.size = 0.3,
      segment.color = "black",
      segment.alpha = 0.6,
      force = 5,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = c(
      "Known" = "#1f78b4",
      "Novel" = "#0b3c8c",
      "Other" = "#a6cee3",
      "NonSignificant" = "grey70"
    )) +
    
    theme_classic() +
    labs(
      title = dataset_name,
      x = "dIF",
      y = expression(-log[10](stageR.padj)),
      color = "Category"
    )
  
  return(p)
}

plot_fig_6E("DLPFC_GvG")
plot_fig_6E("DLPFC_NvN")
plot_fig_6E("OFC_GvG")
plot_fig_6E("OFC_NvN")
