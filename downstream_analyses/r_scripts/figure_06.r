# R SCRIPT FOR GENERATING PANELS FOR FIGURE 6 AND THEIR SUPPLEMENTAL FIGURES

library(ggtranscript) # v.1.0.0
library(rtracklayer) # v.1.66.0
library(tidyverse) # v.2.0.0
library(ggrepel) # v0.9.6


# ——— MAIN FIGURES ———

# Location of input files: repo/isoform_annotation/data

# ——————
# FIG.6A - DENSITY OF NOVEL ISOFORMS PER DEGs
# ——————

ref <- readRDS("DLPFC_reference.rds")
DEG_GvG <- readRDS("DLPFC_DEG_GvG_summary.rds") # GvG = GABA vs GLU
DEG_NvN <- readRDS("DLPFC_DEG_NvN_summary.rds") # NvN = Neuron vs Nonneuron (OLIG)

generate_novel_density <- function(DEG, ref) {
  
  # FILTER FOR SIGNIFICANT DEGs
  DEG <- DEG %>% filter(significance == "Significant")

  # ASSIGN NEW STRUCTURAL CATEGORIES
  ref <- ref %>%
    mutate(group = case_when(
      structural_category %in% c("full-splice_match", "incomplete-splice_match") ~ "Known",
      structural_category %in% c("novel_in_catalog", "novel_not_in_catalog") ~ "Novel",
      TRUE ~ "Other"
    ))
  
  # SUMMARIZE NUMBER OF ISOFORMS PER GENE
  iso_summary <- ref %>%
    group_by(gene.id) %>%
    summarise(
      total_isoforms = n(),
      novel_isoforms = sum(group == "Novel"),
      .groups = "drop"
    ) %>%
    mutate(novel_ratio = novel_isoforms / total_isoforms)
  
  # MERGE WITH DEG TABLE
  DEG <- DEG %>%
    left_join(iso_summary, by = "gene.id") %>%
    filter(total_isoforms >= 2)
  
  # DENSITY PLOT
  p <- ggplot(DEG, aes(x = novel_ratio)) +
    geom_density(fill="purple", size = 2) +
    geom_vline(
      xintercept = mean(DEG$novel_ratio, na.rm = TRUE),
      linetype = "dashed",
      size = 1.25
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    theme_bw()
  
  return(p)
}

plot_fig6A_GvG <- generate_novel_density(DEG_GvG, ref)
plot_fig6A_NvN <- generate_novel_density(DEG_NvN, ref)

# ——————
# FIG.6B - PROPORTION OF KNOWN AND NOVEL TRANSCRIPTS IN EACH DTE CELL TYPE COMPARISON
# ——————

DET_GvG <- readRDS("DLPFC_DET_GvG_summary.rds")
DET_NvN <- readRDS("DLPFC_DET_NvN_summary.rds")

generate_novel_DET_prop <- function(DET) {
  
  # ASSIGN NEW STRUCTURAL CATEGORIES
  DET <- DET %>%
    mutate(group = case_when(
      structural_category %in% c("full-splice_match", "incomplete-splice_match") ~ "Known",
      structural_category %in% c("novel_in_catalog", "novel_not_in_catalog")     ~ "Novel",
      TRUE                                                                       ~ "Other"
    ))
  
  # FILTER FOR SIGNIFICANT DETs
  DET <- DET %>% filter(group != "Other", significance == "Significant")
  
  # SUMMARIZE KNOWN AND NOVEL PROPORTIONS
  plot_df <- DET %>%
    count(group) %>%
    mutate(prop = n / sum(n) * 100)
  
  # STACKED BARPLOT 
  p <- ggplot(plot_df, aes(x = "Isoforms", y = prop, fill = group)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = c(
      "Known" = "#140e37ff",
      "Novel" = "#bab0f4ff"
    )) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 100)
    ) +
    theme_bw() +
    labs(x = NULL, y = NULL)
  
  return(p)
}

plot_fig6B_GvG <- generate_novel_DET_prop(DET_GvG)
plot_fig6B_NvN <- generate_novel_DET_prop(DET_NvN)


# ——————
# FIG.6C - NUMBER OF DEGs WITH/WITHOUT DTE ACROSS CELL TYPE COMPARIONS
# ——————

overlaps <- readRDS("DEG_DET_overlaps.rds")

generate_DEG_DET_overlap <- function(overlaps, conditions) {
  
  # FILTER FOR DLPFC AND CELL TYPE CONDITIONS
  overlaps <- overlaps %>% filter(region == "DLPFC", condition %in% conditions)

  overlaps$group <- factor(overlaps$group, levels = c("DGE_only", "DTE_only", "Both"))
  
  # STACKED BARCHART
  p <- ggplot(overlaps, aes(x = condition, y = n_genes, fill = group)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = c(
      "Both"     = "#d7d7d7ff",
      "DTE_only" = "#ae9cf7ff",
      "DGE_only" = "#5b22ddff"
    )) +
    scale_y_continuous(limits = c(0, 3000)) +
    labs(
      x = "Condition",
      y = "Number of Genes",
      fill = "Group"
    ) +
    theme_bw()
  
  return(p)
}

plot_fig6C_GvG <- generate_DEG_DET_overlap(overlaps, c("GABA", "GLU"))
plot_fig6C_NvN <- generate_DEG_DET_overlap(overlaps, c("Neuron", "OLIG"))



# ——————
# FIG.6D - MAGNITUDE OF ISOFORM FRACTION ACROSS CELL TYPE COMPARIONS (DTU ANALYSIS)
# ——————

DTU_GvG <- readRDS("DLPFC_DTU_GvG_summary.rds")
DTU_NvN <- readRDS("DLPFC_DTU_NvN_summary.rds")

generate_volcano_plot <- function(DTU) {
  
  # ASSIGN NOVEL AND KNOWN ISOFORM CLASS
  DTU <- DTU %>%
    mutate(isoform_class = case_when(
      structural_category %in% c("FSM", "ISM") ~ "Known",
      structural_category %in% c("NIC", "NNC") ~ "Novel",
      TRUE ~ "Other"
    ))
  
  # ASSIGN COLOR GROUPS FOR PLOTTING
  DTU <- DTU %>%
    mutate(color_group = case_when(
      significance != "Significant" ~ "Nonsignificant",
      isoform_class == "Known" ~ "Known",
      isoform_class == "Novel" ~ "Novel",
      TRUE ~ "Other"
    ))
  
  # SET LEVELS AND CLEAN UP DATA
  DTU$color_group <- factor(DTU$color_group, levels = c("Nonsignificant", "Known", "Novel"))
  DTU <- DTU %>% arrange(color_group) %>% na.omit()
  
  # IDENTIFY TOP SIGNIFICANT DTU EVENTS
  sig_points <- DTU %>% filter(significance == "Significant")
  top_right <- sig_points %>% top_n(15,  DTU_dIF)
  top_left  <- sig_points %>% top_n(-15, DTU_dIF)
  top_points <- bind_rows(top_right, top_left)
  
  # VOLCANO PLOT
  p <- ggplot(DTU, aes(x = DTU_dIF, y = -log10(DTU_qval))) +
    geom_point(aes(color = color_group), size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(
      data = top_points,
      aes(label = gene_symbol),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5
    ) +
    labs(
      x = "Δ Isoform Fraction (dIF)",
      y = "-log10(q-value)",
      color = "Structural Category"
    ) +
    theme_bw()
  
  return(p)
}

plot_fig6D_GvG <- generate_volcano_plot(DTU_GvG)
plot_fig6D_NvN <- generate_volcano_plot(DTU_NvN)


# ——————
# FIG.6E - VISUALIZE ISOFORMS AND ITS DGE, DTE AND DTU RESULTS
# ——————

# EXAMPLE: TWF1 (4 steps)

load("TWF1.RData") # CONTENTS: gene_expr, iso_expr, iso_prop, isoforms

# 1. PLOT GENE EXPRESSION BETWEEN NEURONS AND OLIG FOR TWF1
TWF1_gene_expression <- data.frame(
  group = c("Neurons", "OLIG"),
  mean = c(gene_expr$AVG.neuron[1], gene_expr$AVG.olig[1]),
  sem = c(gene_expr$SEM.neuron[1], gene_expr$SEM.olig[1])
)

plot_fig6E_geneExpr <- ggplot(TWF1_gene_expression, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15, linewidth = 0.6) +
  labs(x = "", y = "Expression") +
  theme_bw()

# 2. PLOT ISOFORM EXPRESSION 
iso_expr_long <- iso_expr %>%
  select(transcript_id, AVG.neuron, AVG.olig, SEM.neuron, SEM.olig) %>%
  pivot_longer(cols = c(AVG.neuron, AVG.olig, SEM.neuron, SEM.olig), names_to = c(".value", "group"), names_sep = "\\.")

plot_fig6E_isoExpr <- ggplot(iso_expr_long, aes(x = transcript_id, y = AVG, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = AVG - SEM, ymax = AVG + SEM), width = 0.2, position = position_dodge(width = 0.7)) +
  labs(x = "Isoform", y = "Expression", fill = "") +
  theme_bw()

# 3. PLOT ISOFORM PROPORTIONS
iso_prop_long <- iso_prop %>%
  pivot_longer(cols = c(Neuron.prop, Nonneuron.prop), names_to = "group", values_to = "prop") %>%
  mutate(transcript.id = fct_reorder(transcript.id, prop * ifelse(group == "Neuron.prop", 1, NA), .fun = max, na.rm = TRUE))

plot_fig6E_isoProp <- ggplot(iso_prop_long, aes(x = transcript.id, y = prop * 100, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("Neuron.prop" = "blue", "Nonneuron.prop" = "red"),
    labels = c("Neuron.prop" = "Neuron", "Nonneuron.prop" = "Non-neuron")
  ) +
  labs(x = "Isoform", y = "Proportion (%)", fill = "") +
  theme_bw()


# 4. PLOT ISOFORMS WITH GGTRANSCRIPT
gene_of_interest <- "TWF1"

gene_annotation_from_gtf <- isoforms %>% filter(!is.na(gene_id), gene_id==gene_of_interest)

gene_annotation_from_gtf <- gene_annotation_from_gtf %>% 
    select(seqnames, start, end, strand, type, gene_id, transcript_id, transcript_biotype)

gene_of_interest_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")

exons_rescaled <- shorten_gaps(
  exons = gene_of_interest_exons, 
  introns = to_intron(gene_of_interest_exons, "transcript_id"), 
  group_var = "transcript_id"
)

rescaled_genes <- exons_rescaled %>% dplyr::filter(type == "exon") 
rescaled_genes_introns <- exons_rescaled %>% dplyr::filter(type == "intron") 

plot_fig6E_isoforms <- rescaled_genes %>%
    ggplot(aes(xstart = start, xend = end, y = transcript_id)) +
    geom_range(aes(fill = transcript_biotype, color = transcript_biotype), height=0.8, size=0.5, alpha=0.5) +
    geom_intron(data = rescaled_genes_introns, aes(strand = strand), arrow.min.intron.length = 10000) + 
    theme_bw()