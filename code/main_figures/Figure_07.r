# ———————————
# FIGURE_07.R — R SCRIPT FOR GENERATING FIGURE 7 PANELS
# ———————————

library(tidyverse) # v.2.0.0
library(ggrepel) # v.0.9.8



# GvG = GABA vs GLU
# NvN = Neurons (GABA+GLU) vs Nonneurons (OLIG)

iso_switch_results <- list(
  DLPFC_GvG = readRDS("DLPFC_Switch_GvG_Summary.rds"),
  DLPFC_NvN = readRDS("DLPFC_Switch_NvN_Summary.rds"),
  OFC_GvG = readRDS("OFC_Switch_GvG_Summary.rds"),
  OFC_NvN = readRDS("OFC_Switch_NvN_Summary.rds")  
)


# ——————
# FIG.7A
# ——————
plot_fig_7A <- function(comparison, data_list=iso_switch_results) {
  
    region_df <- data_list[[comparison]]

    top_genes <- region_df %>%
        filter(color.group=="Significant") %>%
        group_by(direction = ifelse(delta.prop<0,"neg","pos")) %>%
        arrange(t.test.p.value) %>% slice_head(n=10) %>% ungroup()

    p <- ggplot(DF, aes(x = delta.prop, y = -log10(t.test.p.value), color=color.group)) +
        geom_point(alpha=0.75, size=3) +
        geom_text_repel(
            data = top_genes,
            aes(label = gene.symbol),
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
        theme_bw() + theme(aspect.ratio=1) + 
        labs(title=comparison, x="DELTA ISOFORM PROPORTION", 
            y = expression(-log[10](t.test.p.value)), color="SIGNIFICANCE")

    return(p)

}

plot_fig_7A("DLPFC_GvG")
plot_fig_7A("DLPFC_NvN")
plot_fig_7A("OFC_GvG")
plot_fig_7A("OFC_NvN")
