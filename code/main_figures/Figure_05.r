# ———————————
# FIGURE_05.R — R SCRIPT FOR GENERATING FIGURE 5 PANELS
# ———————————

library(tidyverse) # v.2.0.0
library(UpSetR) # v.1.4.0


suppa_results <- list(

  DLPFC = readRDS("DLPFC_SUPPA_Summary.rds") %>%
    mutate(
      cell_type = factor(cell_type, levels=c("GABA","GLU","OLIG")),
      splice.event = factor(splice.event, levels=c("RI", "SE", "A3", "A5", "AF", "MX", "AL"))
    ),
  
  OFC = readRDS("OFC_SUPPA_Summary.rds") %>%
    mutate(
      cell_type = factor(cell_type, levels=c("GABA","GLU","OLIG","AST","MG")),
      splice.event = factor(splice.event, levels=c("RI", "SE", "A3", "A5", "AF", "MX", "AL"))
    )

)


# —————————
# FIG.5B-5C 
# —————————
plot_fig_5B_5C <- function(region) {
  
  region_df <- suppa_results[[region]]
  
  splice_event_counts <- region_df %>%
    group_by(splice.event, cell_type) %>%
    summarise(
      n_isoforms = n_distinct(isoform),
      .groups = "drop"
    )
  
  ggplot(splice_event_counts, aes(x=splice.event, y=n_isoforms, fill=cell_type)) +
    geom_col(position = position_dodge(width = 0.8),width = 0.7) +
    labs(x="SUPPA SPLICE EVENTS", y="NUMBER OF ISOFORMS", fill="Cell type", title=region) +
    theme_bw() +
    theme(
      axis.text.x=element_text(angle=45, hjust=1),
      plot.title=element_text(face="bold", hjust=0.5)
    )

}

plot_fig_5B_5C("DLPFC")
plot_fig_5B_5C("OFC")


# ——————
# FIG.5D 
# ——————
plot_fig_5D <- function(region) {
  
  region_df <- suppa_results[[region]]
  
  tpm_props <- region_df %>%
    filter(!is.na(TPM_category)) %>%
    group_by(splice.event, TPM_category) %>%
    summarise(n=n_distinct(isoform), .groups="drop") %>%
    group_by(splice.event) %>% mutate(pct = n / sum(n) * 100)

  tpm_props$TPM_category <- factor(tpm_props$TPM_category, levels=c("No expression","Low expression","High expression"))

  ggplot(tpm_props, aes(x=splice.event, y=pct, fill=TPM_category)) +    
    geom_col(position="stack", width=0.7) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x="SPLICE EVENT", y="PROPORTION OF ISOFORMS (%)", fill="TPM", title=region) +
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      plot.title = element_text(face="bold", hjust=0.5)
    )
}

plot_fig_5D("DLPFC")
plot_fig_5D("OFC")


# ——————
# FIG.5E 
# ——————
plot_fig_5E <- function(region) {
  
  region_df <- suppa_results[[region]]
  
  structural_category_props <- region_df %>%
    group_by(splice.event, structural_category) %>%
    summarise(n=n_distinct(isoform), .groups="drop") %>%
    group_by(splice.event) %>% mutate(pct = n / sum(n) * 100)

  structural_category_props$structural_category <- factor(
    structural_category_props$structural_category,
    levels=c("Other","NNC","NIC","ISM","FSM")
  )
  
  ggplot(structural_category_props,aes(x=splice.event, y=pct, fill=structural_category)) +
    geom_col(position="stack", width=0.7) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x="SPLICE EVENT", y="PROPORTION OF ISOFORMS (%)", fill="STRUCTURAL CATEGORY", title=region) +
    theme_bw()
}

plot_fig_5E("DLPFC")
plot_fig_5E("OFC")


# ——————
# FIG.5F 
# ——————
plot_fig_5F <- function(region) {
  
  region_df <- suppa_results[[region]]
  
  nmd_props <- region_df %>%
    group_by(splice.event, NMD) %>%
    summarise(n=n_distinct(isoform), .groups="drop") %>%
    group_by(splice.event) %>% mutate(pct = n / sum(n) * 100)

  nmd_props$NMD <- factor(nmd_props$NMD, levels=c("Noncoding","TRUE","FALSE"))

  ggplot(nmd_props, aes(x=splice.event, y=pct, fill=NMD)) +    
    geom_col(position="stack", width=0.7) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x="SPLICE EVENT", y="PROPORTION OF ISOFORMS (%)", fill="NMD", title=region) +
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      plot.title = element_text(face="bold", hjust=0.5)
    )
}

plot_fig_5F("DLPFC")
plot_fig_5F("OFC")
