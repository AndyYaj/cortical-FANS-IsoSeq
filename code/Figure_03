# ———————————
# FIGURE_03.R — R SCRIPT FOR GENERATING FIGURE 2 PANELS
# ———————————


library(tidyverse) # v.2.0.0


# LOAD IN FL ISOFORM SUMMARY
fl_isoforms <- list(

    DLPFC = readRDS("DLPFC_FL_Transcript_summary.rds") %>%
        mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG"))),
  
    OFC = readRDS("OFC_FL_Transcript_summary.rds") %>%
        mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG","AST","MG")))

)

# LOAD IN PROTEIN-CODING PREDICTIONS
coding_potential <- list(

    DLPFC = readRDS("DLPFC_coding_potential.rds") %>%
        mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG"))),
  
    OFC = readRDS("OFC_coding_potential.rds") %>%
        mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG","AST","MG")))

)


# ——————
# FIG.3B - PROPORTION OF STRUCTURAL CATEGORIES ACROSS CELL TYPES
# ——————
plot_fig3B <- function(region){

    structural_props <- fl_isoforms[[region]] %>%
        group_by(cell_type, structural_category) %>%
        summarize(n=n(), .groups="drop") %>%
        group_by(cell_type) %>%
        mutate(prop = n / sum(n) * 100) %>%
        ungroup()

    structural_props$cell_type <- factor(
        structural_props$cell_type, levels=c("MG","AST","OLIG","GLU","GABA")
    )

    structural_props$structural_category <- factor(
        structural_props$structural_category, levels=c("Other","NNC","NIC","ISM","FSM")
    )

    ggplot(structural_props, aes(x=prop, y=cell_type, fill=structural_category)) +
        geom_bar(stat="identity", color="black") + theme_bw() + theme(aspect.ratio=1)
}

plot_fig3B("DLPFC")
plot_fig3B("OFC")


# ——————
# FIG.3C - NUMBER OF KNOWN ISOFORMS PER CELL TYPE
# ——————
plot_fig3C <- function(region){

    isoform_counts <- fl_isoforms[[region]] %>%
        filter(isoform_class=="Known") %>%
        group_by(cell_type) %>%
        summarise(n_isoforms = n(), .groups="drop")

    ggplot(isoform_counts, aes(x=cell_type, y=n_isoforms, fill=cell_type)) +
        geom_col(color="black") +
        geom_text(aes(label=n_isoforms), vjust=-0.3, size=3.5) +
        labs(x="CELL TYPE", y="TOTAL ISOFORMS RECOVERED", title=region) + 
        theme_bw() + theme(legend.position = "none", aspect.ratio=1)
}

plot_fig3C("DLPFC")
plot_fig3C("OFC")


# ——————
# FIG.3D - NUMBER OF NOVEL ISOFORMS PER CELL TYPE
# ——————
plot_fig3D <- function(region){

    novel_counts <- fl_isoforms[[region]] %>%
        filter(isoform_class=="Novel") %>%
        group_by(cell_type) %>%
        summarise(n_isoforms = n(), .groups="drop")

    ggplot(novel_counts, aes(x=cell_type, y=n_isoforms, fill=cell_type)) +
        geom_col(color="black") +
        geom_text(aes(label=n_isoforms), vjust=-0.3, size=3.5) +
        labs(x="CELL TYPE", y="TOTAL ISOFORMS RECOVERED", title=region) + 
        theme_bw() + theme(legend.position = "none", aspect.ratio=1)
}

plot_fig3D("DLPFC")
plot_fig3D("OFC")


# ——————
# FIG.3E - DISTRIBUTION OF TRANSCRIPT LENGTH IN KNOWN AND NOVEL ISOFORMS ACROSS CELL TYPES
# ——————
plot_fig3E <- function(region){

    known_and_novel_iso <- fl_isoforms[[region]] %>% filter(isoform_class!="Other")

    ggplot(known_and_novel_iso, aes(x = isoform_class, y = length, fill = cell_type)) +
        geom_boxplot(position = position_dodge(width = 0.8)) +
        labs(x = "", y = "Transcript length (bp)", title=region) + 
        theme_bw() + theme(aspect.ratio=1)
}

plot_fig3E("DLPFC")
plot_fig3E("OFC")


# ——————
# FIG.3F - DISTRIBUTION OF EXON COUNT IN KNOWN AND NOVEL ISOFORMS ACROSS CELL TYPES
# ——————
plot_fig3E <- function(region){

    known_and_novel_iso <- fl_isoforms[[region]] %>% filter(isoform_class!="Other")

    ggplot(known_and_novel_iso, aes(x = isoform_class, y = exons, fill = cell_type)) +
        geom_boxplot(position = position_dodge(width = 0.8)) +
        labs(x = "", y = "Transcript length (bp)", title=region) + 
        theme_bw() + theme(aspect.ratio=1)
}

plot_fig3E("DLPFC")
plot_fig3E("OFC")


# ——————
# FIG.3G - DISTRIBUTION OF EXON COUNT IN NOVEL TRANSCRIPTS BELONING TO ANNOTATED GENES
# ——————
plot_fig3G <- function(region){

    annotated_genes <- fl_isoforms[[region]] %>% filter(gene_class=="Annotated")
    cell_order <- annotated_genes %>% count(cell_type, name="n") %>% arrange(n) %>% pull(cell_type)
    annotated_genes <- annotated_genes %>% mutate(cell_type=factor(cell_type, levels=rev(cell_order)))

    ggplot(annotated_genes, aes(x = exons, fill = cell_type)) +
        geom_histogram(binwidth = 1, color = "black", position = "identity") +
        labs(x = "NUMBER OF EXONS", y = "NUMBER OF NOVEL TRANSCRIPTS FROM ANNOTATED GENES", title=region) +
        theme_bw()

}

plot_fig3G("DLPFC")
plot_fig3G("OFC")


# ——————
# FIG.3H - DISTRIBUTION OF ISOFORM EXPRESSION IN KNOWN AND NOVEL TRANSCRIPTS ACROSS CELL TYPES
# ——————
plot_fig3H <- function(region){

  TPM_props <- fl_isoforms[[region]] %>%
    filter(
      cell_type %in% c("GABA","GLU","OLIG"),
      isoform_class %in% c("Known","Novel"),
      !is.na(tpm_category)
    ) %>%
    count(cell_type, isoform_class, tpm_category) %>%
    group_by(cell_type, isoform_class) %>%
    mutate(prop = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(
      tpm_category=factor(tpm_category, levels=c("undetected","low","high")),
      isoform_class=factor(isoform_class, levels=c("Known","Novel"))
    )

  ggplot(TPM_props, aes(x=cell_type, y=prop, fill=tpm_category)) +
    geom_col(color="black") +
    facet_wrap(~ isoform_class) +
    scale_y_continuous(labels=scales::percent_format(scale=1)) +
    labs(x="", y="ISOFORMS (%)", fill="TPM CATEGORY", title=region) +
    theme_bw() + theme(aspect.ratio=1)

}

plot_fig3H("DLPFC")
plot_fig3H("OFC")


# ————————
# FIG.3I-J - PROTEIN-CODING POTENTIAL OF NOVEL ISOFORMS IN DLPFC AND OFC
# ————————
plot_fig_3I_3J <- function(region_name) {

    region_df <- coding_potential[[region_name]]
    region_df <- region_df %>% filter(structural_category %in% c("NIC", "NNC")) 
    region_df$coding <- factor(region_df$coding, levels=c("noncoding","coding"))
    region_df$tool <- factor(region_df$tool, levels=c("GMST","CPAT","CPC2","PFAM"))
    
    ggplot(region_df, aes(x=tool, fill=coding)) +
        geom_bar(position="fill", color="black") +
        facet_wrap(~ cell_type) +
        scale_y_continuous(labels=scales::percent) +
        labs(x="PROTEIN-CODING PREDICTION TOOL", y="% NOVEL ISOFORMS", fill="CODING STATUS", title=region_name) +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))

}

plot_fig_3I_3J("DLPFC")
plot_fig_3I_3J("OFC")


# ——————
# FIG.3K - PROTEIN-CODING POTENTIAL ACROSS ALL NOVEL ISOFORMS IN THE DLPFC AND OFC
# ——————
combined_coding_potential <- bind_rows(coding_potential, .id="region") %>%
  filter(structural_category %in% c("NIC","NNC")) %>%
  count(region, coding) %>%
  group_by(region) %>%
  mutate(percent = n / sum(n))

combined_coding_potential$coding <- factor(combined_coding_potential$coding, levels=c("noncoding","coding"))
  
plot_fig_3K <- ggplot(combined_coding_potential, aes(x=region, y=percent, fill=coding)) +
  geom_col(color="black") +
  scale_y_continuous(labels = scales::percent) +
  labs(x="Region", y="% NOVEL ISOFORMS", fill="CODING STATUS") + 
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))

plot_fig_3K
