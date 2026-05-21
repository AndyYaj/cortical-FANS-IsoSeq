# ———————————
# FIGURE_02.R — R SCRIPT FOR GENERATING FIGURE 2 PANELS
# ———————————


# PATH TO INPUT FILES: repo/isoform_annotation/data

library(tidyverse) # v.2.0.0
library(UpSetR) # v.1.4.0


# LOAD IN FULL-LENGTH (FL) ISOFORMS RECOVERED ACROSS CORTICAL CELL TYPES IN THE DLPFC AND OFC
fl_isoforms <- list(

  DLPFC = readRDS("Figure_02_DLPFC_FL_Transcripts.rds") %>%
    mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG"))),
  
  OFC = readRDS("Figure_02_OFC_FL_Transcripts.rds") %>%
    mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG","AST","MG")))

)

# ——————
# FIG.2A - TOTAL FL TRANSCRIPTS RECOVERED ACROSS CELL TYPES IN THE DLPFC AND OFC
# ——————
plot_fig_2A <- function(region) {    

    isoform_counts <- fl_isoforms[[region]] %>%
        group_by(cell_type) %>%
        summarize(n_isoforms=n(), .groups="drop")
    
    ggplot(isoform_counts, ggplot2::aes(x=cell_type, y=n_isoforms)) +
        geom_col(fill="gray", color="black") +
        geom_text(aes(label=n_isoforms), vjust=-0.3, size=3.5) +
        labs(x="CELL TYPE", y="TOTAL ISOFORMS RECOVERED", title=region) +
        theme_bw() + theme(legend.position = "none", aspect.ratio = 1)

}

plot_fig_2A("DLPFC")
plot_fig_2A("OFC")


# ——————
# FIG.2B - TOTAL UNIQUE GENES ACROSS CELL TYPES IN THE DLPFC AND OFC
# ——————
plot_fig_2B <- function(region) {

    gene_counts <- fl_isoforms[[region]] %>%
        group_by(cell_type) %>%
        summarize(n_genes=n_distinct(gene_id), .groups="drop")
        
    ggplot(gene_counts, aes(x=cell_type, y=n_genes)) +
        geom_col(fill="gray", color="black") +
        geom_text(aes(label=n_genes), vjust=-0.3, size=3.5) +
        labs(x="CELL TYPE", y="TOTAL UNIQUE GENES", title=region) +
        theme_bw() + theme(legend.position="none", aspect.ratio=1)

}

plot_fig_2B("DLPFC")
plot_fig_2B("OFC")


# ————————
# FIG.2C-D - OVERLAP OF FL TRANSCRIPTS AND UNIQUE GENES
# ————————
plot_fig_2C_2D <- function(region, feature, file = NULL, sets = NULL) {
    
    upset_df <- fl_isoforms[[region]] %>% select(cell_type, {{ feature }})
    
    if (is.null(sets)) {
        if (is.factor(upset_df$cell_type)) {
            sets <- intersect(levels(upset_df$cell_type), unique(upset_df$cell_type))
        } else {
            sets <- unique(upset_df$cell_type)
        }
    }

    cell_type_sets <- upset_df %>%
        distinct(cell_type, {{ feature }}) %>%
        split(.$cell_type) %>%
        lapply(function(x) x[[2]])  
    
    upset_data <- UpSetR::fromList(cell_type_sets)
    UpSetR::upset(upset_data, sets=sets, order.by="freq", decreasing=TRUE)

}

# FIG.2C UPSET - NUMBER OF TRANSCRIPTS
plot_fig_2C_2D("DLPFC", MSTRG_transcript)
plot_fig_2C_2D("OFC", MSTRG_transcript)

# FIG.2D UPSET - NUMBER OF UNIQUE GENES
plot_fig_2C_2D("DLPFC", MSTRG_gene)
plot_fig_2C_2D("OFC", MSTRG_gene)


# ————————
# FIG.2E-F — DISTRIBUTION OF TRANSCRIPT LENGTH AND EXON COUNT
# ————————
plot_fig_2E_2F <- function(region, x_var, binwidth) {
    
    region_df <- fl_isoforms[[region]]
    cell_order <- region_df %>% count(cell_type, name="n") %>% arrange(n) %>% pull(cell_type)
    region_df <- region_df %>% mutate(cell_type=factor(cell_type, levels=rev(cell_order)))
    
    ggplot(region_df, aes(x=.data[[x_var]], fill=cell_type)) +
        geom_histogram(binwidth=binwidth, color="black", position="identity") +
        labs(x=toupper(x_var), y="NUMBER OF TRANSCRIPTS", title=region) +
        theme_bw()
}

# FIG.2E — DISTRIBUTION OF TRANSCRIPT LENGTH
plot_fig_2E_2F("DLPFC", "length", 150)
plot_fig_2E_2F("OFC", "length", 150)

# FIG.2F — DISTRIBUTION OF EXON COUNT
plot_fig_2E_2F("DLPFC", "exons", 1)
plot_fig_2E_2F("OFC", "exons", 1)


# ——————
# FIG.2G - NUMBER OF ISOFORMS PER GENE
# ——————
plot_fig_2G <- function(region){

    region_df <- fl_isoforms[[region]]

    iso_per_gene <- region_df %>%
        group_by(cell_type, gene_id) %>%
        summarize(n_isoforms=n_distinct(isoform_id), .groups="drop") %>%
        mutate(iso_category=case_when(
            n_isoforms == 1 ~ "1", n_isoforms %in% 2:3 ~ "2-3",
            n_isoforms %in% 4:5 ~ "4-5", n_isoforms %in% 6:7 ~ "6-7",
            n_isoforms %in% 8:9 ~ "8-9", TRUE ~ "≥10"
        )
    )

    iso_category_summary <- iso_per_gene %>%
        count(iso_category, cell_type) %>%
        group_by(cell_type) %>%
        mutate(percent = 100 *n / sum(n)) %>%
        ungroup()

    iso_category_summary$iso_category <- factor(
        iso_category_summary$iso_category, levels=c("1","2-3","4-5","6-7","8-9","≥10")
    )  

    ggplot(iso_category_summary, aes(x=iso_category, y=percent, fill=cell_type, color=cell_type)) +
        geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.8, color="black") +
        labs(x="NUMBER OF TRANSCRIPTS", y="NUMBER OF ISOFORM PER GENE", title=region) +
        theme_bw() + theme(aspect.ratio=1)

}

plot_fig_2G("DLPFC")
plot_fig_2G("OFC")
