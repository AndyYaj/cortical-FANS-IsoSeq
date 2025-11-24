# R SCRIPT FOR GENERATING PANELS FOR FIGURE 3

library(ggVennDiagram) # v.1.5.4
library(modelsummary) # v.2.5.0
library(tidyverse) # v.2.0.0
library(biomaRt) # v.2.62.1


# ——— HELPER FUNCTION ———

# GENERATE SUMMARY STATISTICS ACROSS CELL TYPES AND REGIONS
summarize_by_group <- function(df, numeric_var, group_var = "cell_type") {
  
  numeric_var <- rlang::ensym(numeric_var)
  group_var <- rlang::ensym(group_var)
  
  group_summary <- df %>%
    group_by(!!group_var) %>%
    summarise(
      Mean = mean(!!numeric_var, na.rm = TRUE),
      SD = sd(!!numeric_var, na.rm = TRUE),
      Min = min(!!numeric_var, na.rm = TRUE),
      Max = max(!!numeric_var, na.rm = TRUE),
      Median = median(!!numeric_var, na.rm = TRUE),
      .groups = "drop"
    )

  combined_summary <- df %>%
    summarise(
      Mean = mean(!!numeric_var, na.rm = TRUE),
      SD = sd(!!numeric_var, na.rm = TRUE),
      Min = min(!!numeric_var, na.rm = TRUE),
      Max = max(!!numeric_var, na.rm = TRUE),
      Median = median(!!numeric_var, na.rm = TRUE)
    ) %>%
    mutate(!!group_var := "Combined") %>%
    select(!!group_var, everything())

  final_summary <- bind_rows(group_summary, combined_summary)
  datasummary_df(final_summary)
}


# ——— MAIN FIGURES ———

# Location of input files: repo/isoform_annotation/data

# ——————
# FIG.3A - NUMBER OF ISOFORMS PER CELL TYPE
# ——————

# READ AND COMBINE INPUT TABLES
DLPFC <- readRDS("DLPFC_cell_types.rds")

# GET HGNC GENE SYMBOLS
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
DLPFC$cleaned_ensembl_id <- sub("\\..*", "", DLPFC$ensembl_id)
unique_ids <- unique(DLPFC$cleaned_ensembl_id)

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique_ids,
  mart = ensembl
)

DLPFC <- DLPFC %>%
  left_join(gene_mapping, by = c("cleaned_ensembl_id" = "ensembl_gene_id")) %>%
  mutate(
    gene_symbol = coalesce(na_if(hgnc_symbol, ""), cleaned_ensembl_id)
  ) %>%
  dplyr::select(-cleaned_ensembl_id, -hgnc_symbol)

# SUMMARIZE COUNTS PER CELL TYPE
isoform_counts <- DLPFC %>%
  group_by(cell_type) %>%
  summarise(n_isoforms = n())

# PLOT NUMBER OF ISOFORMS PER CELL TYPE
plot_fig3A <- ggplot(isoform_counts, aes(x = cell_type, y = n_isoforms, fill = cell_type)) +
  geom_col() +
  theme_classic() + 
  labs(x = "Cell Type", y = "Number of recovered isoforms") +
  theme(legend.position = "none")


# ——————
# FIG.3A - NUMBER OF UNIQUE GENES PER CELL TYPE
# ——————

# SUMMARIZE COUNTS PER CELL TYPE
gene_counts <- DLPFC %>%
  group_by(cell_type) %>%
  summarise(n_genes = n_distinct(gene_id))

# PLOT NUMBER OF ISOFORMS PER CELL TYPE
plot_fig3B <- ggplot(gene_counts, aes(x = cell_type, y = n_genes, fill = cell_type)) +
  geom_col() +
  theme_classic() + 
  labs(x = "Cell Type", y = "Number of unique genes") +
  theme(legend.position = "none")


# ——————
# FIG.3C - ISOFORM OVERLAP
# ——————

# CREATE TRANSCRIPT SETS FOR EACH CELL TYPE
transcript_sets <- DLPFC %>%
  group_by(cell_type) %>%
  summarise(transcripts = list(unique(MSTRG_transcript[MSTRG_transcript != ""]))) %>%
  deframe() 

# PLOT OVERLAP OF ISOFORMS ACROSS CELL TYPES
plot_fig3C <- ggVennDiagram(transcript_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "blue") 


# ——————
# FIG.3D - GENE OVERLAP
# ——————

# CREATE GENE SETS FOR EACH CELL TYPE
gene_sets <- DLPFC %>%
  group_by(cell_type) %>%
  summarise(transcripts = list(unique(MSTRG_gene[MSTRG_gene != ""]))) %>%
  deframe() 

# PLOT OVERLAP OF GENES ACROSS CELL TYPES
plot_fig3D <- ggVennDiagram(gene_sets, label_alpha = 0) + 
    scale_fill_gradient(low = "white", high = "blue")


# ——————
# FIG.3E - DISTRIBUTION OF TRANSCRIPT LENGTH
# ——————

DLPFC <- DLPFC %>% mutate(length = as.numeric(length))
DLPFC$cell_type <- factor(DLPFC$cell_type, levels=c("OLIG","GABA","GLU"))

# PLOT TRANSCRIPT LENGTH DISTRIBUTION BY CELL TYPE
plot_fig3E <- ggplot(DLPFC, aes(x = length, fill = cell_type)) +
  geom_histogram(binwidth = 150, color="black", position = "identity") +
  labs(x = "Transcript length (bp)", y = "Number of transcripts") +
  theme_classic()

# TRANSCRIPT LENGTH SUMMARY STATISTICS
summarize_by_group(DLPFC, length)


# ——————
# FIG.3F - DISTRIBUTION OF EXON COUNT
# ——————

DLPFC <- DLPFC %>% mutate(exons = as.numeric(exons))

# PLOT EXON COUNT DISTRIBUTION BY CELL TYPE
plot_fig3F <- ggplot(DLPFC, aes(x = exons, fill = cell_type)) +
  geom_histogram(binwidth = 1, color="black", position = "identity") +
  labs(x = "Exon count", y = "Number of transcripts") +
  theme_classic()

# EXON COUNT SUMMARY STATISTICS
summarize_by_group(DLPFC, exons)


# ——————
# FIG.3G - NUMBER OF ISOFORMS PER GENE
# ——————

# DETERMINE NUMBER OF ISOFORMS PER GENE
iso_per_gene <- DLPFC %>%
  group_by(cell_type, gene_id) %>%            
  summarise(n_isoforms = n_distinct(isoform_id), .groups = "drop") %>%  
  # CREATE ISO COUNT CATEGORIES
  mutate(
    iso_category = case_when(
      n_isoforms == 1 ~ "one_isoform",
      n_isoforms %in% 2:3 ~ "two_to_three",
      n_isoforms %in% 4:5 ~ "four_to_five",
      n_isoforms %in% 6:7 ~ "six_to_seven",
      n_isoforms %in% 8:9 ~ "eight_to_nine",
      TRUE ~ "ten_or_more"
    )
  )

# COMPUTE % GENES PER ISOFORM CATEGORY x CELL TYPE
iso_category_summary <- iso_per_gene %>%
  count(iso_category, cell_type) %>%
  group_by(cell_type) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

iso_category_summary$iso_category <- factor(iso_category_summary$iso_category, levels=c("one_isoform","two_to_three","four_to_five","six_to_seven","eight_to_nine","ten_or_more"))
iso_category_summary$cell_type <- factor(iso_category_summary$cell_type, levels=c("GABA","GLU","OLIG"))

# PLOT ISOFORM PER GENE ACROSS CELL TYPES
plot_fig3G <- ggplot(iso_category_summary, aes(x = iso_category, y = percent, fill = cell_type, color=cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width=0.8, alpha=0.75, linewidth=1.5) +
  theme_classic()

# ISO PER GENE SUMMARY STATISTICS
summarize_by_group(isoforms_per_gene, n_isoforms)
