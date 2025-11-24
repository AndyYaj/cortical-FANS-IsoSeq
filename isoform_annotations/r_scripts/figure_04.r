# R SCRIPT FOR GENERATING PANELS FOR FIGURE 3 AND THEIR SUPPLEMENTAL FIGURES

library(modelsummary) # v.2.5.0
library(tidyverse) # v.2.0.0

# ——— MAIN FIGURES ———

# Location of input files: repo/isoform_annotation/data

# ——————
# FIG.4B - FRACTION OF STRUCTURAL CATEGORIES ACROSS CELL TYPES
# ——————


# READ AND COMBINE INPUT TABLES
DLPFC <- readRDS("DLPFC_CT_summary.rds")

# FOR OFC, REPLACE DLPFC WITH OFC:
    # OFC <- readRDS("OFC_CT_summary.rds")

DLPFC$structural_category <- dplyr::case_when(
  DLPFC$structural_category == "full-splice_match" ~ "FSM",
  DLPFC$structural_category == "incomplete-splice_match" ~ "ISM",
  DLPFC$structural_category == "novel_in_catalog" ~ "NIC",
  DLPFC$structural_category == "novel_not_in_catalog" ~ "NNC",
  TRUE ~ "Other"
)

structural_props <- DLPFC %>%
  group_by(cell_type, structural_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(prop = n / sum(n) * 100)


structural_props$cell_type <- factor(structural_props$cell_type, levels=c("OLIG","GLU","GABA"))
structural_props$structural_category <- factor(structural_props$structural_category, levels=c("Other","NNC","NIC","ISM","FSM"))

plot_fig4B <- ggplot(structural_props, aes(x = prop, y = cell_type, fill = structural_category)) +
  geom_bar(stat = "identity") +
  theme_classic()


# ——————
# FIG.4C - NUMBER OF KNOWN ISOFORMS PER CELL TYPE
# ——————

DLPFC <- DLPFC %>%
  mutate(isoform_class = case_when(
    structural_category %in% c("FSM", "ISM") ~ "Known",
    structural_category %in% c("NIC", "NNC") ~ "Novel",
    TRUE ~ "Other"
  ))

# SUMMARIZE KNOWN ISOFORM COUNTS ACROSS CELL TYPES
known_iso_counts <- DLPFC %>%
  filter(isoform_class == "Known") %>%      
  group_by(cell_type) %>%
  summarise(n_isoforms = n())

# PLOT NUMBER OF KNOWN ISOFORMS ACROSS CELL TYPES
plot_fig4C <- ggplot(known_iso_counts, aes(x = cell_type, y = n_isoforms, fill = cell_type)) +
  geom_col() +
  theme_classic() + 
  labs(x = "Cell Type", y = "Number of recovered isoforms") +
  theme(legend.position = "none")


# ——————
# FIG.4D - NUMBER OF NOVEL ISOFORMS PER CELL TYPE
# ——————

# SUMMARIZE NOVEL ISOFORM COUNTS ACROSS CELL TYPES
novel_iso_counts <- DLPFC %>%
  filter(isoform_class == "Novel") %>%      
  group_by(cell_type) %>%
  summarise(n_isoforms = n())

# PLOT NUMBER OF NOVEL ISOFORMS ACROSS CELL TYPES
plot_fig4D <- ggplot(novel_iso_counts, aes(x = cell_type, y = n_isoforms, fill = cell_type)) +
  geom_col() +
  theme_classic() + 
  labs(x = "Cell Type", y = "Number of recovered isoforms") +
  theme(legend.position = "none")


# ——————
# FIG.4E - DISTRIBUTION OF TRANSCRIPT LENGTH IN KNOWN AND NOVEL READS ACROSS CELL TYPES
# ——————

DLPFC <- DLPFC %>% mutate(length = as.numeric(length))
known_and_novel_iso <- DLPFC %>% filter(class != "Other")
summarize_by_group(known_and_novel_iso, length)

plot_fig4E <- ggplot(known_and_novel_iso, aes(x = class, y = length, fill = cell_type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "", y = "Transcript length (bp)") + 
  theme_bw()


# ——————
# FIG.4F - DISTRIBUTION OF EXON COUNT IN KNOWN AND NOVEL READS ACROSS CELL TYPES
# ——————

DLPFC <- DLPFC %>% mutate(exons = as.numeric(exons))
known_and_novel_iso <- DLPFC %>% filter(class != "Other")
summarize_by_group(known_and_novel_iso, length)

plot_fig4E <- ggplot(known_and_novel_iso, aes(x = class, y = exons, fill = cell_type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "", y = "Exon count") + 
  theme_bw()


# ——————
# FIG.4G - DISTRIBUTION OF EXON COUNT IN NOVEL TRANSCRIPTS BELONING TO ANNOTATED GENES
# ——————

DLPFC <- DLPFC %>%
  mutate(gene_class = ifelse(str_starts(ensembl_id, "novelGene_"), "Novel", "Annotated"))

annotated_genes <- DLPFC %>% filter(gene_class == "Annotated")
annotated_genes$cell_type <- factor(annotated_genes$cell_type, levels=c("OLIG","GABA","GLU"))

plot_fig4G <- annotated_genes %>%
  ggplot(aes(x = exons, fill = cell_type)) +
  geom_histogram(binwidth = 1, color = "black", position = "identity") +
  labs(x = "Exon count", y = "Number of novel transcripts from annotated genes") +
  theme_bw()

plot_fig4G


# ——————
# FIG.4H - DISTRIBUTION OF ISOFORM EXPRESSION IN KNOWN AND NOVEL TRANSCRIPTS ACROSS CELL TYPES
# ——————

# LOAD IN TPM COUNTS
DLPFC_counts <- readRDS("DLPFC_CT_TPM.rds")

counts_list <- list(
  GABA = DLPFC_counts$DLPFC_GABA_TPM,
  GLU  = DLPFC_counts$DLPFC_GLU_TPM,
  OLIG = DLPFC_counts$DLPFC_OLIG_TPM
)

# ADD TPM BY CELL TYPE
DLPFC <- DLPFC %>%
  mutate(
    avg_TPM = map2_dbl(cell_type, isoform_id, ~ {
      df <- counts_list[[.x]]
      df$avg_TPM[df$isoform_id == .y]
    }),
    # TRANSFORM TPM VALUES 
    logTPM = log2(avg_TPM + 1),
    # CREATE TPM CATEGORIES
    TPM_category = case_when(
      logTPM == 0 ~ "Undetected",
      logTPM > 0 & logTPM < 2.5 ~ "Low",
      logTPM >= 2.5 ~ "High",
      TRUE ~ NA_character_
    )
  )

# SUMMARIZE TPM CATEGORY PROPORTIONS 
TPM_props <- DLPFC %>%
  group_by(cell_type, TPM_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(prop = n / sum(n) * 100)

summarize_by_group(DLPFC, logTPM) # SUMMARY STATS FOR TRANSFORMED TPM VALUES ACROSS CELL TYPES

TPM_props$TPM_category <- factor(TPM_props$TPM_category, levels=c("Undetected","Low","High"))

plot_fig4H <- ggplot(TPM_props, aes(x = cell_type, y = prop, fill = TPM_category)) +
  geom_col(color = "black") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "", y = "Proportion of Isoforms (%)", fill = "TPM Category") +
  theme_bw()

plot_fig4H


# ——————
# FIG.4I - PROTEIN CODING POTENTIAL ACROSS DIFFERENT TOOLS
# ——————

DLPFC_coding <- readRDS("DLPFC_coding_potential.rds")

# SUMMARIZE PROTEIN CODING POTENTIAL PROPORTIONS (OUTPUT FROM SQANTI3 - GMST PREDICTIONS)
coding_prop <- DLPFC_coding %>%
  group_by(tool, coding) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tool) %>%
  mutate(prop = n / sum(n))

coding_prop$coding <- factor(coding_prop$coding, levels=c("noncoding","coding"))
coding_prop$tool <- factor(coding_prop$tool, levels=c("GMST","CPAT","CPC2","PFAM"))

plot_fig4I <- ggplot(coding_prop, aes(x = tool, y = prop, fill = coding)) +
  geom_col(color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Protein-coding prediction tool", y = "Isoform proportion (%)", fill = "Coding Status") +
  theme_bw()


# ——————
# FIG.4J - PROTEIN CODING POTENTIAL OF NOVEL TRANSCRIPTS ACROSS DIFFERENT TOOLS
# ——————

DLPFC_novel_coding <- DLPFC_coding %>% filter(structural_category %in% c("novel_in_catalog","novel_not_in_catalog"))

novel_coding_prop <- DLPFC_novel_coding %>%
  group_by(tool, coding) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tool) %>%
  mutate(prop = n / sum(n))

novel_coding_prop$coding <- factor(novel_coding_prop$coding, levels=c("noncoding","coding"))
novel_coding_prop$tool <- factor(novel_coding_prop$tool, levels=c("GMST","CPAT","CPC2","PFAM"))

plot_fig4J <- ggplot(novel_coding_prop, aes(x = tool, y = prop, fill = coding)) +
  geom_col(color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Protein-coding prediction tool", y = "Isoform proportion (%)", fill = "Coding Status") +
  theme_bw()


# ——————
# FIG.4K - PROPORTION OF ISOFORMS WITH FULL LENGTH SUPPORT (5' and 3' END SUPPORT)
# ——————

# CONVERT COLUMNS INTO LOGICAL VALUES
DLPFC <- DLPFC %>%
  mutate(
    within_CAGE = within_CAGE == "TRUE",
    within_polyA = within_polyA == "TRUE",
    # DEFINE FL SUPPORT: WHEN ISOFORM HAS BOTH CAGE AND POLYA SUPPORT
    FL_support = within_CAGE & within_polyA
  )

fl_props <- DLPFC %>%
  count(FL_support) %>%
  mutate(prop = n / sum(n))

# PLOT PROPORTION OF ISOFORMS THAT HAVE BOTH 5' AND 3' END SUPPORT
plot_fig4K <- ggplot(fl_props, aes(x = "FL_support", y = prop, fill = FL_support)) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "FL support", y = "Proportion", fill = "FL support") +
  theme_bw()

# ——————
# FIG.4L - PROPORTION OF ISOFORMS WITH CAGE AND POLYA SUPPORT ACROSS STRUCTURAL CATEGORIES
# ——————

# SUMMARIZE END SUPPORT PROPORTIONS
end_support_props <- DLPFC %>%
  pivot_longer(
    cols = c(within_CAGE, within_polyA),
    names_to = "support_type",
    values_to = "supported"
  ) %>%
  group_by(structural_category, support_type, supported) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(structural_category, support_type) %>%
  mutate(prop = n / sum(n))

# PLOT PROPORTION OF ISOFORMS WITH CAGE AND POLYA END SUPPORT ACROSS STRUCTURAL CATEGORIES
plot_fig4L <- ggplot(end_support_props %>% filter(supported == TRUE), aes(x = prop, y = structural_category)) +
  geom_col(fill = "steelblue", color = "black") +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = "Structural Category", y = "Proportion Supported") +
  facet_wrap(~support_type, ncol = 2) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))