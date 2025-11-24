# R SCRIPT FOR GENERATING PANELS FOR FIGURE 3 AND THEIR SUPPLEMENTAL FIGURES

library(tidyverse) # v.2.0.0

# ——— MAIN FIGURES ———

# Location of input files: repo/isoform_annotation/data

# ——————
# FIG.5B - AS EVENTS ACROSS CELL TYPES
# ——————
SUPPA_events <- readRDS("SUPPA_events.rds")

DLPFC_events <- SUPPA_events %>% filter(region=="DLPFC") # REPLACE "DLPFC" WITH "OFC" TO ANALYZE OFC SUPPA EVENTS

# COUNT ISOFORMS PER AS EVENT PER CELL TYPE
as_counts <- DLPFC_events %>% count(AS_event, cell_type)

as_counts$AS_event <- factor(as_counts$AS_event, levels=c("RI","SE","A3","A5","AF","MX","AL"))

plot_fig5B <- ggplot(as_counts, aes(x = AS_event, y = n, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  labs(x = "AS Event", y = "Isoform Count", fill = "Cell Type") +
  theme_bw()

# ——————
# FIG.5D - TPM CATEGORIES ACROSS CELL TYPES
# ——————

# LOAD IN TPM COUNTS
DLPFC_counts <- readRDS("DLPFC_CT_TPM.rds")

counts_list <- list(
  GABA = DLPFC_counts$DLPFC_GABA_TPM,
  GLU  = DLPFC_counts$DLPFC_GLU_TPM,
  OLIG = DLPFC_counts$DLPFC_OLIG_TPM
)

# ADD TPM BY CELL TYPE
DLPFC_events <- DLPFC_events %>%
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

# SUMMARIZE ISOFORM COUNT IN EACH AS EVENT x TPM CATEGORY
tpm_props <- DLPFC_events %>%
  count(AS_event, TPM_category) %>%
  group_by(AS_event) %>%
  mutate(prop = n / sum(n)) %>%     
  ungroup()


tpm_props$AS_event <- factor(tpm_props$AS_event, levels=c("RI","SE","A3","A5","AF","MX","AL"))
tpm_props$TPM_category <- factor(tpm_props$TPM_category, levels=c("Undetected","Low","High"))

# PLOT PROPORTION OF TPM CATEGORY ACROSS AS EVENT
plot_fig5D <- ggplot(tpm_props, aes(x = AS_event, y = prop, fill = TPM_category)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "AS Event", y = "Isoform proportion (%)", fill = "TPM Category") +
  theme_bw()

# ——————
# FIG.5E - STRUCTURAL CATEGORIES
# ——————

# RECODE STRUCTURAL CATEGORIES
DLPFC_events$structural_category <- dplyr::case_when(
  DLPFC_events$structural_category == "full-splice_match" ~ "FSM",
  DLPFC_events$structural_category == "incomplete-splice_match" ~ "ISM",
  DLPFC_events$structural_category == "novel_in_catalog" ~ "NIC",
  DLPFC_events$structural_category == "novel_not_in_catalog" ~ "NNC",
  TRUE ~ "Other"
)

# SUMMARIZE STRUCTURAL CATEGORIES ACROSS AS EVENTS
structural_props <- DLPFC_events %>%
  count(AS_event, structural_category) %>%
  group_by(AS_event) %>%
  mutate(prop = n / sum(n)) %>%   
  ungroup()

structural_props$AS_event <- factor(structural_props$AS_event, levels=c("RI","SE","A3","A5","AF","MX","AL"))
structural_props$structural_category <- factor(structural_props$structural_category, levels=c("Other","NNC","NIC","ISM","FSM"))

# PLOT PROPORTION OF STRUCTURAL CATEGORIES ACROSS AS EVENTS
plot_fig5E <- ggplot(structural_props, aes(x = AS_event, y = prop, fill = structural_category)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "AS Event", y = "Isoform proportion (%)", fill = "Structural Category") +
  theme_bw()

# ——————
# FIG.5F - NMD
# ——————

# SQANTI3 NMD: TRUE IF PREDICTED ORF & CDS ENDS ARE AT LEAST 50BP BEFORE LAST JUNCTION (NA = NONCODING)

# DETERMINE WHICH NMD VALUES ARE NONCODING
DLPFC_events <- DLPFC_events %>%
  mutate(
    NMD = as.character(NMD), 
    NMD = if_else(coding == "non_coding", "Noncoding", NMD)
  )

# SUMMARIZE PROPORTION OF NMD BY AS EVENT
nmd_props <- DLPFC_events %>%
  count(AS_event, NMD) %>%
  group_by(AS_event) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

nmd_props$AS_event <- factor(nmd_props$AS_event, levels=c("RI","SE","A3","A5","AF","MX","AL"))
nmd_props$NMD <- factor(nmd_props$NMD, levels=c("Noncoding","FALSE","TRUE"))

# PLOT PROPORTION OF NDM BY AS EVENT
ggplot(nmd_props, aes(x = AS_event, y = prop, fill = NMD)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "AS Event", y = "Isoform proportion (%)", fill = "NMD Status") +
  theme_bw()


# PANELS G AND H WERE CREATED VIA THE UCSC GENOME BROWSER
# IR LOCATION ACROSS ALL TRANSCRIPT MODELS CAN BE FOUND IN THE SUPPLEMENTAL TABLES