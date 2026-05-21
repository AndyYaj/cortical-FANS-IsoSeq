# ———————————
# FIGURE_04.R — R SCRIPT FOR GENERATING FIGURE 4 PANELS
# ———————————


library(tidyverse) # v.2.0.0

end_support <- list(
  DLPFC = readRDS("DLPFC_Final_End_Support.rds") %>%
    mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG"))),
  OFC = readRDS("OFC_Final_End_Support.rds") %>%
    mutate(cell_type=factor(cell_type, levels=c("GABA","GLU","OLIG","AST","MG")))
)

end_support_heatmap <- readRDS("End_Support_Heatmap.rds")
long_read_comparison <- readRDS("Matched_Isoforms_Summary.rds")
region_struct_cell <- readRDS("levels.rds")

# ——————
# FIG.4B 
# ——————
# PLOT 3'-END SUPPORT (polyA site AND polyA motif)
plot_fig_4B_3end <- function(region){

    region_df <- end_support[[region]]

    three_end_support <- region_df %>% 
        select(structural_category, within_polyA_site, polyA_motif_found) %>%
        pivot_longer(
            cols = c(within_polyA_site, polyA_motif_found),
            names_to = "feature",
            values_to = "value"
        ) %>%
        group_by(structural_category, feature) %>%
        summarise(
            pct_true = mean(value==TRUE) * 100,
            .groups = "drop"
        )

    ggplot(three_end_support, aes(x=structural_category, y=pct_true, fill=feature)) +
        geom_col(position=position_dodge(width=0.9), width=0.9, color="black") +
        labs(x="STRUCTURAL CATEGORY", y="% 3' END SUPPORT", fill="END SUPPORT", title=region) + 
        theme_bw() + theme(aspect.ratio=1)

}

plot_fig_4B_3end("DLPFC")
plot_fig_4B_3end("OFC")


# PLOT 5'-END SUPPORT (CAGE AND ATAC)
plot_fig_4B_5end <- function(region){

    region_df <- end_support[[region]]

    # CAGE PEAK
    cage_df <- region_df %>%
        group_by(structural_category) %>%
        summarise(
            value=mean(within_CAGE_peak==TRUE) * 100,
            .groups="drop"
        ) %>%
        mutate(type="CAGE")
    
    # ATAC (STACKED BY NUMBER OF STUDIES)
    atac_df <- region_df %>%
        mutate(atac_bin=factor(num_ATAC_studies)) %>%
        group_by(structural_category, atac_bin) %>%
        summarise(n=n(), .groups="drop") %>%
        group_by(structural_category) %>%
        mutate(value = n / sum(n) * 100) %>%
        ungroup() %>%
        mutate(type="ATAC")

    # COMBINE CAGE AND ATAC STRUCTURE POSITIONS
    all_cats <- unique(region_df$structural_category)
    x_map <- setNames(seq_along(all_cats), all_cats)
    cage_df$x <- x_map[cage_df$structural_category] - 0.2
    atac_df$x <- x_map[atac_df$structural_category] + 0.2

    ggplot() +
        geom_col(data=cage_df, aes(x=x, y=value), width=0.4, fill="gray", color="black") +
        geom_col(data=atac_df, aes(x=x, y=value, fill=atac_bin), width=0.4, position="stack", color="black") +
        scale_x_continuous(breaks=seq_along(all_cats), labels=all_cats) +
        labs(x="STRUCTURAL CATEGORY", y="5' END SUPPORT (%)", fill="ATAC studies", title=region) +
        theme_bw() + theme(aspect.ratio=1)
}

plot_fig_4B_5end("DLPFC")
plot_fig_4B_5end("OFC")


# ————————
# FIG.4C-D
# ————————
plot_fig_4C_4D <- function(region_input) {

  celltype_order <- list(
    DLPFC=c("GABA","GLU","OLIG"),
    OFC=c("GABA","GLU","OLIG","AST","MG")
  )

  category_order <- c("FSM","ISM","NIC","NNC","Other")
  ct_order <- celltype_order[[region_input]]

  region_heatmap <- end_support_heatmap %>%
    filter(grepl(paste0("^", region_input, "_"), label)) %>%
    separate(label, into=c("region","celltype","category"), sep="_") %>%
    mutate(
      celltype=factor(celltype, levels=ct_order),
      category=factor(category, levels=category_order)
    ) %>%
    arrange(category, celltype) %>%
    mutate(
      label=paste(region, celltype, category, sep="_"),
      label=factor(label, levels=rev(unique(label)))
    ) %>%
    pivot_longer(
      cols=starts_with("tier"),
      names_to="tier",
      values_to="value"
    )

  ggplot(region_heatmap, aes(x=tier, y=label, fill=value)) +
    geom_tile(color="black") +
    geom_text(aes(label=sprintf("%.1f", value)), size=4) +
    scale_fill_gradient(low="white", high="midnightblue", name="Support") +
    theme_minimal(base_size=14) +
    labs(x="Tier", y="Label", title=paste(region_input)) +
    theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())

}

plot_fig_4C_4D("DLPFC")
plot_fig_4C_4D("OFC")


# ————————
# FIG.4E-F
# ————————
plot_fig_4E_4F <- function(region){

    region_df <- end_support[[region]]

    fl_support_df <- region_df %>%
        group_by(structural_category, Tier2) %>% summarise(n=n(), .groups="drop") %>%
        group_by(structural_category) %>%mutate(prop = n / sum(n))

    ggplot(fl_support_df, aes(x=structural_category, y=prop, fill=Tier2)) +
        geom_col(position="stack", color="black") + 
        scale_y_continuous(labels=scales::percent_format()) +
        labs(x="", y="% FL SUPPORT (TIER 2)", title=region) +
        theme_bw() + theme(aspect.ratio=1)

}

plot_fig_4E_4F("DLPFC")
plot_fig_4E_4F("OFC")


# ————————
# FIG.4H-I
# ————————
plot_fig_4H_4I <- function(region){

    region_df <- long_read_comparison %>% 
        filter(grepl(paste0("^", region), as.character(region_struct_cell)))

    ggplot(region_df, aes(x=region_struct_cell, y=proportion, fill=shared_count)) +
        geom_col(width=0.9, color="black", linewidth=0.4) +
        scale_fill_brewer(palette="Blues") +
        scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), aspect.ratio=1)           
}

plot_fig_4H_4I("DLPFC") 
plot_fig_4H_4I("OFC") 


# ————————
# FIG.4J-K
# ————————
plot_fig_4J_4K <- function(region) {
  
  files <- c(
    paste0(region, "_FSM_phastCons.gz"),
    paste0(region, "_ISM_phastCons.gz"),
    paste0(region, "_NIC_phastCons.gz"),
    paste0(region, "_NNC_phastCons.gz"),
    paste0(region, "_OTHER_phastCons.gz"),
    paste0(region, "_Shuffled_phastCons")
  )
  
  labels <- c("FSM","ISM","NIC","NNC","Other","Shuff")
  
  process_phastCons_matrix <- function(file, label) {
    mat <- fread(file, header = FALSE)
    
    long <- mat %>%
      rename(region=V1) %>%
      pivot_longer(
        cols = -region,
        names_to = "bin",
        values_to = "score",
        values_transform = list(score = as.numeric)
      ) %>%
      mutate(bin = as.numeric(gsub("V", "", bin)),
             type = label)
    
    long %>%
      group_by(bin, type) %>%
      summarise(
        mean_score = mean(score, na.rm=TRUE),
        sd_score = sd(score,  na.r =TRUE),
        .groups="drop"
      )
  }
  
  all_profiles <- map2_dfr(files, labels, process_phastCons_matrix) %>% mutate(smoothed = mean_score)
  
  LEFT_BIN <- 5
  TSS_BIN <- 125
  TTS_BIN <- 328
  RIGHT_BIN <- 448
  
  all_profiles$type <- factor(all_profiles$type, levels=c("Shuff","Other","ISM","NNC","NIC","FSM"))
  
  plot <- ggplot(all_profiles, aes(x=bin, y=smoothed, color=type, group=type)) + 
    geom_line(linewidth=2.5) +
    scale_y_continuous(limits=c(0, 0.9)) +
    scale_x_continuous(limits=c(LEFT_BIN, RIGHT_BIN), breaks=c(LEFT_BIN, TSS_BIN, TTS_BIN, RIGHT_BIN), labels=c("-3 kb", "TSS", "TTS", "+3 kb")) +
    geom_hline(yintercept=0, color="grey30", linewidth=0.5) +
    scale_color_manual(values = c(
      "FSM"="#374e58ff",
      "ISM"="#568282ff",
      "NNC"="#7e47efff",
      "NIC"="#e200ecff",
      "Other"="#bdd7e7ff",
      "Shuff"="#bec3d1ff")) +
      labs(x="", y="Conservation Score (PhastCons)", title=region) +
    theme_bw() + theme(aspect.ratio=1)

  return(plot)
}

# NOTE: PLOTS TAKE ~3 MINUTES TO GENERATE
plot_fig_4J_4K("DLPFC")
plot_fig_4J_4K("OFC")
