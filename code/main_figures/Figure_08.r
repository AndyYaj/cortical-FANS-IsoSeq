# ———————————
# FIGURE_08.R — R SCRIPT FOR GENERATING FIGURE 8 PANELS
# ———————————


library(tidyverse) # v.2.0.0
library(ggrepel) # v.0.9.8


fig_8A_df <- read.delim("Figure_8A.tsv", header=T, sep="\t")
fig_8C_df <- read.delim("Figure_8C.tsv", header=T, sep="\t")
fig_8D_df <- read.delim("Figure_8D.tsv", header=T, sep="\t")
fig_8E_df <- read.delim("Figure_8E.tsv", header=T, sep="\t")
fig_8F_8G_df <- read.delim("Figure_8F_8G.tsv", header=T, sep="\t")
fig_8F_8G_inset_df <- read.delim("Figure_8F_8G_inset.tsv", header=T, sep="\t")


# ——————
# FIG.8A
# ——————
fig_8A_df$structural_category <- factor(fig_8A_df$structural_category, levels=c("antisense","intergenic","genic","fusion","novel_not_in_catalog","novel_in_catalog","incomplete-splice_match","full-splice_match"))
fig_8A_df$region <- factor(fig_8A_df$region, levels=c("OFC","DLPFC"))

plot_fig_8A <- ggplot(fig_8A_df, aes(x=unique_variants, y=structural_category, fill=region)) +
  geom_col(position=position_dodge(width=0.9)) +
  theme_bw() + theme(aspect.ratio=1)


# ——————
# FIG.8C
# ——————
fig_8C_df$region <- factor(fig_8C_df$region, levels=c("OFC","DLPFC"))
fig_8C_df$class <- factor(fig_8C_df$class, levels=c("novel_exon_intergenic","novel_exon_intronic","internal_exon_extension","novel_donor","novel_acceptor"))

plot_fig_8C <- ggplot(fig_8C_df, aes(x=count, y=class, fill=region)) +
  geom_col(position=position_dodge(width=0.9)) +
  theme_bw() + theme(aspect.ratio=1)


# ——————
# FIG.8D
# ——————
fig_8D_df$region <- factor(fig_8D_df$region, levels=c("OFC","DLPFC"))
fig_8D_df$class <- factor(fig_8D_df$class, levels=c("novel_exon_intergenic","novel_exon_intronic","internal_exon_extension","novel_donor","novel_acceptor"))

plot_fig_8C <- ggplot(fig_8D_df, aes(x=unique_variants, y=class, fill=region)) +
  geom_col(position=position_dodge(width=0.9)) +
  theme_bw() + theme(aspect.ratio=1)


# ——————
# FIG.8E
# ——————
fig_8E_df$region <- factor(fig_8E_df$region, levels=c("DLPFC","OFC"))
fig_8E_df$class <- factor(fig_8E_df$class, levels=c("novel_exon_intergenic","novel_exon_intronic","internal_exon_extension","novel_donor","novel_acceptor"))

plot_fig_8C <- ggplot(fig_8E_df, aes(x=odds_ratio, y=class, color=region)) +
  geom_point(position=position_dodge(width=0.4), size=3) +
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper), position=position_dodge(width=0.4), width=0.2) +
  theme_bw() + theme(aspect.ratio = 1/0.5)


# ———————
# FIG. 8F
# ———————
DLPFC_df <- fig_8F_8G_df %>% filter(region=="DLPFC")
top_genes <- DLPFC_df %>% arrange(desc(score)) %>% slice(1:25) 

plot_fig_8F <- ggplot(DLPFC_df, aes(x=novel_splice_events, y=novel_event_overlap_variants)) +
  geom_point(size=7, shape=21, color="black", stroke=1) +
  geom_point(size=7, color="#e2e5eaff") +   
  geom_point(data=top_genes, color="#445dafff", size=8) +
  theme_bw() + theme(aspect.ratio=1/2)

# INSET PLOT 
DLPFC_inset_df <- fig_8F_8G_inset_df %>% filter(region=="DLPFC")

plot_fig_8F_inset <- ggplot(DLPFC_inset_df, aes(x=log10(splice_density), y=log10(overlap_variant_density))) +
  geom_point(size=10, shape=21, color="black", stroke=1.25) +
  geom_point(size=10) + 
  geom_smooth(method="lm", se=FALSE, size=5) +
  theme_bw() + theme(aspect.ratio = 1/0.75)


# ———————
# FIG. 8G
# ———————
OFC_df <- fig_8F_8G_df %>% filter(region=="OFC")
top_genes <- OFC_df %>% arrange(desc(score)) %>% slice(1:25) 

plot_fig_8F <- ggplot(OFC_df, aes(x = novel_splice_events, y=novel_event_overlap_variants)) +
  geom_point(size=7, shape=21, color="black", stroke=1) +
  geom_point(size=7, color="#e2e5eaff") +   
  geom_point(data=top_genes, color="#445dafff", size=8) +
  theme_bw() + theme(aspect.ratio=1/2)


# INSET PLOT 
OFC_inset_df <- fig_8F_8G_inset_df %>% filter(region=="OFC")

plot_fig_8G_inset <- ggplot(OFC_inset_df, aes(x=log10(splice_density), y=log10(overlap_variant_density))) +
  geom_point(size=10, shape=21, color="black", stroke=1.25) +
  geom_point(size=10) + 
  geom_smooth(method="lm", se = FALSE, size=5) +
  theme_bw() + theme(aspect.ratio=1/0.75)
