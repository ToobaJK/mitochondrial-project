# ================================
# Figure 3: PRISM drug response across mitochondrial metabolic states
# ================================

# -------- SETTINGS --------
base_dir <- "D:/UAEU/Figure 3"
output_dir <- file.path(base_dir, "outputR1")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
expr_file      <- file.path(base_dir, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
cell_meta_file <- file.path(base_dir, "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv")
drug_meta_file <- file.path(base_dir, "Repurposing_Public_24Q2_Treatment_Meta_Data.csv")
lfc_file       <- file.path(base_dir, "Repurposing_Public_24Q2_LFC_COLLAPSED.csv")

# -------- LIBRARIES --------
libs <- c("data.table", "dplyr", "tidyr", "ggplot2", "ggrepel", "readr", "stringr", "forcats")
to_install <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(stringr)
library(forcats)

# -------- SAVE PLOT HELPER --------
save_both <- function(plot_obj, filename, width = 10, height = 7, dpi = 300) {
  ggsave(file.path(output_dir, paste0(filename, ".pdf")),
         plot = plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(filename, ".jpg")),
         plot = plot_obj, width = width, height = height, dpi = dpi)
}
# -------- CONSISTENT STATE ORDER + COLORS --------
state_levels <- c(
  "Mito-High / Glyco-Low" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High" = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High" = "Both-high\n(Hybrid)",
  "Both-Low" = "Both-low\n(Low)"
)

state_labels <- c(
  "Mito-High / Glyco-Low" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High" = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High" = "Both-high\n(Hybrid)",
  "Both-Low" = "Both-low\n(Low)"
)

state_colors <- c(
  "Mito-High / Glyco-Low" = "darkgreen",
  "Mito-Low / Glyco-High" = "orange",
  "Both-High" = "purple",
  "Both-Low" = "maroon"
)
# -------- GENE SETS --------
mito_genes <- c(
  "VDAC1","VDAC2","VDAC3",
  "SLC25A4","SLC25A5","SLC25A6","SLC25A54",
  "UCP2","UCP3",
  "PPIF","PPARGC1A","PPARGC1B",
  "TFAM","NRF1",
  "MFN1","MFN2","DNM1L","FIS1","MFF","OPA1",
  "BAX","BAK1","BCL2","BCL2L1","BBC3","BID","PMAIP1","MCL1","AIFM2",
  "GPX4","GSR","PRDX3","PRDX5","TXN2",
  "ACSL4","SLC7A11","FTL","FTH1"
)

glyco_genes <- c(
  "SLC2A1","HK1","HK2","GPI","PFKL","PFKM","PFKP",
  "ALDOA","ALDOB","ALDOC","TPI1","GAPDH","PGK1",
  "PGAM1","PGAM4","ENO1","ENO2","ENO3","PKM","LDHA","LDHB"
)

# -------- READ EXPRESSION --------
message("Reading CCLE expression...")

library(data.table)
library(dplyr)
library(tidyr)

expr <- fread(expr_file, data.table = FALSE)

# Extract cell line ID
expr$cell_line <- expr$ModelID

# Remove metadata columns
meta_cols <- c(
  "SequencingID","ModelID","IsDefault","ModelCondition",
  "IsDefaultEntryForModel","OncotreeLineage","OncotreePrimaryDisease"
)

meta_cols <- intersect(meta_cols, colnames(expr))
expr <- expr[, !(colnames(expr) %in% meta_cols)]

# Keep only numeric gene columns
gene_cols <- sapply(expr, is.numeric)
expr <- expr[, gene_cols | colnames(expr) == "cell_line"]

# Convert to long format
expr_long <- expr %>%
  pivot_longer(
    cols = -cell_line,
    names_to = "gene",
    values_to = "expr"
  ) %>%
  mutate(expr = as.numeric(expr))

# Remove NA
expr_long <- expr %>%
  pivot_longer(
    cols = -cell_line,
    names_to = "gene",
    values_to = "expr"
  ) %>%
  mutate(
    expr = as.numeric(expr),
    gene = toupper(gene),
    gene = sub(" \\(.*\\)$", "", gene)   # removes " (2268)" part
  ) %>%
  filter(!is.na(expr))

expr_long


# Make gene lists uppercase
mito_genes  <- toupper(mito_genes)
glyco_genes <- toupper(glyco_genes)

###Quick check before scoring

head(unique(expr_long$gene), 20)
intersect(mito_genes, unique(expr_long$gene))
intersect(glyco_genes, unique(expr_long$gene))
length(intersect(mito_genes, unique(expr_long$gene)))
length(intersect(glyco_genes, unique(expr_long$gene)))
# -------- SCORE MITO AND GLYCOLYSIS --------
message("Scoring metabolic programs...")

score_program <- function(df_long, genes, score_name) {
  df_long %>%
    filter(gene %in% genes) %>%
    group_by(cell_line) %>%
    summarise(!!score_name := mean(expr, na.rm = TRUE), .groups = "drop")
}

mito_score  <- score_program(expr_long, mito_genes, "mito_score")
glyco_score <- score_program(expr_long, glyco_genes, "glyco_score")

scores <- mito_score %>%
  inner_join(glyco_score, by = "cell_line") %>%
  mutate(
    mito_z  = as.numeric(scale(mito_score)),
    glyco_z = as.numeric(scale(glyco_score))
  )

# -------- ASSIGN STATES --------
scores <- scores %>%
  mutate(
    metabolic_state = case_when(
      mito_z >= 0 & glyco_z >= 0 ~ "Both-High",
      mito_z <  0 & glyco_z <  0 ~ "Both-Low",
      mito_z >= 0 & glyco_z <  0 ~ "Mito-High / Glyco-Low",
      mito_z <  0 & glyco_z >= 0 ~ "Mito-Low / Glyco-High",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    metabolic_state = factor(
      metabolic_state,
      levels = c("Both-High", "Both-Low", "Mito-High / Glyco-Low", "Mito-Low / Glyco-High")
    )
  )

write.csv(scores, file.path(output_dir, "CellLine_Metabolic_States.csv"), row.names = FALSE)

# -------- READ CELL META --------
message("Reading cell line metadata...")
cell_meta <- fread(cell_meta_file, data.table = FALSE)
colnames(cell_meta) <- make.names(colnames(cell_meta))
message("Cell meta columns: ", paste(colnames(cell_meta), collapse = ", "))

# -------- READ DRUG META --------
message("Reading treatment metadata...")
drug_meta <- fread(drug_meta_file, data.table = FALSE)
colnames(drug_meta) <- make.names(colnames(drug_meta))
message("Drug meta columns: ", paste(colnames(drug_meta), collapse = ", "))

# -------- READ PRISM LFC --------
message("Reading PRISM LFC...")
lfc <- fread(lfc_file, data.table = FALSE)
colnames(lfc) <- make.names(colnames(lfc))
message("LFC columns: ", paste(colnames(lfc), collapse = ", "))

# This PRISM file is already LONG format
# Expected columns include: row_id, broad_id, dose, compound, screen, culture, LFC

# standardize names
lfc_long <- lfc %>%
  dplyr::rename(
    cell_line_raw = row_id,
    treatment_id = broad_id,
    lfc = LFC
  ) %>%
  mutate(
    lfc = as.numeric(lfc)
  ) %>%
  filter(!is.na(lfc))

# extract DepMap ACH ID from cell_line_raw
# example: ACH-000824::P107::PR500A::REP1M -> ACH-000824
lfc_long <- lfc_long %>%
  mutate(
    cell_line = sub("::.*$", "", cell_line_raw)
  )

# save quick check
write.csv(lfc_long, file.path(output_dir, "PRISM_LFC_Long_Cleaned.csv"), row.names = FALSE)

# -------- MERGE STATES WITH DRUG RESPONSE --------
drug_state <- lfc_long %>%
  inner_join(scores, by = "cell_line")

write.csv(
  drug_state,
  file.path(output_dir, "Drug_Response_with_Metabolic_States.csv"),
  row.names = FALSE
)

# -------- MAP DRUG NAMES --------
# usually treatment metadata contains broad_id and compound/name columns

drug_meta_sub <- drug_meta

if ("broad_id" %in% colnames(drug_meta_sub)) {
  drug_meta_sub$treatment_id <- drug_meta_sub$broad_id
}

if ("name" %in% colnames(drug_meta_sub)) {
  drug_meta_sub$drug_name <- drug_meta_sub$name
} else if ("compound" %in% colnames(drug_meta_sub)) {
  drug_meta_sub$drug_name <- drug_meta_sub$compound
} else if ("pert_iname" %in% colnames(drug_meta_sub)) {
  drug_meta_sub$drug_name <- drug_meta_sub$pert_iname
} else {
  drug_meta_sub$drug_name <- drug_meta_sub$treatment_id
}

drug_state <- drug_state %>%
  left_join(
    drug_meta_sub %>% dplyr::select(treatment_id, drug_name) %>% distinct(),
    by = "treatment_id"
  )

write.csv(
  drug_state,
  file.path(output_dir, "Drug_Response_with_Metabolic_States_Annotated.csv"),
  row.names = FALSE
)

# quick checks
cat("Rows in lfc_long:", nrow(lfc_long), "\n")
cat("Rows in drug_state:", nrow(drug_state), "\n")
cat("Unique drugs:", length(unique(drug_state$treatment_id)), "\n")
cat("Unique cell lines matched:", length(unique(drug_state$cell_line)), "\n")

# ============================================================
# 3) SUMMARY + ANOVA ACROSS ALL 4 STATES
# ============================================================

# -------- DRUG-LEVEL SUMMARY --------
drug_summary <- drug_state %>%
  filter(!is.na(metabolic_state), !is.na(lfc)) %>%
  filter(!grepl("QC Failure", drug_name, ignore.case = TRUE)) %>%
  group_by(treatment_id, drug_name, metabolic_state) %>%
  summarise(
    mean_lfc = mean(lfc, na.rm = TRUE),
    median_lfc = median(lfc, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )
write.csv(
  drug_summary,
  file.path(output_dir, "Drug_State_Summary.csv"),
  row.names = FALSE
)

# -------- ANOVA PER DRUG --------
# -------- ANOVA PER DRUG --------
drug_anova <- drug_state %>%
  filter(!is.na(metabolic_state), !is.na(lfc)) %>%
  filter(!grepl("QC Failure", drug_name, ignore.case = TRUE)) %>%
  group_by(treatment_id, drug_name) %>%
  group_modify(~{
    dat <- .x
    if (n_distinct(dat$metabolic_state) < 2 || nrow(dat) < 10) {
      return(data.frame(p_value = NA_real_))
    }
    fit <- tryCatch(aov(lfc ~ metabolic_state, data = dat), error = function(e) NULL)
    if (is.null(fit)) return(data.frame(p_value = NA_real_))
    p <- tryCatch(summary(fit)[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    data.frame(p_value = p)
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj, p_value)

top_drugs <- drug_anova %>%
  filter(!is.na(p_adj)) %>%
  filter(!grepl("QC Failure", drug_name, ignore.case = TRUE)) %>%
  slice_head(n = 20)


write.csv(
  top_drugs,
  file.path(output_dir, "Top20_Differential_Drugs.csv"),
  row.names = FALSE
)

top_drugs50 <- drug_anova %>%
  filter(!is.na(p_adj)) %>%
  filter(!grepl("QC Failure", drug_name, ignore.case = TRUE)) %>%
  slice_head(n = 50)


write.csv(
  top_drugs50,
  file.path(output_dir, "Top50_Differential_Drugs.csv"),
  row.names = FALSE
)
# ============================================================
# FIGURE 3A: HEATMAP OF TOP STATE-DIFFERENTIAL DRUGS
# ============================================================

state_levels <- c(
  "Mito-High / Glyco-Low",
  "Mito-Low / Glyco-High",
  "Both-High",
  "Both-Low"
)

state_labels <- c(
  "Mito-High / Glyco-Low" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High" = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High" = "Both-high\n(Hybrid)",
  "Both-Low" = "Both-low\n(Low)"
)
top_heat <- drug_summary %>%
  semi_join(top_drugs50, by = c("treatment_id", "drug_name")) %>%
  mutate(
    drug_name = stringr::str_to_sentence(drug_name),
    metabolic_state = factor(metabolic_state, levels = state_levels)
  ) %>%
  filter(!is.na(metabolic_state))

p3A <- ggplot(
  top_heat,
  aes(
    x = metabolic_state,
    y = fct_reorder(drug_name, mean_lfc),
    fill = mean_lfc
  )
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = median(top_heat$mean_lfc, na.rm = TRUE),
    oob = scales::squish,
    name = "Mean LFC"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Top 50 State-Differential Drug Sensitivities",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_both(p3A, "Figure3AASupplementary_TopDifferential_Drugs_Heatmap", width = 11, height = 9)

# ============================================================
# FIGURE 3B: ANOVA EFFECT SIZE PLOT ACROSS ALL 4 STATES
# ============================================================

drug_anova <- drug_state %>%
  filter(!is.na(metabolic_state), !is.na(lfc)) %>%
  filter(!grepl("QC Failure", drug_name, ignore.case = TRUE)) %>%
  group_by(treatment_id, drug_name) %>%
  group_modify(~{
    dat <- .x
    
    if (n_distinct(dat$metabolic_state) < 4 || nrow(dat) < 20) {
      return(data.frame(
        p_value = NA_real_,
        effect_range = NA_real_
      ))
    }
    
    fit <- tryCatch(aov(lfc ~ metabolic_state, data = dat), error = function(e) NULL)
    
    p <- if (is.null(fit)) {
      NA_real_
    } else {
      summary(fit)[[1]][["Pr(>F)"]][1]
    }
    
    means <- dat %>%
      group_by(metabolic_state) %>%
      summarise(mean_lfc = mean(lfc, na.rm = TRUE), .groups = "drop")
    
    data.frame(
      p_value = p,
      effect_range = max(means$mean_lfc, na.rm = TRUE) - min(means$mean_lfc, na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = ifelse(p_adj < 0.05, "FDR < 0.05", "Not significant")
  ) %>%
  arrange(p_adj, desc(effect_range))

colnames(drug_anova)

write.csv(
  drug_anova %>% arrange(p_adj),
  file.path(output_dir, "Drug_ANOVA_All4States_sorted.csv"),
  row.names = FALSE
)
#===========================================================

anova_plot_df <- drug_anova %>%
  filter(!is.na(effect_range), !is.na(p_adj)) %>%
  arrange(p_adj, desc(effect_range)) %>%
  slice_head(n = 40) %>%
  mutate(
    drug_name = stringr::str_to_sentence(drug_name),   
    drug_label = ifelse(p_adj < 0.05, paste0(drug_name, "*"), drug_name),
    drug_label = factor(drug_label, levels = rev(drug_label))
  )

p3B <- ggplot(
  anova_plot_df,
  aes(x = effect_range, y = drug_label, fill = significant)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c(
    "FDR < 0.05" = "orangered",
    "Not significant" = "gray50"
  )) +
  theme_classic(base_size = 12) +
  labs(
    title = "All-State Differential Drug Sensitivity",
    subtitle = "ANOVA across four metabolic states; * FDR < 0.05",
    x = "Effect size (range of mean PRISM LFC across states)",
    y = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

save_both(p3B, "Figure3B_ANOVA_All4States_EffectSize", width = 9, height = 10)

# ============================================================
# FIGURE 3C: VOLCANO OXPHOS-DOMINANT VS GLYCOLYTIC 
# ============================================================

contrast_df <- drug_state %>%
  filter(metabolic_state %in% c("Mito-High / Glyco-Low", 
                                "Mito-Low / Glyco-High")) %>%
  filter(!grepl("QC Failure", drug_name, ignore.case = TRUE)) %>%
  mutate(drug_name = stringr::str_to_sentence(drug_name)) %>%
  group_by(treatment_id, drug_name) %>%
  group_modify(~{
    dat <- .x
    
    g1 <- dat %>%
      filter(metabolic_state == "Mito-High / Glyco-Low") %>%
      pull(lfc)
    
    g2 <- dat %>%
      filter(metabolic_state == "Mito-Low / Glyco-High") %>%
      pull(lfc)
    
    if (length(g1) < 3 || length(g2) < 3) {
      return(data.frame(effect = NA_real_, p_value = NA_real_))
    }
    
    tt <- tryCatch(t.test(g1, g2), error = function(e) NULL)
    if (is.null(tt)) {
      return(data.frame(effect = NA_real_, p_value = NA_real_))
    }
    
    data.frame(
      effect = mean(g1, na.rm = TRUE) - mean(g2, na.rm = TRUE),
      p_value = tt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neglog10_p = -log10(p_value),
    
    sensitivity_group = case_when(
      effect < 0 ~ "Mito-High / Glyco-Low",
      effect > 0 ~ "Mito-Low / Glyco-High",
      TRUE ~ "Neutral"
    ),
    
    #  FIRST create label
    sensitivity_group_label = case_when(
      sensitivity_group == "Mito-High / Glyco-Low" ~ 
        "Mito-high\nGlyco-low\n(OXPHOS-dependent)",
      sensitivity_group == "Mito-Low / Glyco-High" ~ 
        "Mito-low\nGlyco-high\n(Glycolytic)",
      TRUE ~ "Neutral"
    ),
    
    # THEN convert to factor
    sensitivity_group_label = factor(
      sensitivity_group_label,
      levels = c(
        "Mito-high\nGlyco-low\n(OXPHOS-dependent)",
        "Mito-low\nGlyco-high\n(Glycolytic)",
        "Neutral"
      )
    ),
    
    side = sensitivity_group,
    
    label_text = ifelse(
      !is.na(p_adj) & p_adj < 0.05,
      paste0(drug_name, "*"),
      drug_name
    )
  )

write.csv(
  contrast_df,
  file.path(output_dir, "Drug_Contrast_OXPHOS_vs_Glycolytic.csv"),
  row.names = FALSE
)

label_df <- contrast_df %>%
  filter(!is.na(effect), !is.na(p_adj)) %>%
  group_by(side) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 12) %>%
  ungroup()

count_label <- contrast_df %>%
  filter(sensitivity_group %in% c("Mito-High / Glyco-Low", 
                                  "Mito-Low / Glyco-High")) %>%
  group_by(sensitivity_group_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(text = paste0(
    ifelse(grepl("OXPHOS", sensitivity_group_label),
           "Mito-high (OXPHOS-dependent)",
           "Mito-low (Glycolytic)"),
    ": ", n
  )) %>%
  pull(text) %>%
  paste(collapse = "\n")

p3C <- ggplot(
  contrast_df,
  aes(x = effect, y = neglog10_p, color = sensitivity_group)
) +
  geom_point(alpha = 0.75, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_df,
    aes(label = label_text),
    size = 3,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(state_colors, "Neutral" = "grey70"),
    labels = c(state_labels, "Neutral" = "Neutral")
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Differential Drug Sensitivity: OXPHOS-dominant vs Glycolytic States",
    x = "Mean LFC difference",
    y = "-log10(P value)",
    color = "More sensitive state"
  )

save_both(
  p3C,
  "Figure3C_Volcano_OXPHOSdominant_vs_Glycolytic",
  width = 10,
  height = 8
)


####+======================================================
####+SUpplementary 3
####+ -=====================================================

label_df <- contrast_df %>%
  filter(!is.na(effect), !is.na(p_adj)) %>%
  group_by(side) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 20) %>%
  ungroup()

count_label <- contrast_df %>%
  filter(sensitivity_group %in% c("Mito-High / Glyco-Low", 
                                  "Mito-Low / Glyco-High")) %>%
  group_by(sensitivity_group_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(text = paste0(
    ifelse(grepl("OXPHOS", sensitivity_group_label),
           "Mito-high (OXPHOS-dependent)",
           "Mito-low (Glycolytic)"),
    ": ", n
  )) %>%
  pull(text) %>%
  paste(collapse = "\n")

p3C <- ggplot(
  contrast_df,
  aes(x = effect, y = neglog10_p, color = sensitivity_group)
) +
  geom_point(alpha = 0.75, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_df,
    aes(label = label_text),
    size = 3,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(state_colors, "Neutral" = "grey70"),
    labels = c(state_labels, "Neutral" = "Neutral")
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Differential Drug Sensitivity: OXPHOS-dominant vs Glycolytic States",
    x = "Mean LFC difference",
    y = "-log10(P value)",
    color = "More sensitive state"
  )

save_both(
  p3C,
  "Figure3Bsupplementary_Volcano_OXPHOSdominant_vs_Glycolytic",
  width = 10,
  height = 8
)

# ============================================================
# FIGURE 3D: REPRESENTATIVE BOXPLOTS, ALL 4 STATES
# ============================================================

# For main Figure 3D: balanced top 6
top_oxphos <- contrast_df %>%
  filter(effect < 0, !is.na(p_adj)) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 3)

top_glyco <- contrast_df %>%
  filter(effect > 0, !is.na(p_adj)) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 3)

top_combined <- bind_rows(top_oxphos, top_glyco) %>%
  mutate(
    drug_name = stringr::str_to_sentence(drug_name),
    sensitive_state = ifelse(effect < 0, "OXPHOS-sensitive", "Glycolytic-sensitive"),
    drug_label = ifelse(
      p_adj < 0.05,
      paste0(drug_name, "*\n[", sensitive_state, "]"),
      paste0(drug_name, "\n[", sensitive_state, "]")
    )
  )

plot_box_dat <- drug_state %>%
  mutate(drug_name = stringr::str_to_sentence(drug_name)) %>%
  filter(drug_name %in% top_combined$drug_name) %>%
  left_join(top_combined %>% select(drug_name, drug_label), by = "drug_name") %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels),
    drug_label = factor(drug_label, levels = top_combined$drug_label)
  )

p3D <- ggplot(
  plot_box_dat,
  aes(x = metabolic_state, y = lfc, fill = metabolic_state)
) +
  geom_boxplot(outlier.shape = NA, width = 0.65, color = "black") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 0.55, color = "black") +
  facet_wrap(~ drug_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Representative State-Specific Drug Responses",
    x = NULL,
    y = "PRISM log2 fold-change viability (LFC)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 8),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_both(p3D, "Figure3D_Representative_Drug_Boxplots", width = 12, height = 8)

# ============================================================
# FIGURE 3supplementary 3c: REPRESENTATIVE BOXPLOTS, ALL 4 STATES
# ============================================================


top_oxphos_supp <- contrast_df %>%
  filter(effect < 0, !is.na(p_adj)) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 10)

top_glyco_supp <- contrast_df %>%
  filter(effect > 0, !is.na(p_adj)) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 10)

top_combined_supp <- bind_rows(top_oxphos_supp, top_glyco_supp) %>%
  mutate(
    drug_name = stringr::str_to_sentence(drug_name),
    sensitive_state = ifelse(effect < 0, "OXPHOS-sensitive", "Glycolytic-sensitive"),
    drug_label = ifelse(
      p_adj < 0.05,
      paste0(drug_name, "*\n[", sensitive_state, "]"),
      paste0(drug_name, "\n[", sensitive_state, "]")
    )
  )

plot_box_dat_supp <- drug_state %>%
  mutate(drug_name = stringr::str_to_sentence(drug_name)) %>%
  filter(drug_name %in% top_combined_supp$drug_name) %>%
  left_join(
    top_combined_supp %>% select(drug_name, drug_label),
    by = "drug_name"
  ) %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels),
    drug_label = factor(drug_label, levels = top_combined_supp$drug_label)
  )

p3D_supp <- ggplot(
  plot_box_dat_supp,
  aes(x = metabolic_state, y = lfc, fill = metabolic_state)
) +
  geom_boxplot(outlier.shape = NA, width = 0.65, color = "black") +
  geom_jitter(width = 0.18, alpha = 0.25, size = 0.5, color = "black") +
  
  facet_wrap(~ drug_label, scales = "free_y", ncol = 4) +  
  
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  
  theme_bw(base_size = 11) +
  labs(
    title = "Supplementary: State-Specific Drug Responses (Top 20 Balanced)",
    x = NULL,
    y = "PRISM log2 fold-change viability (LFC)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 7),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_both(
  p3D_supp,
  "FigureS3Csupplementary_Top20_Balanced_Drug_Boxplots",
  width = 16,
  height = 14
)

# ============================================================
# FIGURE 3E: BALANCED TOP DIFFERENTIAL DRUGS
# Top 15 OXPHOS-sensitive + Top 15 Glycolytic-sensitive
# ============================================================

library(dplyr)
library(stringr)
library(ggplot2)

# ---- Select top OXPHOS-sensitive drugs ----
top_oxphos_E <- contrast_df %>%
  filter(!is.na(effect), !is.na(p_adj), effect < 0) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 15)

# ---- Select top Glycolytic-sensitive drugs ----
top_glyco_E <- contrast_df %>%
  filter(!is.na(effect), !is.na(p_adj), effect > 0) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 15)

# ---- Combine and format ----
p3E_df <- bind_rows(top_oxphos_E, top_glyco_E) %>%
  mutate(
    drug_name = stringr::str_to_sentence(drug_name),
    sensitivity_class = ifelse(effect < 0, 
                               "OXPHOS-sensitive", 
                               "Glycolytic-sensitive"),
    drug_label = ifelse(
      p_adj < 0.05,
      paste0(drug_name, "*"),
      drug_name
    ),
    drug_label = factor(drug_label, levels = drug_label[order(effect)])
  )

# ---- Save source data ----
write.csv(
  p3E_df,
  file.path(output_dir, "Figure3ER1_Balanced_Top_Differential_Drugs_Data.csv"),
  row.names = FALSE
)

# ---- Plot ----
p3E <- ggplot(
  p3E_df,
  aes(x = effect, y = drug_label, fill = sensitivity_class)
) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(
    "OXPHOS-sensitive" = "darkgreen",
    "Glycolytic-sensitive" = "orange"
  )) +
  theme_classic(base_size = 12) +
  labs(
    title = "Top OXPHOS- and Glycolytic-Sensitive Drugs Across Metabolic States",
    subtitle = "Top 15 OXPHOS-sensitive and top 15 glycolytic-sensitive drugs; * FDR < 0.05",
    x = "Mean LFC difference (OXPHOS − Glycolytic)",
    y = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p3E)

# ---- Save ----
save_both(
  p3E,
  "Figure3ER1_Balanced_Top15_OXPHOS_Top15_Glycolytic",
  width = 9,
  height = 8
)

####====================================================
##Supplementary Figure 3c: Correlation volcano plot
####====================================================

label_drugs <- corr_summary_all %>%
  filter(padj < 0.05, R >= 0.15)

# ---- Correlation across ALL drugs ----
corr_summary_all <- drug_state %>%
  filter(!is.na(mito_z), !is.na(lfc), !is.na(drug_name)) %>%
  group_by(drug_name) %>%
  summarise(
    n = n(),
    R = suppressWarnings(cor(mito_z, lfc, method = "spearman")),
    p = suppressWarnings(cor.test(mito_z, lfc, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  filter(n >= 5) %>%
  mutate(
    padj = p.adjust(p, method = "BH"),
    neglog10_padj = -log10(padj),
    direction = case_when(
      abs(R) < 0.15 ~ "Weak association",
      R >= 0.15 & padj < 0.05 ~ "Positive (meaningful)",
      TRUE ~ "Not significant"
    ),
    drug_name = str_to_sentence(drug_name),
    drug_name = str_remove(drug_name, "\\s*-\\s*qc failure")
  )

corr_colors <- c(
  "Positive (meaningful)" = "orangered",
  "Weak association" = "grey70",
  "Not significant" = "black"
)

# ---- Volcano plot ----
pS_corr_volcano <- ggplot(corr_summary_all, aes(x = R, y = neglog10_padj)) +
  geom_point(aes(color = direction), alpha = 0.75, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0.15, linetype = "dashed", color = "red") +   
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = label_drugs,
    aes(label = drug_name),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5
  ) +
  scale_color_manual(values = corr_colors) +
  theme_bw(base_size = 12) +
  labs(
    title = "Global association between mitochondrial score and PRISM drug sensitivity",
    x = "Spearman correlation with mitochondrial score (R)",
    y = expression(-log[10]~"(BH-adjusted p-value)"),
    color = "Association"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

save_both(
  pS_corr_volcano,
  "Supplementary_Global_MitoScore_PRISM_Correlation_Volcano",
  width = 10,
  height = 6
)

# =========================================================
# Figure 3F: Pathway-level drug sensitivity
# Top 15 OXPHOS-sensitive + top 15 glycolytic-sensitive
# =========================================================
# ---- 1. Select top 15 from each state ----
top15_each <- top_combined %>%
  mutate(
    drug_name = str_to_sentence(drug_name),
    sensitive_state = case_when(
      effect < 0 ~ "OXPHOS-sensitive",
      effect > 0 ~ "Glycolytic-sensitive",
      TRUE ~ NA_character_
    ),
    sensitive_state = factor(
      sensitive_state,
      levels = c("OXPHOS-sensitive", "Glycolytic-sensitive")
    )
  ) %>%
  filter(!is.na(effect), !is.na(p_adj), !is.na(sensitive_state)) %>%
  group_by(sensitive_state) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 15) %>%
  ungroup()

print(top15_each %>% select(drug_name, sensitive_state, effect, p_adj))

write.csv(
  top15_each %>% select(drug_name, sensitive_state, effect, p_adj),
  "Top15_each_state_drugs_for_pathway_annotation.csv",
  row.names = FALSE
)

# ---- 2. Curated pathway annotation ----
drug_pathway_map <- data.frame(
  drug_name = c(
    "Cep-40783", "Bms-986158", "Azaguanine-8", "Verdinexor",
    "Senaparib", "Ilorasertib", "Encorafenib", "Pd-161570",
    "Eravacycline", "Mps1-in-5", "Guadecitabine", "Ai-10-49",
    "Pbt-1033", "Pemigatinib", "Baloxavir-marboxil",
    "Elexacaftor", "Decloxizine", "Nbi-98782", "Efonidipine",
    "Pbi-4050", "Uprifosbuvir", "Sar439859", "Jnj-54175446",
    "Furamidine", "Setmelanotide", "Cloperastine-fendizoate",
    "Ml213", "Nsc-636819", "Lylamine", "Spautin-1"
  ),
  pathway = c(
    "Multi-kinase signaling",
    "BET / epigenetic regulation",
    "Purine metabolism",
    "Nuclear export",
    "PARP / DNA damage response",
    "Aurora kinase / mitosis",
    "BRAF / MAPK signaling",
    "Receptor tyrosine kinase signaling",
    "Protein synthesis inhibition",
    "MPS1 / mitotic checkpoint",
    "DNA methylation",
    "Transcriptional regulation",
    "Metabolic stress response",
    "FGFR signaling",
    "RNA metabolism / antiviral",
    "Ion transport / CFTR modulation",
    "Antihistamine / ion channel activity",
    "VMAT2 / vesicular transport",
    "Calcium channel signaling",
    "Fibrosis / TGF-beta signaling",
    "Viral replication / protease",
    "Estrogen receptor signaling",
    "P2X7 receptor signaling",
    "DNA minor-groove binding",
    "Melanocortin receptor signaling",
    "Histamine receptor signaling",
    "Potassium channel modulation",
    "Unannotated / experimental",
    "Lipid / amine metabolism",
    "Autophagy"
  )
)

# ---- 3. Merge pathway annotation ----
pathway_dat <- top15_each %>%
  left_join(drug_pathway_map, by = "drug_name") %>%
  mutate(
    pathway = ifelse(is.na(pathway), "Unannotated / experimental", pathway)
  )

print(
  pathway_dat %>%
    select(drug_name, sensitive_state, effect, p_adj, pathway)
)

# ---- 4. Summarize mean effect by pathway ----
pathway_summary <- pathway_dat %>%
  group_by(sensitive_state, pathway) %>%
  summarise(
    mean_effect = mean(effect, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    pathway_label = paste0(pathway, " (n=", n, ")"),
    pathway_label = fct_reorder(pathway_label, mean_effect)
  )

# For main Figure 3D: balanced top 6
top_oxphos <- contrast_df %>%
  filter(effect < 0, !is.na(p_adj)) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 15)

top_glyco <- contrast_df %>%
  filter(effect > 0, !is.na(p_adj)) %>%
  arrange(p_adj, desc(abs(effect))) %>%
  slice_head(n = 15)

top_combined <- bind_rows(top_oxphos, top_glyco) %>%
  mutate(
    drug_name = stringr::str_to_sentence(drug_name),
    sensitive_state = ifelse(effect < 0, "OXPHOS-sensitive", "Glycolytic-sensitive"),
    drug_label = ifelse(
      p_adj < 0.05,
      paste0(drug_name, "*\n[", sensitive_state, "]"),
      paste0(drug_name, "\n[", sensitive_state, "]")
    )
  )
# ---- 5. Plot ----
state_colors <- c(
  "OXPHOS-sensitive" = "darkgreen",
  "Glycolytic-sensitive" = "orange"
)

p3F_pathway <- ggplot(
  pathway_dat,
  aes(x = effect, y = reorder(drug_name, effect), fill = pathway)
) +
  geom_col(color = "black") +
  facet_wrap(~sensitive_state, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Pathway-associated drug sensitivities across metabolic states",
    x = "LFC difference (OXPHOS − Glycolytic)",
    y = "Drug",
    fill = "Pathway"
  ) +
  guides(
    fill = guide_legend(
      ncol = 1   
    )
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

save_both(
  p3F_pathway,
  "Figure3F_Pathway_Level_Top15_Each_State",
  width = 10,
  height = 8
)
