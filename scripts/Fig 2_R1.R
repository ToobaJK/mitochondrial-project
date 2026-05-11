# ============================================================
# FIGURE 2 (DepMap): State-specific CRISPR dependencies
CREATED BY: TOOBA
# ============================================================

# ---- Packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
})

# ---- Paths  ----
base_dir  <- "D:/UAEU/Depmap data"
expr_file <- file.path(base_dir, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
crispr_file <- file.path(base_dir, "CRISPRGeneEffect.csv")
model_file  <- file.path(base_dir, "Model.csv")

# ---- Output folder for Figure 2 ----
out_dir <- file.path(base_dir, "Figure2_outputsR2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Gene sets ----
gatekeepers <- c(
  "VDAC1","VDAC2","VDAC3",
  "SLC25A4","SLC25A5","SLC25A6",
  "PPIF","UCP2","UCP3",
  "BAX","BAK1","BID","MCL1","BCL2","BCL2L1","BBC3",
  "DNM1L","MFN1","MFN2","FIS1","MFF","OPA1","TFAM","NRF1",
  "PPARGC1A","PPARGC1B",
  "SLC7A11","GPX4","AIFM2","PRDX3","PRDX5","TXN2","GSR","ACSL4",
  "FTH1","FTL","PMAIP1"
)

# Optional glycolysis set (only if you want 4-state like Fig1D)
glyco_genes <- c("SLC2A1","HK2","PFKP","ALDOA","GAPDH","PGK1","ENO1","PKM","LDHA")

# ---- Helpers ----
safe_keep_cols <- function(dt, cols) {
  cols <- cols[cols %in% colnames(dt)]
  dt[, ..cols]
}

zscore_vec <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

# ============================================================
# 1) Load data (FAST)
# ============================================================
message("Loading DepMap files (these are large)...")
expr   <- fread(expr_file)
crispr <- fread(crispr_file)
model  <- fread(model_file)

# ============================================================
# 1b) FIX  DepMap column naming issues
# ============================================================

# --- Standardize ID column to ModelID ---
# CRISPR: V1 is the cell line ID (ACH-xxxx)
if ("V1" %in% names(crispr) && !"ModelID" %in% names(crispr)) {
  setnames(crispr, "V1", "ModelID")
}

# Expression: in your file ModelID already exists, but keep this as safety
if ("V1" %in% names(expr) && !"ModelID" %in% names(expr)) {
  setnames(expr, "V1", "ModelID")
}

# Model: ModelID already exists, but keep safety
if ("V1" %in% names(model) && !"ModelID" %in% names(model)) {
  setnames(model, "V1", "ModelID")
}

# --- Clean gene column names: "TP53 (7157)" -> "TP53" ---
clean_gene_names <- function(nms) gsub(" \\(.*\\)$", "", nms)

# Clean CRISPR gene columns (keep ModelID unchanged)
setnames(crispr,
         old = names(crispr),
         new = ifelse(names(crispr) == "ModelID", "ModelID", clean_gene_names(names(crispr))))

# Clean expression gene columns too (optional but safe)
setnames(expr,
         old = names(expr),
         new = ifelse(names(expr) == "ModelID", "ModelID", clean_gene_names(names(expr))))

# ============================================================
# 1c) Check IDs exist now
# ============================================================
stopifnot("ModelID" %in% names(expr),
          "ModelID" %in% names(crispr),
          "ModelID" %in% names(model))

message("OK: ModelID is present in expr, crispr, and model. Proceeding...")

# ============================================================
# 2) Build mitochondrial score + state assignment
# ============================================================

# Keep only needed columns
expr_sub <- safe_keep_cols(expr, c("ModelID", gatekeepers, glyco_genes))

# Ensure numeric
gene_cols <- setdiff(colnames(expr_sub), "ModelID")
expr_sub[, (gene_cols) := lapply(.SD, as.numeric), .SDcols = gene_cols]

# ---- Calculate Mito score ----
present_gatekeepers <- gatekeepers[gatekeepers %in% colnames(expr_sub)]

expr_sub[, MitoScore := rowMeans(.SD, na.rm = TRUE),
         .SDcols = present_gatekeepers]

# ---- Calculate Glycolysis score (if available) ----
present_glyco <- glyco_genes[glyco_genes %in% colnames(expr_sub)]

if (length(present_glyco) > 0) {
  expr_sub[, GlycoScore := rowMeans(.SD, na.rm = TRUE),
           .SDcols = present_glyco]
} else {
  expr_sub[, GlycoScore := NA_real_]
}

# ---- Define metabolic states (median split) ----
# ---- Z-score normalization of metabolic scores ----
expr_sub[, `:=`(
  MitoScore_Z  = as.numeric(scale(MitoScore)),
  GlycoScore_Z = as.numeric(scale(GlycoScore))
)]

# ---- Define metabolic states (full descriptive + biological labels) ----
expr_sub[, State := fifelse(
  MitoScore_Z >= 0 & GlycoScore_Z < 0,  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  fifelse(
    MitoScore_Z < 0 & GlycoScore_Z >= 0, "Mito-Low / Glyco-High (Glycolytic)",
    fifelse(
      MitoScore_Z >= 0 & GlycoScore_Z >= 0, "Both-High (Hybrid)",
      "Both-Low (Low metabolic)"
    )
  )
)]

# ---- Set factor levels (consistent order across figures) ----
expr_sub[, State := factor(State,
                           levels = c(
                             "Mito-High / Glyco-Low (OXPHOS-dominant)",
                             "Mito-Low / Glyco-High (Glycolytic)",
                             "Both-High (Hybrid)",
                             "Both-Low (Low metabolic)"
                           ))]

# ---- Check distribution ----
table(expr_sub$State)

# ============================================================
# 3) Merge with Model metadata
# ============================================================

meta_cols <- c(
  "ModelID",
  "CellLineName",
  "StrippedCellLineName",
  "OncotreeLineage",
  "OncotreePrimaryDisease",
  "PrimaryDisease",
  "Lineage"
)

meta <- safe_keep_cols(model, meta_cols)

lineage_col <- intersect(c("OncotreeLineage", "Lineage"), colnames(meta))[1]
disease_col <- intersect(c("OncotreePrimaryDisease", "PrimaryDisease"), colnames(meta))[1]
name_col    <- intersect(c("StrippedCellLineName", "CellLineName", "ModelID"), colnames(meta))[1]

# IMPORTANT: keep Z scores also
expr_meta <- merge(
  expr_sub[, .(
    ModelID,
    MitoScore,
    GlycoScore,
    MitoScore_Z,
    GlycoScore_Z
  )],
  meta,
  by = "ModelID",
  all.x = TRUE
)

# ============================================================
# 4) Merge with CRISPR dependency
# ============================================================

crispr_sub <- safe_keep_cols(crispr, c("ModelID", gatekeepers))

dat <- merge(expr_meta, crispr_sub, by = "ModelID", all = FALSE)

# ============================================================
# 5) Re-create State AFTER merge
# ============================================================

dat[, State := fifelse(
  MitoScore_Z >= 0 & GlycoScore_Z < 0,
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  fifelse(
    MitoScore_Z < 0 & GlycoScore_Z >= 0,
    "Mito-Low / Glyco-High (Glycolytic)",
    fifelse(
      MitoScore_Z >= 0 & GlycoScore_Z >= 0,
      "Both-High (Hybrid)",
      "Both-Low (Low metabolic)"
    )
  )
)]

state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

dat[, State := factor(State, levels = state_levels)]

# Check distribution
table(dat$State, useNA = "ifany")

# Save fixed merged table
fwrite(dat, file.path(out_dir, "Figure2_DepMap_Merged_Table_FIXED.csv"))



# ============================================================
# FIGURE 2 – DepMap CRISPR dependency across ALL gatekeeper genes
# ============================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Settings
# ---------------------------
# ---- Define consistent state levels ----
state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

# ---- Ensure State column uses these levels ----
dat$State <- factor(dat$State, levels = state_levels)

# ---- Genes to analyze ----
genes_to_run <- gatekeepers[gatekeepers %in% colnames(dat)]

# ---- Pairwise comparisons (same naming) ----
my_comparisons <- list(
  c(state_levels[1], state_levels[2]),
  c(state_levels[1], state_levels[3]),
  c(state_levels[1], state_levels[4]),
  c(state_levels[2], state_levels[3]),
  c(state_levels[2], state_levels[4]),
  c(state_levels[3], state_levels[4])
)

# ---- Plot settings ----
show_pairwise   <- FALSE
label_fontsize  <- 2.8
plot_w <- 9
plot_h <- 5
jpg_dpi <- 300

# ---------------------------
# Prepare base dataframe once
# ---------------------------
df_base <- dat %>%
  filter(!is.na(State)) %>%
  mutate(State = factor(State, levels = state_levels))

# ---------------------------
# Storage for master outputs
# ---------------------------
master_overall  <- list()
master_pairwise <- list()

write.csv(dat, file.path(out_dir, "DepMap_with_states1.csv"), row.names = FALSE)
# ============================================================
# Loop over genes
# ============================================================
for (g in genes_to_run) {
  
  # Gene-specific dataframe
  df <- df_base %>%
    filter(!is.na(.data[[g]])) %>%
    transmute(ModelID, State, value = .data[[g]])
  
  if (nrow(df) < 30) next
  
  # ---- Summary stats for labels (mean±SD + n) ----
  sumtab <- df %>%
    group_by(State) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      ymax = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("n=", n, "\nmean=", sprintf("%.2f", mean), "\nSD=", sprintf("%.2f", sd))
    ) %>%
    as.data.frame()
  
  # Use ONE common label height per gene (prevents some groups getting cut)
  top_y <- max(df$value, na.rm = TRUE)
  sumtab$y <- top_y + 0.25
  
  # Save per-gene summary stats CSV
  write.csv(
    sumtab[, c("State","n","mean","sd")],
    file.path(out_dir, paste0("Fig2_", g, "_SummaryStats_meanSD1.csv")),
    row.names = FALSE
  )
  
  # ---- Overall p-value (Kruskal–Wallis) ----
  p_overall <- tryCatch(
    kruskal.test(value ~ State, data = df)$p.value,
    error = function(e) NA_real_
  )
  
  master_overall[[g]] <- data.frame(
    Gene = g,
    Test = "Kruskal",
    Pvalue = p_overall,
    N = nrow(df)
  )
  
  # ---- Pairwise Wilcoxon (BH) ----
  pw <- tryCatch(
    pairwise.wilcox.test(df$value, df$State, p.adjust.method = "BH"),
    error = function(e) NULL
  )
  
  if (!is.null(pw)) {
    pw_df <- as.data.frame(as.table(pw$p.value))
    colnames(pw_df) <- c("Group1","Group2","p_adj")
    pw_df <- pw_df[!is.na(pw_df$p_adj), ]
    if (nrow(pw_df) > 0) {
      pw_df$Gene <- g
      master_pairwise[[g]] <- pw_df
      write.csv(
        pw_df,
        file.path(out_dir, paste0("Fig2_", g, "_PairwiseP_BH.csv")),
        row.names = FALSE
      )
    }
  }
  
  # ---- Build plot ----
  p <- ggplot(df, aes(x = State, y = value)) +
    geom_boxplot(outlier.size = 0.35) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.6) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(g, " Dependency Across Mitochondrial Metabolic States"),
      subtitle = paste0("Overall Kruskal–Wallis p = ", signif(p_overall, 3)),
      y = "CERES gene effect (more negative = more essential)",
      x = NULL
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    geom_text(
      data = sumtab,
      aes(x = State, y = y, label = label),
      inherit.aes = FALSE,
      size = label_fontsize
    ) +
    # Prevent label clipping (adds headroom + margins)
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(12, 12, 12, 12))
  
  if (show_pairwise) {
    p <- p + stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      p.adjust.method = "BH",
      label = "p.adj",
      hide.ns = TRUE
    )
  }
  
  print(p)
  
  # ---- Save figure (PDF + JPG) ----
  ggsave(
    filename = file.path(out_dir, paste0("Fig2_", g, "_Dependency.pdf")),
    plot = p,
    width = plot_w, height = plot_h
  )
  
  ggsave(
    filename = file.path(out_dir, paste0("Fig2_", g, "_Dependency.jpg")),
    plot = p,
    width = plot_w, height = plot_h, dpi = jpg_dpi
  )
}

# ============================================================
# Save MASTER tables (after loop)
# ============================================================

overall_tbl <- do.call(rbind, master_overall)
overall_tbl$FDR_allGenes <- p.adjust(overall_tbl$Pvalue, method = "BH")
overall_tbl <- overall_tbl[order(overall_tbl$FDR_allGenes), ]

write.csv(
  overall_tbl,
  file.path(out_dir, "Fig2_AllGenes_OverallP_Kruskal_BH.csv"),
  row.names = FALSE
)

if (length(master_pairwise) > 0) {
  pairwise_tbl <- do.call(rbind, master_pairwise)
  pairwise_tbl <- pairwise_tbl[, c("Gene","Group1","Group2","p_adj")]
  write.csv(
    pairwise_tbl,
    file.path(out_dir, "Fig2_AllGenes_PairwiseP_Wilcox_BH.csv"),
    row.names = FALSE
  )
}

message("DONE. All plots (PDF+JPG) + per-gene CSV + master CSVs saved to: ", out_dir)
message("Tip: If you want pairwise p-values on the plot, set show_pairwise <- TRUE and rerun.")


# ============================================================
# 6) supplementarey 2A
# ============================================================

state_levels <- c( "Mito-High / Glyco-Low (OXPHOS-dominant)",
                   "Mito-Low / Glyco-High (Glycolytic)",
                   "Both-High (Hybrid)",
                   "Both-Low (Low metabolic)")

genes_to_run <- gatekeepers[gatekeepers %in% colnames(dat)]

# IMPORTANT: convert dat (data.table) -> data.frame so dplyr works
dat_df <- as.data.frame(dat)

# Long format
df_long <- dat_df %>%
  filter(!is.na(State)) %>%
  mutate(State = factor(State, levels = state_levels)) %>%
  dplyr::select(ModelID, State, all_of(genes_to_run)) %>%
  pivot_longer(
    cols = -c(ModelID, State),
    names_to = "Gene",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# Per-gene overall p-values (Kruskal)
pvals <- df_long %>%
  group_by(Gene) %>%
  summarise(p = kruskal.test(value ~ State)$p.value, .groups = "drop") %>%
  mutate(p_lab = paste0("Kruskal–Wallis p=", signif(p, 3)))

df_long <- left_join(df_long, pvals, by = "Gene")

# Put p-value text inside each facet (top-right)
pos <- df_long %>%
  group_by(Gene) %>%
  summarise(
    y = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pvals, by = "Gene") %>%
  mutate(x = "Mito-High / Glyco-Low (OXPHOS-dominant)")  # choose a position column

state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)
# Remove rows where State did not match the defined levels
df_long <- df_long %>%
  filter(!is.na(State))
p_all <- ggplot(df_long, aes(x = State, y = value, fill = State)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_jitter(width = 0.12, alpha = 0.2, size = 0.35) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  
  geom_text(
    data = pos,
    aes(x = 2.5, y = y, label = p_lab),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = -0.3,
    size = 2.5
  ) +
  
  #scale_fill_manual(values = state_colors, drop = FALSE) +
  
  scale_x_discrete(
    drop = FALSE,
    labels = c(
      "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
      "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
      "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
      "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
    )
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 10) +
  labs(
    title = "DepMap CRISPR Dependency Across Mitochondrial Metabolic States",
    y = "CERES gene effect (more negative = more essential)",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(15, 15, 25, 15)
  )

print(p_all)

ggsave(file.path(out_dir, "Figure2Supplementary_AllGenes_Multipanel.pdf"),
       p_all, width = 19, height = 14)

ggsave(file.path(out_dir, "Figure2Supplementary_AllGenes_Multipanel.jpg"),
       p_all, width = 19, height = 14, dpi = 300)

message("Saved supplementary all-gene multipanel figure (PDF + JPG) in: ", out_dir)
# ============================================================
# 6) FIGURE 2A: State-specific dependency for key genes 
# ============================================================
key_genes <- c("GPX4","SLC7A11","DNM1L","MFN2","SLC25A4","VDAC1","BAX","BAK1")
key_genes <- key_genes[key_genes %in% colnames(dat)]

make_dep_box <- function(g) {
  
  df <- as.data.frame(dat)
  
  # Keep only needed columns (and lineage if you want it later)
  df <- dplyr::select(df, ModelID, State, dplyr::all_of(g))
  
  # Create a common column name "Dependency" (no rename() needed)
  df$Dependency <- df[[g]]
  
  # Filter
  df <- dplyr::filter(df, !is.na(State) & !is.na(Dependency))
  
  # Optional: remove extreme outliers (comment out if you want full range)
  df <- dplyr::filter(df, Dependency > -3 & Dependency < 1)
  
  state_colors <- c(
    "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
    "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
    "Both-High (Hybrid)"                      = "purple",
    "Both-Low (Low metabolic)"                = "maroon"
  )
  
  ggplot(df, aes(x = State, y = Dependency, fill = State)) +
    geom_boxplot(outlier.size = 0.7) +
    geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
    
    scale_fill_manual(values = state_colors, drop = FALSE) +
    
    scale_x_discrete(
      drop = FALSE,
      labels = c(
        "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dependent)",
        "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
        "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
        "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
      )
    ) +
    
    theme_classic(base_size = 12) +
    labs(
      title = paste0(g, " Dependency (CRISPR gene effect)"),
      x = NULL,
      y = "CERES gene effect (more negative = more essential)"
    ) +
    
    theme(
      axis.text.x = element_text(size = 8),
      legend.position = "none"
    )
}

if (length(key_genes) > 0) {
  for (g in key_genes) {
    p <- make_dep_box(g)
    ggsave(file.path(out_dir, paste0("Fig2B_Dependency_Box_", g, ".pdf")),
           p, width = 7.5, height = 5)
    ggsave(file.path(out_dir, paste0("Fig2B_Dependency_Box_", g, ".png")),
           p, width = 7.5, height = 5, dpi = 300)
  }
}

message("Saved key gene dependency plots (PDF + PNG) in: ", out_dir)



# ============================================================
# FIGURE 2A: Key gene dependencies....1 plot
# ============================================================

key_genes <- c("GPX4","SLC7A11","DNM1L","MFN2",
               "SLC25A4","VDAC1","BAX","BAK1")

key_genes <- key_genes[key_genes %in% colnames(dat)]

state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

df_long <- as.data.frame(dat) %>%
  dplyr::select(ModelID, State, dplyr::all_of(key_genes)) %>%
  mutate(State = factor(State, levels = state_levels)) %>%
  pivot_longer(
    cols = -c(ModelID, State),
    names_to = "Gene",
    values_to = "Dependency"
  ) %>%
  filter(!is.na(State), !is.na(Dependency)) %>%
  filter(Dependency > -3 & Dependency < 1)

pvals_key <- df_long %>%
  group_by(Gene) %>%
  summarise(
    p = kruskal.test(Dependency ~ State)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_lab = paste0("Kruskal–Wallis p = ", signif(p, 3)))

pos_key <- df_long %>%
  group_by(Gene) %>%
  summarise(
    y = max(Dependency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pvals_key, by = "Gene") %>%
  mutate(x = "Mito-High / Glyco-Low (OXPHOS-dominant)")

state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)

p_multi <- ggplot(df_long, aes(x = State, y = Dependency, fill = State)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  
  geom_text(
    data = pos_key,
    aes(x = 2.5, y = y, label = p_lab),
    inherit.aes = FALSE,
    size = 2.8,
    hjust = 0.5,
    vjust = -0.3
  ) +
  
  scale_fill_manual(values = state_colors) +
  
  scale_x_discrete(
    labels = c(
      "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
      "Mito-Low / Glyco-High (Glycolytic)" = "Mito-low\nGlyco-high\n(Glycolytic)",
      "Both-High (Hybrid)" = "Both-high\n(Hybrid)",
      "Both-Low (Low metabolic)" = "Both-low\n(Low)"
    )
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 11) +
  labs(
    title = "State-Specific CRISPR Dependency of Mitochondrial Metabolic Genes",
    x = NULL,
    y = "CERES gene effect (more negative = more essential)"
  ) +
  theme(
    axis.text.x = element_text(size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(15, 15, 25, 15)
  )

print(p_multi)

ggsave(file.path(out_dir, "Figure2AA_KeyGene_Dependency_Multipanel.pdf"),
       p_multi, width = 14, height = 8)

ggsave(file.path(out_dir, "Figure2AA_KeyGene_Dependency_Multipanel.jpg"),
       p_multi, width = 14, height = 8, dpi = 300)

message("Saved multipanel Figure 2A (PDF + JPG)")

# ============================================================
# Fig 2B: Top state-differential dependencies (heatmap) 
# ============================================================

# Gatekeeper genes that exist in dat
dep_genes_all <- gatekeepers[gatekeepers %in% colnames(dat)]

if (length(dep_genes_all) >= 10) {
  
  # Make sure dat is a plain data.frame
  dat_df <- as.data.frame(dat)
  
  # Long format: ModelID / State / gene / dep
  df_long <- dat_df %>%
    dplyr::select(ModelID, State, dplyr::all_of(dep_genes_all)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(dep_genes_all),
      names_to = "gene",
      values_to = "dep"
    ) %>%
    dplyr::filter(!is.na(State), !is.na(dep))
  
  # ANOVA per gene (dep ~ State)
  ptab <- df_long %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      p = tryCatch(stats::anova(stats::lm(dep ~ State))$`Pr(>F)`[1],
                   error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj = p.adjust(p, method = "BH")) %>%
    dplyr::arrange(p_adj)
  
  # Save stats
  write.csv(ptab, file.path(out_dir, "Fig2_TopDifferential_DependencyStats.csv"), row.names = FALSE)
  
  # Top N genes (edit N if needed)
  topN <- ptab %>%
    dplyr::filter(!is.na(p_adj)) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::pull(gene)
  
  # Mean dependency by state (matrix)
  mat <- df_long %>%
    dplyr::filter(gene %in% topN) %>%
    dplyr::group_by(gene, State) %>%
    dplyr::summarise(mean_dep = mean(dep, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = State, values_from = mean_dep) %>%
    as.data.frame()
  
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  
  # Row z-score
  zscore_vec <- function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (is.na(sx) || sx == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sx
  }
  mat_z <- t(apply(as.matrix(mat), 1, zscore_vec))
  
  # Heatmap df for ggplot
  heat_df <- as.data.frame(mat_z)
  heat_df$gene <- rownames(heat_df)
  
  heat_df <- heat_df %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "State",
      values_to = "z"
    )
  
  # Keep your state order (prevents random column order)
  state_levels <- c(
    "Mito-High / Glyco-Low (OXPHOS-dominant)",
    "Mito-Low / Glyco-High (Glycolytic)",
    "Both-High (Hybrid)",
    "Both-Low (Low metabolic)"
  )
  heat_df$State <- factor(heat_df$State, levels = state_levels)
  
  # Plot (blue/white/red like your ComplexHeatmap palette)
  # Plot (blue/white/red like your ComplexHeatmap palette)
  pD <- ggplot2::ggplot(heat_df, ggplot2::aes(x = State, y = gene, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#2b6cb0",
      mid = "white",
      high = "#c53030",
      midpoint = 0,
      limits = c(-2, 2),
      name = "Row-wise\nZ-score"
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        "Mito-High / Glyco-Low (OXPHOS-dominant)" =
          "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
        "Mito-Low / Glyco-High (Glycolytic)" =
          "Mito-low\nGlyco-high\n(Glycolytic)",
        "Both-High (Hybrid)" =
          "Both-high\n(Hybrid)",
        "Both-Low (Low metabolic)" =
          "Both-low\n(Low)"
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "Top State-Differential Dependencies of Mitochondrial-Associated Genes",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = 9,
        lineheight = 0.9
      ),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin = ggplot2::margin(12, 20, 55, 20)
    )
  print(pD)
  
  ggsave(file.path(out_dir, "Figure2BB_TopDifferential_Dependency_Heatmap.pdf"),
         pD, width = 10, height = 9)
  
  ggsave(file.path(out_dir, "Figure2BB_TopDifferential_Dependency_Heatmap.png"),
         pD, width = 10, height = 9, dpi = 300)
  
  message("Fig2B saved + stats CSV saved in: ", out_dir)
  
} else {
  message("Not enough gatekeeper genes found in dat to run Fig2B.")
}

# ============================================================
# Lineage-adjusted State effect (ALL genes → ONE file) =====not needed
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

genes_all <- gatekeepers[gatekeepers %in% colnames(dat)]

if (length(genes_all) == 0) {
  stop("No gatekeeper genes found in dat.")
}

if (is.na(lineage_col) || length(lineage_col) == 0) {
  stop("lineage_col not defined.")
}

dat_df <- as.data.frame(dat)
dat_df$Lineage <- dat_df[[lineage_col]]

results <- lapply(genes_all, function(g) {
  
  df <- dat_df %>%
    dplyr::select(State, Lineage, all_of(g)) %>%
    dplyr::rename(dep = all_of(g)) %>%
    dplyr::filter(!is.na(State), !is.na(Lineage), !is.na(dep))
  
  if (nrow(df) < 50) {
    return(data.frame(
      Gene = g,
      N = nrow(df),
      P_state = NA_real_
    ))
  }
  
  df$State   <- as.factor(df$State)
  df$Lineage <- as.factor(df$Lineage)
  
  fit <- tryCatch(lm(dep ~ State + Lineage, data = df),
                  error = function(e) NULL)
  
  if (is.null(fit)) {
    return(data.frame(
      Gene = g,
      N = nrow(df),
      P_state = NA_real_
    ))
  }
  
  a <- anova(fit)
  p_state <- a$`Pr(>F)`[match("State", rownames(a))]
  
  data.frame(
    Gene = g,
    N = nrow(df),
    P_state = p_state
  )
})

res <- bind_rows(results)

# BH correction across ALL genes
res$FDR_state_BH <- p.adjust(res$P_state, method = "BH")

# Order by FDR
res <- res[order(res$FDR_state_BH), ]

# Save ONE file
write.csv(res,
          file.path(out_dir, "Fig2E_LineageAdjusted_StateEffect_ALL_Genes.csv"),
          row.names = FALSE)

message("Saved ONE file with lineage-adjusted state effects for all genes.")



library(ggplot2)

res <- read.csv(file.path(out_dir,
                          "Fig2E_LineageAdjusted_StateEffect_ALL_Genes.csv"))

res$logFDR <- -log10(res$FDR_state_BH)

p_volcano <- ggplot(res, aes(x = P_state, y = logFDR)) +
  geom_point(alpha = 0.7) +
  theme_classic(base_size = 12) +
  labs(
    title = "Lineage-adjusted State Effect on Dependency",
    x = "Raw p-value (State term)",
    y = "-log10(FDR)"
  )

print(p_volcano)

ggsave(file.path(out_dir, "Fig2E_LineageAdjusted_Volcano.pdf"),
       p_volcano, width = 7, height = 5)

ggsave(file.path(out_dir, "Fig2E_LineageAdjusted_Volcano.png"),
       p_volcano, width = 7, height = 5, dpi = 300)


# ============================================================
# Figure 2C: TRUE volcano (Effect size vs -log10(FDR))
# Lineage-adjusted: dep ~ State + Lineage
# Effect size = (max mean across states) - (min mean across states)
# ============================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

dep_genes_all <- gatekeepers[gatekeepers %in% colnames(dat)]
stopifnot(length(dep_genes_all) >= 5)

df0 <- as.data.frame(dat) %>%
  filter(!is.na(State)) %>%
  mutate(
    State = factor(State, levels = state_levels),
    Lineage = .data[[lineage_col]]
  ) %>%
  filter(!is.na(Lineage))

df_long <- df0 %>%
  dplyr::select(ModelID, State, Lineage, dplyr::all_of(dep_genes_all)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(dep_genes_all),
    names_to = "Gene",
    values_to = "dep"
  ) %>%
  filter(!is.na(dep))

res <- df_long %>%
  group_by(Gene) %>%
  summarise(
    N = n(),
    p_state = tryCatch({
      fit <- lm(dep ~ State + Lineage)
      a <- anova(fit)
      a$`Pr(>F)`[which(rownames(a) == "State")]
    }, error = function(e) NA_real_),
    effect_range = {
      m <- tapply(dep, State, mean, na.rm = TRUE)
      as.numeric(max(m, na.rm = TRUE) - min(m, na.rm = TRUE))
    },
    .groups = "drop"
  ) %>%
  mutate(
    FDR_state_BH = p.adjust(p_state, method = "BH"),
    neglog10FDR  = -log10(FDR_state_BH),
    Significance = ifelse(FDR_state_BH < 0.05,
                          "Significant (FDR < 0.05)",
                          "Not significant")
  ) %>%
  arrange(FDR_state_BH)

write.csv(res,
          file.path(out_dir, "Figure2C_LineageAdjusted_VolcanoStats.csv"),
          row.names = FALSE)

p_volcano <- ggplot(res, aes(x = effect_range, y = neglog10FDR)) +
  geom_point(aes(color = Significance), size = 2.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    aes(label = Gene, color = Significance),
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.15,
    segment.color = "gray50",
    min.segment.length = 0,
    force = 1,
    force_pull = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Significant (FDR < 0.05)" = "orangered",
    "Not significant" = "gray50"
  )) +
  theme(
    plot.margin = margin(15, 30, 20, 15),  # increase RIGHT margin
    axis.title.x = element_text(margin = margin(t = 10))
    ) +
  labs(
    title = "Lineage-Adjusted State-Specific Dependency Effects of Mitochondrial-Associated Genes",
    subtitle = "Effect size = range of state mean CERES (max − min); dashed line = FDR 0.05",
    x = "Effect size (range of mean CERES dependency across metabolic states)",
    y = "-log10(FDR, BH)",
    color = NULL
  ) +
  theme(legend.position = "top")



ggsave(file.path(out_dir, "Figure2C_LineageAdjusted_Volcano.pdf"),
       p_volcano, width = 9, height = 6)

ggsave(file.path(out_dir, "Figure2C_LineageAdjusted_Volcano.png"),
       p_volcano, width = 9, height = 6, dpi = 300)
# ============================================================
# Figure 2D: Ranked state-specific dependency effect size
# ============================================================

res_rank <- res %>%
  filter(!is.na(effect_range), !is.na(FDR_state_BH)) %>%
  arrange(effect_range) %>%
  mutate(
    Gene = factor(Gene, levels = Gene),
    Significant = ifelse(FDR_state_BH < 0.05, "FDR < 0.05", "Not significant")
  )

p_rank <- ggplot(res_rank,
                 aes(x = Gene,
                     y = effect_range,
                     fill = Significant)) +
  geom_col(width = 0.75) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(values = c(
    "FDR < 0.05" = "orangered",
    "Not significant" = "gray50"
  )) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 35, 20, 10)
  ) +
  labs(
    title = "Ranked State-specific Dependency Effects",
    x = NULL,
    y = "Effect size (range of mean CERES dependency across metabolic states)",
    fill = NULL
  )
print(p_rank)

ggsave(file.path(out_dir, "Figure2D_Ranked_Dependency_EffectSize.pdf"),
       p_rank, width = 8.8, height = 9)

ggsave(file.path(out_dir, "Figure2D_Ranked_Dependency_EffectSize.png"),
       p_rank, width = 8.8, height = 9, dpi = 300)

# ============================================================
# FIGURE 2E: Global mitochondrial gatekeeper dependency burden
# Across metabolic states
# ============================================================

# ---- State order ----
state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

# ---- Genes present in DepMap data ----
dep_genes <- gatekeepers[gatekeepers %in% colnames(dat)]

# ---- Calculate mean mitochondrial dependency per cell line ----
dat_2E <- as.data.frame(dat) %>%
  mutate(State = factor(State, levels = state_levels)) %>%
  filter(!is.na(State)) %>%
  rowwise() %>%
  mutate(
    Global_Mito_Dependency = mean(c_across(all_of(dep_genes)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(!is.na(Global_Mito_Dependency))

# ---- Overall Kruskal-Wallis test ----
kw_p <- kruskal.test(Global_Mito_Dependency ~ State, data = dat_2E)$p.value

# ---- Colors ----
state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)

# ---- Plot Figure 2E ----
p_2E <- ggplot(dat_2E,
               aes(x = State,
                   y = Global_Mito_Dependency,
                   fill = State)) +
  geom_violin(trim = TRUE, alpha = 0.75, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.16,
               outlier.size = 0.4,
               alpha = 0.9,
               color = "black") +
  geom_jitter(width = 0.12,
              size = 0.45,
              alpha = 0.25,
              color = "black") +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(
    labels = c(
      "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
      "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
      "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
      "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
    )
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Global Mitochondrial-Associated Gene Dependency Across Metabolic States",
    subtitle = paste0("Kruskal–Wallis p = ", signif(kw_p, 3)),
    x = NULL,
    y = "Mean CERES gene effect\n(more negative = stronger dependency)"
  ) +
  theme(
    axis.text.x = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

print(p_2E)

# ---- Save ----
ggsave(file.path(out_dir, "Figure2E_Global_Mitochondrial_Dependency.pdf"),
       p_2E, width = 8, height = 5)

ggsave(file.path(out_dir, "Figure2E_Global_Mitochondrial_Dependency.png"),
       p_2E, width = 8, height = 5, dpi = 300)

# ---- Save source data ----
write.csv(dat_2E,
          file.path(out_dir, "Figure2E_Global_Mitochondrial_Dependency_Data.csv"),
          row.names = FALSE)

message("Figure 2E saved successfully in: ", out_dir)
