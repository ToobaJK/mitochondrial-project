# =========================================================
# FIGURE 5: Pan-cancer miRNA landscape across metabolic states
# Panels A-F
# =========================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ggpubr)
library(janitor)
library(RColorBrewer)

# -----------------------------
# 1. Paths
# -----------------------------
setwd("D:/UAEU/Dr. Ajaz/New paper_Dr. Ajaz/Fig 7")

output_dir <- "outputR3"
if (!dir.exists(output_dir)) dir.create(output_dir)

mirna_file <- "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena"
state_file <- "Figure1E_TCGA_MetabolicStates_Table.tsv"
pheno_file <- "TCGA_phenotype_denseDataOnlyDownload.tsv"
mrna_file <- "tcga_RSEM_gene_tpm"

tcga_expr <- fread(mrna_file)

# Ensure header is correct
if (!"sample" %in% colnames(tcga_expr)[1]) {
  cat("Warning: header might be wrong\n")
}

# Force column names
setnames(tcga_expr, make.names(colnames(tcga_expr)))



# -----------------------------
# 2. Helper functions
# -----------------------------
clean_tcga_id <- function(x) {
  x <- gsub("\\.", "-", x)
  substr(x, 1, 15)
}

save_plot <- function(plot, name, width = 8, height = 6) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")),
         plot, width = width, height = height)
  ggsave(file.path(output_dir, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 300)
}

# -----------------------------
# 3. Load miRNA expression
# -----------------------------
# Reload miRNA matrix from original raw object/file
mirna_raw <- read.delim(
  mirna_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(mirna_raw)[1] <- "miRNA"

mirna_mat <- mirna_raw
rownames(mirna_mat) <- mirna_mat$miRNA
mirna_mat$miRNA <- NULL

colnames(mirna_mat) <- gsub("\\.", "-", colnames(mirna_mat))
colnames(mirna_mat) <- substr(colnames(mirna_mat), 1, 15)

mirna_mat <- as.data.frame(
  lapply(mirna_mat, function(x) as.numeric(as.character(x)))
)

rownames(mirna_mat) <- mirna_raw$miRNA

mirna_mat <- mirna_mat[, !duplicated(colnames(mirna_mat))]

dim(mirna_mat)
head(colnames(mirna_mat), 10)

cat("miRNA matrix:\n")
print(dim(mirna_mat))
print(colnames(mirna_mat)[1:5])
print(rownames(mirna_mat)[1:5])

# -----------------------------
# 4. Load phenotype file
# -----------------------------
pheno <- read.delim(
  pheno_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pheno$sample <- clean_tcga_id(pheno$sample)

# -----------------------------
# 5. Disease to TCGA code mapping
# -----------------------------
disease_to_code <- c(
  "adrenocortical cancer" = "ACC",
  "bladder urothelial carcinoma" = "BLCA",
  "breast invasive carcinoma" = "BRCA",
  "cervical & endocervical cancer" = "CESC",
  "cholangiocarcinoma" = "CHOL",
  "colon adenocarcinoma" = "COAD",
  "diffuse large B-cell lymphoma" = "DLBC",
  "esophageal carcinoma" = "ESCA",
  "glioblastoma multiforme" = "GBM",
  "head & neck squamous cell carcinoma" = "HNSC",
  "kidney chromophobe" = "KICH",
  "kidney clear cell carcinoma" = "KIRC",
  "kidney papillary cell carcinoma" = "KIRP",
  "acute myeloid leukemia" = "LAML",
  "brain lower grade glioma" = "LGG",
  "liver hepatocellular carcinoma" = "LIHC",
  "lung adenocarcinoma" = "LUAD",
  "lung squamous cell carcinoma" = "LUSC",
  "mesothelioma" = "MESO",
  "ovarian serous cystadenocarcinoma" = "OV",
  "pancreatic adenocarcinoma" = "PAAD",
  "pheochromocytoma & paraganglioma" = "PCPG",
  "prostate adenocarcinoma" = "PRAD",
  "rectum adenocarcinoma" = "READ",
  "sarcoma" = "SARC",
  "skin cutaneous melanoma" = "SKCM",
  "stomach adenocarcinoma" = "STAD",
  "testicular germ cell tumor" = "TGCT",
  "thyroid carcinoma" = "THCA",
  "thymoma" = "THYM",
  "uterine corpus endometrioid carcinoma" = "UCEC",
  "uterine carcinosarcoma" = "UCS",
  "uveal melanoma" = "UVM"
)
# -----------------------------
# Build phenotype map correctly
# -----------------------------
pheno_map <- pheno %>%
  mutate(
    sample = clean_tcga_id(sample),
    cancer_code = disease_to_code[`_primary_disease`]
  ) %>%
  select(sample, cancer_code, sample_type, `_primary_disease`) %>%
  filter(!is.na(cancer_code)) %>%
  distinct(sample, .keep_all = TRUE)

dim(pheno_map)
head(pheno_map)
table(pheno_map$cancer_code)

# -----------------------------
# 6. Load metabolic state file
# -----------------------------
tcga_states <- read.delim(
  state_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tcga_states <- tcga_states %>%
  mutate(
    cancer_code = toupper(cancer_code),
    metabolic_state = state
  )


# Clean miRNA sample names
colnames(mirna_mat) <- gsub("\\.", "-", colnames(mirna_mat))
colnames(mirna_mat) <- substr(colnames(mirna_mat), 1, 15)

# Clean phenotype sample names
pheno_map$sample <- gsub("\\.", "-", pheno_map$sample)
pheno_map$sample <- substr(pheno_map$sample, 1, 15)

# Check overlap
length(intersect(colnames(mirna_mat), pheno_map$sample))
head(intersect(colnames(mirna_mat), pheno_map$sample))
# -----------------------------
# 7. Merge miRNA sample info with states
# -----------------------------

# clean again
colnames(mirna_mat) <- gsub("\\.", "-", colnames(mirna_mat))
colnames(mirna_mat) <- substr(colnames(mirna_mat), 1, 15)

pheno_map$sample <- gsub("\\.", "-", pheno_map$sample)
pheno_map$sample <- substr(pheno_map$sample, 1, 15)

# remove duplicate miRNA columns
mirna_mat <- mirna_mat[, !duplicated(colnames(mirna_mat))]

sample_info <- data.frame(
  sample = colnames(mirna_mat),
  stringsAsFactors = FALSE
) %>%
  left_join(pheno_map, by = "sample") %>%
  left_join(
    tcga_states %>%
      select(cancer_code, metabolic_state, mito_gatekeeper_score, glycolysis_score),
    by = "cancer_code"
  ) %>%
  filter(!is.na(metabolic_state)) %>%
  distinct(sample, .keep_all = TRUE)

cat("Matched miRNA samples:\n")
print(nrow(sample_info))
print(table(sample_info$metabolic_state))

# keep only samples that really exist in miRNA matrix
common_samples <- intersect(colnames(mirna_mat), sample_info$sample)

mirna_mat <- mirna_mat[, common_samples, drop = FALSE]

sample_info <- sample_info %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, colnames(mirna_mat)))

# check
print(dim(mirna_mat))
print(nrow(sample_info))
print(head(sample_info$sample))
print(head(colnames(mirna_mat)))

stopifnot(identical(sample_info$sample, colnames(mirna_mat)))
# -----------------------------
# 8. Long format
# -----------------------------
mirna_long <- as.data.frame(mirna_mat)
mirna_long$miRNA <- rownames(mirna_mat)

mirna_long <- mirna_long %>%
  pivot_longer(
    cols = -miRNA,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  left_join(sample_info, by = "sample") %>%
  filter(!is.na(expression), !is.na(metabolic_state))

mirna_long$metabolic_state <- factor(
  mirna_long$metabolic_state,
  levels = unique(sample_info$metabolic_state)
)

mirna_long$metabolic_state <- as.character(mirna_long$metabolic_state)

mirna_long$metabolic_state <- dplyr::case_when(
  grepl("OXPHOS", mirna_long$metabolic_state) ~ "Mito-high Glyco-low (OXPHOS-dominant)",
  grepl("Glyco-High", mirna_long$metabolic_state) ~ "Mito-low Glyco-high (Glycolytic)",
  grepl("Both-High", mirna_long$metabolic_state) ~ "Both-high (Hybrid)",
  grepl("Both-Low", mirna_long$metabolic_state) ~ "Both-low (Low metabolic)",
  TRUE ~ NA_character_
)

dim(mirna_long)
table(mirna_long$metabolic_state)
sum(is.na(mirna_long$metabolic_state))

state_order <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)",
  "Mito-low Glyco-high (Glycolytic)",
  "Both-high (Hybrid)",
  "Both-low (Low metabolic)"
)

state_colors <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)" = "darkgreen",
  "Mito-low Glyco-high (Glycolytic)"      = "orange",
  "Both-high (Hybrid)"                    = "purple",
  "Both-low (Low metabolic)"              = "maroon"
)
# =========================================================
# FIGURE 5A: PCA
# =========================================================

sample_info %>%
  count(metabolic_state)

state_order <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)

# Remove zero-variance miRNAs
mirna_var <- apply(mirna_mat, 1, var, na.rm = TRUE)
mirna_mat_pca <- mirna_mat[mirna_var > 0, ]

# Top 500 variable miRNAs
mirna_var <- apply(mirna_mat_pca, 1, var, na.rm = TRUE)
top_var_mirna <- names(sort(mirna_var, decreasing = TRUE))[1:500]

pca_mat <- t(mirna_mat_pca[top_var_mirna, ])
pca_mat[is.na(pca_mat)] <- 0

pca_res <- prcomp(pca_mat, scale. = TRUE)

# Check original labels from sample_info
unique(sample_info$metabolic_state)

# Rebuild pca_df fresh
pca_df <- data.frame(
  sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
)

pca_df$sample <- gsub("\\.", "-", pca_df$sample)
pca_df$sample <- substr(pca_df$sample, 1, 15)

sample_info$sample <- gsub("\\.", "-", sample_info$sample)
sample_info$sample <- substr(sample_info$sample, 1, 15)

pca_df <- pca_df %>%
  left_join(sample_info[, c("sample", "metabolic_state")], by = "sample")

# Rename any OXPHOS label variation
pca_df$metabolic_state <- as.character(pca_df$metabolic_state)
pca_df$metabolic_state[grepl("OXPHOS", pca_df$metabolic_state)] <- 
  "Mito-High / Glyco-Low (OXPHOS-dominant)"

# Now factor again
pca_df$metabolic_state <- factor(pca_df$metabolic_state, levels = state_order)

table(pca_df$metabolic_state, useNA = "ifany")
sum(is.na(pca_df$metabolic_state))

p5A <- ggplot(pca_df, aes(PC1, PC2, color = metabolic_state)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = state_colors, breaks = state_order) +
  theme_bw(base_size = 12) +
  labs(
    title = "Pan-cancer miRNA Expression Landscape",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"),
    color = "Metabolic state"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p5A)
save_plot(p5A, "Figure5A_miRNA_PCA", 7, 6)


# =========================================================
# FIGURE 5B: Differential miRNA volcano plot
# OXPHOS-dominant vs Glycolytic
# =========================================================
ox_state <- "Mito-high Glyco-low (OXPHOS-dominant)"
gly_state <- "Mito-low Glyco-high (Glycolytic)"


diff_mirna <- mirna_long %>%
  filter(metabolic_state %in% c(ox_state, gly_state)) %>%
  group_by(miRNA) %>%
  summarise(
    mean_oxphos = mean(expression[metabolic_state == ox_state], na.rm = TRUE),
    mean_glyco  = mean(expression[metabolic_state == gly_state], na.rm = TRUE),
    effect = mean_oxphos - mean_glyco,
    p_value = tryCatch(
      wilcox.test(
        expression[metabolic_state == ox_state],
        expression[metabolic_state == gly_state]
      )$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    neglog10 = -log10(p_adj),
    regulation = case_when(
      p_adj < 0.05 & effect > 0 ~ "Mito-High / Glyco-Low (OXPHOS-dominant)",
      p_adj < 0.05 & effect < 0 ~ "Mito-Low / Glyco-High (Glycolytic)",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  diff_mirna,
  file.path(output_dir, "Figure5Ba_Differential_miRNA_OXPHOS_vs_Glycolytic.csv"),
  row.names = FALSE
)

cat("\nMetabolic state counts:\n\n")
print(table(diff_mirna$regulation))

# Cap extreme values (recommended)
diff_mirna$neglog10 <- pmin(diff_mirna$neglog10, 200)

top_labels <- diff_mirna %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  slice_head(n = 15)

p5B <- ggplot(diff_mirna, aes(x = effect, y = neglog10, color = regulation)) +
  geom_point(alpha = 0.75, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_labels,
    aes(label = miRNA),
    size = 3,
    max.overlaps = 50
  ) +
  scale_color_manual(
    values = c(
      "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
      "Mito-Low / Glyco-High (Glycolytic)" = "orange",
      "Not significant" = "grey50"
    )
  ) +
  coord_cartesian(ylim = c(0, 200)) +  # cleaner y-axis
  theme_bw(base_size = 12) +
  labs(
    title = "Metabolic State-Associated miRNAs",
    x = "Mean difference: OXPHOS-dominant − Glycolytic",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )



ggsave(
  file.path(output_dir, "Figure5Ba_miRNA_Volcano_OXPHOS_vs_Glycolytic.pdf"),
  plot = p5B,
  width = 10,
  height = 7
)

ggsave(
  file.path(output_dir, "Figure5Ba_miRNA_Volcano_OXPHOS_vs_Glycolytic.png"),
  plot = p5B,
  width = 10,
  height = 7,
  dpi = 300
)


# =========================================================
# FIGURE 5C: Pan-cancer miRNA expression heatmap
# Row-wise Z-score across metabolic states
# Saves PDF + JPG
# =========================================================

library(dplyr)
library(ggplot2)
library(scales)

# -----------------------------
# 1. State order, labels, colors
# -----------------------------
state_order <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)",
  "Mito-low Glyco-high (Glycolytic)",
  "Both-high (Hybrid)",
  "Both-low (Low metabolic)"
)

state_labels <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-low Glyco-high (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-high (Hybrid)"                    = "Both-high\n(Hybrid)",
  "Both-low (Low metabolic)"              = "Both-low\n(Low)"
)

# -----------------------------
# 2. Select top variable miRNAs
# -----------------------------
top_mirnas <- names(
  sort(apply(mirna_mat, 1, var, na.rm = TRUE), decreasing = TRUE)
)[1:40]

# -----------------------------
# 3. Mean miRNA expression per metabolic state
# -----------------------------
top_heat <- mirna_long %>%
  filter(miRNA %in% top_mirnas) %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_order)
  ) %>%
  group_by(miRNA, metabolic_state) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 4. Row-wise Z-score
# -----------------------------
top_heat_z <- top_heat %>%
  group_by(miRNA) %>%
  mutate(
    z_expr = as.numeric(scale(mean_expr)),
    mean_global = mean(mean_expr, na.rm = TRUE)
  ) %>%
  ungroup()

# -----------------------------
# 5. Plot heatmap
# -----------------------------
p5C <- ggplot(
  top_heat_z,
  aes(
    x = metabolic_state,
    y = reorder(miRNA, mean_global),
    fill = z_expr
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Pan-cancer miRNA Expression Across Metabolic States",
    x = NULL,
    y = NULL,
    fill = "Row-wise\nZ-score"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

print(p5C)

# -----------------------------
# 6. Save PDF + JPG
# -----------------------------
ggsave(
  file.path(output_dir, "Figure5Ca_PanCancer_miRNA_Heatmap_RowZscore.pdf"),
  plot = p5C,
  width = 8,
  height = 10
)

ggsave(
  file.path(output_dir, "Figure5Ca_PanCancer_miRNA_Heatmap_RowZscore.jpg"),
  plot = p5C,
  width = 8,
  height = 10,
  dpi = 300
)


# =========================================================
# FIGURE 5D: Pan-cancer miRNA heatmap annotated by metabolic state
# miRNA x cancer type + metabolic-state annotation
# =========================================================

library(dplyr)
library(ggplot2)
library(scales)

# -----------------------------
# 1. Select miRNAs
# Option A: top variable miRNAs
# -----------------------------
top_mirnas <- names(
  sort(apply(mirna_mat, 1, var, na.rm = TRUE), decreasing = TRUE)
)[1:40]

# -----------------------------
# 2. Mean miRNA expression per cancer type
# -----------------------------
heat_cancer <- mirna_long %>%
  filter(miRNA %in% top_mirnas) %>%
  group_by(miRNA, cancer_code, metabolic_state) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 3. Row-wise Z-score across cancer types
# -----------------------------
heat_cancer <- heat_cancer %>%
  group_by(miRNA) %>%
  mutate(z_expr = as.numeric(scale(mean_expr))) %>%
  ungroup()

# -----------------------------
# 4. Order cancer types by metabolic state
# -----------------------------
state_order <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)",
  "Mito-low Glyco-high (Glycolytic)",
  "Both-high (Hybrid)",
  "Both-low (Low metabolic)"
)

state_labels <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-low Glyco-high (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-high (Hybrid)"                    = "Both-high\n(Hybrid)",
  "Both-low (Low metabolic)"              = "Both-low\n(Low)"
)

state_colors <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)" = "darkgreen",
  "Mito-low Glyco-high (Glycolytic)"      = "orange",
  "Both-high (Hybrid)"                    = "purple",
  "Both-low (Low metabolic)"              = "maroon"
)

cancer_order <- heat_cancer %>%
  distinct(cancer_code, metabolic_state) %>%
  mutate(metabolic_state = factor(metabolic_state, levels = state_order)) %>%
  arrange(metabolic_state, cancer_code) %>%
  pull(cancer_code)

heat_cancer$cancer_code <- factor(heat_cancer$cancer_code, levels = cancer_order)

# -----------------------------
# 5. Main heatmap
# -----------------------------
p5D_main <- ggplot(
  heat_cancer,
  aes(
    x = cancer_code,
    y = miRNA,
    fill = z_expr
  )
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Pan-cancer miRNA Expression Landscape by Metabolic State",
    x = NULL,
    y = NULL,
    fill = "Row Z-score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p5D_main)

ggsave(
  file.path(output_dir, "Figure5Da_PanCancer_miRNA_ByCancer_MetabolicState.pdf"),
  plot = p5D_main,
  width = 12,
  height = 10
)

ggsave(
  file.path(output_dir, "Figure5Da_PanCancer_miRNA_ByCancer_MetabolicState.jpg"),
  plot = p5D_main,
  width = 12,
  height = 10,
  dpi = 300
)



# =========================================================
# FIGURE 5E: Integrated pan-cancer miRNA heatmap
# X-axis = TCGA cancers grouped by metabolic state
# Top bar = metabolic state
# Color = row-wise Z-score
# Saves PDF + JPG
# =========================================================

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

# -----------------------------
# 1. State order and colors
# -----------------------------
state_order <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)",
  "Mito-low Glyco-high (Glycolytic)",
  "Both-high (Hybrid)",
  "Both-low (Low metabolic)"
)

state_colors <- c(
  "Mito-high Glyco-low (OXPHOS-dominant)" = "darkgreen",
  "Mito-low Glyco-high (Glycolytic)"      = "orange",
  "Both-high (Hybrid)"                    = "purple",
  "Both-low (Low metabolic)"              = "maroon"
)

# -----------------------------
# 2. Select miRNAs
# -----------------------------
selected_mirnas <- unique(top_mirnas)

# -----------------------------
# 3. Mean miRNA expression per cancer type
# -----------------------------
heat_cancer <- mirna_long %>%
  filter(miRNA %in% selected_mirnas) %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_order)
  ) %>%
  group_by(miRNA, cancer_code, metabolic_state) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 4. Row-wise Z-score across cancers
# -----------------------------
heat_cancer <- heat_cancer %>%
  group_by(miRNA) %>%
  mutate(
    z_expr = as.numeric(scale(mean_expr))
  ) %>%
  ungroup()

# -----------------------------
# 5. Order cancer types by metabolic state
# -----------------------------
cancer_order <- heat_cancer %>%
  distinct(cancer_code, metabolic_state) %>%
  arrange(metabolic_state, cancer_code) %>%
  pull(cancer_code)

heat_cancer$cancer_code <- factor(
  heat_cancer$cancer_code,
  levels = cancer_order
)

# -----------------------------
# 6. Annotation bar data
# -----------------------------
ann_df <- heat_cancer %>%
  distinct(cancer_code, metabolic_state) %>%
  mutate(
    cancer_code = factor(cancer_code, levels = cancer_order),
    annotation = "Metabolic\nstate"
  )

# -----------------------------
# 7. Top metabolic-state annotation bar
# -----------------------------
p_top <- ggplot(
  ann_df,
  aes(x = cancer_code, y = annotation, fill = metabolic_state)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = state_colors,
    name = "Metabolic state"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 0, r = 5, b = -8, l = 5)
  )

# -----------------------------
# 8. Main heatmap
# -----------------------------
p_main <- ggplot(
  heat_cancer,
  aes(x = cancer_code, y = miRNA, fill = z_expr)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Row-wise\nZ-score"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Integrated miRNA Landscape Across TCGA Cancers and Metabolic States",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.4, "cm"),
    plot.margin = margin(t = -5, r = 5, b = 5, l = 5)
  )

# -----------------------------
# 9. Combine plots and collect legends
# -----------------------------
p5E <- (p_top / p_main) +
  plot_layout(
    heights = c(0.45, 7),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 5)
  )

print(p5E)

# -----------------------------
# 10. Save with smaller width
# -----------------------------
ggsave(
  file.path(output_dir, "Figure5E_miRNA_TCGA_Cancers_Grouped_by_MetabolicState.pdf"),
  plot = p5E,
  width = 11.5,
  height = 10
)

ggsave(
  file.path(output_dir, "Figure5E_miRNA_TCGA_Cancers_Grouped_by_MetabolicState.jpg"),
  plot = p5E,
  width = 11.5,
  height = 10,
  dpi = 300
)


# =========================================================
# FIGURE 5F: miRNA–metabolic gene correlation heatmap
# TCGA miRNA vs TCGA mRNA expression
# =========================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(AnnotationDbi)
library(org.Hs.eg.db)


# =========================================================
# Load TCGA mRNA expression exactly like Figure 1
# =========================================================

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

base_dir <- "D:/UAEU/Dr. Ajaz/New paper_Dr. Ajaz/Pancancer data"
mrna_file <- file.path(base_dir, "tcga_RSEM_gene_tpm")

tcga_expr <- fread(mrna_file)

gene_col <- colnames(tcga_expr)[1]

tcga_mat <- as.matrix(tcga_expr[, -1, with = FALSE])
rownames(tcga_mat) <- tcga_expr[[gene_col]]

tcga_mat <- log2(tcga_mat + 0.001)

# Same as Figure 1
colnames(tcga_mat) <- substr(colnames(tcga_mat), 1, 15)

tcga_mat <- tcga_mat[, !duplicated(colnames(tcga_mat)), drop = FALSE]

cat("mRNA matrix loaded:\n")
print(dim(tcga_mat))
print(head(colnames(tcga_mat)))
print(head(rownames(tcga_mat)))

# -----------------------------
# 2. Convert Ensembl IDs to gene symbols
# -----------------------------

ens <- rownames(tcga_mat)
ens_clean <- sub("\\..*$", "", ens)

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

keep <- !is.na(symbol_map) & symbol_map != ""

tcga_mat2 <- tcga_mat[keep, , drop = FALSE]
rownames(tcga_mat2) <- symbol_map[keep]

# Collapse duplicate gene symbols by mean
tcga_dt2 <- as.data.table(tcga_mat2, keep.rownames = "gene")
tcga_dt2 <- tcga_dt2[, lapply(.SD, mean, na.rm = TRUE), by = gene]

tcga_mat_symbol <- as.matrix(tcga_dt2[, -1, with = FALSE])
rownames(tcga_mat_symbol) <- tcga_dt2$gene

cat("mRNA converted to gene symbols:\n")
print(dim(tcga_mat_symbol))
print(head(rownames(tcga_mat_symbol)))


# -----------------------------
# 3. Clean miRNA sample IDs
# -----------------------------

mirna_mat_corr <- mirna_mat

colnames(mirna_mat_corr) <- gsub("\\.", "-", colnames(mirna_mat_corr))
colnames(mirna_mat_corr) <- substr(colnames(mirna_mat_corr), 1, 15)

mirna_mat_corr <- mirna_mat_corr[, !duplicated(colnames(mirna_mat_corr)), drop = FALSE]

cat("miRNA matrix for correlation:\n")
print(dim(mirna_mat_corr))
print(head(colnames(mirna_mat_corr)))


# -----------------------------
# 4. Match miRNA and mRNA samples
# -----------------------------

# Re-create miRNA correlation matrix safely
mirna_mat_corr <- mirna_mat

# Make sure it is a matrix/data frame with colnames
mirna_mat_corr <- as.data.frame(mirna_mat_corr)

colnames(mirna_mat_corr) <- as.character(colnames(mirna_mat_corr))
colnames(mirna_mat_corr) <- gsub("\\.", "-", colnames(mirna_mat_corr))
colnames(mirna_mat_corr) <- substr(colnames(mirna_mat_corr), 1, 15)

mirna_mat_corr <- mirna_mat_corr[, !duplicated(colnames(mirna_mat_corr)), drop = FALSE]

# Clean mRNA symbol matrix column names safely
colnames(tcga_mat_symbol) <- as.character(colnames(tcga_mat_symbol))
colnames(tcga_mat_symbol) <- substr(colnames(tcga_mat_symbol), 1, 15)

tcga_mat_symbol <- tcga_mat_symbol[, !duplicated(colnames(tcga_mat_symbol)), drop = FALSE]

# Check both
cat("mRNA samples:\n")
print(head(colnames(tcga_mat_symbol)))

cat("miRNA samples:\n")
print(head(colnames(mirna_mat_corr)))

common_samples_corr <- intersect(
  colnames(tcga_mat_symbol),
  colnames(mirna_mat_corr)
)

cat("Matched miRNA-mRNA samples:", length(common_samples_corr), "\n")

tcga_mat_symbol <- tcga_mat_symbol[, common_samples_corr, drop = FALSE]
mirna_mat_corr  <- mirna_mat_corr[, common_samples_corr, drop = FALSE]

stopifnot(identical(colnames(tcga_mat_symbol), colnames(mirna_mat_corr)))

# -----------------------------
# 5. Select metabolic genes
# -----------------------------

metabolic_genes <- c(
  "GPX4", "SLC7A11", "VDAC1", "VDAC2", "VDAC3",
  "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1",
  "MFN1", "MFN2", "OPA1", "DNM1L", "FIS1",
  "SOD2", "PRDX3", "PRDX5",
  "HK2", "PKM", "LDHA", "SLC2A1", "PDK1"
)

genes <- metabolic_genes[metabolic_genes %in% rownames(tcga_mat_symbol)]

cat("Metabolic genes found:\n")
print(genes)


# -----------------------------
# 6. Select top miRNAs from Figure 5B
# -----------------------------

top_mirnas_corr <- diff_mirna %>%
  filter(!is.na(p_adj), p_adj < 0.05) %>%
  arrange(p_adj) %>%
  slice_head(n = 20) %>%
  pull(miRNA)

top_mirnas_corr <- top_mirnas_corr[
  top_mirnas_corr %in% rownames(mirna_mat_corr)
]

cat("Top miRNAs found:\n")
print(top_mirnas_corr)

if (length(top_mirnas_corr) == 0) {
  stop("No miRNAs found. Check rownames(mirna_mat_corr) and diff_mirna$miRNA.")
}


# -----------------------------
# 7. Subset matrices
# -----------------------------

mrna <- tcga_mat_symbol[genes, , drop = FALSE]
mir  <- mirna_mat_corr[top_mirnas_corr, , drop = FALSE]


# -----------------------------
# 8. Spearman correlation
# -----------------------------

cor_results <- expand.grid(
  miRNA = rownames(mir),
  Gene = rownames(mrna),
  stringsAsFactors = FALSE
)

cor_results <- cor_results %>%
  rowwise() %>%
  mutate(
    Correlation = cor(
      as.numeric(mir[miRNA, ]),
      as.numeric(mrna[Gene, ]),
      method = "spearman",
      use = "pairwise.complete.obs"
    ),
    p_value = cor.test(
      as.numeric(mir[miRNA, ]),
      as.numeric(mrna[Gene, ]),
      method = "spearman"
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

write.csv(
  cor_results,
  file.path(output_dir, "Figure5F_miRNA_MetabolicGene_Correlation_Table.csv"),
  row.names = FALSE
)


# -----------------------------
# 9. Plot heatmap
# -----------------------------

p5F <- ggplot(cor_results, aes(x = Gene, y = miRNA, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 3) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "miRNA–Metabolic Gene Correlation Landscape",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p5F)


# -----------------------------
# 10. Save PDF + JPG
# -----------------------------

ggsave(
  file.path(output_dir, "Figure5F_miRNA_MetabolicGene_Correlation.pdf"),
  plot = p5F,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_dir, "Figure5F_miRNA_MetabolicGene_Correlation.jpg"),
  plot = p5F,
  width = 10,
  height = 8,
  dpi = 300
)


# =========================================================
# FIGURE 5G: miRNA–Metabolic Gene Regulatory Network
# Uses significant miRNA–gene correlations from Figure 5F
# =========================================================

library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)

# -----------------------------
# 1. Use correct correlation table
# -----------------------------

cor_net <- cor_results

# Check required columns
stopifnot(all(c("miRNA", "Gene", "Correlation", "p_adj") %in% colnames(cor_net)))


# -----------------------------
# 2. Select significant correlations
# -----------------------------

sig_pairs <- cor_net %>%
  filter(!is.na(Correlation)) %>%
  filter(!is.na(p_adj)) %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  slice_head(n = 80)

cat("Network edges selected:", nrow(sig_pairs), "\n")


# -----------------------------
# 3. Build edge table
# -----------------------------

edges <- sig_pairs %>%
  mutate(
    from = miRNA,
    to = Gene,
    regulation = ifelse(Correlation > 0, "Positive", "Negative"),
    weight = abs(Correlation)
  ) %>%
  dplyr::select(from, to, Correlation, p_adj, regulation, weight)


# -----------------------------
# 4. Build node table
# -----------------------------

nodes <- data.frame(
  name = unique(c(edges$from, edges$to)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    type = ifelse(grepl("^hsa-miR", name), "miRNA", "Metabolic gene")
  )


# -----------------------------
# 5. Build network
# -----------------------------

g_net <- igraph::graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = FALSE
)


# -----------------------------
# 6. Plot network
# -----------------------------

p5G <- ggraph(g_net, layout = "fr") +
  geom_edge_link(
    aes(width = weight, color = regulation),
    alpha = 0.7
  ) +
  geom_node_point(
    aes(shape = type),
    size = 4,
    color = "black"
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_edge_color_manual(
    values = c(
      "Positive" = "#c53030",
      "Negative" = "#2b6cb0"
    )
  ) +
  scale_edge_width(range = c(0.3, 1.8)) +
  theme_void(base_size = 12) +
  labs(
    title = "miRNA–Metabolic Gene Regulatory Network",
    edge_color = "Correlation",
    edge_width = "|Spearman rho|",
    shape = "Node type"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p5G)


# -----------------------------
# 7. Save edge and node tables
# -----------------------------

write.csv(
  edges,
  file.path(output_dir, "Figure5G_miRNA_Gene_Network_Edges.csv"),
  row.names = FALSE
)

write.csv(
  nodes,
  file.path(output_dir, "Figure5G_miRNA_Gene_Network_Nodes.csv"),
  row.names = FALSE
)


# -----------------------------
# 8. Save figure
# -----------------------------

ggsave(
  file.path(output_dir, "Figure5G_miRNA_MetabolicGene_Network.pdf"),
  plot = p5G,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_dir, "Figure5G_miRNA_MetabolicGene_Network.jpg"),
  plot = p5G,
  width = 10,
  height = 8,
  dpi = 300
)

# =========================================================
# FIGURE 5H: State-specific miRNA–gene regulation bubble plot
# Correlations computed within each metabolic state
# x = metabolic gene
# y = miRNA
# size = |Spearman rho|
# color = correlation direction/strength
# facet = metabolic state
# =========================================================

library(dplyr)
library(ggplot2)
library(scales)

# -----------------------------
# 1. Prepare sample metadata
# -----------------------------

meta_state <- sample_info %>%
  mutate(sample15 = substr(sample, 1, 15)) %>%
  dplyr::select(sample15, cancer_code, metabolic_state) %>%
  distinct(sample15, .keep_all = TRUE)

meta_state <- meta_state %>%
  dplyr::filter(
    sample15 %in% colnames(mirna_mat_corr),
    sample15 %in% colnames(tcga_mat_symbol)
  )

cat("Matched samples with metabolic state:", nrow(meta_state), "\n")
print(table(meta_state$metabolic_state))


# -----------------------------
# 2. Settings
# -----------------------------

state_order <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

min_samples <- 20
corr_cutoff <- 0.25

genes_use <- genes
mirnas_use <- top_mirnas_corr


# -----------------------------
# 3. Compute correlations within each metabolic state
# -----------------------------

state_cor_list <- list()

for (st in unique(meta_state$metabolic_state)) {
  
  samples_st <- meta_state %>%
    dplyr::filter(metabolic_state == st) %>%
    dplyr::pull(sample15)
  
  samples_st <- intersect(samples_st, colnames(tcga_mat_symbol))
  samples_st <- intersect(samples_st, colnames(mirna_mat_corr))
  
  if (length(samples_st) < min_samples) next
  
  mrna_st <- tcga_mat_symbol[genes_use, samples_st, drop = FALSE]
  mir_st  <- mirna_mat_corr[mirnas_use, samples_st, drop = FALSE]
  
  tmp <- expand.grid(
    miRNA = rownames(mir_st),
    Gene = rownames(mrna_st),
    stringsAsFactors = FALSE
  )
  
  tmp$Correlation <- NA_real_
  
  for (k in seq_len(nrow(tmp))) {
    tmp$Correlation[k] <- cor(
      as.numeric(mir_st[tmp$miRNA[k], ]),
      as.numeric(mrna_st[tmp$Gene[k], ]),
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  }
  
  tmp$metabolic_state <- st
  tmp$n_samples <- length(samples_st)
  
  state_cor_list[[st]] <- tmp
}

state_cor_df <- dplyr::bind_rows(state_cor_list)


# -----------------------------
# 4. Clean labels + filter and format
# -----------------------------

state_order <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)

# See exact labels before cleaning
cat("Original labels:\n")
print(unique(state_cor_df$metabolic_state))

# Robustly standardize labels
state_cor_df <- state_cor_df %>%
  dplyr::mutate(
    metabolic_state = trimws(as.character(metabolic_state)),
    metabolic_state = dplyr::case_when(
      grepl("OXPHOS", metabolic_state, ignore.case = TRUE) ~
        "Mito-High / Glyco-Low (OXPHOS-dominant)",
      grepl("Glycolytic", metabolic_state, ignore.case = TRUE) ~
        "Mito-Low / Glyco-High (Glycolytic)",
      grepl("Both-High|Hybrid", metabolic_state, ignore.case = TRUE) ~
        "Both-High (Hybrid)",
      grepl("Both-Low|Low metabolic", metabolic_state, ignore.case = TRUE) ~
        "Both-Low (Low metabolic)",
      TRUE ~ NA_character_
    )
  )

cat("Cleaned labels:\n")
print(table(state_cor_df$metabolic_state, useNA = "ifany"))

bubble_state_df <- state_cor_df %>%
  dplyr::filter(!is.na(Correlation)) %>%
  dplyr::filter(!is.na(metabolic_state)) %>%
  dplyr::mutate(
    abs_cor = abs(Correlation),
    direction = ifelse(Correlation > 0, "Positive", "Negative"),
    metabolic_state = factor(metabolic_state, levels = state_order)
  ) %>%
  dplyr::group_by(metabolic_state) %>%
  dplyr::arrange(desc(abs_cor), .by_group = TRUE) %>%
  dplyr::slice_head(n = 80) %>%
  dplyr::ungroup()

cat("Final plotted counts:\n")
print(table(bubble_state_df$metabolic_state))
write.csv(
  bubble_state_df,
  file.path(output_dir, "Figure5H_StateSpecific_miRNA_Gene_Correlations.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. Bubble plot
# -----------------------------

p5H <- ggplot(
  bubble_state_df,
  aes(
    x = Gene,
    y = miRNA,
    size = abs_cor,
    fill = Correlation
  )
) +
  geom_point(shape = 21, color = "black", alpha = 0.85, stroke = 0.25) +
  facet_wrap(~ metabolic_state, ncol = 2) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(
    range = c(1.5, 6),
    name = "|Spearman rho|"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "State-specific miRNA–Metabolic Gene Regulation",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p5H)


# -----------------------------
# 6. Save
# -----------------------------

ggsave(
  file.path(output_dir, "Figure5H_StateSpecific_miRNA_Gene_BubblePlot.pdf"),
  plot = p5H,
  width = 13,
  height = 10
)

ggsave(
  file.path(output_dir, "Figure5H_StateSpecific_miRNA_Gene_BubblePlot.jpg"),
  plot = p5H,
  width = 13,
  height = 10,
  dpi = 300
)

# =========================================================
# FIGURE 5I: Clinical relevance of key metabolic-state miRNAs
# Overall survival analysis
# =========================================================

library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(patchwork)

# -----------------------------
# 1. Add survival file path near your other file paths
# -----------------------------
survival_file <- "Survival_SupplementalTable_S1_20171025_xena_sp"

# -----------------------------
# 2. Read survival file
# -----------------------------
surv <- fread(survival_file, data.table = FALSE)
colnames(surv) <- make.names(colnames(surv))

sample_col_surv <- grep(
  "sample|bcr_patient_barcode|sampleID",
  colnames(surv),
  ignore.case = TRUE,
  value = TRUE
)[1]

surv$sample <- clean_tcga_id(surv[[sample_col_surv]])

surv$OS.time <- as.numeric(surv$OS.time)

if (is.character(surv$OS)) {
  surv$OS <- ifelse(grepl("DECEASED|DEAD|1", toupper(surv$OS)), 1, 0)
} else {
  surv$OS <- ifelse(as.numeric(surv$OS) > 0, 1, 0)
}

surv2 <- surv[, c("sample", "OS.time", "OS"), drop = FALSE]
surv2 <- surv2[!duplicated(surv2$sample), , drop = FALSE]

cat("Survival rows:\n")
print(nrow(surv2))

# -----------------------------
# 3. Prepare miRNA matrix for survival
# -----------------------------
mirna_surv_mat <- mirna_mat

colnames(mirna_surv_mat) <- gsub("\\.", "-", colnames(mirna_surv_mat))
colnames(mirna_surv_mat) <- substr(colnames(mirna_surv_mat), 1, 15)

mirna_surv_mat <- mirna_surv_mat[, !duplicated(colnames(mirna_surv_mat)), drop = FALSE]

# -----------------------------
# 4. Select key miRNAs
# -----------------------------
key_mirnas <- diff_mirna %>%
  filter(p_adj < 0.01) %>%
  arrange(p_adj) %>%
  slice_head(n = 10) %>%
  pull(miRNA)

key_mirnas <- key_mirnas[key_mirnas %in% rownames(mirna_surv_mat)]

cat("Key miRNAs used for survival:\n")
print(key_mirnas)

# -----------------------------
# 5. Kaplan-Meier survival plots
# -----------------------------
plot_list <- list()
cox_results <- list()

for (mir in key_mirnas) {
  
  expr_vec <- mirna_surv_mat[mir, ]
  
  dat <- data.frame(
    sample = names(expr_vec),
    miRNA_expr = as.numeric(expr_vec),
    stringsAsFactors = FALSE
  )
  
  dat$sample <- clean_tcga_id(dat$sample)
  
  dat <- merge(dat, surv2, by = "sample")
  
  dat <- dat %>%
    filter(!is.na(OS.time), !is.na(OS), OS.time > 0, !is.na(miRNA_expr))
  
  dat$group <- ifelse(
    dat$miRNA_expr >= median(dat$miRNA_expr, na.rm = TRUE),
    "High",
    "Low"
  )
  
  dat$group <- factor(dat$group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS.time, OS) ~ group, data = dat)
  cox_fit <- coxph(Surv(OS.time, OS) ~ group, data = dat)
  cox_sum <- summary(cox_fit)
  
  hr <- round(cox_sum$conf.int[1, "exp(coef)"], 2)
  lower <- round(cox_sum$conf.int[1, "lower .95"], 2)
  upper <- round(cox_sum$conf.int[1, "upper .95"], 2)
  pval <- signif(cox_sum$coefficients[1, "Pr(>|z|)"], 3)
  
  cox_results[[mir]] <- data.frame(
    miRNA = mir,
    HR = hr,
    lower95 = lower,
    upper95 = upper,
    p_value = pval,
    n = nrow(dat)
  )
  
  p <- ggsurvplot(
    fit,
    data = dat,
    pval = TRUE,
    risk.table = FALSE,
    conf.int = FALSE,
    size = 1.1,
    palette = c("#2b6cb0", "#c53030"),
    legend.title = mir,
    legend.labs = c("Low", "High"),
    ggtheme = theme_bw(base_size = 12),
    xlab = "Time (days)",
    ylab = "Overall survival probability"
  )
  
  p_plot <- p$plot +
    annotate(
      "label",
      x = max(dat$OS.time, na.rm = TRUE) * 0.45,
      y = 0.18,
      label = paste0(
        "HR = ", hr,
        " (", lower, "-", upper, ")\n",
        "Cox p = ", pval
      ),
      size = 3,
      label.size = 0.25,
      fill = "white"
    ) +
    labs(
      title = paste0("Overall survival by ", mir, " expression")
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  plot_list[[mir]] <- p_plot
}

# -----------------------------
# 6. Combine plots
# -----------------------------
p5I <- wrap_plots(plot_list, ncol = length(plot_list))

print(p5I)

# -----------------------------
# 7. Save Cox table
# -----------------------------
cox_table_miRNA <- bind_rows(cox_results)

write.csv(
  cox_table_miRNA,
  file.path(output_dir, "Figure5I_miRNA_Survival_Cox_Table.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Save figure
# -----------------------------
ggsave(
  file.path(output_dir, "Figure5I_miRNA_Survival.pdf"),
  plot = p5I,
  width = 12,
  height = 5.5
)

ggsave(
  file.path(output_dir, "Figure5I_miRNA_Survival.jpg"),
  plot = p5I,
  width = 12,
  height = 5.5,
  dpi = 300
)


# =========================================================
# FIGURE 5I + SUPPLEMENTARY S5I
# Cox forest plots for metabolic-state-associated miRNAs
# Main: Top 10 miRNAs by differential adjusted p-value
# Supplementary: All significant miRNAs
# Cox model stratified by cancer type
# =========================================================

library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(forcats)

# -----------------------------
# 1. Survival file
# -----------------------------
survival_file <- "Survival_SupplementalTable_S1_20171025_xena_sp"

surv <- fread(survival_file, data.table = FALSE)
colnames(surv) <- make.names(colnames(surv))

sample_col_surv <- grep(
  "sample|bcr_patient_barcode|sampleID",
  colnames(surv),
  ignore.case = TRUE,
  value = TRUE
)[1]

surv$sample <- clean_tcga_id(surv[[sample_col_surv]])
surv$OS.time <- as.numeric(surv$OS.time)

if (is.character(surv$OS)) {
  surv$OS <- ifelse(grepl("DECEASED|DEAD|1", toupper(surv$OS)), 1, 0)
} else {
  surv$OS <- ifelse(as.numeric(surv$OS) > 0, 1, 0)
}

surv2 <- surv[, c("sample", "OS.time", "OS"), drop = FALSE]
surv2 <- surv2[!duplicated(surv2$sample), ]

# -----------------------------
# 2. Cancer type metadata
# -----------------------------
cancer_meta <- sample_info %>%
  dplyr::select(sample, cancer_code) %>%
  distinct(sample, .keep_all = TRUE)

# -----------------------------
# 3. miRNA matrix for survival
# -----------------------------
mirna_surv_mat <- mirna_mat

colnames(mirna_surv_mat) <- gsub("\\.", "-", colnames(mirna_surv_mat))
colnames(mirna_surv_mat) <- substr(colnames(mirna_surv_mat), 1, 15)

mirna_surv_mat <- mirna_surv_mat[, !duplicated(colnames(mirna_surv_mat)), drop = FALSE]

# -----------------------------
# 4. Select miRNAs
# -----------------------------
all_sig_mirnas <- diff_mirna %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  pull(miRNA)

all_sig_mirnas <- all_sig_mirnas[all_sig_mirnas %in% rownames(mirna_surv_mat)]

top10_mirnas <- diff_mirna %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  slice_head(n = 10) %>%
  pull(miRNA)

top10_mirnas <- top10_mirnas[top10_mirnas %in% rownames(mirna_surv_mat)]

cat("Top 10 miRNAs used in main figure:\n")
print(top10_mirnas)

cat("All significant miRNAs used in supplementary:\n")
cat("n =", length(all_sig_mirnas), "\n")

# -----------------------------
# 5. Function: Cox model per miRNA
# -----------------------------
run_mirna_cox <- function(mirna_vector) {
  
  cox_results <- list()
  
  for (mir in mirna_vector) {
    
    expr_vec <- mirna_surv_mat[mir, ]
    
    dat <- data.frame(
      sample = names(expr_vec),
      miRNA_expr = as.numeric(expr_vec),
      stringsAsFactors = FALSE
    )
    
    dat$sample <- clean_tcga_id(dat$sample)
    
    dat <- dat %>%
      left_join(surv2, by = "sample") %>%
      left_join(cancer_meta, by = "sample") %>%
      filter(
        !is.na(OS.time),
        !is.na(OS),
        OS.time > 0,
        !is.na(miRNA_expr),
        !is.na(cancer_code)
      )
    
    dat$group <- ifelse(
      dat$miRNA_expr >= median(dat$miRNA_expr, na.rm = TRUE),
      "High",
      "Low"
    )
    
    dat$group <- factor(dat$group, levels = c("Low", "High"))
    dat$cancer_code <- factor(dat$cancer_code)
    
    # Need enough events and both groups
    if (nrow(dat) < 50 || length(unique(dat$group)) < 2 || sum(dat$OS == 1, na.rm = TRUE) < 10) next
    
    cox_fit <- tryCatch(
      coxph(Surv(OS.time, OS) ~ group + strata(cancer_code), data = dat),
      error = function(e) NULL
    )
    
    if (is.null(cox_fit)) next
    
    cox_sum <- summary(cox_fit)
    
    cox_results[[mir]] <- data.frame(
      miRNA = mir,
      HR = cox_sum$conf.int[1, "exp(coef)"],
      lower95 = cox_sum$conf.int[1, "lower .95"],
      upper95 = cox_sum$conf.int[1, "upper .95"],
      p_value = cox_sum$coefficients[1, "Pr(>|z|)"],
      n = nrow(dat),
      events = sum(dat$OS == 1, na.rm = TRUE),
      model = "Cox model stratified by cancer type",
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(cox_results)
}

# -----------------------------
# 6. Run Cox models
# -----------------------------
cox_top10 <- run_mirna_cox(top10_mirnas)
cox_all   <- run_mirna_cox(all_sig_mirnas)

cox_all <- cox_all %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    direction = ifelse(HR > 1, "Risk", "Protective"),
    HR_label = paste0(
      round(HR, 2), " (",
      round(lower95, 2), "-",
      round(upper95, 2), ")"
    )
  ) %>%
  arrange(p_value)

cox_top10 <- cox_all %>%
  filter(miRNA %in% top10_mirnas) %>%
  arrange(match(miRNA, top10_mirnas))

# -----------------------------
# 7. Save supplementary Cox table
# -----------------------------
write.csv(
  cox_all,
  file.path(output_dir, "Supplementary_Table_miRNA_Cox_AllSignificant.csv"),
  row.names = FALSE
)

write.csv(
  cox_top10,
  file.path(output_dir, "Figure5I_Top10_miRNA_Cox_Table.csv"),
  row.names = FALSE
)

cat("\nMain Figure 5I Cox results:\n")
print(cox_top10)

cat("\nSupplementary Cox results summary:\n")
print(table(cox_all$direction))

# -----------------------------
# 8. Main forest plot: top 10
# -----------------------------
p5I <- ggplot(
  cox_top10,
  aes(
    x = HR,
    y = fct_reorder(miRNA, HR),
    color = direction
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = lower95, xmax = upper95),
    height = 0.25,
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Risk" = "#c53030",
      "Protective" = "#2b6cb0"
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Prognostic Impact of Metabolic-State-Associated miRNAs",
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    color = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p5I)

ggsave(
  file.path(output_dir, "Figure5I_Top10_miRNA_Cox_ForestPlot.pdf"),
  plot = p5I,
  width = 10,
  height = 8
)

ggsave(
  file.path(output_dir, "Figure5I_Top10_miRNA_Cox_ForestPlot.jpg"),
  plot = p5I,
  width = 10,
  height = 8,
  dpi = 300
)

# =========================================================
# SUPPLEMENTARY FIGURE S5J:
# Cancer-type-specific Cox forest plot for key miRNAs
# =========================================================

library(dplyr)
library(survival)
library(ggplot2)
library(forcats)

# -----------------------------
# 1. Cancer-type-specific Cox
# -----------------------------

cox_by_cancer <- list()

for (mir in all_sig_mirnas) {
  
  expr_vec <- mirna_surv_mat[mir, ]
  
  dat <- data.frame(
    sample = names(expr_vec),
    miRNA_expr = as.numeric(expr_vec),
    stringsAsFactors = FALSE
  )
  
  dat$sample <- clean_tcga_id(dat$sample)
  
  dat <- dat %>%
    left_join(surv2, by = "sample") %>%
    left_join(cancer_meta, by = "sample") %>%
    filter(
      !is.na(OS.time),
      !is.na(OS),
      OS.time > 0,
      !is.na(miRNA_expr),
      !is.na(cancer_code)
    )
  
  # Median cutoff within each cancer type
  dat <- dat %>%
    group_by(cancer_code) %>%
    mutate(
      group = ifelse(
        miRNA_expr >= median(miRNA_expr, na.rm = TRUE),
        "High", "Low"
      )
    ) %>%
    ungroup()
  
  dat$group <- factor(dat$group, levels = c("Low", "High"))
  
  for (ct in unique(dat$cancer_code)) {
    
    dct <- dat %>% filter(cancer_code == ct)
    
    if (nrow(dct) < 40) next
    if (sum(dct$OS == 1, na.rm = TRUE) < 10) next
    if (length(unique(dct$group)) < 2) next
    
    fit <- tryCatch(
      coxph(Surv(OS.time, OS) ~ group, data = dct),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    s <- summary(fit)
    
    cox_by_cancer[[paste(mir, ct, sep = "_")]] <- data.frame(
      miRNA = mir,
      cancer_type = ct,
      HR = s$conf.int[1, "exp(coef)"],
      lower95 = s$conf.int[1, "lower .95"],
      upper95 = s$conf.int[1, "upper .95"],
      p_value = s$coefficients[1, "Pr(>|z|)"],
      n = nrow(dct),
      events = sum(dct$OS == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

cox_by_cancer_table <- bind_rows(cox_by_cancer) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    direction = ifelse(HR > 1, "Risk", "Protective"),
    HR_label = paste0(
      round(HR, 2), " (",
      round(lower95, 2), "-",
      round(upper95, 2), ")"
    )
  ) %>%
  arrange(FDR, p_value)

# -----------------------------
# 2. Save supplementary table
# -----------------------------

write.csv(
  cox_by_cancer_table,
  file.path(output_dir, "Supplementary_Table_miRNA_Cox_ByCancerType.csv"),
  row.names = FALSE
)

# -----------------------------
# 3. Select significant cancer-specific results for plotting
# -----------------------------

plot_cancer_cox <- cox_by_cancer_table %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice_head(n = 30) %>%
  mutate(
    label = paste0(miRNA, " | ", cancer_type),
    label = fct_reorder(label, HR)
  )

# -----------------------------
# 4. Supplementary forest plot
# -----------------------------

pS5J <- ggplot(
  plot_cancer_cox,
  aes(
    x = HR,
    y = label,
    color = direction
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = lower95, xmax = upper95),
    height = 0.25,
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Risk" = "#c53030",
      "Protective" = "#2b6cb0"
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Cancer-Type-Specific Prognostic Associations of Metabolic-State miRNAs",
    x = "Hazard Ratio for High vs Low miRNA Expression (95% CI)",
    y = NULL,
    color = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(pS5J)

ggsave(
  file.path(output_dir, "Supplementary_Figure_S5J_miRNA_Cox_ByCancerType.pdf"),
  plot = pS5J,
  width = 11,
  height = 9
)

ggsave(
  file.path(output_dir, "Supplementary_Figure_S5J_miRNA_Cox_ByCancerType.jpg"),
  plot = pS5J,
  width = 11,
  height = 9,
  dpi = 300
)

# =========================================================
# Supplementary Figure S5J
# Data-driven selected cancer-type-specific miRNA Cox forest plot
# =========================================================

library(dplyr)
library(ggplot2)
library(forcats)

# -----------------------------
# 1. Select cancer types objectively
# -----------------------------

cancer_summary <- cox_by_cancer_table %>%
  group_by(cancer_type) %>%
  summarise(
    n_total = max(n, na.rm = TRUE),
    events_total = max(events, na.rm = TRUE),
    n_models = n(),
    n_significant = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_significant), desc(events_total), desc(n_total))

write.csv(
  cancer_summary,
  file.path(output_dir, "Supplementary_Table_CancerType_Selection_Criteria.csv"),
  row.names = FALSE
)

selected_cancers <- cancer_summary %>%
  filter(
    n_total >= 80,
    events_total >= 20,
    n_models >= 5
  ) %>%
  arrange(desc(n_significant), desc(events_total), desc(n_total)) %>%
  slice_head(n = 10) %>%
  pull(cancer_type)

cat("Selected cancer types for S5J:\n")
print(selected_cancers)

# -----------------------------
# 2. Select top 3 miRNA associations per cancer type
# -----------------------------

plot_selected_cancers_all <- cox_by_cancer_table %>%
  filter(
    cancer_type %in% selected_cancers,
    !is.na(HR),
    !is.na(lower95),
    !is.na(upper95),
    !is.na(p_value),
    !is.na(FDR)
  ) %>%
  group_by(cancer_type) %>%
  arrange(FDR, p_value, .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(
    direction = ifelse(HR > 1, "Risk", "Protective"),
    label = paste0(miRNA, " | ", cancer_type),
    HR_label = paste0(
      round(HR, 2), " (",
      round(lower95, 2), "-",
      round(upper95, 2), ")"
    )
  )

write.csv(
  plot_selected_cancers_all,
  file.path(output_dir, "Supplementary_Table_S5J_SelectedCancerType_miRNA_Cox_AllSelected.csv"),
  row.names = FALSE
)

# Remove unstable very wide CI associations from plot only
plot_selected_cancers_plot <- plot_selected_cancers_all %>%
  filter(upper95 <= 8) %>%
  mutate(
    label = fct_reorder(label, HR)
  )

write.csv(
  plot_selected_cancers_plot,
  file.path(output_dir, "Supplementary_Table_S5J_SelectedCancerType_miRNA_Cox_Plotted.csv"),
  row.names = FALSE
)

# -----------------------------
# 3. Forest plot
# -----------------------------

pS5J_selected <- ggplot(
  plot_selected_cancers_plot,
  aes(
    x = HR,
    y = label,
    color = direction
  )
) +
  geom_point(size = 2.8) +
  geom_errorbarh(
    aes(xmin = lower95, xmax = upper95),
    height = 0.22,
    linewidth = 0.65
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Risk" = "#c53030",
      "Protective" = "#2b6cb0"
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Cancer-Type-Specific Prognostic Associations of Selected Metabolic-State miRNAs",
    x = "Hazard Ratio for High vs Low miRNA Expression (95% CI)",
    y = NULL,
    color = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(pS5J_selected)

ggsave(
  file.path(output_dir, "Supplementary_Figure_S5J_SelectedCancerTypes_miRNA_Cox_ForestPlot.pdf"),
  plot = pS5J_selected,
  width = 11,
  height = 8
)

ggsave(
  file.path(output_dir, "Supplementary_Figure_S5J_SelectedCancerTypes_miRNA_Cox_ForestPlot.jpg"),
  plot = pS5J_selected,
  width = 11,
  height = 8,
  dpi = 300
)