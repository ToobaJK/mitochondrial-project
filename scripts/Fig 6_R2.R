# =========================================================
# Check metabolic state table level + create sample-level table
# Folder: D:/UAEU/Fig 6
# Created By: TOOBA
# =========================================================

library(data.table)
library(dplyr)

# -----------------------------
# 1. Set folder paths
# -----------------------------

base_dir <- "D:/UAEU/Fig 6"

out_dir <- file.path(base_dir, "Output_lncRNA_Figure6")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output folder created at:\n", out_dir, "\n\n")

# -----------------------------
# 2. Load metabolic state table
# -----------------------------

state_file <- file.path(base_dir, "Figure1E_TCGA_MetabolicStates_Table.tsv")

meta_state <- fread(state_file)

cat("Metabolic state table dimensions:\n")
print(dim(meta_state))

cat("\nColumn names:\n")
print(colnames(meta_state))

cat("\nFirst few rows:\n")
print(head(meta_state))

# -----------------------------
# 3. Check whether cancer-level or sample-level
# -----------------------------

cat("\nChecking table level...\n")

has_sample_col <- any(grepl("sample|barcode|patient", colnames(meta_state), ignore.case = TRUE))

if ("cancer_code" %in% colnames(meta_state)) {
  n_unique_cancer <- length(unique(meta_state$cancer_code))
} else {
  n_unique_cancer <- NA
}

cat("Number of rows:", nrow(meta_state), "\n")
cat("Unique cancer codes:", n_unique_cancer, "\n")
cat("Has sample-like column:", has_sample_col, "\n")

if (!has_sample_col && nrow(meta_state) <= 40) {
  cat("\nConclusion: This is CANCER-LEVEL metabolic state table.\n")
} else {
  cat("\nConclusion: This may already be SAMPLE-LEVEL table.\n")
}

# -----------------------------
# 4. Load phenotype file
# -----------------------------

pheno_file <- list.files(
  base_dir,
  pattern = "^TCGA_phenotype_denseDataOnlyDownload",
  full.names = TRUE
)[1]

pheno <- fread(pheno_file)

# first column is sample ID
setnames(pheno, old = colnames(pheno)[1], new = "sample")

pheno <- pheno %>%
  mutate(
    sample15 = substr(sample, 1, 15)
  )

cat("\nPhenotype dimensions:\n")
print(dim(pheno))

cat("\nPhenotype column names:\n")
print(colnames(pheno))

# -----------------------------
# 5. Find cancer type column
# -----------------------------

possible_cancer_cols <- c(
  "_primary_disease",
  "primary_disease",
  "cancer_type",
  "disease"
)

cancer_col <- possible_cancer_cols[possible_cancer_cols %in% colnames(pheno)][1]

if (is.na(cancer_col)) {
  stop("Could not find cancer type column in phenotype file.")
}

cat("\nUsing cancer column:", cancer_col, "\n")

# -----------------------------
# 6. TCGA cancer name to TCGA code mapping
# -----------------------------

tcga_abbrev_map <- c(
  "Acute Myeloid Leukemia" = "LAML",
  "Adrenocortical Cancer" = "ACC",
  "Bladder Urothelial Carcinoma" = "BLCA",
  "Brain Lower Grade Glioma" = "LGG",
  "Breast Invasive Carcinoma" = "BRCA",
  "Cervical & Endocervical Cancer" = "CESC",
  "Cholangiocarcinoma" = "CHOL",
  "Colon Adenocarcinoma" = "COAD",
  "Diffuse Large B-Cell Lymphoma" = "DLBC",
  "Esophageal Carcinoma" = "ESCA",
  "Glioblastoma Multiforme" = "GBM",
  "Head & Neck Squamous Cell Carcinoma" = "HNSC",
  "Kidney Chromophobe" = "KICH",
  "Kidney Clear Cell Carcinoma" = "KIRC",
  "Kidney Papillary Cell Carcinoma" = "KIRP",
  "Liver Hepatocellular Carcinoma" = "LIHC",
  "Lung Adenocarcinoma" = "LUAD",
  "Lung Squamous Cell Carcinoma" = "LUSC",
  "Mesothelioma" = "MESO",
  "Ovarian Serous Cystadenocarcinoma" = "OV",
  "Pancreatic Adenocarcinoma" = "PAAD",
  "Pheochromocytoma & Paraganglioma" = "PCPG",
  "Prostate Adenocarcinoma" = "PRAD",
  "Rectum Adenocarcinoma" = "READ",
  "Sarcoma" = "SARC",
  "Skin Cutaneous Melanoma" = "SKCM",
  "Stomach Adenocarcinoma" = "STAD",
  "Testicular Germ Cell Tumor" = "TGCT",
  "Thymoma" = "THYM",
  "Thyroid Carcinoma" = "THCA",
  "Uterine Corpus Endometrioid Carcinoma" = "UCEC",
  "Uterine Carcinosarcoma" = "UCS",
  "Uveal Melanoma" = "UVM"
)

to_tcga_code <- function(x) {
  x2 <- tools::toTitleCase(as.character(x))
  out <- unname(tcga_abbrev_map[x2])
  out[is.na(out)] <- x2[is.na(out)]
  out
}

# -----------------------------
# 7. Create sample-level phenotype table
# -----------------------------

sample_info <- pheno %>%
  transmute(
    sample,
    sample15,
    cancer_type = .data[[cancer_col]],
    cancer_code = to_tcga_code(.data[[cancer_col]])
  ) %>%
  distinct(sample15, .keep_all = TRUE)

cat("\nSample info preview:\n")
print(head(sample_info))

# -----------------------------
# 8. Join cancer-level metabolic state to samples
# -----------------------------

meta_state_sample <- sample_info %>%
  left_join(meta_state, by = "cancer_code")

cat("\nSample-level metabolic state table dimensions:\n")
print(dim(meta_state_sample))

cat("\nNumber of samples with assigned state:\n")
print(sum(!is.na(meta_state_sample$state)))

cat("\nState distribution across samples:\n")
print(table(meta_state_sample$state, useNA = "ifany"))

# -----------------------------
# 9. Save output files
# -----------------------------

write.table(
  meta_state_sample,
  file.path(out_dir, "Sample_Level_TCGA_MetabolicStates_Table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  meta_state_sample,
  file.path(out_dir, "Sample_Level_TCGA_MetabolicStates_Table.csv"),
  row.names = FALSE
)

cat("\nDone. Files saved in:\n", out_dir, "\n")

# =========================================================
# Run this first to check lncRNA
# =========================================================

library(data.table)

gtf_file <- "D:/UAEU/Fig 6/gencode.v22.annotation.gtf.gz"

# Step 1: read lines properly
gtf_lines <- readLines(gzfile(gtf_file))

# Step 2: remove comments
gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

# Step 3: split by TAB
gtf_split <- strsplit(gtf_lines, "\t", fixed = TRUE)

# Step 4: keep only correct rows (9 columns)
gtf_split <- gtf_split[lengths(gtf_split) == 9]

# Step 5: convert to data.frame
gtf <- as.data.frame(do.call(rbind, gtf_split), stringsAsFactors = FALSE)
colnames(gtf) <- paste0("V", 1:9)

# Step 6: check structure
cat("GTF structure:\n")
print(dim(gtf))
print(head(gtf[,1:3]))

# Step 7: NOW this should work
table(gtf$V3)[1:10]

###Extract genes
gtf_gene <- gtf[gtf$V3 == "gene", ]

cat("Gene rows:\n")
print(dim(gtf_gene))

head(gtf_gene$V9, 3)

##Now check lncRNA
grep("lincRNA", gtf_gene$V9, value = TRUE)[1:5]


extract_attr <- function(x, key) {
  pattern <- paste0(key, " \"[^\"]+\"")
  hit <- regmatches(x, regexpr(pattern, x))
  hit <- gsub(paste0(key, " \"|\""), "", hit)
  hit[hit == ""] <- NA
  hit
}

annot <- data.frame(
  gene_id = extract_attr(gtf_gene$V9, "gene_id"),
  gene_name = extract_attr(gtf_gene$V9, "gene_name"),
  gene_type = extract_attr(gtf_gene$V9, "gene_type"),
  stringsAsFactors = FALSE
)

annot$gene_id_clean <- sub("\\..*$", "", annot$gene_id)

lnc_types <- c(
  "lincRNA",
  "antisense",
  "processed_transcript",
  "sense_intronic",
  "sense_overlapping",
  "3prime_overlapping_ncRNA",
  "macro_lncRNA",
  "non_coding"
)

lnc_annot <- annot %>%
  dplyr::filter(gene_type %in% lnc_types) %>%
  dplyr::distinct(gene_id_clean, gene_name, gene_type)

cat("lncRNA genes in annotation:", nrow(lnc_annot), "\n")
head(lnc_annot)

# =========================================================
# RAM-OPTIMIZED FIGURE 6: lncRNA analysis
# =========================================================

# -----------------------------
# 0. Paths
# ----------------------------
expr_file  <- file.path(base_dir, "tcga_RSEM_gene_tpm")
state_file <- file.path(base_dir, "Figure1E_TCGA_MetabolicStates_Table.tsv")
pheno_file <- list.files(base_dir, pattern = "^TCGA_phenotype_denseDataOnlyDownload", full.names = TRUE)[1]
surv_file  <- list.files(base_dir, pattern = "^Survival_SupplementalTable", full.names = TRUE)[1]

to_tcga15 <- function(x) substr(x, 1, 15)

zscore_rows <- function(m) {
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  z
}

# -----------------------------
# 1. Check lncRNA annotation
# -----------------------------
stopifnot(exists("lnc_annot"))
cat("lncRNA genes in annotation:", nrow(lnc_annot), "\n")
print(head(lnc_annot))

# -----------------------------
# 2. Chunk extract lncRNAs only - SAFE VERSION
# -----------------------------

header <- names(fread(expr_file, nrows = 0))
sample_names <- to_tcga15(header[-1])

chunk_size <- 1000
skip_line <- 1
lnc_chunks <- list()
chunk_id <- 1

repeat {
  
  dt <- tryCatch(
    fread(
      expr_file,
      skip = skip_line,
      nrows = chunk_size,
      header = FALSE,
      data.table = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dt) || nrow(dt) == 0) {
    cat("Reached end of file at line:", skip_line, "\n")
    break
  }
  
  cat("Reading chunk starting at line:", skip_line, "\n")
  
  gene_ids <- sub("\\..*$", "", dt[[1]])
  keep <- gene_ids %in% lnc_annot$gene_id_clean
  
  if (any(keep)) {
    sub_dt <- dt[keep]
    sub_gene_ids <- sub("\\..*$", "", sub_dt[[1]])
    
    gene_names <- lnc_annot$gene_name[
      match(sub_gene_ids, lnc_annot$gene_id_clean)
    ]
    
    mat <- as.matrix(sub_dt[, -1, with = FALSE])
    rownames(mat) <- gene_names
    colnames(mat) <- sample_names
    
    lnc_chunks[[chunk_id]] <- mat
    chunk_id <- chunk_id + 1
  }
  
  rm(dt)
  gc()
  
  skip_line <- skip_line + chunk_size
}

lnc_mat_raw <- do.call(rbind, lnc_chunks)

cat("Raw lncRNA matrix:", dim(lnc_mat_raw), "\n")
# -----------------------------
# 3. Remove duplicated samples and duplicated lncRNA names
# -----------------------------

keep_cols <- !duplicated(colnames(lnc_mat_raw))
lnc_mat_raw <- lnc_mat_raw[, keep_cols, drop = FALSE]

lnc_dt <- as.data.table(lnc_mat_raw, keep.rownames = "lncRNA")
lnc_dt <- lnc_dt[, lapply(.SD, mean, na.rm = TRUE), by = lncRNA]

lnc_mat <- as.matrix(lnc_dt[, -1])
rownames(lnc_mat) <- lnc_dt$lncRNA

rm(lnc_dt, lnc_mat_raw)
gc()

# Log transform
lnc_mat <- log2(lnc_mat + 0.001)

# Filter low-expression lncRNAs
lnc_mat <- lnc_mat[
  rowMeans(lnc_mat > log2(0.1 + 0.001), na.rm = TRUE) >= 0.10,
  ,
  drop = FALSE
]

cat("Final lncRNA matrix:", dim(lnc_mat), "\n")

write.csv(
  data.frame(lncRNA = rownames(lnc_mat)),
  file.path(out_dir, "Figure6_lncRNAs_detected.csv"),
  row.names = FALSE
)

# -----------------------------
# 4. Sample-level metabolic state table
# -----------------------------

meta_state <- fread(state_file)

pheno <- fread(pheno_file)
setnames(pheno, old = colnames(pheno)[1], new = "sample")

tcga_abbrev_map <- c(
  "Acute Myeloid Leukemia"="LAML",
  "Adrenocortical Cancer"="ACC",
  "Bladder Urothelial Carcinoma"="BLCA",
  "Brain Lower Grade Glioma"="LGG",
  "Breast Invasive Carcinoma"="BRCA",
  "Cervical & Endocervical Cancer"="CESC",
  "Cholangiocarcinoma"="CHOL",
  "Colon Adenocarcinoma"="COAD",
  "Diffuse Large B-Cell Lymphoma"="DLBC",
  "Esophageal Carcinoma"="ESCA",
  "Glioblastoma Multiforme"="GBM",
  "Head & Neck Squamous Cell Carcinoma"="HNSC",
  "Kidney Chromophobe"="KICH",
  "Kidney Clear Cell Carcinoma"="KIRC",
  "Kidney Papillary Cell Carcinoma"="KIRP",
  "Liver Hepatocellular Carcinoma"="LIHC",
  "Lung Adenocarcinoma"="LUAD",
  "Lung Squamous Cell Carcinoma"="LUSC",
  "Mesothelioma"="MESO",
  "Ovarian Serous Cystadenocarcinoma"="OV",
  "Pancreatic Adenocarcinoma"="PAAD",
  "Pheochromocytoma & Paraganglioma"="PCPG",
  "Prostate Adenocarcinoma"="PRAD",
  "Rectum Adenocarcinoma"="READ",
  "Sarcoma"="SARC",
  "Skin Cutaneous Melanoma"="SKCM",
  "Stomach Adenocarcinoma"="STAD",
  "Testicular Germ Cell Tumor"="TGCT",
  "Thymoma"="THYM",
  "Thyroid Carcinoma"="THCA",
  "Uterine Corpus Endometrioid Carcinoma"="UCEC",
  "Uterine Carcinosarcoma"="UCS",
  "Uveal Melanoma"="UVM"
)

to_tcga_code <- function(x) {
  x2 <- tools::toTitleCase(as.character(x))
  out <- unname(tcga_abbrev_map[x2])
  out[is.na(out)] <- x2[is.na(out)]
  out
}

sample_info <- pheno %>%
  mutate(sample15 = to_tcga15(sample)) %>%
  transmute(
    sample,
    sample15,
    cancer_type = .data[["_primary_disease"]],
    cancer_code = to_tcga_code(.data[["_primary_disease"]])
  ) %>%
  distinct(sample15, .keep_all = TRUE)

meta_state_sample <- sample_info %>%
  left_join(meta_state, by = "cancer_code")

write.table(
  meta_state_sample,
  file.path(out_dir, "Sample_Level_TCGA_MetabolicStates_Table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 5. Match lncRNA samples with metabolic states
# -----------------------------


common_samples <- intersect(colnames(lnc_mat), meta_state_sample$sample15)

lnc_mat2 <- lnc_mat[, common_samples, drop = FALSE]

meta2 <- meta_state_sample %>%
  filter(sample15 %in% common_samples) %>%
  arrange(match(sample15, common_samples))

lnc_mat2 <- lnc_mat2[, meta2$sample15, drop = FALSE]

state_levels <- unique(meta_state_sample$state)

meta2$state <- factor(meta2$state, levels = state_levels)

table(meta2$state, useNA = "ifany")

state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

meta2$state <- factor(meta2$state, levels = state_levels)

cat("Matched lncRNA samples:", ncol(lnc_mat2), "\n")
print(table(meta2$state, useNA = "ifany"))

# =========================================================
# FIGURE 6A: PCA of lncRNA expression
# =========================================================

# -----------------------------
# Fix state levels
# -----------------------------
meta2$state <- trimws(as.character(meta2$state))

state_levels <- c(
  "Mito-High / Glyco-Low (OXPHOS dominant)",
  "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)",
  "Both-Low (Low metabolic)"
)

meta2$state <- factor(meta2$state, levels = state_levels)

print(table(meta2$state, useNA = "ifany"))

# -----------------------------
# Clean lncRNA matrix before PCA - safer version
# -----------------------------

lnc_mat2[is.infinite(lnc_mat2)] <- NA

# Keep lncRNAs with at least 80% non-missing samples and non-zero variance
keep_lnc_pca <- rowMeans(!is.na(lnc_mat2)) >= 0.80
keep_lnc_pca <- keep_lnc_pca & apply(lnc_mat2, 1, var, na.rm = TRUE) > 0

lnc_mat_pca <- lnc_mat2[keep_lnc_pca, , drop = FALSE]

cat("lncRNAs kept for PCA:", nrow(lnc_mat_pca), "\n")

# Replace remaining NA with row median
for (i in seq_len(nrow(lnc_mat_pca))) {
  row_median <- median(lnc_mat_pca[i, ], na.rm = TRUE)
  lnc_mat_pca[i, is.na(lnc_mat_pca[i, ])] <- row_median
}

# Select top variable lncRNAs
n_top <- min(1000, nrow(lnc_mat_pca))

top_var_lnc <- names(sort(
  apply(lnc_mat_pca, 1, var, na.rm = TRUE),
  decreasing = TRUE
))[1:n_top]

pca <- prcomp(t(lnc_mat_pca[top_var_lnc, , drop = FALSE]), scale. = TRUE)

# Percent variance explained
var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pc1_lab <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2], 1), "%)")

pca_df <- data.frame(
  sample15 = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  dplyr::left_join(meta2, by = "sample15")

# Make display label with OXPHOS-dominant
pca_df$state_label <- dplyr::recode(
  as.character(pca_df$state),
  "Mito-High / Glyco-Low (OXPHOS dominant)" = "Mito-High / Glyco-Low (OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)" = "Mito-Low / Glyco-High (Glycolytic)",
  "Both-High (Hybrid)" = "Both-High (Hybrid)",
  "Both-Low (Low metabolic)" = "Both-Low (Low)"
)

pca_df$state_label <- factor(
  pca_df$state_label,
  levels = c(
    "Mito-High / Glyco-Low (OXPHOS-dominant)",
    "Mito-Low / Glyco-High (Glycolytic)",
    "Both-High (Hybrid)",
    "Both-Low (Low)"
  )
)
p6A <- ggplot(pca_df, aes(PC1, PC2, color = state_label)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_color_manual(values = c(
    "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
    "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
    "Both-High (Hybrid)"                      = "purple",
    "Both-Low (Low)"                          = "maroon"
  )) +
  theme_classic(base_size = 12) +
  labs(
    title = "lncRNA Expression Landscape Across Metabolic States",
    x = pc1_lab,
    y = pc2_lab,
    color = "Metabolic state"
  )

ggsave(file.path(out_dir, "Figure6A_lncRNA_PCA.pdf"), p6A, width = 7, height = 5)
ggsave(file.path(out_dir, "Figure6A_lncRNA_PCA.jpg"), p6A, width = 7, height = 5, dpi = 400)



# =========================================================
# FIGURE 6B: Differential lncRNA expression
# =========================================================

stopifnot(ncol(lnc_mat2) == nrow(meta2))

ox_idx <- which(as.character(meta2$state) == "Mito-High / Glyco-Low (OXPHOS dominant)")
gly_idx <- which(as.character(meta2$state) == "Mito-Low / Glyco-High (Glycolytic)")

cat("OXPHOS samples:", length(ox_idx), "\n")
cat("Glycolytic samples:", length(gly_idx), "\n")

diff_lnc <- lapply(seq_len(nrow(lnc_mat2)), function(i) {
  
  x <- as.numeric(lnc_mat2[i, ox_idx])
  y <- as.numeric(lnc_mat2[i, gly_idx])
  
  x <- x[!is.na(x) & is.finite(x)]
  y <- y[!is.na(y) & is.finite(y)]
  
  p <- tryCatch(
    wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA
  )
  
  data.frame(
    lncRNA = rownames(lnc_mat2)[i],
    mean_OXPHOS = mean(x, na.rm = TRUE),
    mean_Glycolytic = mean(y, na.rm = TRUE),
    effect = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    p_value = p
  )
}) %>% dplyr::bind_rows()

diff_lnc <- diff_lnc %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    neglog10 = -log10(FDR),
    regulation = dplyr::case_when(
      FDR < 0.05 & effect > 0 ~ "OXPHOS-enriched",
      FDR < 0.05 & effect < 0 ~ "Glycolytic-enriched",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  diff_lnc,
  file.path(out_dir, "Supplementary_Table_lncRNA_Differential_Expression.csv"),
  row.names = FALSE
)

diff_lnc <- diff_lnc %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    neglog10 = -log10(FDR),
    regulation = dplyr::case_when(
      FDR < 0.05 & effect > 0 ~ "Mito-High / Glyco-Low (OXPHOS-dominant)",
      FDR < 0.05 & effect < 0 ~ "Mito-Low / Glyco-High (Glycolytic)",
      TRUE ~ "Not significant"
    )
  )

top_labels <- diff_lnc %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::group_by(regulation) %>%
  dplyr::arrange(FDR, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

p6B <- ggplot(diff_lnc, aes(x = effect, y = neglog10, color = regulation)) +
  geom_point(alpha = 0.7, size = 1.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = lncRNA, color = regulation),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
    "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
    "Not significant"                         = "gray50"
  )) +
  theme_classic(base_size = 12) +
  labs(
    title = "Differential lncRNA Expression Between OXPHOS-dominant and Glycolytic States",
    x = "Mean difference: OXPHOS-dominant - Glycolytic",
    y = "-log10(FDR)",
    color = NULL
  )

ggsave(file.path(out_dir, "Figure6B_lncRNA_Volcano.pdf"), p6B, width = 8, height = 6)
ggsave(file.path(out_dir, "Figure6B_lncRNA_Volcano.jpg"), p6B, width = 8, height = 6, dpi = 400)
# =========================================================
# Figure 6C: ggplot heatmap of top lncRNAs by metabolic state
# =========================================================
top_lnc <- diff_lnc %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice_head(n = 50) %>%
  pull(lncRNA)

lnc_state_avg <- sapply(state_levels, function(st) {
  sam <- meta2$sample15[as.character(meta2$state) == st]
  rowMeans(lnc_mat2[top_lnc, sam, drop = FALSE], na.rm = TRUE)
})

rownames(lnc_state_avg) <- top_lnc

colnames(lnc_state_avg) <- c(
  "Mito-High\nGlyco-Low\n(OXPHOS-dominant)",
  "Mito-Low\nGlyco-High\n(Glycolytic)",
  "Both-High\n(Hybrid)",
  "Both-Low\n(Low)"
)

lnc_state_z <- zscore_rows(lnc_state_avg)

lnc_df <- as.data.frame(lnc_state_z) %>%
  rownames_to_column("lncRNA") %>%
  pivot_longer(
    cols = -lncRNA,
    names_to = "Metabolic_state",
    values_to = "Row_z"
  )

lnc_df$Metabolic_state <- factor(
  lnc_df$Metabolic_state,
  levels = colnames(lnc_state_z)
)

lnc_df$lncRNA <- factor(
  lnc_df$lncRNA,
  levels = rev(rownames(lnc_state_z))
)

p6C <- ggplot(lnc_df, aes(x = Metabolic_state, y = lncRNA, fill = Row_z)) +
  geom_tile(color = NA) +
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
    title = "Top Metabolic-State-Associated lncRNAs",
    x = NULL,
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, lineheight = 0.85),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  file.path(out_dir, "Figure6C_lncRNA_State_Heatmap.pdf"),
  p6C, width = 7.2, height = 7.8
)

ggsave(
  file.path(out_dir, "Figure6C_lncRNA_State_Heatmap.jpg"),
  p6C, width = 7.2, height = 7.8, dpi = 400
)
# =========================================================
# FIGURE 6D: lncRNA Cox survival analysis
# =========================================================

surv <- fread(surv_file)
setnames(surv, old = colnames(surv)[1], new = "sample")
surv <- surv %>% mutate(sample15 = to_tcga15(sample))

cox_meta <- meta2 %>%
  left_join(surv %>% dplyr::select(sample15, OS, OS.time), by = "sample15") %>%
  filter(!is.na(OS), !is.na(OS.time), OS.time > 0)

cox_lnc <- diff_lnc %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice_head(n = 100) %>%
  pull(lncRNA)

cox_lnc <- cox_lnc[cox_lnc %in% rownames(lnc_mat2)]

cox_results <- lapply(cox_lnc, function(g) {
  dat <- cox_meta %>%
    filter(sample15 %in% colnames(lnc_mat2))
  
  dat$expr <- as.numeric(lnc_mat2[g, dat$sample15])
  
  fit <- tryCatch(
    coxph(Surv(OS.time, OS) ~ scale(expr) + strata(cancer_code), data = dat),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  s <- summary(fit)
  
  data.frame(
    lncRNA = g,
    HR = s$coefficients[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"],
    p_value = s$coefficients[1, "Pr(>|z|)"],
    n = nrow(dat),
    events = sum(dat$OS == 1, na.rm = TRUE)
  )
}) %>% bind_rows()

cox_results <- cox_results %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    direction = ifelse(HR > 1, "Risk", "Protective"),
    HR_label = paste0(round(HR, 2), " (", round(lower95, 2), "-", round(upper95, 2), ")")
  ) %>%
  arrange(FDR)

write.csv(
  cox_results,
  file.path(out_dir, "Supplementary_Table_lncRNA_Cox_Survival.csv"),
  row.names = FALSE
)

plot_cox <- cox_results %>%
  filter(!is.na(FDR)) %>%
  arrange(FDR) %>%
  slice_head(n = 15) %>%
  mutate(lncRNA = factor(lncRNA, levels = rev(lncRNA)))

p6D <- ggplot(plot_cox, aes(x = HR, y = lncRNA, color = direction)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_classic(base_size = 12) +
  labs(
    title = "Prognostic lncRNAs Associated With Metabolic States",
    x = "Hazard ratio",
    y = NULL,
    color = NULL
  )

ggsave(file.path(out_dir, "Figure6D_lncRNA_Cox_Forest.pdf"), p6D, width = 7, height = 5)
ggsave(file.path(out_dir, "Figure6D_lncRNA_Cox_Forest.jpg"), p6D, width = 7, height = 5, dpi = 400)

cat("\nDone. lncRNA Figure 6 saved in:\n", out_dir, "\n")

# =========================================================
# FIGURE 6E-F + SUPPLEMENT:
# =========================================================

# -----------------------------
# 0. Required objects
# -----------------------------
# Required from previous analysis:
# lnc_mat2        = lncRNA expression matrix, rows = lncRNAs, cols = sample15
# mirna_mat_corr  = miRNA expression matrix, rows = miRNAs, cols = sample15
# tcga_mat_symbol = mRNA expression matrix, rows = gene symbols, cols = sample15
# meta2           = sample metadata with sample15, state, cancer_code
# diff_lnc        = lncRNA differential results
# diff_mirna      = miRNA differential results from Figure 5

# -----------------------------
# 1. Harmonize samples
# -----------------------------

common_samples <- Reduce(intersect, list(
  colnames(lnc_mat2),
  colnames(mirna_mat_corr),
  colnames(tcga_mat_symbol),
  meta2$sample15
))

cat("Common matched samples:", length(common_samples), "\n")

lnc_net_mat   <- lnc_mat2[, common_samples, drop = FALSE]
mirna_net_mat <- mirna_mat_corr[, common_samples, drop = FALSE]
mrna_net_mat  <- tcga_mat_symbol[, common_samples, drop = FALSE]

meta_net <- meta2 %>%
  filter(sample15 %in% common_samples) %>%
  arrange(match(sample15, common_samples))

lnc_net_mat   <- lnc_net_mat[, meta_net$sample15, drop = FALSE]
mirna_net_mat <- mirna_net_mat[, meta_net$sample15, drop = FALSE]
mrna_net_mat  <- mrna_net_mat[, meta_net$sample15, drop = FALSE]

# -----------------------------
# 2. Select features
# -----------------------------

metabolic_genes <- c(
  "GPX4", "SLC7A11", "ACSL4", "AIFM2",
  "VDAC1", "VDAC2", "SLC25A4",
  "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1",
  "MFN1", "MFN2", "DNM1L", "FIS1",
  "HK2", "LDHA", "PKM", "SLC2A1"
)

metabolic_genes <- intersect(metabolic_genes, rownames(mrna_net_mat))

top_lnc_net <- diff_lnc %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  slice_head(n = 20) %>%
  pull(lncRNA)

top_lnc_net <- intersect(top_lnc_net, rownames(lnc_net_mat))

# Fix miRNA adjusted p-value column
if ("FDR" %in% colnames(diff_mirna)) {
  diff_mirna$miRNA_FDR <- diff_mirna$FDR
} else if ("p_adj" %in% colnames(diff_mirna)) {
  diff_mirna$miRNA_FDR <- diff_mirna$p_adj
} else if ("padj" %in% colnames(diff_mirna)) {
  diff_mirna$miRNA_FDR <- diff_mirna$padj
} else {
  stop("No adjusted p-value column found in diff_mirna")
}

top_mirna_net <- diff_mirna %>%
  filter(miRNA_FDR < 0.05) %>%
  arrange(miRNA_FDR) %>%
  slice_head(n = 20) %>%
  pull(miRNA)

top_mirna_net <- intersect(top_mirna_net, rownames(mirna_net_mat))

cat("lncRNAs selected:", length(top_lnc_net), "\n")
cat("miRNAs selected:", length(top_mirna_net), "\n")
cat("mRNAs selected:", length(metabolic_genes), "\n")

# -----------------------------
# 3. Correlation helper
# -----------------------------

cor_pairs <- function(mat1, mat2, type1, type2, samples, min_abs_r = 0.25, p_cut = 0.05) {
  
  res <- list()
  k <- 1
  
  for (i in rownames(mat1)) {
    for (j in rownames(mat2)) {
      
      x <- as.numeric(mat1[i, samples])
      y <- as.numeric(mat2[j, samples])
      
      keep <- is.finite(x) & is.finite(y)
      
      if (sum(keep) < 30) next
      
      ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
      
      res[[k]] <- data.frame(
        source = i,
        target = j,
        source_type = type1,
        target_type = type2,
        rho = as.numeric(ct$estimate),
        p_value = ct$p.value,
        n = sum(keep)
      )
      
      k <- k + 1
    }
  }
  
  out <- bind_rows(res)
  
  if (nrow(out) == 0) return(out)
  
  out <- out %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      direction = ifelse(rho > 0, "Positive", "Negative"),
      abs_rho = abs(rho)
    ) %>%
    filter(abs_rho >= min_abs_r, FDR < p_cut)
  
  out
}

# =========================================================
# FIGURE 6E: Integrated miRNA-lncRNA-mRNA network
# =========================================================

edges_mirna_mrna <- cor_pairs(
  mirna_net_mat[top_mirna_net, , drop = FALSE],
  mrna_net_mat[metabolic_genes, , drop = FALSE],
  "miRNA", "mRNA",
  common_samples,
  min_abs_r = 0.25,
  p_cut = 0.05
) %>%
  filter(rho < 0)   # keep negative miRNA-mRNA links

edges_lnc_mrna <- cor_pairs(
  lnc_net_mat[top_lnc_net, , drop = FALSE],
  mrna_net_mat[metabolic_genes, , drop = FALSE],
  "lncRNA", "mRNA",
  common_samples,
  min_abs_r = 0.25,
  p_cut = 0.05
)

edges_mirna_lnc <- cor_pairs(
  mirna_net_mat[top_mirna_net, , drop = FALSE],
  lnc_net_mat[top_lnc_net, , drop = FALSE],
  "miRNA", "lncRNA",
  common_samples,
  min_abs_r = 0.25,
  p_cut = 0.05
) %>%
  filter(rho < 0)   # possible competing/sponging-like pattern

edges_all <- bind_rows(
  edges_mirna_mrna,
  edges_lnc_mrna,
  edges_mirna_lnc
) %>%
  arrange(desc(abs_rho)) %>%
  slice_head(n = 120)

write.csv(
  edges_all,
  file.path(out_dir, "Figure6E_Integrated_miRNA_lncRNA_mRNA_Edges.csv"),
  row.names = FALSE
)

nodes <- data.frame(
  name = unique(c(edges_all$source, edges_all$target))
) %>%
  mutate(
    type = case_when(
      name %in% top_mirna_net ~ "miRNA",
      name %in% top_lnc_net ~ "lncRNA",
      name %in% metabolic_genes ~ "mRNA",
      TRUE ~ "Other"
    )
  )

g <- graph_from_data_frame(
  d = edges_all %>% dplyr::select(source, target, rho, direction, abs_rho),
  vertices = nodes,
  directed = FALSE
)
p6E <- ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(width = abs_rho, color = direction),
    alpha = 0.65
  ) +
  geom_node_point(
    aes(color = type),
    size = 4
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_edge_color_manual(values = c(
    "Positive" = "#c53030",
    "Negative" = "#2b6cb0"
  )) +
  scale_color_manual(values = c(
    "lncRNA" = "#6a3d9a",
    "miRNA" = "#1b9e77",
    "mRNA" = "#e31a1c",
    "Other" = "gray50"
  )) +
  scale_edge_width(range = c(0.3, 1.8)) +
  theme_void(base_size = 12) +
  labs(
    title = "Integrated miRNA-lncRNA-mRNA Regulatory Network",
    color = "Node type",
    edge_color = "Correlation",
    edge_width = "|Spearman rho|"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "Figure6E_Integrated_Regulatory_Network.pdf"),
       p6E, width = 9, height = 7)

ggsave(file.path(out_dir, "Figure6E_Integrated_Regulatory_Network.jpg"),
       p6E, width = 9, height = 7, dpi = 400)

# =========================================================
# FIGURE 6F: State-specific regulatory rewiring bubble plot
# =========================================================

selected_regulators <- unique(c(
  edges_all$source[edges_all$source_type %in% c("miRNA", "lncRNA")],
  edges_all$target[edges_all$target_type %in% c("miRNA", "lncRNA")]
))

selected_regulators <- selected_regulators[
  selected_regulators %in% c(rownames(mirna_net_mat), rownames(lnc_net_mat))
]

selected_genes <- unique(c(edges_all$source, edges_all$target))
selected_genes <- intersect(selected_genes, metabolic_genes)

state_corr_list <- list()
k <- 1

for (st in unique(as.character(meta_net$state))) {
  
  sam_st <- meta_net$sample15[as.character(meta_net$state) == st]
  
  if (length(sam_st) < 30) next
  
  for (reg in selected_regulators) {
    
    if (reg %in% rownames(mirna_net_mat)) {
      reg_vec <- as.numeric(mirna_net_mat[reg, sam_st])
      reg_type <- "miRNA"
    } else if (reg %in% rownames(lnc_net_mat)) {
      reg_vec <- as.numeric(lnc_net_mat[reg, sam_st])
      reg_type <- "lncRNA"
    } else {
      next
    }
    
    for (gene in selected_genes) {
      
      gene_vec <- as.numeric(mrna_net_mat[gene, sam_st])
      keep <- is.finite(reg_vec) & is.finite(gene_vec)
      
      if (sum(keep) < 20) next
      
      ct <- suppressWarnings(cor.test(reg_vec[keep], gene_vec[keep], method = "spearman"))
      
      state_corr_list[[k]] <- data.frame(
        state = st,
        regulator = reg,
        regulator_type = reg_type,
        gene = gene,
        rho = as.numeric(ct$estimate),
        p_value = ct$p.value,
        n = sum(keep)
      )
      
      k <- k + 1
    }
  }
}

state_corr <- bind_rows(state_corr_list) %>%
  group_by(state) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    direction = ifelse(rho > 0, "Positive", "Negative"),
    abs_rho = abs(rho)
  ) %>%
  filter(abs_rho >= 0.25, FDR < 0.05)

write.csv(
  state_corr,
  file.path(out_dir, "Figure6F_State_Specific_Regulatory_Rewiring.csv"),
  row.names = FALSE
)

# keep only most frequent regulators to avoid overcrowding
top_regs_plot <- state_corr %>%
  count(regulator, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(regulator)

plot_6F <- state_corr %>%
  filter(regulator %in% top_regs_plot) %>%
  mutate(
    state = factor(state, levels = state_levels),
    regulator = factor(regulator, levels = rev(top_regs_plot)),
    gene = factor(gene, levels = selected_genes)
  )

p6F <- ggplot(plot_6F, aes(x = gene, y = regulator)) +
  geom_point(aes(size = abs_rho, color = rho), alpha = 0.85) +
  facet_wrap(~ state, ncol = 2) +
  scale_color_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(range = c(1.5, 5), name = "|rho|") +
  theme_bw(base_size = 11) +
  labs(
    title = "State-Specific Rewiring of Noncoding RNA-Metabolic Gene Associations",
    x = "Metabolic gene",
    y = "miRNA / lncRNA regulator"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "Figure6F_State_Specific_Rewiring_BubblePlot.pdf"),
       p6F, width = 12, height = 9)

ggsave(file.path(out_dir, "Figure6F_State_Specific_Rewiring_BubblePlot.jpg"),
       p6F, width = 12, height = 9, dpi = 400)

# =========================================================
# SUPPLEMENT: Pan-cancer validation heatmap
# Cancer-type-specific conservation of top interactions
# =========================================================

top_edges_for_pan <- edges_all %>%
  arrange(desc(abs_rho)) %>%
  slice_head(n = 40) %>%
  mutate(edge_id = paste(source, target, sep = " ~ "))

pan_corr_list <- list()
k <- 1

for (cc in unique(meta_net$cancer_code)) {
  
  sam_cc <- meta_net$sample15[meta_net$cancer_code == cc]
  
  if (length(sam_cc) < 30) next
  
  for (i in seq_len(nrow(top_edges_for_pan))) {
    
    s <- top_edges_for_pan$source[i]
    t <- top_edges_for_pan$target[i]
    
    get_vec <- function(feature, samples) {
      if (feature %in% rownames(mirna_net_mat)) return(as.numeric(mirna_net_mat[feature, samples]))
      if (feature %in% rownames(lnc_net_mat)) return(as.numeric(lnc_net_mat[feature, samples]))
      if (feature %in% rownames(mrna_net_mat)) return(as.numeric(mrna_net_mat[feature, samples]))
      return(rep(NA_real_, length(samples)))
    }
    
    x <- get_vec(s, sam_cc)
    y <- get_vec(t, sam_cc)
    
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 25) next
    
    ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
    
    pan_corr_list[[k]] <- data.frame(
      cancer_code = cc,
      edge_id = top_edges_for_pan$edge_id[i],
      rho = as.numeric(ct$estimate),
      p_value = ct$p.value,
      n = sum(keep)
    )
    
    k <- k + 1
  }
}

pan_corr <- bind_rows(pan_corr_list) %>%
  group_by(edge_id) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(
  pan_corr,
  file.path(out_dir, "Supplementary_PanCancer_Regulatory_Edge_Validation.csv"),
  row.names = FALSE
)

p_pan <- ggplot(pan_corr, aes(x = cancer_code, y = edge_id, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    name = "Spearman\nrho"
  ) +
  theme_bw(base_size = 10) +
  labs(
    title = "Pan-Cancer Validation of Integrated Regulatory Interactions",
    x = "TCGA cancer type",
    y = "Regulatory interaction"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "Supplementary_PanCancer_Regulatory_Validation_Heatmap.pdf"),
       p_pan, width = 11, height = 10)

ggsave(file.path(out_dir, "Supplementary_PanCancer_Regulatory_Validation_Heatmap.jpg"),
       p_pan, width = 11, height = 10, dpi = 400)

# =========================================================
# FIGURE 6F: ceRNA lncRNA-miRNA-mRNA triplet analysis
# =========================================================
# -----------------------------
# 1. Install/load multiMiR
# -----------------------------

if (!requireNamespace("multiMiR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("multiMiR")
}

library(multiMiR)

# -----------------------------
# 2. First create the file
# -----------------------------

# load your edges file
edges_all <- fread(file.path(out_dir, "Figure6E_Integrated_miRNA_lncRNA_mRNA_Edges.csv"))

# check columns
colnames(edges_all)

lnc_mirna_db <- edges_all %>%
  dplyr::filter(
    source_type == "miRNA",
    target_type == "lncRNA",
    rho < 0   # negative correlation (IMPORTANT)
  ) %>%
  dplyr::select(
    lncRNA = target,
    miRNA = source
  ) %>%
  dplyr::distinct()

write.csv(
  lnc_mirna_db,
  file.path(base_dir, "lnc_mirna_interactions.csv"),
  row.names = FALSE
)

head(lnc_mirna_db)
nrow(lnc_mirna_db)

table(edges_all$source_type, edges_all$target_type)


# -----------------------------
# 2. lncRNA-miRNA interaction file
# -----------------------------
# You still need this file from starBase / miRcode:
# columns should be: lncRNA, miRNA

lnc_mirna_file <- file.path(base_dir, "lnc_mirna_interactions.csv")

lnc_mirna_db <- data.table::fread(lnc_mirna_file)

colnames(lnc_mirna_db)[1:2] <- c("lncRNA", "miRNA")

lnc_mirna_db <- lnc_mirna_db %>%
  dplyr::select(lncRNA, miRNA) %>%
  dplyr::distinct()

cat("lncRNA-miRNA pairs loaded:", nrow(lnc_mirna_db), "\n")

# -----------------------------
# 3. Create miRNA-mRNA targets from your existing edge file
# -----------------------------

mirna_mrna_db <- edges_all %>%
  dplyr::filter(
    source_type == "miRNA",
    target_type == "mRNA",
    rho < 0   # negative miRNA-mRNA correlation
  ) %>%
  dplyr::select(
    miRNA = source,
    mRNA = target
  ) %>%
  dplyr::distinct()

data.table::fwrite(
  mirna_mrna_db,
  file.path(base_dir, "mirna_mrna_targets_from_edges.csv")
)

cat("miRNA-mRNA pairs from existing edges:", nrow(mirna_mrna_db), "\n")
head(mirna_mrna_db)


# -----------------------------
# 4. Match samples across lncRNA, miRNA, mRNA
# Memory-safe version
# -----------------------------

# RNAs actually needed for ceRNA analysis
needed_lncRNAs <- unique(lnc_mirna_db$lncRNA)
needed_miRNAs  <- unique(c(lnc_mirna_db$miRNA, mirna_mrna_db$miRNA))
needed_mRNAs   <- unique(mirna_mrna_db$mRNA)

needed_lncRNAs <- intersect(needed_lncRNAs, rownames(lnc_mat2))
needed_miRNAs  <- intersect(needed_miRNAs, rownames(mirna_mat_corr))
needed_mRNAs   <- intersect(needed_mRNAs, rownames(tcga_mat_symbol))

cat("Needed lncRNAs:", length(needed_lncRNAs), "\n")
cat("Needed miRNAs:", length(needed_miRNAs), "\n")
cat("Needed mRNAs:", length(needed_mRNAs), "\n")

common_samples <- Reduce(intersect, list(
  colnames(lnc_mat2),
  colnames(mirna_mat_corr),
  colnames(tcga_mat_symbol),
  meta2$sample15
))

# memory-safe subsetting
lnc_net_mat   <- lnc_mat2[needed_lncRNAs, common_samples, drop = FALSE]
mirna_net_mat <- mirna_mat_corr[needed_miRNAs, common_samples, drop = FALSE]
mrna_net_mat  <- tcga_mat_symbol[needed_mRNAs, common_samples, drop = FALSE]

meta_net <- meta2 %>%
  dplyr::filter(sample15 %in% common_samples)

cat("Matched samples:", length(common_samples), "\n")

# -----------------------------
# 5. Select candidate lncRNAs
# -----------------------------

if (!"FDR" %in% colnames(diff_lnc)) {
  stop("diff_lnc must contain column named FDR")
}

top_lnc_net <- diff_lnc %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  dplyr::slice_head(n = 100) %>%
  dplyr::pull(lncRNA)

top_lnc_net <- intersect(top_lnc_net, rownames(lnc_net_mat))

cat("lncRNAs selected:", length(top_lnc_net), "\n")

# -----------------------------
# 6. Select candidate miRNAs
# -----------------------------

if ("FDR" %in% colnames(diff_mirna)) {
  diff_mirna$miRNA_FDR <- diff_mirna$FDR
} else if ("p_adj" %in% colnames(diff_mirna)) {
  diff_mirna$miRNA_FDR <- diff_mirna$p_adj
} else if ("padj" %in% colnames(diff_mirna)) {
  diff_mirna$miRNA_FDR <- diff_mirna$padj
} else {
  stop("No adjusted p-value column found in diff_mirna")
}

top_mirna_net <- diff_mirna %>%
  dplyr::filter(miRNA_FDR < 0.05) %>%
  dplyr::arrange(miRNA_FDR) %>%
  dplyr::slice_head(n = 100) %>%
  dplyr::pull(miRNA)

top_mirna_net <- intersect(top_mirna_net, rownames(mirna_net_mat))

cat("miRNAs selected:", length(top_mirna_net), "\n")

# -----------------------------
# 7. Select metabolic mRNAs
# -----------------------------

metabolic_genes <- c(
  "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1",
  "MFN1", "MFN2", "FIS1", "DNM1L", "TFAM",
  "GPX4", "SLC7A11", "ACSL4", "PRDX3", "PRDX5",
  "HK2", "LDHA", "PDK1", "PKM", "SLC2A1"
)

metabolic_genes <- intersect(metabolic_genes, rownames(mrna_net_mat))

cat("mRNAs selected:", length(metabolic_genes), "\n")

# -----------------------------
# 8. Harmonize miRNA names
# -----------------------------

lnc_mirna_db$miRNA <- gsub("HSA-", "hsa-", lnc_mirna_db$miRNA)
lnc_mirna_db$miRNA <- gsub("MIR", "miR", lnc_mirna_db$miRNA)

mirna_mrna_db$miRNA <- gsub("HSA-", "hsa-", mirna_mrna_db$miRNA)
mirna_mrna_db$miRNA <- gsub("MIR", "miR", mirna_mrna_db$miRNA)

# -----------------------------
# 9. Build candidate triplets
# -----------------------------

cat("lnc_mirna_db miRNAs:", length(unique(lnc_mirna_db$miRNA)), "\n")
cat("mirna_mrna_db miRNAs:", length(unique(mirna_mrna_db$miRNA)), "\n")
cat("Overlap miRNAs:", length(intersect(lnc_mirna_db$miRNA, mirna_mrna_db$miRNA)), "\n")

candidate_triplets <- lnc_mirna_db %>%
  dplyr::inner_join(mirna_mrna_db, by = "miRNA") %>%
  dplyr::filter(
    lncRNA %in% rownames(lnc_net_mat),
    miRNA %in% rownames(mirna_net_mat),
    mRNA %in% rownames(mrna_net_mat)
  ) %>%
  dplyr::distinct(lncRNA, miRNA, mRNA)

cat("Candidate triplets:", nrow(candidate_triplets), "\n")

if (nrow(candidate_triplets) == 0) {
  cat("No overlap between lncRNA-miRNA and miRNA-mRNA pairs.\n")
  cat("Check these:\n")
  print(intersect(lnc_mirna_db$miRNA, mirna_mrna_db$miRNA))
}

unique(mirna_mrna_db$miRNA)
unique(lnc_mirna_db$miRNA)
# -----------------------------
# 10. Correlation helper
# -----------------------------

get_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  
  if (sum(keep) < 30) {
    return(c(rho = NA, p = NA, n = sum(keep)))
  }
  
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
  
  c(
    rho = as.numeric(ct$estimate),
    p = ct$p.value,
    n = sum(keep)
  )
}

# -----------------------------
# 11. Test ceRNA correlation pattern
# -----------------------------

triplet_results <- candidate_triplets %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    lnc_mrna_rho = get_cor(
      as.numeric(lnc_net_mat[lncRNA, common_samples]),
      as.numeric(mrna_net_mat[mRNA, common_samples])
    )[1],
    lnc_mrna_p = get_cor(
      as.numeric(lnc_net_mat[lncRNA, common_samples]),
      as.numeric(mrna_net_mat[mRNA, common_samples])
    )[2],
    
    mir_lnc_rho = get_cor(
      as.numeric(mirna_net_mat[miRNA, common_samples]),
      as.numeric(lnc_net_mat[lncRNA, common_samples])
    )[1],
    mir_lnc_p = get_cor(
      as.numeric(mirna_net_mat[miRNA, common_samples]),
      as.numeric(lnc_net_mat[lncRNA, common_samples])
    )[2],
    
    mir_mrna_rho = get_cor(
      as.numeric(mirna_net_mat[miRNA, common_samples]),
      as.numeric(mrna_net_mat[mRNA, common_samples])
    )[1],
    mir_mrna_p = get_cor(
      as.numeric(mirna_net_mat[miRNA, common_samples]),
      as.numeric(mrna_net_mat[mRNA, common_samples])
    )[2]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    lnc_mrna_FDR = p.adjust(lnc_mrna_p, method = "BH"),
    mir_lnc_FDR  = p.adjust(mir_lnc_p, method = "BH"),
    mir_mrna_FDR = p.adjust(mir_mrna_p, method = "BH"),
    
    ceRNA_supported =
      lnc_mrna_rho > 0 &
      mir_lnc_rho < 0 &
      mir_mrna_rho < 0 &
      lnc_mrna_FDR < 0.05 &
      mir_lnc_FDR < 0.05 &
      mir_mrna_FDR < 0.05
  )

ceRNA_triplets <- triplet_results %>%
  dplyr::filter(ceRNA_supported == TRUE) %>%
  dplyr::arrange(lnc_mrna_FDR, mir_mrna_FDR, mir_lnc_FDR)

write.csv(
  triplet_results,
  file.path(out_dir, "Supplementary_Table_All_ceRNA_Triplets.csv"),
  row.names = FALSE
)

write.csv(
  ceRNA_triplets,
  file.path(out_dir, "Supplementary_Table_Significant_ceRNA_Triplets.csv"),
  row.names = FALSE
)

cat("Significant ceRNA triplets:", nrow(ceRNA_triplets), "\n")

if (nrow(ceRNA_triplets) == 0) {
  stop("No significant ceRNA triplets found under FDR < 0.05. Try using FDR < 0.10 or increasing top RNAs.")
}

# -----------------------------
# 12. Build network for top triplets
# -----------------------------

plot_triplets <- ceRNA_triplets %>%
  dplyr::slice_head(n = 50)

edges_ceRNA <- dplyr::bind_rows(
  plot_triplets %>%
    dplyr::transmute(
      source = lncRNA,
      target = miRNA,
      interaction = "lncRNA-miRNA",
      correlation = mir_lnc_rho
    ),
  plot_triplets %>%
    dplyr::transmute(
      source = miRNA,
      target = mRNA,
      interaction = "miRNA-mRNA",
      correlation = mir_mrna_rho
    ),
  plot_triplets %>%
    dplyr::transmute(
      source = lncRNA,
      target = mRNA,
      interaction = "lncRNA-mRNA",
      correlation = lnc_mrna_rho
    )
) %>%
  dplyr::distinct()

nodes_ceRNA <- data.frame(
  name = unique(c(edges_ceRNA$source, edges_ceRNA$target))
) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      name %in% rownames(lnc_net_mat) ~ "lncRNA",
      name %in% rownames(mirna_net_mat) ~ "miRNA",
      name %in% rownames(mrna_net_mat) ~ "mRNA",
      TRUE ~ "Other"
    )
  )

g_ceRNA <- igraph::graph_from_data_frame(
  d = edges_ceRNA %>%
    dplyr::select(source, target, interaction, correlation),
  vertices = nodes_ceRNA,
  directed = FALSE
)

p6E <- ggraph::ggraph(g_ceRNA, layout = "fr") +
  ggraph::geom_edge_link(
    aes(width = abs(correlation), color = interaction),
    alpha = 0.7
  ) +
  ggraph::geom_node_point(
    aes(color = type),
    size = 4
  ) +
  ggraph::geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_color_manual(values = c(
    "lncRNA" = "#6a3d9a",
    "miRNA" = "#ff9900",
    "mRNA" = "#1b9e77",
    "Other" = "gray60"
  )) +
  scale_edge_color_manual(values = c(
    "lncRNA-miRNA" = "#2b6cb0",
    "miRNA-mRNA" = "magenta",
    "lncRNA-mRNA" = "#c53030"
  )) +
  scale_edge_width(range = c(0.3, 1.6)) +
  theme_void() +
  labs(
    title = "ceRNA lncRNA-miRNA-mRNA Network",
    color = "Node type",
    edge_color = "Interaction",
    edge_width = "|Correlation|"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

ggsave(
  file.path(out_dir, "Figure6E_ceRNA_Network.pdf"),
  p6E,
  width = 8,
  height = 6
)

ggsave(
  file.path(out_dir, "Figure6E_ceRNA_Network.jpg"),
  p6E,
  width = 8,
  height = 6,
  dpi = 400
)

cat("Figure 6E ceRNA network saved.\n")

# =========================================================
# SUPPLEMENTARY ceRNA NETWORK: Expanded lncRNA-miRNA-mRNA network
# =========================================================

# -----------------------------
# 1. Expanded lncRNA-miRNA edges
# -----------------------------

lnc_mirna_db_ext <- edges_all %>%
  dplyr::filter(
    source_type == "miRNA",
    target_type == "lncRNA",
    rho < 0,
    FDR < 0.10
  ) %>%
  dplyr::select(
    lncRNA = target,
    miRNA = source
  ) %>%
  dplyr::distinct()

# -----------------------------
# 2. Expanded miRNA-mRNA edges
# -----------------------------

mirna_mrna_db_ext <- edges_all %>%
  dplyr::filter(
    source_type == "miRNA",
    target_type == "mRNA",
    rho < 0,
    FDR < 0.10
  ) %>%
  dplyr::select(
    miRNA = source,
    mRNA = target
  ) %>%
  dplyr::distinct()

# -----------------------------
# 3. Expanded lncRNA-mRNA positive edges
# -----------------------------

lnc_mrna_db_ext <- edges_all %>%
  dplyr::filter(
    source_type == "lncRNA",
    target_type == "mRNA",
    rho > 0,
    FDR < 0.10
  ) %>%
  dplyr::select(
    lncRNA = source,
    mRNA = target,
    lnc_mrna_rho = rho,
    lnc_mrna_FDR = FDR
  ) %>%
  dplyr::distinct()

# -----------------------------
# 4. Build expanded triplets
# -----------------------------

candidate_triplets_ext <- lnc_mirna_db_ext %>%
  dplyr::inner_join(mirna_mrna_db_ext, by = "miRNA") %>%
  dplyr::inner_join(lnc_mrna_db_ext, by = c("lncRNA", "mRNA")) %>%
  dplyr::filter(
    lncRNA %in% rownames(lnc_mat2),
    miRNA %in% rownames(mirna_mat_corr),
    mRNA %in% rownames(tcga_mat_symbol)
  ) %>%
  dplyr::distinct(lncRNA, miRNA, mRNA, lnc_mrna_rho, lnc_mrna_FDR)

cat("Expanded triplets:", nrow(candidate_triplets_ext), "\n")
cat("Expanded lncRNAs:", length(unique(candidate_triplets_ext$lncRNA)), "\n")
cat("Expanded miRNAs:", length(unique(candidate_triplets_ext$miRNA)), "\n")
cat("Expanded mRNAs:", length(unique(candidate_triplets_ext$mRNA)), "\n")

write.csv(
  candidate_triplets_ext,
  file.path(out_dir, "Supplementary_Table_Expanded_ceRNA_Triplets.csv"),
  row.names = FALSE
)

# -----------------------------
# 5. Limit network size for plotting
# -----------------------------

plot_triplets_ext <- candidate_triplets_ext %>%
  dplyr::arrange(lnc_mrna_FDR, dplyr::desc(abs(lnc_mrna_rho))) %>%
  dplyr::slice_head(n = 100)

# -----------------------------
# 6. Build edges
# -----------------------------

edges_ext <- dplyr::bind_rows(
  plot_triplets_ext %>%
    dplyr::transmute(
      source = lncRNA,
      target = miRNA,
      interaction = "lncRNA-miRNA"
    ),
  plot_triplets_ext %>%
    dplyr::transmute(
      source = miRNA,
      target = mRNA,
      interaction = "miRNA-mRNA"
    ),
  plot_triplets_ext %>%
    dplyr::transmute(
      source = lncRNA,
      target = mRNA,
      interaction = "lncRNA-mRNA"
    )
) %>%
  dplyr::distinct()

# add correlation values from edges_all
edges_ext <- edges_ext %>%
  dplyr::left_join(
    edges_all %>%
      dplyr::select(source, target, rho, FDR),
    by = c("source", "target")
  ) %>%
  dplyr::mutate(
    abs_rho = abs(rho)
  )

# -----------------------------
# 7. Build nodes
# -----------------------------

nodes_ext <- data.frame(
  name = unique(c(edges_ext$source, edges_ext$target))
) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      name %in% rownames(lnc_mat2) ~ "lncRNA",
      name %in% rownames(mirna_mat_corr) ~ "miRNA",
      name %in% rownames(tcga_mat_symbol) ~ "mRNA",
      TRUE ~ "Other"
    )
  )

# -----------------------------
# 8. Plot expanded supplementary network
# -----------------------------

edges_ext <- edges_ext %>%
  dplyr::filter(!is.na(abs_rho))

p_ext <- ggraph::ggraph(g_ext, layout = "fr") +
  ggraph::geom_edge_link(
    aes(width = abs_rho, color = interaction),
    alpha = 0.6
  ) +
  ggraph::geom_node_point(
    aes(color = type),
    size = 3.5
  ) +
  ggraph::geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 2.8
  ) +
  scale_color_manual(values = c(
    "lncRNA" = "#6a3d9a",
    "miRNA" = "#ff9900",
    "mRNA" = "#1b9e77",
    "Other" = "gray60"
  )) +
  scale_edge_color_manual(values = c(
    "lncRNA-miRNA" = "#2b6cb0",
    "miRNA-mRNA" = "magenta",
    "lncRNA-mRNA" = "#c53030"
  )) +
  scale_edge_width(range = c(0.2, 1.4)) +
  theme_void() +
  labs(
    title = "Expanded ceRNA lncRNA-miRNA-mRNA Network",
    color = "Node type",
    edge_color = "Interaction",
    edge_width = "|Correlation|"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

ggsave(
  file.path(out_dir, "Supplementary_Figure_Expanded_ceRNA_Network.pdf"),
  p_ext,
  width = 10,
  height = 8
)

ggsave(
  file.path(out_dir, "Supplementary_Figure_Expanded_ceRNA_Network.jpg"),
  p_ext,
  width = 10,
  height = 8,
  dpi = 400
)

cat("Expanded supplementary ceRNA network saved.\n")
