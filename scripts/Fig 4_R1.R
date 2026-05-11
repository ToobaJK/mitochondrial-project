# =========================================================
# Figure 4. Clinical and molecular landscape of
# mitochondrial metabolic states
# Created By: TOOBA
# =========================================================

# -----------------------------
# SETTINGS
# -----------------------------
base_dir   <- "D:/UAEU/Figure 4"
output_dir <- file.path(base_dir, "outputR1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# LIBRARIES
# -----------------------------
libs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "stringr",
  "survival", "survminer", "ComplexHeatmap", "circlize",
  "tibble", "forcats", "grid"
)
to_install <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(forcats)
library(grid)

# -----------------------------
# HELPERS
# -----------------------------
save_both <- function(plot_obj, filename, width = 10, height = 7, dpi = 300) {
  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot = plot_obj, width = width, height = height, device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, paste0(filename, ".jpg")),
    plot = plot_obj, width = width, height = height, dpi = dpi
  )
}

save_heatmap_both <- function(ht, filename, width = 9, height = 6, res = 300) {
  pdf(file.path(output_dir, paste0(filename, ".pdf")), width = width, height = height)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  jpeg(file.path(output_dir, paste0(filename, ".jpg")),
       width = width * res, height = height * res, res = res)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}

find_file <- function(pattern, base_dir) {
  x <- list.files(base_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(x) == 0) stop(paste("File not found for pattern:", pattern))
  x[1]
}

clean_tcga_id <- function(x) {
  x <- gsub("\\.", "-", x)
  x <- substr(x, 1, 15)
  x
}

# -----------------------------
# FILES
# -----------------------------
expr_file      <- find_file("^tcga_RSEM_gene_tpm", base_dir)
pheno_file     <- find_file("^TCGA_phenotype_denseDataOnlyDownload.*\\.tsv$", base_dir)
immune_file    <- find_file("^Subtype_Immune_Model_Based", base_dir)
subtype_file   <- find_file("^TCGASubtype\\.20170308\\.tsv$", base_dir)
survival_file  <- find_file("^Survival_SupplementalTable_S1_20171025", base_dir)

cat("Expression file: ", expr_file, "\n")
cat("Phenotype file:  ", pheno_file, "\n")
cat("Immune file:     ", immune_file, "\n")
cat("Subtype file:    ", subtype_file, "\n")
cat("Survival file:   ", survival_file, "\n")

# -----------------------------
# GENE SETS
# -----------------------------

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

oxphos_genes <- c(
  "NDUFA1","NDUFB8","NDUFS1","SDHA","SDHB",
  "UQCRC1","UQCRC2","COX4I1","COX5A",
  "ATP5F1A","ATP5F1B"
)

ros_genes <- c(
  "TXN","TXN2","PRDX1","PRDX3","PRDX5",
  "GSR","GPX4","SOD2","CAT","NFE2L2"
)

hypoxia_genes <- c(
  "HIF1A","VEGFA","CA9","SLC2A1","LDHA","PGK1","ENO1","ALDOA"
)

fa_genes <- c(
  "CPT1A","CPT2","ACADM","ACADVL",
  "HADHA","HADHB","ACOX1"
)

redox_genes <- c(
  "SLC7A11","GPX4","AIFM2","FSP1",
  "GSR","TXN","TXN2","PRDX3","PRDX5",
  "NFE2L2","ACSL4"
)

# Cyt genes
cyt_genes <- c("GZMA", "PRF1")

# Combine ALL genes needed in pipeline
all_needed_genes <- unique(toupper(c(
  mito_genes,
  glyco_genes,
  oxphos_genes,
  ros_genes,
  hypoxia_genes,
  fa_genes,
  redox_genes,
  cyt_genes
)))

target_symbols <- all_needed_genes

# -----------------------------
# READ EXPRESSION (MEMORY-SAFE)
# -----------------------------
cat("Reading expression in memory-safe mode...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(AnnotationDbi)
library(org.Hs.eg.db)

# gene symbols you need
target_symbols <- all_needed_genes

# map SYMBOL -> ENSEMBL
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = target_symbols,
  keytype = "SYMBOL",
  columns = c("ENSEMBL", "SYMBOL")
)

# map SYMBOL -> ENSEMBL
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = target_symbols,
  keytype = "SYMBOL",
  columns = c("ENSEMBL", "SYMBOL")
)

symbol_map <- symbol_map %>%
  filter(!is.na(ENSEMBL), !is.na(SYMBOL)) %>%
  mutate(
    SYMBOL = toupper(SYMBOL),
    ENSEMBL_FULL = ENSEMBL
  ) %>%
  distinct(ENSEMBL_FULL, SYMBOL)

target_ensembl <- unique(symbol_map$ENSEMBL_FULL)

cat("Target symbols:", length(target_symbols), "\n")
cat("Mapped Ensembl IDs:", length(target_ensembl), "\n")

# -----------------------------
# READ HEADER ONLY
# -----------------------------
con_header <- gzfile(expr_file, open = "rt")
header_line <- readLines(con_header, n = 1)
close(con_header)

header <- strsplit(header_line, "\t")[[1]]
header <- strsplit(header_line, "\t")[[1]]

gene_col <- header[1]
sample_names <- header[-1]

cat("Total sample columns:", length(sample_names), "\n")
cat("Example sample names:\n")
print(head(sample_names, 10))

# keep primary tumor sample columns
primary_cols <- grepl("-01", sample_names)
sample_names_primary <- sample_names[primary_cols]

cat("Primary tumor sample columns:", length(sample_names_primary), "\n")

# positions in the original file line
keep_idx <- c(1, which(primary_cols) + 1)

# -----------------------------
# STREAM THE FILE IN CHUNKS
# -----------------------------
con <- file(expr_file, open = "r")
on.exit(close(con))

# skip header
readLines(con, n = 1)

chunk_size <- 50000
matched_list <- list()
chunk_id <- 1
total_matched <- 0

repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0) break
  
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  
  # extract Ensembl IDs without version
  ens_ids <- vapply(split_lines, function(x) sub("\\..*$", "", x[1]), character(1))
  keep_rows <- ens_ids %in% target_ensembl
  
  if (any(keep_rows)) {
    kept <- split_lines[keep_rows]
    
    chunk_df <- do.call(rbind, lapply(kept, function(x) {
      x <- x[keep_idx]
      length(x) <- length(keep_idx)
      x
    }))
    
    chunk_df <- as.data.frame(chunk_df, stringsAsFactors = FALSE)
    colnames(chunk_df) <- c("ensembl_id_raw", sample_names_primary)
    
    matched_list[[chunk_id]] <- chunk_df
    total_matched <- total_matched + nrow(chunk_df)
    chunk_id <- chunk_id + 1
  }
  
  if (chunk_id %% 5 == 0) {
    cat("Processed chunks:", chunk_id - 1, 
        "| matched rows so far:", total_matched, "\n")
  }
}

expr_small <- bind_rows(matched_list)

cat("Matched expression rows:", nrow(expr_small), "\n")

if (nrow(expr_small) == 0) {
  stop("No target genes were found in the expression file.")
}

# -----------------------------
# CLEAN GENE IDS AND MAP BACK TO SYMBOLS
# -----------------------------
expr_small <- expr_small %>%
  mutate(
    ENSEMBL_FULL = sub("\\..*$", "", ensembl_id_raw)
  ) %>%
  left_join(symbol_map, by = "ENSEMBL_FULL") %>%
  filter(!is.na(SYMBOL))

write.csv(
  expr_small,
  file.path(output_dir, "Figure4_TargetGene_Expression_RawExtract.csv"),
  row.names = FALSE
)

# -----------------------------
# BUILD EXPRESSION MATRIX
# -----------------------------
expr_num <- expr_small[, c("SYMBOL", sample_names_primary), drop = FALSE]

expr_num_mat <- as.matrix(expr_num[, -1, drop = FALSE])
mode(expr_num_mat) <- "numeric"

rownames(expr_num_mat) <- expr_num$SYMBOL
colnames(expr_num_mat) <- sample_names_primary

# collapse duplicated symbols by mean
expr_mat <- rowsum(expr_num_mat, group = rownames(expr_num_mat), reorder = FALSE) /
  as.vector(table(rownames(expr_num_mat)))

clean_tcga_id <- function(x) {
  x <- gsub("\\.", "-", x)
  substr(x, 1, 15)
}

colnames(expr_mat) <- clean_tcga_id(colnames(expr_mat))
expr_mat <- expr_mat[, !duplicated(colnames(expr_mat)), drop = FALSE]

cat("Final expression matrix dimensions:", dim(expr_mat)[1], "genes x", dim(expr_mat)[2], "samples\n")
print(head(rownames(expr_mat), 20))
print(head(colnames(expr_mat), 10))

cat("Final expression matrix dimensions:", dim(expr_mat)[1], "genes x", dim(expr_mat)[2], "samples\n")
# -----------------------------
# CALCULATE TCGA METABOLIC STATES

state_levels <- c(
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

score_program <- function(mat, genes) {
  genes <- intersect(toupper(genes), rownames(mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
}

tcga_states <- data.frame(
  sample = colnames(expr_mat),
  mito_score  = score_program(expr_mat, mito_genes),
  glyco_score = score_program(expr_mat, glyco_genes),
  stringsAsFactors = FALSE
) %>%
  mutate(
    mito_z  = as.numeric(scale(mito_score)),
    glyco_z = as.numeric(scale(glyco_score)),
    metabolic_state = case_when(
      mito_z >= 0 & glyco_z <  0 ~ "Mito-High / Glyco-Low (OXPHOS-dominant)",
      mito_z <  0 & glyco_z >= 0 ~ "Mito-Low / Glyco-High (Glycolytic)",
      mito_z >= 0 & glyco_z >= 0 ~ "Both-High (Hybrid)",
      mito_z <  0 & glyco_z <  0 ~ "Both-Low (Low metabolic)",
      TRUE ~ NA_character_
    ),
    metabolic_state = factor(metabolic_state, levels = state_levels)
  ) %>%
  filter(!is.na(metabolic_state))

write.csv(
  tcga_states,
  file.path(output_dir, "TCGA_Metabolic_States.csv"),
  row.names = FALSE
)

table(tcga_states$metabolic_state, useNA = "ifany")
# -----------------------------
# READ PHENOTYPE
# -----------------------------
cat("Reading phenotype...\n")
pheno <- fread(pheno_file, data.table = FALSE)
colnames(pheno) <- make.names(colnames(pheno))

# sample column
sample_col_pheno <- grep("^sample$|sample", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
pheno$sample <- clean_tcga_id(pheno[[sample_col_pheno]])

# cancer type column from your file
pheno$cancer_type <- pheno$X_primary_disease
if (!"cancer_type" %in% colnames(pheno) || all(is.na(pheno$cancer_type))) {
  pheno$cancer_type <- pheno$._primary_disease
}
if (all(is.na(pheno$cancer_type))) {
  pheno$cancer_type <- pheno$`_primary_disease`
}

# keep primary tumor where possible
if ("sample_type_id" %in% colnames(pheno)) {
  pheno <- pheno %>% filter(sample_type_id == 1)
} else {
  pheno <- pheno %>% filter(grepl("-01$", sample))
}

pheno2 <- pheno %>%
  dplyr::select(sample, cancer_type, dplyr::everything()) %>%
  distinct(sample, .keep_all = TRUE)

cat("Phenotype rows after filtering:", nrow(pheno2), "\n")
cat("Example cancer types:\n")
print(head(unique(pheno2$cancer_type), 20))

# -----------------------------
# READ IMMUNE SUBTYPE
# -----------------------------
cat("Reading immune subtype...\n")
immune <- fread(immune_file, data.table = FALSE)
colnames(immune) <- make.names(colnames(immune))

sample_col_immune <- grep("^sample$|sample", colnames(immune), ignore.case = TRUE, value = TRUE)[1]
immune$sample <- clean_tcga_id(immune[[sample_col_immune]])

immune_sub_col <- grep("Immune|Subtype", colnames(immune), ignore.case = TRUE, value = TRUE)
immune_sub_col <- immune_sub_col[immune_sub_col != "sample"][1]
immune$immune_subtype <- immune[[immune_sub_col]]

immune2 <- immune[, c("sample", "immune_subtype"), drop = FALSE]
immune2 <- immune2[!duplicated(immune2$sample), , drop = FALSE]

cat("Immune rows:", nrow(immune2), "\n")
print(head(immune2))

# -----------------------------
# READ SURVIVAL
# -----------------------------
cat("Reading survival...\n")
surv <- fread(survival_file, data.table = FALSE)
colnames(surv) <- make.names(colnames(surv))

sample_col_surv <- grep("sample|bcr_patient_barcode|sampleID", colnames(surv), ignore.case = TRUE, value = TRUE)[1]
surv$sample <- clean_tcga_id(surv[[sample_col_surv]])

# your file already has these columns
os_time_col <- "OS.time"
os_col <- "OS"

surv2 <- surv
surv2$OS.time <- as.numeric(surv2[[os_time_col]])
surv2$OS <- surv2[[os_col]]

if (is.character(surv2$OS)) {
  surv2$OS <- ifelse(grepl("DECEASED|DEAD|1", toupper(surv2$OS)), 1, 0)
} else {
  surv2$OS <- as.numeric(surv2$OS)
  surv2$OS <- ifelse(surv2$OS > 0, 1, 0)
}

surv2 <- surv2[, c("sample", "OS.time", "OS"), drop = FALSE]
surv2 <- surv2[!duplicated(surv2$sample), , drop = FALSE]

# merge cancer type from phenotype
pheno_tmp <- pheno2[, c("sample", "cancer_type"), drop = FALSE]
pheno_tmp <- pheno_tmp[!duplicated(pheno_tmp$sample), , drop = FALSE]

surv2 <- merge(surv2, pheno_tmp, by = "sample", all.x = TRUE)

write.csv(
  surv2,
  file.path(output_dir, "Figure4A_Survival_Input.csv"),
  row.names = FALSE
)

cat("Survival rows:", nrow(surv2), "\n")
print(head(surv2))

# -----------------------------
# MERGE MASTER TABLE
# -----------------------------
master <- tcga_states

pheno_tmp  <- pheno2[, c("sample", "cancer_type"), drop = FALSE]
immune_tmp <- immune2[, c("sample", "immune_subtype"), drop = FALSE]
surv_tmp   <- surv2[, c("sample", "OS.time", "OS"), drop = FALSE]

pheno_tmp  <- pheno_tmp[!duplicated(pheno_tmp$sample), , drop = FALSE]
immune_tmp <- immune_tmp[!duplicated(immune_tmp$sample), , drop = FALSE]
surv_tmp   <- surv_tmp[!duplicated(surv_tmp$sample), , drop = FALSE]

master <- merge(master, pheno_tmp,  by = "sample", all.x = TRUE)
master <- merge(master, immune_tmp, by = "sample", all.x = TRUE)
master <- merge(master, surv_tmp,   by = "sample", all.x = TRUE)

write.csv(
  master,
  file.path(output_dir, "Figure4_Master_Merged_Table.csv"),
  row.names = FALSE
)

cat("Master rows:", nrow(master), "\n")
print(head(master))



# apply everywhere
master$metabolic_state <- factor(master$metabolic_state, levels = state_levels)
tcga_states$metabolic_state <- factor(tcga_states$metabolic_state, levels = state_levels)


# =========================================================
# FIGURE 4A
# Overall survival by metabolic state with HR + 95% CI
# =========================================================

fig4a_dat <- master %>%
  filter(!is.na(OS.time), !is.na(OS), OS.time > 0, !is.na(metabolic_state)) %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels),
    metabolic_state_cox = relevel(metabolic_state, ref = "Both-Low (Low metabolic)")
  )

print(table(fig4a_dat$metabolic_state, useNA = "ifany"))

fit_os <- survfit(Surv(OS.time, OS) ~ metabolic_state, data = fig4a_dat)

cox_os <- coxph(Surv(OS.time, OS) ~ metabolic_state_cox, data = fig4a_dat)
cox_sum <- summary(cox_os)

hr_df <- data.frame(
  State = gsub("metabolic_state_cox", "", rownames(cox_sum$coefficients)),
  HR = round(cox_sum$conf.int[, "exp(coef)"], 2),
  lower = round(cox_sum$conf.int[, "lower .95"], 2),
  upper = round(cox_sum$conf.int[, "upper .95"], 2),
  p = signif(cox_sum$coefficients[, "Pr(>|z|)"], 3)
)

hr_df <- rbind(
  data.frame(
    State = "Both-Low (Low metabolic)",
    HR = 1.00,
    lower = NA,
    upper = NA,
    p = NA
  ),
  hr_df
)

hr_df$State <- factor(hr_df$State, levels = state_levels)
hr_df <- hr_df %>% arrange(State)
hr_df$State <- as.character(hr_df$State)

hr_df$State_short <- c(
  "OXPHOS-dominant",
  "Glycolytic",
  "Hybrid",
  "Low metabolic"
)

hr_df$HR_95CI <- ifelse(
  hr_df$State == "Both-Low (Low metabolic)",
  "1.00 (reference)",
  paste0(hr_df$HR, " (", hr_df$lower, "-", hr_df$upper, ")")
)

hr_text_inside <- paste0(
  "Cox Hazard Ratio (reference: Low metabolic)\n",
  paste0(hr_df$State_short, ": ", hr_df$HR_95CI, collapse = "\n")
)

p4a <- ggsurvplot(
  fit_os,
  data = fig4a_dat,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  size = 1.15,
  palette = unname(state_colors[state_levels]),
  legend.labs = state_levels,
  ggtheme = theme_bw(base_size = 14),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "Metabolic State"
)

p4a$plot <- p4a$plot +
  annotate(
    "label",
    x = max(fig4a_dat$OS.time, na.rm = TRUE) * 0.35,
    y = 0.93,
    label = hr_text_inside,
    hjust = 0,
    vjust = 1,
    size = 3.1,
    label.size = 0.25,
    fill = "white",
    alpha = 0.9
  ) +
  labs(title = "Overall Survival Across Mitochondrial Metabolic States") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(output_dir, "Figure4A_OS_KM_HR.pdf"),
  plot = p4a$plot,
  width = 11.5,
  height = 6,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "Figure4A_OS_KM_HR.jpg"),
  plot = p4a$plot,
  width = 11.5,
  height = 6,
  dpi = 300
)

write.csv(
  hr_df,
  file.path(output_dir, "Figure4A_Cox_HR_Table.csv"),
  row.names = FALSE
)

print(hr_df)
# =========================================================
# FIGURE Supplementary 4A
# Immune subtype distribution across metabolic states
# =========================================================
fig4b_dat <- master %>%
  filter(!is.na(immune_subtype), !is.na(metabolic_state)) %>%
  count(metabolic_state, immune_subtype) %>%
  group_by(metabolic_state) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

write.csv(fig4b_dat,
          file.path(output_dir, "Figure4B_Immune_Subtype_Composition.csv"),
          row.names = FALSE)


# Custom x-axis labels (multi-line like your image)
state_labels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
  "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
)

p4b <- ggplot(fig4b_dat, aes(x = metabolic_state, y = frac, fill = immune_subtype)) +
  geom_bar(stat = "identity", position = "fill", color = "white") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = state_labels) +
  theme_bw(base_size = 12) +
  labs(
    title = "Immune Subtype Composition Across Mitochondrial Metabolic States",
    x = NULL,
    y = "Fraction of samples",
    fill = "Immune subtype"
  ) +
  theme(
    axis.text.x = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_both(p4b, "FigureS4A_ImmuneSubtype_StackedBar", width = 11, height = 7)

# =========================================================
# FIGURE 4B
# Hypoxia program across mitochondrial metabolic states
# =========================================================

library(ggpubr)

hypoxia_sig <- c("HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "PGK1", "ENO1", "ALDOA")

fig4b_dat <- data.frame(
  sample = colnames(expr_mat),
  HypoxiaScore = score_program(expr_mat, hypoxia_sig),
  stringsAsFactors = FALSE
)

fig4b_dat <- merge(
  fig4b_dat,
  tcga_states[, c("sample", "metabolic_state")],
  by = "sample",
  all.x = FALSE
)

fig4b_dat <- fig4b_dat %>%
  filter(!is.na(HypoxiaScore), !is.na(metabolic_state)) %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels)
  )

write.csv(
  fig4b_dat,
  file.path(output_dir, "Figure4B_Hypoxia_Score.csv"),
  row.names = FALSE
)

# Overall Kruskal-Wallis test
kw_hypoxia <- kruskal.test(HypoxiaScore ~ metabolic_state, data = fig4b_dat)
print(kw_hypoxia)

# Format global p-value label
kw_label <- ifelse(
  kw_hypoxia$p.value < 2.2e-16,
  "Kruskal-Wallis p < 2e-16",
  paste0("Kruskal-Wallis p = ", signif(kw_hypoxia$p.value, 3))
)

# Pairwise Wilcoxon test, BH corrected
pairwise_hypoxia <- pairwise.wilcox.test(
  fig4b_dat$HypoxiaScore,
  fig4b_dat$metabolic_state,
  p.adjust.method = "BH"
)

print(pairwise_hypoxia)

write.csv(
  as.data.frame(pairwise_hypoxia$p.value),
  file.path(output_dir, "Figure4B_Hypoxia_Pairwise_Wilcox_BH.csv")
)

# Pairwise comparisons using Both-Low as reference
my_comparisons_hypoxia <- list(
  c("Both-Low (Low metabolic)", "Mito-High / Glyco-Low (OXPHOS-dominant)"),
  c("Both-Low (Low metabolic)", "Mito-Low / Glyco-High (Glycolytic)"),
  c("Both-Low (Low metabolic)", "Both-High (Hybrid)")
)

state_labels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
  "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
)

p4b <- ggplot(fig4b_dat, aes(x = metabolic_state, y = HypoxiaScore, fill = metabolic_state)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85, linewidth = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.18, size = 0.45) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  stat_compare_means(
    comparisons = my_comparisons_hypoxia,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.08,
    size = 4
  ) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    title = "Hypoxia Program Across Mitochondrial Metabolic States",
    subtitle = kw_label,
    x = NULL,
    y = "Hypoxia signature score"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p4b)

save_both(
  p4b,
  "Figure4B_HypoxiaScore_Boxplot_with_BothLow_Pvalues",
  width = 8,
  height = 6
)
# =========================================================
# FIGURE 4C
# Cytolytic activity across metabolic states
# =========================================================

library(ggpubr)

cyt_genes <- c("GZMA", "PRF1")

score_program_safe <- function(mat, genes) {
  mat_genes <- toupper(rownames(mat))
  genes <- toupper(genes)
  matched_rows <- which(mat_genes %in% genes)
  
  cat("Matched genes:", rownames(mat)[matched_rows], "\n")
  
  if (length(matched_rows) == 0) {
    warning("No CYT genes found in expression matrix.")
    return(rep(NA_real_, ncol(mat)))
  }
  
  colMeans(mat[matched_rows, , drop = FALSE], na.rm = TRUE)
}

fig4c_dat <- data.frame(
  sample = colnames(expr_mat),
  CYT = score_program_safe(expr_mat, cyt_genes),
  stringsAsFactors = FALSE
)

fig4c_dat$sample <- clean_tcga_id(fig4c_dat$sample)
tcga_states$sample <- clean_tcga_id(tcga_states$sample)

fig4c_dat <- merge(
  fig4c_dat,
  tcga_states[, c("sample", "metabolic_state")],
  by = "sample",
  all.x = FALSE
)

fig4c_dat <- fig4c_dat %>%
  mutate(
    metabolic_state = factor(as.character(metabolic_state), levels = state_levels)
  ) %>%
  filter(!is.na(CYT), !is.na(metabolic_state))

cat("Rows in Figure 4C:", nrow(fig4c_dat), "\n")
print(table(fig4c_dat$metabolic_state, useNA = "ifany"))
summary(fig4c_dat$CYT)

# Overall statistics
kw_cyt <- kruskal.test(CYT ~ metabolic_state, data = fig4c_dat)

# Pairwise Wilcoxon, BH corrected
pairwise_cyt <- pairwise.wilcox.test(
  fig4c_dat$CYT,
  fig4c_dat$metabolic_state,
  p.adjust.method = "BH"
)

print(kw_cyt)
print(pairwise_cyt)

write.csv(
  fig4c_dat,
  file.path(output_dir, "Figure4C_CYT_Score.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(pairwise_cyt$p.value),
  file.path(output_dir, "Figure4C_CYT_Pairwise_Wilcox_BH.csv")
)

state_labels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
  "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
)

# Pairwise comparisons using Both-Low as reference
my_comparisons_cyt <- list(
  c("Both-Low (Low metabolic)", "Mito-High / Glyco-Low (OXPHOS-dominant)"),
  c("Both-Low (Low metabolic)", "Mito-Low / Glyco-High (Glycolytic)"),
  c("Both-Low (Low metabolic)", "Both-High (Hybrid)")
)

p4c <- ggplot(fig4c_dat, aes(x = metabolic_state, y = CYT, fill = metabolic_state)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85, linewidth = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.3) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  stat_compare_means(
    comparisons = my_comparisons_cyt,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.08,
    size = 4
  ) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    title = "Cytolytic Activity Across Mitochondrial Metabolic States",
    subtitle = paste0("Kruskal-Wallis p = ", signif(kw_cyt$p.value, 3)),
    x = NULL,
    y = "Cytolytic activity score"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p4c)

save_both(
  p4c,
  "Figure4C_CYT_Boxplot_with_BothLow_Pvalues",
  width = 8,
  height = 6
)
# =========================================================
# FIGURE Supplementary 4B
# Pathway activity heatmap across metabolic states
# =========================================================

fig4c_corr <- merge(
  fig4c_dat,
  tcga_states[, c("sample", "glyco_score", "mito_score")],
  by = "sample"
)

cor_glyco <- cor.test(fig4c_corr$CYT, fig4c_corr$glyco_score)
cor_mito  <- cor.test(fig4c_corr$CYT, fig4c_corr$mito_score)

print(cor_glyco)
print(cor_mito)

##PLot CYT vs glycolysis
p_corr1 <- ggplot(fig4c_corr, aes(x = glyco_score, y = CYT)) +
  geom_point(alpha = 0.15, size = 0.5, color = "darkorange") +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "Cytolytic Activity vs Glycolysis",
    subtitle = paste0(
      "r = ", round(cor_glyco$estimate, 2),
      ", p = ", signif(cor_glyco$p.value, 3)
    ),
    x = "Glycolysis score",
    y = "Cytolytic activity score"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

##Plot CYt vs OXPHOS
p_corr2 <- ggplot(fig4c_corr, aes(x = mito_score, y = CYT)) +
  geom_point(alpha = 0.15, size = 0.5, color = "darkgreen") +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "Cytolytic Activity vs Mitochondrial Program",
    subtitle = paste0(
      "r = ", round(cor_mito$estimate, 2),
      ", p = ", signif(cor_mito$p.value, 3)
    ),
    x = "Mitochondrial score",
    y = "Cytolytic activity score"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

save_both(p_corr1, "FigureS4B_CYT_vs_Glycolysis", 6, 5)
save_both(p_corr2, "FigureS4B_CYT_vs_Mito", 6, 5)

# =========================================================
# SUPPLEMENTARY FIGURE S4C
# Cytolytic activity vs total metabolic activity
# =========================================================

fig4c_corr <- merge(
  fig4c_dat,
  tcga_states[, c("sample", "glyco_score", "mito_score", "metabolic_state")],
  by = "sample",
  all.x = FALSE
)

fig4c_corr <- fig4c_corr %>%
  mutate(
    total_metabolism = glyco_score + mito_score,
    metabolic_state = factor(as.character(metabolic_state.x), levels = state_levels)
  ) %>%
  filter(!is.na(CYT), !is.na(total_metabolism), !is.na(metabolic_state))

cor_total <- cor.test(
  fig4c_corr$CYT,
  fig4c_corr$total_metabolism,
  method = "pearson"
)

print(cor_total)

write.csv(
  fig4c_corr,
  file.path(output_dir, "FigureS4D_CYT_vs_TotalMetabolism_Data.csv"),
  row.names = FALSE
)

p_corr_total <- ggplot(fig4c_corr, aes(x = total_metabolism, y = CYT)) +
  geom_point(alpha = 0.15, size = 0.5, color = "purple") +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "Cytolytic Activity vs Total Metabolic Activity",
    subtitle = paste0(
      "r = ", round(cor_total$estimate, 2),
      ", p = ", signif(cor_total$p.value, 3)
    ),
    x = "Total metabolic activity score (Glycolysis + Mitochondrial)",
    y = "Cytolytic activity score"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

save_both(
  p_corr_total,
  "FigureS4C_CYT_vs_TotalMetabolism",
  width = 6,
  height = 5
)


# =========================================================
# FIGURE 4D
# Pathway activity heatmap across metabolic states
# =========================================================

score_df <- data.frame(
  sample = colnames(expr_mat),
  OXPHOS = score_program_safe(expr_mat, oxphos_genes),
  Glycolysis = score_program_safe(expr_mat, glyco_genes),
  ROS_Pathway = score_program_safe(expr_mat, ros_genes),
  Hypoxia = score_program_safe(expr_mat, hypoxia_genes),
  Fatty_Acid_Metabolism = score_program_safe(expr_mat, fa_genes),
  stringsAsFactors = FALSE
)

tcga_tmp <- tcga_states[, c("sample", "metabolic_state"), drop = FALSE]
tcga_tmp <- tcga_tmp[!duplicated(tcga_tmp$sample), , drop = FALSE]

score_df <- merge(score_df, tcga_tmp, by = "sample", all.x = FALSE)

write.csv(
  score_df,
  file.path(output_dir, "Figure4D_Pathway_Scores_Per_Sample.csv"),
  row.names = FALSE
)

pathway_summary <- aggregate(
  cbind(OXPHOS, Glycolysis, ROS_Pathway, Hypoxia, Fatty_Acid_Metabolism) ~ metabolic_state,
  data = score_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

pathway_summary <- pathway_summary %>%
  mutate(metabolic_state = factor(metabolic_state, levels = state_levels)) %>%
  arrange(metabolic_state)

rownames(pathway_summary) <- pathway_summary$metabolic_state
pathway_summary$metabolic_state <- NULL

pathway_mat <- t(as.matrix(pathway_summary))
pathway_mat <- t(scale(t(pathway_mat)))

write.csv(
  as.data.frame(pathway_mat),
  file.path(output_dir, "Figure4D_Pathway_Heatmap_Matrix.csv")
)

pathway_long <- as.data.frame(pathway_mat) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "metabolic_state",
    values_to = "Zscore"
  )

pathway_long$metabolic_state <- factor(
  pathway_long$metabolic_state,
  levels = state_levels
)

state_labels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
  "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
)

p4d <- ggplot(pathway_long, aes(x = metabolic_state, y = Pathway, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "Z-score"
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Pathway Activity Across Mitochondrial Metabolic States",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      face = "bold",
      size = 11,
      hjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y = element_text(
      #face = "bold",
      size = 12
    ),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

save_both(
  p4d,
  "Figure4D_Pathway_Heatmap",
  width = 8.8,
  height = 5.8
)


# =========================================================
# FIGURE 4E
# Redox and ferroptosis regulators across metabolic states
# =========================================================

# Biologically ordered genes:
# GPX4/SLC7A11 = ferroptosis suppression
# ACSL4 = ferroptosis sensitization
# PRDX3/PRDX5/TXN2 = mitochondrial redox buffering
redox_selected <- c(
  "GPX4",
  "SLC7A11",
  "ACSL4",
  "PRDX3",
  "PRDX5",
  "TXN2"
)

redox_keep <- intersect(redox_selected, rownames(expr_mat))

cat("Matched redox / ferroptosis genes:\n")
print(redox_keep)

expr_redox <- as.data.frame(expr_mat[redox_keep, , drop = FALSE])
expr_redox$gene <- rownames(expr_redox)

expr_redox_long <- expr_redox %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expr"
  )

expr_redox_long$sample <- clean_tcga_id(expr_redox_long$sample)

tcga_tmp <- tcga_states[, c("sample", "metabolic_state"), drop = FALSE]
tcga_tmp <- tcga_tmp[!duplicated(tcga_tmp$sample), , drop = FALSE]

expr_redox_long <- merge(
  expr_redox_long,
  tcga_tmp,
  by = "sample",
  all.x = FALSE
)

expr_redox_long <- expr_redox_long %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels),
    gene = factor(gene, levels = rev(redox_selected))
  )

write.csv(
  expr_redox_long,
  file.path(output_dir, "Figure4E_Redox_Ferroptosis_Long.csv"),
  row.names = FALSE
)

redox_summary <- expr_redox_long %>%
  group_by(gene, metabolic_state) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    z_expr = as.numeric(scale(mean_expr))
  ) %>%
  ungroup()

write.csv(
  redox_summary,
  file.path(output_dir, "Figure4E_Redox_Ferroptosis_Summary.csv"),
  row.names = FALSE
)

state_labels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
  "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
)

p4e <- ggplot(redox_summary, aes(x = metabolic_state, y = gene)) +
  geom_point(aes(size = median_expr, color = z_expr), alpha = 0.95) +
  
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  
  scale_color_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    name = "Mean expression\nZ-score"
  ) +
  
  scale_size_continuous(
    name = "Median\nexpression",
    range = c(2.5, 8)
  ) +
  
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  ) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Redox and Ferroptosis Regulators Are Differentially Activated Across Metabolic States",
    x = NULL,
    y = NULL
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 0,
      face = "bold",
      size = 11,
      hjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold")
  )
save_both(
  p4e,
  "Figure4E_Redox_Ferroptosis_Dotplot",
  width = 12,
  height = 7
)


# =========================================================
# FIGURE Supplementary 4E
# All redox / ferroptosis regulators with p-values
# =========================================================
library(ggpubr)

redox_keep <- intersect(toupper(redox_genes), rownames(expr_mat))

cat("All available redox / ferroptosis genes included:\n")
print(redox_keep)

expr_redox <- as.data.frame(expr_mat[redox_keep, , drop = FALSE])
expr_redox$gene <- rownames(expr_redox)

expr_redox_long <- expr_redox %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expr"
  )

expr_redox_long$sample <- clean_tcga_id(expr_redox_long$sample)

tcga_tmp <- tcga_states[, c("sample", "metabolic_state"), drop = FALSE]
tcga_tmp <- tcga_tmp[!duplicated(tcga_tmp$sample), , drop = FALSE]

expr_redox_long <- merge(
  expr_redox_long,
  tcga_tmp,
  by = "sample",
  all.x = FALSE
)

gene_class <- data.frame(
  gene = c("SLC7A11", "GPX4", "AIFM2", "ACSL4",
           "GSR", "TXN", "TXN2", "PRDX3", "PRDX5", "NFE2L2"),
  class = c("Ferroptosis", "Ferroptosis", "Ferroptosis", "Ferroptosis",
            "Redox", "Redox", "Redox", "Redox", "Redox", "Redox")
)

expr_redox_long <- expr_redox_long %>%
  left_join(gene_class, by = "gene") %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels),
    gene_label = paste0(gene, " (", class, ")")
  )

gene_order <- expr_redox_long %>%
  distinct(gene, class, gene_label) %>%
  arrange(class, gene) %>%
  pull(gene_label)

expr_redox_long$gene_label <- factor(expr_redox_long$gene_label, levels = gene_order)

# Global Kruskal-Wallis p-value per gene
pval_df <- expr_redox_long %>%
  group_by(gene_label) %>%
  summarise(
    p_value = kruskal.test(expr ~ metabolic_state)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = ifelse(
      p_value < 2.2e-16,
      "KW p < 2e-16",
      paste0("KW p = ", signif(p_value, 3))
    )
  )

write.csv(
  expr_redox_long,
  file.path(output_dir, "Figure_Supplementary4E_All_Redox_Ferroptosis_Long.csv"),
  row.names = FALSE
)

write.csv(
  pval_df,
  file.path(output_dir, "Figure_Supplementary4E_All_Redox_Ferroptosis_KW_Pvalues.csv"),
  row.names = FALSE
)

state_labels <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
  "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
  "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
  "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
)

# Pairwise comparisons: Both-Low reference
my_comparisons_redox <- list(
  c("Both-Low (Low metabolic)", "Mito-High / Glyco-Low (OXPHOS-dominant)"),
  c("Both-Low (Low metabolic)", "Mito-Low / Glyco-High (Glycolytic)"),
  c("Both-Low (Low metabolic)", "Both-High (Hybrid)")
)

pS4E <- ggplot(expr_redox_long, aes(x = metabolic_state, y = expr, fill = metabolic_state)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.9,
    color = "black"
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.12,
    size = 0.35,
    color = "black"
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y.npc = 0.98,
    size = 3,
    aes(label = paste0("Kruskal-Wallis ", after_stat(p.format)))
  ) +
  stat_compare_means(
    comparisons = my_comparisons_redox,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    step.increase = 0.10,
    size = 3
  ) +
  facet_wrap(~ gene_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  theme_bw(base_size = 13) +
  labs(
    title = "All Redox and Ferroptosis Regulators Across Mitochondrial Metabolic States",
    x = NULL,
    y = "Expression"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      face = "bold",
      size = 8,
      hjust = 0.5,
      lineheight = 0.85
    ),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

print(pS4E)

save_both(
  pS4E,
  "Figure_Supplementary4E_All_Redox_Ferroptosis_Boxplots_GlobalPairwisePvalues",
  width = 16,
  height = 11
)
# =========================================================
# OPTIONAL: SAVE PANEL QUESTIONS / KEY STATEMENTS
# =========================================================
panel_notes <- data.frame(
  Panel = c("Figure4A", "Figure4B", "Figure4C", "Figure4D"),
  Question = c(
    "Do mitochondrial metabolic states correlate with patient survival across cancers?",
    "Do mitochondrial metabolic states associate with immune subtype distributions across cancers?",
    "Do mitochondrial metabolic states show differential activation of major metabolic pathways?",
    "Are redox and ferroptosis regulators associated with mitochondrial metabolic states?"
  ),
  Key_Statement = c(
    "Distinct mitochondrial metabolic states are associated with different overall survival outcomes.",
    "Mitochondrial metabolic states show different immune subtype compositions across tumors.",
    "Metabolic states correspond to distinct pathway activation signatures reflecting oxidative versus glycolytic metabolism.",
    "Redox and ferroptosis regulators display state-specific expression patterns linked to mitochondrial metabolic programs."
  )
)

write.csv(panel_notes,
          file.path(output_dir, "FigureE_Panel_Questions_and_KeyStatements.csv"),
          row.names = FALSE)

cat("\nDone. Outputs saved in:\n", output_dir, "\n")


# =========================================================
# FIGURE 4F (FINAL)
# Integrated correlation map of metabolic programs
# =========================================================
# -----------------------------
# Program score table
# -----------------------------
fig4f_scores <- data.frame(
  sample = colnames(expr_mat),
  Mitochondrial = score_program_safe(expr_mat, mito_genes),
  Glycolysis = score_program_safe(expr_mat, glyco_genes),
  OXPHOS = score_program_safe(expr_mat, oxphos_genes),
  Hypoxia = score_program_safe(expr_mat, hypoxia_genes),
  ROS_Redox = score_program_safe(expr_mat, ros_genes),
  Fatty_Acid_Metabolism = score_program_safe(expr_mat, fa_genes),
  Cytolytic_Activity = score_program_safe(expr_mat, cyt_genes),
  Ferroptosis_Redox = score_program_safe(expr_mat, redox_genes),
  stringsAsFactors = FALSE
)

fig4f_scores$sample <- clean_tcga_id(fig4f_scores$sample)

fig4f_scores <- merge(
  fig4f_scores,
  tcga_states[, c("sample", "metabolic_state")],
  by = "sample",
  all.x = FALSE
)

fig4f_scores <- fig4f_scores %>%
  filter(!is.na(metabolic_state)) %>%
  mutate(
    metabolic_state = factor(metabolic_state, levels = state_levels)
  )

# -----------------------------
# Correlation matrix
# -----------------------------
score_mat <- fig4f_scores %>%
  dplyr::select(
    Mitochondrial,
    Glycolysis,
    OXPHOS,
    Hypoxia,
    ROS_Redox,
    Fatty_Acid_Metabolism,
    Cytolytic_Activity,
    Ferroptosis_Redox
  )

cor_mat <- cor(score_mat, use = "pairwise.complete.obs", method = "spearman")

# P-values
cor_p_mat <- matrix(NA, ncol(score_mat), ncol(score_mat))
rownames(cor_p_mat) <- colnames(score_mat)
colnames(cor_p_mat) <- colnames(score_mat)

for (i in seq_len(ncol(score_mat))) {
  for (j in seq_len(ncol(score_mat))) {
    cor_p_mat[i, j] <- suppressWarnings(
      cor.test(score_mat[[i]], score_mat[[j]], method = "spearman")$p.value
    )
  }
}

# Long format
cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("Program1") %>%
  pivot_longer(-Program1, names_to = "Program2", values_to = "rho")

p_df <- as.data.frame(cor_p_mat) %>%
  rownames_to_column("Program1") %>%
  pivot_longer(-Program1, names_to = "Program2", values_to = "p_value")

cor_df <- cor_df %>%
  left_join(p_df, by = c("Program1", "Program2")) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Remove diagonal + show only meaningful correlations
cor_df <- cor_df %>%
  mutate(
    label = ifelse(
      Program1 == Program2, "",
      ifelse(abs(rho) >= 0.3, paste0(round(rho, 2), sig), "")
    )
  )

# Order
program_order <- c(
  "Mitochondrial",
  "OXPHOS",
  "Glycolysis",
  "Hypoxia",
  "ROS_Redox",
  "Ferroptosis_Redox",
  "Fatty_Acid_Metabolism",
  "Cytolytic_Activity"
)

program_labels <- c(
  "Mitochondrial" = "Mitochondrial\nprogram",
  "OXPHOS" = "OXPHOS",
  "Glycolysis" = "Glycolysis",
  "Hypoxia" = "Hypoxia",
  "ROS_Redox" = "ROS / redox",
  "Ferroptosis_Redox" = "Ferroptosis /\nredox",
  "Fatty_Acid_Metabolism" = "Fatty acid\nmetabolism",
  "Cytolytic_Activity" = "Cytolytic\nactivity"
)

cor_df <- cor_df %>%
  mutate(
    Program1 = factor(Program1, levels = rev(program_order)),
    Program2 = factor(Program2, levels = program_order)
  )

# Highlight immune axis
cor_df <- cor_df %>%
  mutate(
    highlight = ifelse(
      Program1 == "Cytolytic_Activity" | Program2 == "Cytolytic_Activity",
      "yes", "no"
    )
  )

# -----------------------------
# Plot
# -----------------------------
p4f <- ggplot(cor_df, aes(x = Program2, y = Program1, fill = rho)) +
  geom_tile(aes(color = highlight), linewidth = 0.4) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold") +
  
  scale_color_manual(
    values = c("yes" = "black", "no" = "white"),
    guide = "none"
  ) +
  
  scale_x_discrete(labels = program_labels) +
  scale_y_discrete(labels = program_labels) +
  
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  
  theme_bw(base_size = 13) +
  labs(
    title = "Metabolic Programs Coordinate Hypoxia and Redox but Weakly Associate with Cytolytic Activity",
    subtitle = "Spearman correlations across TCGA tumors (BH-adjusted significance shown)",
    caption = "Black outline: correlations involving cytolytic activity",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 0,face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(
      hjust = 1,
      size = 10,
      face = "italic",
      margin = margin(t = 8)
    ),
    plot.caption.position = "plot",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

print(p4f)

save_both(
  p4f,
  "Figure4F_FINAL_Integrated_Correlation",
  width = 12,
  height = 8
)
