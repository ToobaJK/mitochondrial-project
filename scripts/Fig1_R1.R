# =========================================================
# Figure 1 - Pan-cancer landscape of mitochondrial gatekeeper genes
# Input: UCSC Xena (TCGA TOIL TPM), Phenotype, Survival, Immune Subtype
# Output: PDF + JPG in Output_Figure1
# =========================================================

# -------- Packages --------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(tidyr)
})

# -------- Paths --------
base_dir <- "D:/UAEU/New paper_Dr. Ajaz/Pancancer data"

expr_file <- file.path(base_dir, "tcga_RSEM_gene_tpm")
surv_file <- file.path(base_dir, "Survival_SupplementalTable_S1_20171025_xena_sp")
imm_file  <- file.path(base_dir, "Subtype_Immune_Model_Based")

# Phenotype file has long name; detect automatically:
pheno_file <- list.files(base_dir, pattern = "^TCGA_phenotype_denseDataOnlyDownload", full.names = TRUE)[1]

# Output folder for Figure 1
out_dir <- file.path(base_dir, "Output_Figure1R2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# RDS cache folder
rds_dir <- file.path(base_dir, "RDS_cache")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: standardize TCGA sample IDs to 15 characters
to_tcga15 <- function(x) substr(x, 1, 15)


# =========================================================
# 1) LOAD DATA
# =========================================================

# -------------------------
# 1) READ ONCE + SAVE RDS
# -------------------------

# ---- Expression matrix (big) ----
expr_rds <- file.path(rds_dir, "tcga_expression_matrix_log2TPM_plus0.001.rds")

if (!file.exists(expr_rds)) {
  message("Reading BIG expression file (first time only)...")
  expr_dt <- fread(expr_file)
  
  gene_col <- colnames(expr_dt)[1]
  expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
  rownames(expr_mat) <- expr_dt[[gene_col]]
  # UCSC Xena file = RSEM TPM, not raw counts
  # Log-transform before downstream analysis
  expr_mat <- log2(expr_mat + 0.001)
  
  # Standardize sample names
  colnames(expr_mat) <- to_tcga15(colnames(expr_mat))
  
  saveRDS(expr_mat, expr_rds)
  message("Saved: ", expr_rds)
} else {
  message("Loading cached expression RDS...")
  expr_mat <- readRDS(expr_rds)
}

# ---- Phenotype ----
pheno_rds <- file.path(rds_dir, "tcga_phenotype_primaryDisease_sampleType.rds")

if (!file.exists(pheno_rds)) {
  message("Reading phenotype file (first time only)...")
  pheno_dt <- fread(pheno_file)
  
  pheno_dt <- pheno_dt %>%
    rename(sample = 1) %>%
    mutate(sample15 = to_tcga15(sample))
  
  saveRDS(pheno_dt, pheno_rds)
  message("Saved: ", pheno_rds)
} else {
  message("Loading cached phenotype RDS...")
  pheno_dt <- readRDS(pheno_rds)
}

# ---- Survival (TCGA-CDR) ----
surv_rds <- file.path(rds_dir, "tcga_survival_CDR.rds")

if (!file.exists(surv_rds)) {
  message("Reading survival file (first time only)...")
  surv_dt <- fread(surv_file)
  
  surv_dt <- surv_dt %>%
    rename(sample = 1) %>%
    mutate(sample15 = to_tcga15(sample))
  
  saveRDS(surv_dt, surv_rds)
  message("Saved: ", surv_rds)
} else {
  message("Loading cached survival RDS...")
  surv_dt <- readRDS(surv_rds)
}

# ---- Immune subtype (optional) ----
imm_rds <- file.path(rds_dir, "tcga_immune_subtype.rds")

if (file.exists(imm_file)) {
  if (!file.exists(imm_rds)) {
    message("Reading immune subtype file (first time only)...")
    imm_dt <- fread(imm_file)
    
    imm_dt <- imm_dt %>%
      rename(sample = 1) %>%
      mutate(sample15 = to_tcga15(sample))
    
    saveRDS(imm_dt, imm_rds)
    message("Saved: ", imm_rds)
  } else {
    message("Loading cached immune subtype RDS...")
    imm_dt <- readRDS(imm_rds)
  }
} else {
  imm_dt <- data.frame(sample15 = character(), immune_subtype = character())
}

# -------------------------
# 2) PREPARE METADATA
# -------------------------

# Cancer type column (explicit, based on your file)
cancer_col <- "_primary_disease"

meta <- pheno_dt %>%
  transmute(
    sample15,
    cancer_type = .data[[cancer_col]]
  ) %>%
  left_join(
    surv_dt %>% select(sample15, OS, OS.time, PFI, PFI.time, DSS, DSS.time, DFI, DFI.time),
    by = "sample15"
  )

# Immune subtype column (if present)
if (nrow(imm_dt) > 0) {
  possible_imm_cols <- c("immune_subtype", "ImmuneSubtype", "Subtype", "immune subtype")
  imm_col <- possible_imm_cols[possible_imm_cols %in% colnames(imm_dt)][1]
  if (!is.na(imm_col)) {
    meta <- meta %>%
      left_join(imm_dt %>% transmute(sample15, immune_subtype = .data[[imm_col]]), by = "sample15")
  }
}

# Keep only samples present in expression
meta <- meta %>% filter(sample15 %in% colnames(expr_mat)) %>% filter(!is.na(cancer_type))
expr_mat <- expr_mat[, meta$sample15, drop = FALSE]

# -------------------------
# 3) GENE LIST (starter mitochondrial gatekeepers)
# -------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))

library(AnnotationDbi)
library(org.Hs.eg.db)



# 1) Check if rownames look like Ensembl
head(rownames(expr_mat))

# 2) Remove Ensembl version suffix if present (e.g., ENSG... .1)
ens <- rownames(expr_mat)
ens_clean <- sub("\\..*$", "", ens)

# 3) Map Ensembl -> Symbol
symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ens_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# 4) Keep only rows that mapped
keep <- !is.na(symbol_map) & symbol_map != ""
expr_mat2 <- expr_mat[keep, , drop = FALSE]
symbol_map2 <- symbol_map[keep]

# 5) Replace rownames with gene symbols
rownames(expr_mat2) <- symbol_map2

# 6) If multiple Ensembl map to same symbol, collapse by mean
expr_dt2 <- as.data.table(expr_mat2, keep.rownames = "gene")
expr_dt2 <- expr_dt2[, lapply(.SD, mean, na.rm = TRUE), by = gene]
expr_mat_symbol <- as.matrix(expr_dt2[, -1])
rownames(expr_mat_symbol) <- expr_dt2$gene

# Use this symbol-based matrix from now on
expr_mat <- expr_mat_symbol

cat("✅ Converted to gene symbols. New dimensions:", dim(expr_mat), "\n")
cat("✅ Example rownames:", head(rownames(expr_mat)), "\n")

mito_genes <- c(
  # Channels/transporters
  "VDAC1","VDAC2","VDAC3","SLC25A4","SLC25A5","SLC25A6",
  # Uncouplers
  "UCP2","UCP3",
  # Apoptosis gatekeepers
  "BAX","BAK1","BCL2","BCL2L1","MCL1","BID","BBC3","PMAIP1",
  # Dynamics
  "MFN1","MFN2","OPA1","DNM1L","FIS1","MFF",
  # Biogenesis
  "PPARGC1A","PPARGC1B","TFAM","NRF1",
  # ROS/redox
  "SOD2","TXN2","PRDX3","PRDX5","GPX1","GPX4","GSR",
  # Ferroptosis/redox (optional but relevant)
  "SLC7A11","AIFM2","ACSL4","FTH1","FTL",
  # Pore regulator
  "PPIF"
)

present_genes <- intersect(rownames(expr_mat), mito_genes)
length(present_genes)
present_genes

expr_mito <- expr_mat[present_genes, , drop = FALSE]
dim(expr_mito)

saveRDS(expr_mat, file.path(rds_dir, "tcga_expression_matrix_SYMBOLS_log2TPM.rds"))

### next time use this
expr_mat <- readRDS(file.path(rds_dir, "tcga_expression_matrix_SYMBOLS_log2TPM.rds"))

# -------------------------
# TCGA long name -> TCGA code (single source of truth)
# -------------------------
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
  out[is.na(out)] <- x2[is.na(out)]   # keep original if unmapped
  out
}



# -------------------------
# 4) FIGURE 1A: Heatmap (mean per cancer type, z-score per gene)
# -------------------------
# Merge cancer type
expr_long <- as.data.frame(t(expr_mito))
expr_long$sample15 <- rownames(expr_long)

expr_long <- expr_long %>%
  left_join(dplyr::select(as.data.frame(meta), sample15, cancer_type),
            by = "sample15") %>%
  tidyr::pivot_longer(
    cols = all_of(present_genes),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  group_by(cancer_type, gene) %>%
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")

# Wide format
expr_avg <- expr_long %>%
  pivot_wider(names_from = cancer_type, values_from = expr) %>%
  as.data.frame()

rownames(expr_avg) <- expr_avg$gene
expr_avg$gene <- NULL
expr_avg_mat <- as.matrix(expr_avg)

# Row-wise Z-score per gene across cancer types
expr_avg_z <- t(scale(t(expr_avg_mat)))
expr_avg_z[is.na(expr_avg_z)] <- 0

# Convert cancer names to TCGA codes
colnames(expr_avg_z) <- to_tcga_code(colnames(expr_avg_z))

# TCGA order
tcga_order_fig <- c(
  "COAD","LUSC","READ","GBM","KICH","LGG","KIRP","SKCM","BRCA","LUAD",
  "CESC","ESCA","MESO","DLBC","KIRC","HNSC","THCA","UCS","THYM","BLCA",
  "SARC","PAAD","TGCT","STAD","UCEC","OV","CHOL","LAML","PRAD","ACC",
  "PCPG","UVM","LIHC"
)

# Reorder columns
keep <- intersect(tcga_order_fig, colnames(expr_avg_z))
expr_avg_z <- expr_avg_z[, keep, drop = FALSE]

# Heatmap
ht <- Heatmap(
  expr_avg_z,
  name = "Row-wise\nZ-score",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2b6cb0", "white", "#c53030")
  ),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  column_title = "Pan-Cancer Expression of Mitochondrial-Associated Genes",
  row_title = "Mitochondrial-associated genes",
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Row-wise\nZ-score",
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2")
  )
)

pdf(file.path(out_dir, "Figure1A_Mitochondrial-Associated_Heatmap.pdf"),
    width = 12, height = 9)
draw(ht, heatmap_legend_side = "right")
dev.off()

jpeg(file.path(out_dir, "Figure1A_Mitochondrial-Associated_Heatmap.jpg"),
     width = 3600, height = 2700, res = 400)
draw(ht, heatmap_legend_side = "right")
dev.off()


# -------------------------
# 5) FIGURE 1B: Mito Gatekeeper Score per cancer type
# -------------------------

# Score = mean row-wise Z-score across mitochondrial gatekeeper genes
mito_score <- colMeans(expr_avg_z, na.rm = TRUE)

df_score <- data.frame(
  cancer_type = names(mito_score),
  mito_gatekeeper_score = as.numeric(mito_score)
) %>%
  arrange(desc(mito_gatekeeper_score))

# Keep TCGA codes as labels
df_score$cancer_type <- factor(df_score$cancer_type, levels = df_score$cancer_type)

p_score <- ggplot(df_score,
                  aes(x = cancer_type,
                      y = mito_gatekeeper_score,
                      fill = mito_gatekeeper_score)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40",
             linewidth = 0.4) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#2b6cb0",
    mid = "white",
    high = "#c53030",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = paste0("Mean row-wise Z-score (n = ", length(present_genes), " genes)"),
    title = "Mitochondrial Composite Score Across TCGA Cancer Types"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "Figure1B_MitometabolicScore_byCancer.pdf"),
       p_score, width = 10, height = 8)

ggsave(file.path(out_dir, "Figure1B_MitometabolicScore_byCancer.jpg"),
       p_score, width = 10, height = 8, dpi = 400)
# -------------------------
# 6) FIGURE 1C (Optional): Example genes boxplots across cancers
# -------------------------

# =========================================================
# Figure 1C (FINAL): Representative mitochondrial gatekeeper genes
# - Uses TCGA short codes (LUAD, BRCA, etc.)
# - Facet titles include gene function: e.g., "BAX (Apoptosis gatekeeper)"
# - Replaces missing SOD2 with PRDX3 (mitochondrial redox enzyme)
# - Saves PDF + JPG in out_dir
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# -------------------------
# 0) Choose genes (SOD2 replaced with PRDX3)
# -------------------------
genes_wanted <- c(
  "VDAC1","SLC25A4",         # Transport/ANT
  "BAX","BCL2L1",            # Apoptosis
  "MFN2","DNM1L",            # Dynamics
  "PRDX3","GPX4","SLC7A11"   # Redox/Ferroptosis
)

genes_for_box <- genes_wanted[genes_wanted %in% present_genes]
missing <- setdiff(genes_wanted, genes_for_box)
if (length(missing) > 0) message("Missing genes (not plotted): ", paste(missing, collapse = ", "))

# -------------------------
# 1) Gene -> Function mapping for facet titles
# -------------------------
gene_function <- c(
  VDAC1   = "Transport / OMM pore",
  SLC25A4 = "ANT / Inner membrane exchange",
  BAX     = "Apoptosis regulator",
  BCL2L1  = "Apoptosis regulator",
  MFN2    = "Mitochondrial dynamics",
  DNM1L   = "Mitochondrial dynamics",
  PRDX3   = "Redox / ROS buffering",
  GPX4    = "Redox / Ferroptosis",
  SLC7A11 = "Redox / Ferroptosis"
)

facet_label <- function(g) paste0(g, " (", gene_function[g], ")")


# -------------------------
# 3) Build long table for plotting
# -------------------------
expr_box <- as.data.frame(t(expr_mito[genes_for_box, , drop = FALSE]))
expr_box$sample15 <- rownames(expr_box)

expr_box_long <- expr_box %>%
  left_join(dplyr::select(as.data.frame(meta), sample15, cancer_type),
            by = "sample15") %>%
  tidyr::pivot_longer(
    cols = all_of(genes_for_box),
    names_to = "gene",
    values_to = "expr"
  )

# Use TCGA short codes
expr_box_long$cancer_type <- to_tcga_code(expr_box_long$cancer_type)

# Keep cancer order consistent with Figure 1B (convert df_score too)
df_score$cancer_type <- to_tcga_code(df_score$cancer_type)
ordered_cancers <- df_score$cancer_type
expr_box_long$cancer_type <- factor(expr_box_long$cancer_type, levels = ordered_cancers)

# Keep gene order as in genes_wanted
expr_box_long$gene <- factor(expr_box_long$gene, levels = genes_for_box)

expr_box_long$expr <- pmax(expr_box_long$expr, -5)####will remove extreme value
# -------------------------
# 4) Plot
# -------------------------
p_box <- ggplot(expr_box_long,
                aes(x = cancer_type, y = expr, fill = cancer_type)) +
  geom_boxplot(
    outlier.size = 0.15,
    linewidth = 0.25
  ) +
  facet_wrap(
    ~gene,
    scales = "free_y",
    ncol = 3,
    labeller = labeller(gene = facet_label)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "log2(TPM + 0.001)",
    title = "Representative Mitochondrial-Associated Genes Across TCGA Cancers"
  )
# -------------------------
# 5) Save outputs
# -------------------------
ggsave(file.path(out_dir, "Figure1C_ExampleGenes_Boxplots_FunctionTitles_TCGAcode.pdf"),
       p_box, width = 14, height = 10)

ggsave(file.path(out_dir, "Figure1C_ExampleGenes_Boxplots_FunctionTitles_TCGAcode.jpg"),
       p_box, width = 14, height = 10, dpi = 400)

####----------------------------------------------------
###Figure 1D
###--------------------------------------------------------

# =========================================================
# Figure 1D: "Metabolic rewiring" = Mito gatekeeper score vs Glycolysis score
# - Uses your TCGA PanCan expression matrix (expr_mat; gene x sample, log2(TPM+0.001))
# - Uses meta (sample15, cancer_type) from Xena phenotype/survival merge
# - Converts cancer names to TCGA short codes (LUAD, BRCA, etc.)
# - Saves PDF + JPG to out_dir
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})


# -------------------------
# 1) Define gene sets
#    - Gatekeepers: use your present_genes (from mito list intersection)
#    - Glycolysis: use a compact Hallmark-like core set (edit/expand anytime)
# -------------------------

gatekeeper_genes <- present_genes  # your 38 genes (or set to mito_genes[mito_genes %in% rownames(expr_mat)])

glycolysis_genes <- c(
  "SLC2A1","SLC2A3",          # Glucose transport
  "HK1","HK2","GPI","PFKM","PFKP",
  "ALDOA","ALDOB","TPI1","GAPDH",
  "PGK1","PGAM1","ENO1","ENO2",
  "PKM","LDHA","LDHB",
  "PDK1","PDK3",              # blocks pyruvate entry (rewiring marker)
  "SLC16A1","SLC16A3"         # MCT1/MCT4 lactate transport
)

# keep only genes present in your matrix
gatekeeper_genes_use <- intersect(gatekeeper_genes, rownames(expr_mat))
glycolysis_genes_use <- intersect(glycolysis_genes, rownames(expr_mat))

message("Gatekeeper genes used: ", length(gatekeeper_genes_use))
message("Glycolysis genes used: ", length(glycolysis_genes_use))

if (length(gatekeeper_genes_use) < 10) warning("Low number of gatekeeper genes found in expr_mat. Check rownames(expr_mat).")
if (length(glycolysis_genes_use) < 10) warning("Low number of glycolysis genes found in expr_mat. Consider expanding glycolysis_genes list.")

# -------------------------
# 2) Compute mean expression per cancer type for each gene set
# -------------------------
meta_df <- as.data.frame(meta) %>%
  dplyr::select(sample15, cancer_type) %>%
  mutate(cancer_type = tools::toTitleCase(cancer_type))

# Helper: average a gene set by cancer type (returns gene x cancer matrix)
avg_by_cancer <- function(expr_mat, genes, meta_df) {
  # expr_mat: gene x sample
  mat <- expr_mat[genes, , drop = FALSE]
  df <- as.data.frame(t(mat))
  df$sample15 <- rownames(df)
  
  df_long <- df %>%
    left_join(meta_df, by = "sample15") %>%
    tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(cancer_type, gene) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")
  
  df_wide <- df_long %>%
    tidyr::pivot_wider(names_from = cancer_type, values_from = expr) %>%
    as.data.frame()
  
  rownames(df_wide) <- df_wide$gene
  df_wide$gene <- NULL
  as.matrix(df_wide)
}

gate_avg <- avg_by_cancer(expr_mat, gatekeeper_genes_use, meta_df)
gly_avg  <- avg_by_cancer(expr_mat, glycolysis_genes_use, meta_df)

# -------------------------
# 3) Z-score per gene across cancer types (to compare programs fairly)
# -------------------------
zscore_rows <- function(m) {
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  z
}
gate_z <- zscore_rows(gate_avg)
gly_z  <- zscore_rows(gly_avg)

# Program scores per cancer type (mean z across genes)
mito_gatekeeper_score <- colMeans(gate_z, na.rm = TRUE)
glycolysis_score      <- colMeans(gly_z,  na.rm = TRUE)

df_1d <- data.frame(
  cancer_type = names(mito_gatekeeper_score),
  mito_gatekeeper_score = as.numeric(mito_gatekeeper_score),
  glycolysis_score = as.numeric(glycolysis_score)
)

# Convert to TCGA abbreviations
df_1d$cancer_code <- to_tcga_code(df_1d$cancer_type)

# Optional: classify metabolic state by quadrants (based on 0 since scores are z-based)
# =========================================================
# Figure 1D (UPDATED):
# - Replace "Hybrid" with "Both"
# - Color points by the 4 states (quadrants)
# - Keep dashed 0-lines
# =========================================================

# --- classify metabolic state by quadrants (0-based because scores are z-based)
df_1d$state <- with(df_1d,
                    ifelse(mito_gatekeeper_score >= 0 & glycolysis_score < 0,
                           "Mito-High / Glyco-Low (OXPHOS-dominant)",
                           ifelse(mito_gatekeeper_score < 0 & glycolysis_score >= 0,
                                  "Mito-Low / Glyco-High (Glycolytic)",
                                  ifelse(mito_gatekeeper_score >= 0 & glycolysis_score >= 0,
                                         "Both-High (Hybrid)",
                                         "Both-Low (Low metabolic)"))))

df_1d$state <- factor(df_1d$state,
                      levels = c(
                        "Mito-High / Glyco-Low (OXPHOS-dominant)",
                        "Mito-Low / Glyco-High (Glycolytic)",
                        "Both-High (Hybrid)",
                        "Both-Low (Low metabolic)"
                      ))
# --- plot with colors for states
p_1d <- ggplot(df_1d, aes(x = glycolysis_score, y = mito_gatekeeper_score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = state), size = 3) +
  ggrepel::geom_text_repel(aes(label = cancer_code, color = state),
                           size = 3, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
    "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
    "Both-High (Hybrid)"                      = "purple",
    "Both-Low (Low metabolic)"                = "maroon"
  )) +
  theme_classic(base_size = 12) +
  labs(
    x = "Glycolysis composite score (mean row-wise Z-score)",
    y = "Mitochondrial composite score (mean row-wise Z-score)",
    color = "Metabolic state",
    title = "Metabolic Rewiring Across TCGA Cancers",
    subtitle = "Quadrants highlight glycolysis-high vs mitochondria-high lineages"
  )

ggsave(file.path(out_dir, "Figure1D1_Mito_vs_Glycolysis_Scatter1_TCGAcode.pdf"),
       p_1d, width = 11, height = 8)

ggsave(file.path(out_dir, "Figure1D1_Mito_vs_Glycolysis_Scatter1_TCGAcode.jpg"),
       p_1d, width = 11, height = 8, dpi = 400)


# =========================================================
# Figure 1E: Unsupervised clustering into "mitochondrial metabolic states"
# - Clusters TCGA cancer TYPES (not individual tumors) using expr_avg_z
#   (expr_avg_z = gene x cancer_type z-scored matrix from Figure 1A)
# - Assigns 4 states using k-means (k = 4)
#   Labels: "Mito-High", "Mito-Low", "Both-High", "Both-Low"
#   based on the Mito Gatekeeper score (from expr_avg_z) and Glycolysis score
#   (computed from the same pipeline used for Figure 1D)
# - Produces:
#   (1) Figure1E_Heatmap_Clusters.pdf/.jpg  (heatmap + top annotation for cluster/state)
#   (2) Figure1E_StateComposition_Bar.pdf/.jpg (counts of cancers per state)
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
})

# =========================================================
# 0) TCGA NAME → CODE MAPPING
# =========================================================

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

# =========================================================
# 1) Prepare matrix + metadata
# =========================================================

# =========================================================
# 1) Prepare matrix + metadata (match using TCGA codes)
# =========================================================

# Convert expr_avg_z column names to TCGA codes (if already TCGA, it's fine)
colnames(expr_avg_z) <- to_tcga_code(colnames(expr_avg_z))

# Create a TCGA code column in df_1d
df_1d <- df_1d %>%
  mutate(
    cancer_type = as.character(cancer_type),
    cancer_code = to_tcga_code(cancer_type)
  )

# Match using TCGA codes
common_codes <- intersect(colnames(expr_avg_z), df_1d$cancer_code)

# Subset matrix by codes
expr_avg_z_sub <- expr_avg_z[, common_codes, drop = FALSE]

# Subset df_1d by codes (one row per code)
df_1d_sub <- df_1d %>%
  filter(cancer_code %in% common_codes) %>%
  distinct(cancer_code, .keep_all = TRUE)

# OPTIONAL sanity check
stopifnot(ncol(expr_avg_z_sub) >= 4)  # you need at least 4 for centers=4

# =========================================================
# 2) K-means clustering
# =========================================================
set.seed(123)
km <- kmeans(t(expr_avg_z_sub), centers = 4, nstart = 50)

cluster_df <- data.frame(
  cancer_code = colnames(expr_avg_z_sub),
  cluster = paste0("Cluster-", km$cluster),
  stringsAsFactors = FALSE
)

cluster_df <- dplyr::left_join(
  cluster_df,
  dplyr::select(df_1d_sub, cancer_code, mito_gatekeeper_score, glycolysis_score),
  by = "cancer_code"
)

# =========================================================
# 3) Define metabolic state
# =========================================================
cluster_df$state <- with(cluster_df,
                         ifelse(mito_gatekeeper_score >= 0 & glycolysis_score < 0,
                                "Mito-High / Glyco-Low (OXPHOS-dominant)",
                                ifelse(mito_gatekeeper_score < 0 & glycolysis_score >= 0,
                                       "Mito-Low / Glyco-High (Glycolytic)",
                                       ifelse(mito_gatekeeper_score >= 0 & glycolysis_score >= 0,
                                              "Both-High (Hybrid)",
                                              "Both-Low (Low metabolic)"))))

cluster_df$state <- factor(cluster_df$state,
                           levels = c(
                             "Mito-High / Glyco-Low (OXPHOS-dominant)",
                             "Mito-Low / Glyco-High (Glycolytic)",
                             "Both-High (Hybrid)",
                             "Both-Low (Low metabolic)"
                           ))

# ✅ sanity: make sure state exists
stopifnot(!all(is.na(cluster_df$state)))

# =========================================================
# 4) Order columns by cluster
# =========================================================
ord <- order(km$cluster)
expr_avg_z_ord <- expr_avg_z_sub[, ord, drop = FALSE]
cluster_df <- cluster_df[ord, ]

# column names already TCGA codes (use cancer_code)
colnames(expr_avg_z_ord) <- cluster_df$cancer_code

# =========================================================
# 5) Annotation colors
# =========================================================
state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)
cluster_colors <- setNames(
  c("#111827","#374151","#6b7280","#9ca3af"),
  paste0("Cluster-", 1:4)
)

ha <- HeatmapAnnotation(
  Cluster = cluster_df$cluster,
  State   = cluster_df$state,
  col = list(
    Cluster = cluster_colors,
    State   = state_colors
  )
)
# =========================================================
# 6) Draw heatmap
# =========================================================
#draw(ht1e, heatmap_legend_side = "right", annotation_legend_side = "right")

ht1e <- Heatmap(
  expr_avg_z_ord,
  name = "Row-wise\nZ-score",
  top_annotation = ha,
  col = colorRamp2(c(-2,0,2), c("#2b6cb0","white","#c53030")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_title = "Unsupervised clustering of TCGA Cancers Based on Mitochondrial-Associated Gene Expression",
  row_title = "Mitochondrial-associated genes"
)

pdf(file.path(out_dir, "Figure1E_Heatmap_Clusters.pdf"), width = 12, height = 9)
draw(ht1e, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

jpeg(file.path(out_dir, "Figure1E_Heatmap_Clusters.jpg"),
     width = 3600, height = 2700, res = 300)
draw(ht1e, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
draw(ht1e)
dev.off()

cat("✅ Figure 1E generated without NA labels.\n")

# -------------------------
# 5) Optional: state composition bar plot (how many cancers in each state)
# -------------------------
df_state_counts <- cluster_df %>%
  count(state) %>%
  arrange(desc(n))

p_state <- ggplot(df_state_counts, aes(x = state, y = n, fill = state)) +
  geom_col(width = 0.75) +
  
  # ✅ Add numbers on top of bars
  geom_text(aes(label = n), 
            vjust = -0.3, size = 4, fontface = "bold") +
  
  scale_fill_manual(values = state_colors) +
  scale_x_discrete(
    labels = c(
      "Mito-High / Glyco-Low (OXPHOS-dominant)" = "Mito-high\nGlyco-low\n(OXPHOS-dominant)",
      "Mito-Low / Glyco-High (Glycolytic)"      = "Mito-low\nGlyco-high\n(Glycolytic)",
      "Both-High (Hybrid)"                      = "Both-high\n(Hybrid)",
      "Both-Low (Low metabolic)"                = "Both-low\n(Low)"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(15, 25, 25, 25)
  ) +
  labs(
    x = NULL,
    y = "Number of TCGA cancer types",
    title = "Distribution of TCGA Cancer Types Across Mitochondrial Metabolic States"
  )

ggsave(file.path(out_dir, "Figure1supplementary_StateComposition_Bar.pdf"),
       p_state, width = 10, height = 6)
ggsave(file.path(out_dir, "Figure1supplementary_StateComposition_Bar.jpg"),
       p_state, width = 10, height = 6, dpi = 400)

# Save the cluster assignment table (useful later)
write.table(cluster_df,
            file.path(out_dir, "Figure1E_TCGA_MetabolicStates_Table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =========================================================
# Figure 1F: Barplot of TCGA cancers colored by metabolic state
# =========================================================

library(dplyr)
library(ggplot2)

# Make sure df_1d has state assigned (from Figure 1E section)

df_bar <- df_1d %>%
  mutate(
    cancer_code = to_tcga_code(cancer_type),
    state = factor(state,
                   levels = c(
                     "Mito-High / Glyco-Low (OXPHOS-dominant)",
                     "Mito-Low / Glyco-High (Glycolytic)",
                     "Both-High (Hybrid)",
                     "Both-Low (Low metabolic)"
                   ))
  ) %>%
  arrange(desc(mito_gatekeeper_score))

# Order cancers by mito score
df_bar$cancer_code <- factor(df_bar$cancer_code,
                             levels = df_bar$cancer_code)

state_colors <- c(
  "Mito-High / Glyco-Low (OXPHOS-dominant)" = "darkgreen",
  "Mito-Low / Glyco-High (Glycolytic)"      = "orange",
  "Both-High (Hybrid)"                      = "purple",
  "Both-Low (Low metabolic)"                = "maroon"
)

p_1f <- ggplot(df_bar,
               aes(x = cancer_code,
                   y = mito_gatekeeper_score,
                   fill = state)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = state_colors) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(15, 20, 40, 20)   # 👈 critical for bottom space
  ) +
  labs(
    x = NULL,
    y = "Mitochondrial composite score (mean row-wise Z-score)",
    fill = "Metabolic state",
    title = "TCGA cancers ranked by mitochondrial composite score and annotated by metabolic state"
  )

ggsave(file.path(out_dir, "Figure1F_MetabolicState_Barplot.pdf"),
       p_1f, width = 12, height = 7)

ggsave(file.path(out_dir, "Figure1F_MetabolicState_Barplot.jpg"),
       p_1f, width = 12, height = 7, dpi = 300)

# -------------------------
# 7) Save processed objects used for Figure 1 (fast reload)
# -------------------------
figure1_rds <- file.path(rds_dir, "Figure1_processed_objects.rds")
saveRDS(
  list(
    expr_mito = expr_mito,
    expr_avg_z = expr_avg_z,
    meta = meta,
    mito_genes_present = present_genes
  ),
  figure1_rds
)

# Also export the matrices (optional)
write.table(expr_avg_z, file.path(out_dir, "Figure1A_exprAvgZ_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
fwrite(meta, file.path(out_dir, "Figure1_Metadata_used.tsv"), sep = "\t")

cat("\n✅ Done. Figure 1 saved to: ", out_dir, "\n")
cat("✅ Cached RDS files saved to: ", rds_dir, "\n")

# =========================================================
# Supplementary Table S1:
# Random mitochondrial gene-set control analysis
# Purpose:
# Test whether the curated mitochondrial gene set gives
# stronger metabolic stratification than random gene sets
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(123)

# -----------------------------
# 1. Inputs from your pipeline
# -----------------------------
# expr_mat = gene x sample log2(TPM + 0.001)
# meta = metadata with sample15 and cancer_type
# present_genes = curated mitochondrial gene set used in Figure 1
# glycolysis_genes_use = glycolysis genes present in expr_mat
# out_dir = output folder

meta <- as.data.frame(meta)

curated_genes <- intersect(present_genes, rownames(expr_mat))
n_curated <- length(curated_genes)

cat("Curated genes used:", n_curated, "\n")

# -----------------------------
# 2. Define background gene set
# -----------------------------
# Best option later:
# mito_background <- setdiff(intersect(mitocarta_genes, rownames(expr_mat)), curated_genes)

# Temporary broad background:
mito_background <- setdiff(rownames(expr_mat), curated_genes)

cat("Background genes available:", length(mito_background), "\n")

if (length(mito_background) < n_curated) {
  stop("Background gene set is smaller than curated gene set.")
}

# -----------------------------
# 3. Helper functions
# -----------------------------

avg_by_cancer <- function(expr_mat, genes, meta) {
  
  genes <- intersect(genes, rownames(expr_mat))
  
  if (length(genes) == 0) {
    stop("No genes found in expression matrix.")
  }
  
  mat <- expr_mat[genes, , drop = FALSE]
  df <- as.data.frame(t(mat))
  df$sample15 <- rownames(df)
  
  meta_small <- as.data.frame(meta)[, c("sample15", "cancer_type")]
  
  df_long <- df %>%
    dplyr::left_join(meta_small, by = "sample15") %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(genes),
      names_to = "gene",
      values_to = "expr"
    ) %>%
    dplyr::filter(!is.na(cancer_type)) %>%
    dplyr::group_by(cancer_type, gene) %>%
    dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")
  
  df_wide <- df_long %>%
    tidyr::pivot_wider(
      names_from = cancer_type,
      values_from = expr
    ) %>%
    as.data.frame()
  
  rownames(df_wide) <- df_wide$gene
  df_wide$gene <- NULL
  
  as.matrix(df_wide)
}

zscore_rows <- function(m) {
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  return(z)
}

assign_state <- function(mito, gly) {
  ifelse(mito >= 0 & gly < 0, "OXPHOS-dominant",
         ifelse(mito < 0 & gly >= 0, "Glycolytic",
                ifelse(mito >= 0 & gly >= 0, "Hybrid",
                       "Low metabolic")))
}

get_state_counts <- function(states) {
  tab <- table(states)
  data.frame(
    n_OXPHOS = ifelse("OXPHOS-dominant" %in% names(tab), tab["OXPHOS-dominant"], 0),
    n_Glycolytic = ifelse("Glycolytic" %in% names(tab), tab["Glycolytic"], 0),
    n_Hybrid = ifelse("Hybrid" %in% names(tab), tab["Hybrid"], 0),
    n_Low_metabolic = ifelse("Low metabolic" %in% names(tab), tab["Low metabolic"], 0)
  )
}

# -----------------------------
# 4. Curated mitochondrial score
# -----------------------------

curated_avg <- avg_by_cancer(expr_mat, curated_genes, meta)
curated_z <- zscore_rows(curated_avg)
curated_score <- colMeans(curated_z, na.rm = TRUE)

# Glycolysis score
gly_avg <- avg_by_cancer(expr_mat, glycolysis_genes_use, meta)
gly_z <- zscore_rows(gly_avg)
gly_score <- colMeans(gly_z, na.rm = TRUE)

common_cancers <- intersect(names(curated_score), names(gly_score))

curated_df <- data.frame(
  cancer_type = common_cancers,
  mitochondrial_score = curated_score[common_cancers],
  glycolysis_score = gly_score[common_cancers]
)

curated_df$state <- assign_state(
  curated_df$mitochondrial_score,
  curated_df$glycolysis_score
)

curated_separation <- sd(curated_df$mitochondrial_score, na.rm = TRUE)
curated_state_counts <- get_state_counts(curated_df$state)

cat("Curated separation score:", curated_separation, "\n")
print(table(curated_df$state))

# -----------------------------
# 5. Random gene-set controls
# -----------------------------

n_iter <- 1000
random_results <- vector("list", n_iter)

for (i in seq_len(n_iter)) {
  
  if (i %% 100 == 0) {
    cat("Running iteration:", i, "of", n_iter, "\n")
  }
  
  random_genes <- sample(mito_background, n_curated, replace = FALSE)
  
  random_avg <- avg_by_cancer(expr_mat, random_genes, meta)
  random_z <- zscore_rows(random_avg)
  random_score <- colMeans(random_z, na.rm = TRUE)
  
  common <- intersect(names(random_score), names(gly_score))
  
  random_df <- data.frame(
    cancer_type = common,
    mitochondrial_score = random_score[common],
    glycolysis_score = gly_score[common]
  )
  
  random_df$state <- assign_state(
    random_df$mitochondrial_score,
    random_df$glycolysis_score
  )
  
  state_counts <- get_state_counts(random_df$state)
  
  random_results[[i]] <- data.frame(
    iteration = i,
    separation_score = sd(random_df$mitochondrial_score, na.rm = TRUE),
    n_OXPHOS = state_counts$n_OXPHOS,
    n_Glycolytic = state_counts$n_Glycolytic,
    n_Hybrid = state_counts$n_Hybrid,
    n_Low_metabolic = state_counts$n_Low_metabolic
  )
}

random_results_df <- dplyr::bind_rows(random_results)

# -----------------------------
# 6. Compare curated vs random
# -----------------------------

p_empirical <- mean(random_results_df$separation_score >= curated_separation)

summary_table <- data.frame(
  analysis = c("Curated mitochondrial gene set", "Random gene sets"),
  n_gene_sets = c(1, n_iter),
  genes_per_set = c(n_curated, n_curated),
  separation_score = c(
    curated_separation,
    mean(random_results_df$separation_score, na.rm = TRUE)
  ),
  separation_score_sd = c(
    NA,
    sd(random_results_df$separation_score, na.rm = TRUE)
  ),
  empirical_p_value = c(p_empirical, NA),
  n_OXPHOS = c(curated_state_counts$n_OXPHOS, NA),
  n_Glycolytic = c(curated_state_counts$n_Glycolytic, NA),
  n_Hybrid = c(curated_state_counts$n_Hybrid, NA),
  n_Low_metabolic = c(curated_state_counts$n_Low_metabolic, NA)
)

print(summary_table)

# -----------------------------
# 7. Save supplementary tables
# -----------------------------

write.csv(
  summary_table,
  file.path(out_dir, "Supplementary_Table_S1_Random_Gene_Set_Control_Summary.csv"),
  row.names = FALSE
)

write.csv(
  random_results_df,
  file.path(out_dir, "Supplementary_Table_S1_Random_Gene_Set_Control_AllIterations.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Optional simple supplementary figure
# -----------------------------

plot_df <- random_results_df %>%
  dplyr::mutate(type = "Random gene sets")

p_control <- ggplot(plot_df, aes(x = type, y = separation_score)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(
    yintercept = curated_separation,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Separation score",
    title = "Random Gene-Set Control Analysis",
    subtitle = "Dashed line indicates curated mitochondrial gene set"
  )

ggsave(
  file.path(out_dir, "Supplementary_Figure_S1C_Random_Gene_Set_Control.pdf"),
  p_control,
  width = 5,
  height = 4
)

ggsave(
  file.path(out_dir, "Supplementary_Figure_S1C_Random_Gene_Set_Control.jpg"),
  p_control,
  width = 5,
  height = 4,
  dpi = 400
)

cat("Done. Supplementary control tables and figure saved in:", out_dir, "\n")
