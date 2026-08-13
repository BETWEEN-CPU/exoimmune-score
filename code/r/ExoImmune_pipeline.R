# ============================================================
# ExoImmune Score: A transcriptome-derived signature of exosome-associated
# genes for predicting immunotherapy response
#
# Data: IMvigor210 (training) + TCGA-BLCA / TCGA-SKCM (validation)
# ============================================================
# This script reproduces the complete analytical workflow:
#   Step 1: Load IMvigor210 expression and phenotype data; map gene symbols;
#           align samples; prepare response and survival matrices.
#   Step 2: Intersect ExoCarta and Vesiclepedia human genes with IMvigor210;
#           perform Wilcoxon rank-sum screening (Responder vs NonResponder).
#   Step 3: LASSO logistic regression (10-fold CV, lambda.1se); build
#           ExoImmune Score; evaluate discrimination (ROC/AUC).
#   Step 4: Kaplan-Meier analysis (median dichotomy); univariate and
#           multivariate Cox proportional hazards regression.
#   Step 5: Apply the trained score to TCGA-BLCA and TCGA-SKCM;
#           assess prognostic association and immune module correlations.
#   Step 6: Compare ExoImmune Score correlations with immune vs non-immune
#           control modules; generate summary heatmaps.
# ============================================================

# ============================================================
# IMvigor210 ExoImmune Score: Complete Pipeline
# ============================================================

# Working directory

root_dir <- Sys.getenv("EXOIMMUNE_ROOT", unset = "D:/善意无穷倒退/李新月/EV")

# ============================================================
# Helper functions

# ============================================================

#' 从数据框中按正则表达式模式选择列名
pick_col <- function(df, patterns) {
  cn <- colnames(df)
  hit <- cn[grepl(paste(patterns, collapse = "|"), cn, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

#' 按基因符号合并重复行（取均值）
collapse_by_symbol <- function(expr_mat, gene_symbol) {
  gene_symbol <- as.character(gene_symbol)
  keep <- !is.na(gene_symbol) & gene_symbol != ""
  expr_mat <- expr_mat[keep, , drop = FALSE]
  gene_symbol <- gene_symbol[keep]

  split_idx <- split(seq_along(gene_symbol), gene_symbol)
  collapsed <- t(vapply(split_idx, function(idx) {
    if (length(idx) == 1) {
      expr_mat[idx, ]
    } else {
      colMeans(expr_mat[idx, , drop = FALSE], na.rm = TRUE)
    }
  }, numeric(ncol(expr_mat))))

  rownames(collapsed) <- names(split_idx)
  collapsed
}

#' 标准化基因符号（大写、去空格、去空值）
norm_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "NAN", "NULL")] <- NA
  unique(x[!is.na(x)])
}

#' 判断物种是否为人类
is_human <- function(x) {
  x <- as.character(x)
  grepl("HOMO SAPIENS|HUMAN", toupper(x))
}

#' Wilcoxon筛选：比较Responder vs NonResponder的基因表达
wilcox_screen <- function(expr_mat, group_vec) {
  stopifnot(ncol(expr_mat) == length(group_vec))
  group_vec <- factor(group_vec, levels = c("NonResponder", "Responder"))

  res <- lapply(seq_len(nrow(expr_mat)), function(i) {
    x <- expr_mat[i, ]
    x1 <- x[group_vec == "NonResponder"]
    x2 <- x[group_vec == "Responder"]

    p <- tryCatch(
      wilcox.test(x2, x1)$p.value,
      error = function(e) NA_real_
    )

    c(
      logFC = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
      p.value = p
    )
  })

  res <- as.data.frame(do.call(rbind, res))
  res$gene <- rownames(expr_mat)
  res$padj <- p.adjust(res$p.value, method = "BH")
  res$direction <- ifelse(res$logFC > 0, "Higher_in_Responder", "Higher_in_NonResponder")
  res <- res[order(res$padj, res$p.value), ]
  rownames(res) <- NULL
  res
}

#' 安全数字解析：提取字符串中的第一个数字
to_num <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA
  num <- suppressWarnings(as.numeric(x))
  bad <- is.na(num) & !is.na(x)
  if (any(bad)) {
    extracted <- regmatches(x[bad], regexpr("[0-9]+(\\.[0-9]+)?", x[bad]))
    num[bad] <- suppressWarnings(as.numeric(extracted))
  }
  num
}

#' TCGA样本ID提取前12位作为患者ID
tcga_patient_id <- function(x) substr(as.character(x), 1, 12)

#' 递归查找匹配模式的文件
find_file <- function(root, pattern) {
  x <- list.files(root, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  x[grepl(pattern, x, ignore.case = TRUE)]
}

#' 读取表达矩阵文件（TSV格式，第一列为基因名）
read_expr_table <- function(file) {
  df <- read.delim(file, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
  if (ncol(df) < 2) stop("Expression file has too few columns: ", file)

  first_col <- df[[1]]
  if (is.character(first_col) || is.factor(first_col)) {
    rownames(df) <- as.character(first_col)
    df <- df[, -1, drop = FALSE]
  }

  mat <- as.matrix(df)
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- norm_gene(rownames(mat))
  mat
}

#' 读取临床数据表格
read_clin_table <- function(file) {
  df <- read.delim(file, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
  as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
}

#' 合并表达矩阵中的重复基因名（取均值）
collapse_dups <- function(mat) {
  rn <- rownames(mat)
  keep <- !is.na(rn) & rn != ""
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]

  if (any(duplicated(rn))) {
    idx <- split(seq_along(rn), rn)
    collapsed <- t(vapply(idx, function(ii) {
      if (length(ii) == 1) {
        as.numeric(mat[ii, ])
      } else {
        colMeans(mat[ii, , drop = FALSE], na.rm = TRUE)
      }
    }, numeric(ncol(mat))))
    rownames(collapsed) <- names(idx)
    return(collapsed)
  }

  rownames(mat) <- rn
  mat
}

#' 根据LASSO系数和训练集缩放参数计算ExoImmune Score
make_score <- function(expr, coef_df, scale_tbl) {
  overlap <- intersect(scale_tbl$gene, rownames(expr))
  if (length(overlap) == 0) stop("No overlap between expression and score genes.")

  expr <- expr[overlap, , drop = FALSE]
  scale_tbl <- scale_tbl[match(overlap, scale_tbl$gene), , drop = FALSE]
  stopifnot(all(scale_tbl$gene == rownames(expr)))

  z <- sweep(expr, 1, scale_tbl$mean, "-")
  z <- sweep(z, 1, scale_tbl$sd, "/")
  z[!is.finite(z)] <- 0

  coef_use <- coef_df[coef_df$gene %in% rownames(z) & coef_df$gene != "(Intercept)", , drop = FALSE]
  coef_use <- coef_use[match(rownames(z), coef_use$gene), , drop = FALSE]

  if (any(is.na(coef_use$coef))) {
    stop("Missing coefficients for some genes.")
  }

  intercept <- 0
  if (("(Intercept)") %in% coef_df$gene) {
    intercept <- coef_df$coef[match("(Intercept)", coef_df$gene)]
    if (is.na(intercept)) intercept <- 0
  }

  as.numeric(intercept + t(z) %*% coef_use$coef)
}

#' 中位数分层的KM生存分析 + Cox回归
median_km <- function(df, time_col = "time", event_col = "event", score_col = "score", out_pdf = NULL, title = "") {
  df <- df[!is.na(df[[time_col]]) & !is.na(df[[event_col]]) & !is.na(df[[score_col]]), , drop = FALSE]
  if (nrow(df) < 20 || length(unique(df[[event_col]])) < 2) return(NULL)

  med <- median(df[[score_col]], na.rm = TRUE)
  df$group <- ifelse(df[[score_col]] >= med, "High", "Low")
  df$group <- factor(df$group, levels = c("Low", "High"))

  fit <- survfit(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ group")), data = df)
  sd <- survdiff(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ group")), data = df)
  p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)

  cox <- coxph(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", score_col)), data = df)
  s <- summary(cox)

  out <- list(
    n = nrow(df),
    median = med,
    p_logrank = p,
    hr = unname(s$coefficients[1, "exp(coef)"]),
    ci_low = unname(s$conf.int[1, "lower .95"]),
    ci_high = unname(s$conf.int[1, "upper .95"]),
    p_cox = unname(s$coefficients[1, "Pr(>|z|)"]),
    fit = fit,
    df = df
  )

  if (!is.null(out_pdf)) {
    pdf(out_pdf, width = 7, height = 6)
    plot(fit, col = c("#D95F02", "#1B9E77"), lwd = 2,
         xlab = "Time", ylab = "Survival probability",
         main = paste0(title, "\nLog-rank p = ", signif(p, 3)))
    legend("bottomleft", legend = levels(df$group), col = c("#D95F02", "#1B9E77"), lwd = 2, bty = "n")
    dev.off()
  }

  out
}

#' 计算模块评分（基因集的z-score均值）
module_score <- function(expr, genes) {
  genes <- intersect(norm_gene(genes), rownames(expr))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr)))

  z <- t(scale(t(expr[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  colMeans(z, na.rm = TRUE)
}

#' 计算模块评分（带额外NA保护，10.txt改进版）
module_score_from_expr <- function(expr, genes) {
  genes <- intersect(norm_gene(genes), rownames(expr))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr)))

  z <- t(scale(t(expr[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  out <- colMeans(z, na.rm = TRUE)
  out[!is.finite(out)] <- NA
  out
}

#' 计算score与module score的相关性
cor_one <- function(score, module_score) {
  dd <- data.frame(score = score, module = module_score)
  dd <- dd[complete.cases(dd), , drop = FALSE]
  if (nrow(dd) < 10) {
    return(c(rho = NA_real_, p.value = NA_real_, n = nrow(dd)))
  }
  ct <- suppressWarnings(cor.test(dd$score, dd$module, method = "spearman"))
  c(rho = unname(ct$estimate), p.value = ct$p.value, n = nrow(dd))
}

#' 绘制模块相关性热图
plot_heatmap <- function(df, title, file_out) {
  df$module <- factor(df$module, levels = unique(df$module))
  df$cancer <- factor(df$cancer, levels = c("BLCA", "SKCM"))

  p <- ggplot(df, aes(x = cancer, y = module, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.na(rho), "NA", sprintf("%.2f", rho))), size = 3) +
    scale_fill_gradient2(
      low = "#D73027", mid = "white", high = "#4575B4",
      midpoint = 0, limits = c(-1, 1), na.value = "grey90"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = title, fill = "Spearman rho")

  ggsave(file_out, p, width = 7, height = max(4, 0.35 * length(unique(df$module)) + 1))
}

# ============================================================
# ============================================================
# Step 1: Data loading, cleaning, and modeling input preparation
# ============================================================
# ============================================================
# ============================================================

out_dir_p1 <- file.path(root_dir, "IMvigor210_clean")
dir.create(out_dir_p1, recursive = TRUE, showWarnings = FALSE)

rdata_file <- file.path(root_dir, "IMvigor210", "5ff67-main", "IMvigor210CoreBiologies", "IMvigor210CoreBiologies.Rdata")

# -------------------------
# 1) Load Rdata
# -------------------------
env <- new.env()
load(rdata_file, envir = env)

expreSet <- get("expreSet", envir = env)
annoData <- get("annoData", envir = env)
phenoData <- get("phenoData", envir = env)

# -------------------------
# 2) Build expression matrix (from 2.txt + 3.txt)
# -------------------------
expr_df <- as.data.frame(expreSet, check.names = FALSE)

# If the first column looks like gene IDs, use it as rownames
if (is.character(expr_df[[1]]) || is.factor(expr_df[[1]])) {
  gene_id <- as.character(expr_df[[1]])
  expr_df <- expr_df[, -1, drop = FALSE]
  rownames(expr_df) <- gene_id
} else {
}

expr_mat <- as.matrix(expr_df)
storage.mode(expr_mat) <- "numeric"

# -------------------------
# 3) Clean phenoData (from 3.txt)
# -------------------------
pheno <- phenoData

if (!is.null(rownames(pheno)) && all(nchar(rownames(pheno)) > 0)) {
  pheno$sample_id <- rownames(pheno)
} else if ("ANONPT_ID" %in% colnames(pheno)) {
  pheno$sample_id <- as.character(pheno$ANONPT_ID)
} else {
  stop("Cannot find sample ID in phenoData.")
}

pheno$sample_id <- as.character(pheno$sample_id)

# Build clean response label (from 4.txt)
pheno$response_binary_clean <- NA_character_
pheno$response_binary_clean[pheno$binaryResponse == "CR/PR"] <- "Responder"
pheno$response_binary_clean[pheno$binaryResponse == "SD/PD"] <- "NonResponder"

pheno$response_binary_clean <- factor(
  pheno$response_binary_clean,
  levels = c("NonResponder", "Responder")
)

# -------------------------
# 4) Align samples (from 3.txt)
# -------------------------
common_samples <- intersect(colnames(expr_mat), pheno$sample_id)

if (length(common_samples) == 0) {
  stop("No matching sample IDs between expression matrix and phenoData.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
pheno <- pheno[match(common_samples, pheno$sample_id), , drop = FALSE]

stopifnot(all(colnames(expr_mat) == pheno$sample_id))

# -------------------------
# 5) Map gene symbols (from 4.txt)
# -------------------------

symbol_col <- pick_col(annoData, c("symbol", "gene_name", "gene symbol", "hgnc"))
id_col <- pick_col(annoData, c("gene_id", "ensembl", "id", "probe"))

anno <- as.data.frame(annoData, check.names = FALSE)

# Decide gene symbol source
gene_symbol <- NULL
if (!is.na(symbol_col) && symbol_col %in% colnames(anno)) {
  gene_symbol <- anno[[symbol_col]]
} else if (!is.null(rownames(anno))) {
  gene_symbol <- rownames(anno)
} else if (!is.null(rownames(expr_mat))) {
  gene_symbol <- rownames(expr_mat)
} else {
  stop("Cannot infer gene symbols from annoData or expression matrix.")
}

# Match annotation to expression rows
if (nrow(anno) == nrow(expr_mat)) {
  expr_mat2 <- expr_mat
  gene_symbol2 <- gene_symbol
} else {
  expr_mat2 <- expr_mat
  gene_symbol2 <- gene_symbol[seq_len(min(length(gene_symbol), nrow(expr_mat2)))]
  if (length(gene_symbol2) < nrow(expr_mat2)) {
    gene_symbol2 <- c(gene_symbol2, rep(NA, nrow(expr_mat2) - length(gene_symbol2)))
  }
}

# Set rownames to gene symbol if possible
gene_symbol2 <- as.character(gene_symbol2)
gene_symbol2[gene_symbol2 == ""] <- NA

# Save mapping table
gene_map <- data.frame(
  row_index = seq_len(nrow(expr_mat2)),
  gene_symbol = gene_symbol2,
  stringsAsFactors = FALSE
)

write.table(gene_map, file.path(out_dir_p1, "IMvigor210_gene_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Collapse duplicated symbols
expr_symbol <- collapse_by_symbol(expr_mat2, gene_symbol2)

# -------------------------
# 6) Build analysis tables (from 4.txt)
# -------------------------
analysis_all <- data.frame(
  sample_id = pheno$sample_id,
  response_binary = pheno$binaryResponse,
  response_clean = pheno$response_binary_clean,
  response_best = pheno[["Best Confirmed Overall Response"]],
  os = pheno$os,
  censOS = pheno$censOS,
  stringsAsFactors = FALSE
)

# Response-only subset
resp_df <- analysis_all[!is.na(analysis_all$response_clean), , drop = FALSE]
resp_expr <- expr_symbol[, resp_df$sample_id, drop = FALSE]

# Survival-ready subset
surv_df <- analysis_all[!is.na(analysis_all$os) & !is.na(analysis_all$censOS), , drop = FALSE]
surv_expr <- expr_symbol[, surv_df$sample_id, drop = FALSE]

# -------------------------
# 7) Save Phase 1 outputs
# -------------------------
saveRDS(expr_symbol, file.path(out_dir_p1, "IMvigor210_expr_symbol.rds"))
saveRDS(resp_df, file.path(out_dir_p1, "IMvigor210_response_df.rds"))
saveRDS(surv_df, file.path(out_dir_p1, "IMvigor210_survival_df.rds"))
saveRDS(resp_expr, file.path(out_dir_p1, "IMvigor210_response_expr.rds"))
saveRDS(surv_expr, file.path(out_dir_p1, "IMvigor210_survival_expr.rds"))
saveRDS(pheno, file.path(out_dir_p1, "IMvigor210_pheno_aligned.rds"))

write.table(expr_symbol, file.path(out_dir_p1, "IMvigor210_expr_symbol.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(resp_df, file.path(out_dir_p1, "IMvigor210_response_df.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(surv_df, file.path(out_dir_p1, "IMvigor210_survival_df.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================

# ============================================================
# Step 2: Core exosome gene screening (ExoCarta ∩ Vesiclepedia)
# ============================================================
# ============================================================
# ============================================================

out_dir_p2 <- file.path(root_dir, "IMvigor210_core_screen")
dir.create(out_dir_p2, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load Phase 1 outputs
# -------------------------
clean_dir <- out_dir_p1  # Phase 1 output directory

expr_symbol <- readRDS(file.path(clean_dir, "IMvigor210_expr_symbol.rds"))
resp_df <- readRDS(file.path(clean_dir, "IMvigor210_response_df.rds"))

resp_df <- resp_df[!is.na(resp_df$response_clean), , drop = FALSE]
resp_df$response_clean <- factor(resp_df$response_clean, levels = c("NonResponder", "Responder"))

common_samples <- intersect(colnames(expr_symbol), resp_df$sample_id)
expr_resp <- expr_symbol[, common_samples, drop = FALSE]
resp_df <- resp_df[match(common_samples, resp_df$sample_id), , drop = FALSE]

stopifnot(all(colnames(expr_resp) == resp_df$sample_id))

# Normalize expression gene names for matching
expr_genes_norm <- norm_gene(rownames(expr_resp))

# -------------------------
# Load ExoCarta / Vesiclepedia
# -------------------------
exo_file <- file.path(root_dir, "ExoCarta_protein_mRNA_details_6.txt")
ves_file <- file.path(root_dir, "VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt")

exo <- read.delim(exo_file, sep = "\t", header = TRUE, quote = "", comment.char = "",
                  check.names = FALSE, fill = TRUE)
ves <- read.delim(ves_file, sep = "\t", header = TRUE, quote = "", comment.char = "",
                  check.names = FALSE, fill = TRUE)

# -------------------------
# Extract human gene symbols
# -------------------------
exo_human <- exo[is_human(exo$SPECIES), , drop = FALSE]
ves_human <- ves[is_human(ves$SPECIES), , drop = FALSE]

exo_genes <- norm_gene(exo_human[["GENE SYMBOL"]])
ves_genes <- norm_gene(ves_human[["GENE SYMBOL"]])

# -------------------------
# Core high-confidence gene set (ExoCarta ∩ Vesiclepedia)
# -------------------------
core_genes <- intersect(exo_genes, ves_genes)
core_genes <- core_genes[core_genes != ""]

# Overlap with IMvigor210
core_overlap <- intersect(core_genes, expr_genes_norm)

# Optional broader set for reference
union_genes <- union(exo_genes, ves_genes)
union_overlap <- intersect(union_genes, expr_genes_norm)

# -------------------------
# Save gene lists
# -------------------------
write.table(data.frame(gene = sort(exo_genes)),
            file.path(out_dir_p2, "ExoCarta_human_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(gene = sort(ves_genes)),
            file.path(out_dir_p2, "Vesiclepedia_human_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(gene = sort(core_genes)),
            file.path(out_dir_p2, "core_exosome_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(gene = sort(core_overlap)),
            file.path(out_dir_p2, "core_exosome_genes_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# Build core expression matrix
# -------------------------
core_idx <- match(core_overlap, expr_genes_norm)
core_idx <- core_idx[!is.na(core_idx)]

core_expr <- expr_resp[core_idx, , drop = FALSE]
rownames(core_expr) <- rownames(expr_resp)[core_idx]

saveRDS(core_expr, file.path(out_dir_p2, "core_exosome_expr.rds"))
write.table(core_expr, file.path(out_dir_p2, "core_exosome_expr.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# -------------------------
# Response screening on core genes (Wilcoxon)
# -------------------------
screen_res <- wilcox_screen(core_expr, resp_df$response_clean)

write.table(screen_res, file.path(out_dir_p2, "core_response_screen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(screen_res, file.path(out_dir_p2, "core_response_screen.rds"))

top20 <- head(screen_res, 20)
write.table(top20, file.path(out_dir_p2, "core_response_screen_top20.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# Summary
# -------------------------
summary_file <- file.path(out_dir_p2, "core_screen_summary.txt")

# ============================================================

# ============================================================
# Step 3: LASSO logistic regression and ExoImmune Score construction
# ============================================================
# ============================================================
#         IMvigor210_response_df.rds (Phase 1)
# ============================================================

out_dir_p3 <- file.path(root_dir, "IMvigor210_lasso")
dir.create(out_dir_p3, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

# -------------------------
# Load data
# -------------------------
in_dir_p2 <- out_dir_p2

expr <- readRDS(file.path(in_dir_p2, "core_exosome_expr.rds"))
resp <- readRDS(file.path(out_dir_p1, "IMvigor210_response_df.rds"))

resp <- resp[!is.na(resp$response_clean), , drop = FALSE]
resp$response_clean <- factor(resp$response_clean, levels = c("NonResponder", "Responder"))

common_samples <- intersect(colnames(expr), resp$sample_id)
expr <- expr[, common_samples, drop = FALSE]
resp <- resp[match(common_samples, resp$sample_id), , drop = FALSE]

stopifnot(all(colnames(expr) == resp$sample_id))

# -------------------------
# Load univariate screen results
# -------------------------
screen_file <- file.path(in_dir_p2, "core_response_screen.tsv")
screen <- read.delim(screen_file, sep = "\t", header = TRUE, check.names = FALSE)
screen$gene <- as.character(screen$gene)

# -------------------------
# Candidate filtering
# -------------------------
top_n <- 50

screen2 <- screen[!is.na(screen$p.value), , drop = FALSE]
screen2 <- screen2[order(screen2$padj, screen2$p.value), , drop = FALSE]

sig_genes <- screen2$gene[screen2$padj < 0.1]
top_genes <- head(screen2$gene, top_n)

candidate_genes <- unique(c(sig_genes, top_genes))
candidate_genes <- intersect(candidate_genes, rownames(expr))

# If too few genes, broaden slightly
if (length(candidate_genes) < 10) {
  candidate_genes <- unique(head(screen2$gene, 100))
  candidate_genes <- intersect(candidate_genes, rownames(expr))
}

# -------------------------
# Build modeling matrix
# -------------------------
x <- expr[candidate_genes, , drop = FALSE]

# Z-score each gene across samples
x_z <- t(scale(t(as.matrix(x))))
x_z[is.na(x_z)] <- 0

# Remove zero-variance genes after scaling
keep <- apply(x_z, 1, sd, na.rm = TRUE) > 0
x_z <- x_z[keep, , drop = FALSE]

# Response vector
y <- ifelse(resp$response_clean == "Responder", 1, 0)

# Transpose for glmnet: samples x genes
x_glm <- t(x_z)

# -------------------------
# LASSO Cross-Validation
# -------------------------
set.seed(2026)
cvfit <- cv.glmnet(
  x = x_glm,
  y = y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance",
  standardize = FALSE
)

saveRDS(cvfit, file.path(out_dir_p3, "cvfit.rds"))

# Plot CV curve
pdf(file.path(out_dir_p3, "lasso_cv_curve.pdf"), width = 7, height = 6)
plot(cvfit)
dev.off()

# -------------------------
# Extract coefficients
# -------------------------
coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))
coef_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))

coef_df_min <- data.frame(
  gene = rownames(coef_min),
  coef = as.numeric(coef_min[, 1]),
  stringsAsFactors = FALSE
)
coef_df_min <- coef_df_min[coef_df_min$coef != 0, , drop = FALSE]
coef_df_min <- coef_df_min[order(-abs(coef_df_min$coef)), , drop = FALSE]

coef_df_1se <- data.frame(
  gene = rownames(coef_1se),
  coef = as.numeric(coef_1se[, 1]),
  stringsAsFactors = FALSE
)
coef_df_1se <- coef_df_1se[coef_df_1se$coef != 0, , drop = FALSE]
coef_df_1se <- coef_df_1se[order(-abs(coef_df_1se$coef)), , drop = FALSE]

write.table(coef_df_min, file.path(out_dir_p3, "lasso_coef_lambda_min.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(coef_df_1se, file.path(out_dir_p3, "lasso_coef_lambda_1se.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# Build ExoImmune Score
# -------------------------
# Prefer lambda.1se for parsimony
chosen <- coef_df_1se
if (nrow(chosen) < 2) chosen <- coef_df_min

chosen_genes <- setdiff(chosen$gene, "(Intercept)")
chosen_coefs <- chosen$coef[chosen$gene != "(Intercept)"]

x_score <- t(x_z[chosen_genes, , drop = FALSE])
score <- as.numeric(x_score %*% chosen_coefs)

# Save sample scores
score_df <- data.frame(
  sample_id = resp$sample_id,
  response_clean = resp$response_clean,
  response_binary = resp$response_binary,
  score = score,
  stringsAsFactors = FALSE
)

write.table(score_df, file.path(out_dir_p3, "exoimmune_score.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(score_df, file.path(out_dir_p3, "exoimmune_score.rds"))

# -------------------------
# ROC / AUC
# -------------------------
roc_obj <- pROC::roc(response = score_df$response_clean, predictor = score_df$score,
                     levels = c("NonResponder", "Responder"), direction = "<")

auc_val <- as.numeric(pROC::auc(roc_obj))

pdf(file.path(out_dir_p3, "exoimmune_score_roc.pdf"), width = 6, height = 6)
plot(roc_obj, col = "#2C7FB8", lwd = 2, main = paste0("ExoImmune Score ROC (AUC=", round(auc_val, 3), ")"))
abline(a = 0, b = 1, lty = 2, col = "grey60")
dev.off()

# -------------------------
# Response group comparison
# -------------------------
pdf(file.path(out_dir_p3, "exoimmune_score_boxplot.pdf"), width = 6, height = 5)
boxplot(score ~ response_clean, data = score_df,
        col = c("#D95F02", "#1B9E77"),
        ylab = "ExoImmune Score",
        xlab = "",
        main = "ExoImmune Score by Response")
stripchart(score ~ response_clean, data = score_df,
           vertical = TRUE, method = "jitter", add = TRUE,
           pch = 16, col = rgb(0, 0, 0, 0.4))
dev.off()

# -------------------------
# Save summary
# -------------------------
summary_txt <- file.path(out_dir_p3, "lasso_summary.txt")

# ============================================================

# ============================================================
# Step 4: Survival analysis -- Kaplan-Meier and Cox regression
# ============================================================
# ============================================================
# ============================================================
#       * to_num() 安全数字解析
#       * 年龄/ECOG/性别清洗逻辑
#       * 完整案例分析(complete-case)多变量Cox
# ============================================================

out_dir_p4 <- file.path(root_dir, "IMvigor210_survival")
dir.create(out_dir_p4, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(survival)
})

# -------------------------
# Load score and phenotype
# -------------------------
score_file <- file.path(out_dir_p3, "exoimmune_score.rds")
pheno_file <- file.path(out_dir_p1, "IMvigor210_pheno_aligned.rds")

if (!file.exists(score_file)) stop("Missing file: ", score_file)
if (!file.exists(pheno_file)) stop("Missing file: ", pheno_file)

score_df <- readRDS(score_file)
pheno <- readRDS(pheno_file)

# -------------------------
# Merge data
# -------------------------
pheno$sample_id <- as.character(pheno$sample_id)
score_df$sample_id <- as.character(score_df$sample_id)

df <- merge(score_df, pheno, by = "sample_id", all.x = TRUE)

# -------------------------
# Survival coding (from 7.5.txt)
# -------------------------
# We use event = 1 - censOS so that 1 means death/event occurred.
df$os <- as.numeric(df$os)
df$censOS <- as.numeric(df$censOS)
df$event_os <- ifelse(df$censOS == 0, 1, 0)

# Remove missing survival rows
df <- df[!is.na(df$os) & !is.na(df$event_os) & !is.na(df$score), , drop = FALSE]

# -------------------------
# Group by median score
# -------------------------
median_score <- median(df$score, na.rm = TRUE)
df$score_group <- ifelse(df$score >= median_score, "High", "Low")
df$score_group <- factor(df$score_group, levels = c("Low", "High"))

# -------------------------
# Kaplan-Meier analysis (from 7.5.txt)
# -------------------------
fit_km <- survfit(Surv(os, event_os) ~ score_group, data = df)
logrank <- survdiff(Surv(os, event_os) ~ score_group, data = df)
p_logrank <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

# Save KM plot
pdf(file.path(out_dir_p4, "exoimmune_score_km.pdf"), width = 7, height = 6)
plot(fit_km,
     col = c("#D95F02", "#1B9E77"),
     lwd = 2,
     xlab = "Overall survival time",
     ylab = "Survival probability",
     main = paste0("IMvigor210 KM by ExoImmune Score\nLog-rank p = ", signif(p_logrank, 3)))
legend("bottomleft",
       legend = levels(df$score_group),
       col = c("#D95F02", "#1B9E77"),
       lwd = 2,
       bty = "n")
dev.off()

# -------------------------
# Univariate Cox regression (from 7.5.txt)
# -------------------------
cox_uni <- coxph(Surv(os, event_os) ~ score, data = df)
sum_uni <- summary(cox_uni)

hr <- sum_uni$coefficients[1, "exp(coef)"]
ci_low <- sum_uni$conf.int[1, "lower .95"]
ci_high <- sum_uni$conf.int[1, "upper .95"]
p_uni <- sum_uni$coefficients[1, "Pr(>|z|)"]

# Save univariate Cox results
cox_uni_df <- data.frame(
  model = "univariate",
  variable = "score",
  HR = hr,
  CI_low = ci_low,
  CI_high = ci_high,
  p.value = p_uni,
  stringsAsFactors = FALSE
)

write.table(cox_uni_df,
            file.path(out_dir_p4, "cox_univariate_score.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# Improved multivariate Cox (from 8.txt)
# -------------------------
# Clean core variables
df$score <- as.numeric(df$score)
df$os <- as.numeric(df$os)
df$event_os <- as.numeric(df$event_os)

# Age cleaning with to_num()
df$age_clean <- to_num(df[["Sample age"]])

# ECOG cleaning with to_num()
df$ecog_clean <- to_num(df[["Baseline ECOG Score"]])

# Sex cleaning
sex <- as.character(df[["Sex"]])
sex <- trimws(sex)
sex[sex == ""] <- NA
sex <- toupper(sex)
df$sex_clean <- ifelse(grepl("^M", sex), "Male",
                       ifelse(grepl("^F", sex), "Female", NA))
df$sex_clean <- factor(df$sex_clean, levels = c("Female", "Male"))

# Show cleaned age summary

# Practical multivariable Cox with complete-case analysis
vars <- c("score", "ecog_clean", "sex_clean", "age_clean")
complete_vars <- df[, c("os", "event_os", vars), drop = FALSE]
complete_vars <- complete.cases(complete_vars)

df2 <- df[complete_vars, , drop = FALSE]

cox_multi_df <- NULL
if (nrow(df2) >= 30) {
  # Build formula depending on what survived
  rhs <- c("score")
  if (sum(!is.na(df2$ecog_clean)) > 0) rhs <- c(rhs, "ecog_clean")
  if (sum(!is.na(df2$sex_clean)) > 0) rhs <- c(rhs, "sex_clean")
  if (sum(!is.na(df2$age_clean)) > 0) rhs <- c(rhs, "age_clean")

  fml <- as.formula(paste("Surv(os, event_os) ~", paste(rhs, collapse = " + ")))
  cox_multi <- tryCatch(coxph(fml, data = df2), error = function(e) e)

  if (!inherits(cox_multi, "error")) {
    sum_multi <- summary(cox_multi)
    coef_tab <- as.data.frame(sum_multi$coefficients)
    conf_tab <- as.data.frame(sum_multi$conf.int)

    cox_multi_df <- data.frame(
      variable = rownames(coef_tab),
      HR = conf_tab$`exp(coef)`,
      CI_low = conf_tab$`lower .95`,
      CI_high = conf_tab$`upper .95`,
      p.value = coef_tab$`Pr(>|z|)`,
      stringsAsFactors = FALSE
    )

    write.table(cox_multi_df,
                file.path(out_dir_p4, "cox_multivariate.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

  } else {
  }
} else {
}

# -------------------------
# Score distribution plot (from 7.5.txt)
# -------------------------
pdf(file.path(out_dir_p4, "exoimmune_score_boxplot.pdf"), width = 6, height = 5)
boxplot(score ~ score_group, data = df,
        col = c("#D95F02", "#1B9E77"),
        ylab = "ExoImmune Score",
        xlab = "",
        main = "ExoImmune Score by Survival Group")
stripchart(score ~ score_group, data = df,
           vertical = TRUE, method = "jitter", add = TRUE,
           pch = 16, col = rgb(0, 0, 0, 0.35))
dev.off()

# -------------------------
# Save merged analysis table (from 7.5.txt)
# -------------------------
write.table(df,
            file.path(out_dir_p4, "imvigor210_score_survival_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(df, file.path(out_dir_p4, "imvigor210_score_survival_merged.rds"))

# -------------------------
# Summary file (from 7.5.txt)
# -------------------------
summary_file <- file.path(out_dir_p4, "survival_summary.txt")

if (!is.null(cox_multi_df)) {
}

# ============================================================

# ============================================================
# Step 5: TCGA-BLCA and TCGA-SKCM external validation
# ============================================================
# ============================================================
# ============================================================

out_dir_p5 <- file.path(root_dir, "TCGA_mechanism_prognosis")
dir.create(out_dir_p5, recursive = TRUE, showWarnings = FALSE)

tcga_root <- file.path(root_dir, "EV_immunotherapy_project", "TCGA_processed")

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

# -------------------------
# TCGA bundle reader
# -------------------------
read_tcga_bundle <- function(cancer) {
  expr_file <- find_file(tcga_root, paste0("TCGA-", cancer, ".*log2TPM_symbol_tumor\\.tsv$"))
  if (length(expr_file) == 0) expr_file <- find_file(tcga_root, paste0("TCGA-", cancer, ".*log2TPM_symbol\\.tsv$"))
  clin_file <- find_file(tcga_root, paste0("TCGA-", cancer, ".*clinical\\.tsv$"))

  expr_file <- if (length(expr_file) > 0) expr_file[1] else NA_character_
  clin_file <- if (length(clin_file) > 0) clin_file[1] else NA_character_

  list(expr_file = expr_file, clin_file = clin_file)
}

# -------------------------
# TCGA survival detection
# -------------------------
detect_survival <- function(clin, sample_id_col = NULL) {
  time_col <- pick_col(clin, c("^os$", "overall_survival", "survival_time", "days_to_death", "days_to_last_follow_up", "days_to_last_contact"))
  event_col <- pick_col(clin, c("^event$", "vital_status", "os_status", "death", "censor"))

  time <- rep(NA_real_, nrow(clin))
  event <- rep(NA_integer_, nrow(clin))

  if (!is.na(time_col) && !is.na(event_col)) {
    if (grepl("vital_status", event_col, ignore.case = TRUE)) {
      time <- to_num(clin[[time_col]])
      v <- toupper(trimws(as.character(clin[[event_col]])))
      event <- ifelse(v %in% c("DEAD", "DECEASED", "1"), 1, ifelse(v %in% c("ALIVE", "0"), 0, NA))
    } else {
      time <- to_num(clin[[time_col]])
      event <- to_num(clin[[event_col]])
    }
  }

  # TCGA standard fallback
  if (all(is.na(time)) || all(is.na(event))) {
    dd <- pick_col(clin, c("days_to_death"))
    dfu <- pick_col(clin, c("days_to_last_follow_up", "days_to_last_contact"))
    vital <- pick_col(clin, c("vital_status"))

    if (!is.na(dd) || !is.na(dfu)) {
      death_days <- if (!is.na(dd)) to_num(clin[[dd]]) else rep(NA_real_, nrow(clin))
      follow_days <- if (!is.na(dfu)) to_num(clin[[dfu]]) else rep(NA_real_, nrow(clin))
      v <- if (!is.na(vital)) toupper(trimws(as.character(clin[[vital]]))) else rep(NA_character_, nrow(clin))

      event <- ifelse(!is.na(death_days), 1, ifelse(v %in% c("DEAD", "DECEASED"), 1, 0))
      time <- ifelse(event == 1, ifelse(!is.na(death_days), death_days, follow_days), follow_days)
    }
  }

  data.frame(time = time, event = event, stringsAsFactors = FALSE)
}

# -------------------------
# Single TCGA cancer analysis
# -------------------------
analyze_one_tcga <- function(cancer) {
  bundle <- read_tcga_bundle(cancer)
  if (is.na(bundle$expr_file) || is.na(bundle$clin_file)) {
    message("Missing files for ", cancer)
    return(NULL)
  }

  expr <- read_expr_table(bundle$expr_file)
  expr <- collapse_dups(expr)

  clin <- read_clin_table(bundle$clin_file)

  # Build patient/sample mapping
  expr_ids <- colnames(expr)
  expr_pat <- tcga_patient_id(expr_ids)

  clin_id_col <- pick_col(clin, c("barcode", "sample", "patient", "case", "submitter", "aliquot"))
  if (is.na(clin_id_col)) {
    cand <- clin[[1]]
    if (sum(grepl("^TCGA-", as.character(cand)), na.rm = TRUE) > 0) {
      clin_id_col <- colnames(clin)[1]
    }
  }
  if (is.na(clin_id_col)) {
    stop("Could not detect TCGA ID column in clinical file for ", cancer)
  }

  clin$patient_id <- tcga_patient_id(clin[[clin_id_col]])
  clin <- clin[!is.na(clin$patient_id) & clin$patient_id != "", , drop = FALSE]

  # Score genes
  score <- make_score(expr, coef_df, scale_tbl)
  score_df <- data.frame(
    sample_id = expr_ids,
    patient_id = expr_pat,
    score = score,
    stringsAsFactors = FALSE
  )

  # Merge clinical at patient level
  clin_surv <- detect_survival(clin)
  clin2 <- cbind(clin, clin_surv)

  merged <- merge(score_df, clin2, by = "patient_id", all.x = TRUE)
  merged <- merged[!is.na(merged$score), , drop = FALSE]

  # If multiple clinical rows map to same patient, keep first complete-ish row
  merged <- merged[!duplicated(merged$sample_id), , drop = FALSE]

  # Summary tables
  write.table(score_df, file.path(out_dir_p5, paste0("TCGA-", cancer, "_score.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(merged, file.path(out_dir_p5, paste0("TCGA-", cancer, "_merged.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(merged, file.path(out_dir_p5, paste0("TCGA-", cancer, "_merged.rds")))

  # Survival analysis
  surv_res <- NULL
  if (any(!is.na(merged$time)) && any(!is.na(merged$event))) {
    merged$time <- as.numeric(merged$time)
    merged$event <- as.numeric(merged$event)

    surv_res <- median_km(
      merged,
      time_col = "time",
      event_col = "event",
      score_col = "score",
      out_pdf = file.path(out_dir_p5, paste0("TCGA-", cancer, "_KM.pdf")),
      title = paste0("TCGA-", cancer, " ExoImmune Score")
    )

    if (!is.null(surv_res)) {
      cox_df <- data.frame(
        cancer = cancer,
        n = surv_res$n,
        HR = surv_res$hr,
        CI_low = surv_res$ci_low,
        CI_high = surv_res$ci_high,
        p_cox = surv_res$p_cox,
        p_logrank = surv_res$p_logrank,
        stringsAsFactors = FALSE
      )
      write.table(cox_df, file.path(out_dir_p5, paste0("TCGA-", cancer, "_survival_summary.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  # Mechanism modules
  modules <- list(
    checkpoint = c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "PDCD1LG2", "BTLA"),
    antigen_presentation = c("HLAA", "HLAB", "HLAC", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "CIITA"),
    cytolytic = c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY"),
    exhaustion = c("TOX", "TOX2", "CXCL13", "EOMES", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
    t_cell_inflamed = c("CXCL9", "CXCL10", "CXCL11", "CD8A", "CD8B", "IFNG", "STAT1", "GBP1")
  )

  mech <- data.frame(sample_id = expr_ids, patient_id = expr_pat, score = score, stringsAsFactors = FALSE)
  for (nm in names(modules)) {
    mech[[nm]] <- module_score(expr, modules[[nm]])
  }

  # Correlations
  cor_rows <- lapply(names(modules), function(nm) {
    dd <- mech[, c("score", nm)]
    dd <- dd[complete.cases(dd), , drop = FALSE]
    if (nrow(dd) < 10) return(NULL)
    ct <- suppressWarnings(cor.test(dd$score, dd[[nm]], method = "spearman"))
    data.frame(
      cancer = cancer,
      module = nm,
      rho = unname(ct$estimate),
      p.value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  cor_df <- do.call(rbind, cor_rows)
  if (!is.null(cor_df)) {
    cor_df$padj <- p.adjust(cor_df$p.value, method = "BH")
    write.table(cor_df, file.path(out_dir_p5, paste0("TCGA-", cancer, "_module_correlations.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Save mechanism table
  write.table(mech, file.path(out_dir_p5, paste0("TCGA-", cancer, "_mechanism_scores.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(mech, file.path(out_dir_p5, paste0("TCGA-", cancer, "_mechanism_scores.rds")))

  # Simple scatter plot for best module if available
  if (!is.null(cor_df) && nrow(cor_df) > 0) {
    best_mod <- cor_df$module[which.min(cor_df$p.value)][1]
    p <- ggplot(mech, aes(x = score, y = .data[[best_mod]])) +
      geom_point(size = 1, alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "#D95F02") +
      theme_bw() +
      labs(
        title = paste0("TCGA-", cancer, ": score vs ", best_mod),
        x = "ExoImmune Score",
        y = best_mod
      )
    ggsave(file.path(out_dir_p5, paste0("TCGA-", cancer, "_", best_mod, "_scatter.pdf")), p, width = 6, height = 5)
  }

  list(
    score_df = score_df,
    merged = merged,
    surv = surv_res,
    cor = cor_df,
    mech = mech
  )
}

# -------------------------
# Load model coefficients / scaling (from 9.txt)
# -------------------------
coef_file <- file.path(out_dir_p3, "lasso_coef_lambda_1se.tsv")
if (!file.exists(coef_file)) coef_file <- file.path(out_dir_p3, "lasso_coef_lambda_min.tsv")
if (!file.exists(coef_file)) stop("Cannot find LASSO coefficient file.")

coef_df <- read.delim(coef_file, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
coef_df$gene <- toupper(trimws(as.character(coef_df$gene)))
coef_df$coef <- as.numeric(coef_df$coef)

train_expr_file <- file.path(out_dir_p1, "IMvigor210_expr_symbol.rds")
if (!file.exists(train_expr_file)) stop("Cannot find training expression matrix.")
train_expr <- readRDS(train_expr_file)
train_expr <- as.matrix(train_expr)
suppressWarnings(storage.mode(train_expr) <- "numeric")
rownames(train_expr) <- norm_gene(rownames(train_expr))

selected_genes <- setdiff(coef_df$gene, "(Intercept)")
selected_genes <- intersect(selected_genes, rownames(train_expr))
if (length(selected_genes) == 0) stop("No selected genes overlap with training expression.")

train_sub <- train_expr[selected_genes, , drop = FALSE]
scale_tbl <- data.frame(
  gene = selected_genes,
  mean = rowMeans(train_sub, na.rm = TRUE),
  sd = apply(train_sub, 1, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
scale_tbl$sd[!is.finite(scale_tbl$sd) | scale_tbl$sd == 0] <- 1

write.table(scale_tbl, file.path(out_dir_p5, "training_scaling_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# Run TCGA analyses
# -------------------------
res_blca <- analyze_one_tcga("BLCA")
res_skcm <- analyze_one_tcga("SKCM")

# Combined summary
summary_list <- list()
if (!is.null(res_blca$surv)) {
  summary_list[[length(summary_list) + 1]] <- data.frame(
    cancer = "BLCA",
    n = res_blca$surv$n,
    HR = res_blca$surv$hr,
    CI_low = res_blca$surv$ci_low,
    CI_high = res_blca$surv$ci_high,
    p_cox = res_blca$surv$p_cox,
    p_logrank = res_blca$surv$p_logrank,
    stringsAsFactors = FALSE
  )
}
if (!is.null(res_skcm$surv)) {
  summary_list[[length(summary_list) + 1]] <- data.frame(
    cancer = "SKCM",
    n = res_skcm$surv$n,
    HR = res_skcm$surv$hr,
    CI_low = res_skcm$surv$ci_low,
    CI_high = res_skcm$surv$ci_high,
    p_cox = res_skcm$surv$p_cox,
    p_logrank = res_skcm$surv$p_logrank,
    stringsAsFactors = FALSE
  )
}
if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.table(summary_df, file.path(out_dir_p5, "TCGA_survival_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ============================================================

# ============================================================
# Step 6: Immune mechanism interrogation via module correlation
# ============================================================
# ============================================================
#         TCGA_processed 表达矩阵
# ============================================================

out_dir_p6 <- file.path(root_dir, "TCGA_mechanism_summary")
dir.create(out_dir_p6, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

# -------------------------
# Modules definition
# -------------------------
immune_modules <- list(
  checkpoint = c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "PDCD1LG2", "BTLA"),
  antigen_presentation = c("HLAA", "HLAB", "HLAC", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "CIITA"),
  cytolytic = c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY"),
  exhaustion = c("TOX", "TOX2", "CXCL13", "EOMES", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
  t_cell_inflamed = c("CXCL9", "CXCL10", "CXCL11", "CD8A", "CD8B", "IFNG", "STAT1", "GBP1")
)

control_modules <- list(
  proliferation = c("MKI67", "TOP2A", "CDK1", "CCNB1", "UBE2C", "BIRC5", "PCNA", "MCM2"),
  stroma_ecm = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "SPARC", "ACTA2", "TAGLN"),
  glycolysis = c("HK2", "PFKP", "ALDOA", "GAPDH", "PKM", "LDHA", "ENO1"),
  hypoxia = c("HIF1A", "CA9", "VEGFA", "SLC2A1", "PDK1", "LDHA", "ADM"),
  epithelial_mesenchymal = c("VIM", "ZEB1", "TWIST1", "SNAI1", "SNAI2", "FN1", "ITGA5", "COL1A1")
)

all_modules <- c(immune_modules, control_modules)
module_info <- data.frame(
  module = names(all_modules),
  category = c(rep("immune", length(immune_modules)), rep("control", length(control_modules))),
  stringsAsFactors = FALSE
)

# -------------------------
# Load model score and mapping (from Phase 5)
# -------------------------
mech_dir <- out_dir_p5

blca_mech <- readRDS(file.path(mech_dir, "TCGA-BLCA_mechanism_scores.rds"))
skcm_mech <- readRDS(file.path(mech_dir, "TCGA-SKCM_mechanism_scores.rds"))

# -------------------------
# Load expression
# -------------------------
blca_expr_file <- find_file(tcga_root, "TCGA-BLCA.*log2TPM_symbol_tumor\\.tsv$")
skcm_expr_file <- find_file(tcga_root, "TCGA-SKCM.*log2TPM_symbol_tumor\\.tsv$")

if (length(blca_expr_file) == 0 || length(skcm_expr_file) == 0) {
  stop("Could not find TCGA expression files.")
}

blca_expr <- read_expr_table(blca_expr_file[1])
skcm_expr <- read_expr_table(skcm_expr_file[1])

# -------------------------
# Compute modules and correlations
# -------------------------
run_one <- function(cancer, expr, mech) {
  if (!("sample_id" %in% colnames(mech)) || !("score" %in% colnames(mech))) {
    stop("Mechanism file for ", cancer, " must contain sample_id and score.")
  }

  common <- intersect(colnames(expr), as.character(mech$sample_id))
  expr2 <- expr[, common, drop = FALSE]
  mech2 <- mech[match(common, mech$sample_id), , drop = FALSE]

  stopifnot(all(colnames(expr2) == mech2$sample_id))

  out_list <- lapply(names(all_modules), function(mod) {
    ms <- module_score_from_expr(expr2, all_modules[[mod]])
    cc <- cor_one(mech2$score, ms)
    data.frame(
      cancer = cancer,
      module = mod,
      category = module_info$category[match(mod, module_info$module)],
      rho = unname(cc["rho"]),
      p.value = unname(cc["p.value"]),
      n = unname(cc["n"]),
      stringsAsFactors = FALSE
    )
  })

  mod_df <- do.call(rbind, out_list)
  mod_df$padj <- p.adjust(mod_df$p.value, method = "BH")
  mod_df <- mod_df[order(mod_df$category, mod_df$module), , drop = FALSE]

  # Save per-cancer table
  write.table(mod_df,
              file.path(out_dir_p6, paste0("TCGA-", cancer, "_module_correlations_all.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  mod_df
}

blca_df <- run_one("BLCA", blca_expr, blca_mech)
skcm_df <- run_one("SKCM", skcm_expr, skcm_mech)

all_df <- rbind(blca_df, skcm_df)
write.table(all_df, file.path(out_dir_p6, "module_correlation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# Split and plot heatmaps
# -------------------------
immune_df <- all_df[all_df$category == "immune", , drop = FALSE]
control_df <- all_df[all_df$category == "control", , drop = FALSE]

# Order modules nicely
immune_order <- names(immune_modules)
control_order <- names(control_modules)

immune_df$module <- factor(immune_df$module, levels = immune_order)
control_df$module <- factor(control_df$module, levels = control_order)
all_df$module <- factor(all_df$module, levels = c(immune_order, control_order))

plot_heatmap(immune_df, "ExoImmune Score vs Immune Modules", file.path(out_dir_p6, "immune_module_heatmap.pdf"))
plot_heatmap(control_df, "ExoImmune Score vs Non-immune Control Modules", file.path(out_dir_p6, "control_module_heatmap.pdf"))
plot_heatmap(all_df, "ExoImmune Score vs All Modules", file.path(out_dir_p6, "all_module_heatmap.pdf"))

# -------------------------
# Compare immune vs control
# -------------------------
avg_df <- aggregate(abs(rho) ~ cancer + category, data = all_df, FUN = mean, na.rm = TRUE)
colnames(avg_df)[3] <- "mean_abs_rho"
write.table(avg_df, file.path(out_dir_p6, "immune_vs_control_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p2 <- ggplot(avg_df, aes(x = category, y = mean_abs_rho, fill = category)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  facet_wrap(~ cancer) +
  theme_bw(base_size = 11) +
  ylim(0, 1) +
  labs(title = "Mean absolute correlation: immune vs control", x = NULL, y = "Mean |rho|") +
  scale_fill_manual(values = c("immune" = "#1B9E77", "control" = "#D95F02"))

ggsave(file.path(out_dir_p6, "immune_vs_control_barplot.pdf"), p2, width = 7, height = 4)

# -------------------------
# Compact notes
# -------------------------
