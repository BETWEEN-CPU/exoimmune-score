# ============================================================
# 模块3: 外部免疫治疗队列验证
# ============================================================
# 目的:
#   在独立的 ICI 队列上验证 ExoImmune Score 的疗效预测价值
#   —— 这是原管线最大的缺口 (原"内部验证"与训练集同源)
#
# 本脚本提供两条路径:
#   A) 通用模板: 你已有处理好的外部队列表达矩阵 + 疗效标签时
#   B) 具体示例: 用 GEOquery 自动下载 GSE91061 (Riaz 2017, 黑色素瘤
#      anti-PD-1, n~51 治疗前样本, 有疗效注释) 并验证
#
# 其他候选队列 (可替换):
#   - GSE78220  (Hugo 2016, 黑色素瘤 anti-PD-1, n~28)
#   - GSE135222 (黑色素瘤 anti-PD-1/anti-CTLA4)
#   - PRJEB23709 (Gide 2019, 需申请/整理)
#   - IMmotion150/151 (肾癌 atezo, 见文献补充数据)
#
# 注意:
#   1) 评分系数与缩放参数必须来自训练集 (模块0 资产), 严禁在外部数据上重估
#   2) 外部队列若是 RNA-seq, 表达需 log 转换 (log2(x+1) 或 log-TPM),
#      与 IMvigor210 (log2 normalized) 保持可比
#   3) 疗效标签统一映射: CR/PR -> Responder(1), SD/PD -> NonResponder(0)
#      外部验证时应固定使用训练集中位数作为 High/Low 切点 (若做 KM)
# ============================================================

root_dir <- Sys.getenv("EXOIMMUNE_ROOT", unset = "D:/善意无穷倒退/李新月/EV")

out_dir <- file.path(root_dir, "ExoImmune_external_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(pROC)
})

# -------------------------
# 0) 加载评分资产与通用评分函数
# -------------------------
assets <- readRDS(file.path(root_dir, "ExoImmune_scoring_assets", "exoimmune_scoring_assets.rds"))
score_new_data <- readRDS(file.path(root_dir, "ExoImmune_scoring_assets", "score_new_data.rds"))

# 训练集切点 (用于外部 KM 分组)
train_score <- readRDS(file.path(root_dir, "IMvigor210_lasso", "exoimmune_score.rds"))
train_median <- median(train_score$score, na.rm = TRUE)
cat("训练集中位数切点:", round(train_median, 4), "\n")

# ============================================================
# 通用验证函数
# ============================================================
#' @param expr_mat  基因 x 样本 表达矩阵 (rownames = gene symbol, 已 log 转换)
#' @param response  named vector 或 data.frame(sample_id, y): 1=Responder, 0=NonResponder
#' @param cohort    队列名
validate_external <- function(expr_mat, response, cohort = "external") {

  scores <- score_new_data(expr_mat, assets)
  sc <- data.frame(sample_id = colnames(expr_mat), score = scores, stringsAsFactors = FALSE)

  if (is.data.frame(response)) {
    df <- merge(sc, response, by = "sample_id")
  } else {
    df <- data.frame(sample_id = names(response), y = as.integer(response), stringsAsFactors = FALSE)
    df <- merge(sc, df, by = "sample_id")
  }
  df <- df[!is.na(df$y), , drop = FALSE]

  cat(sprintf("[%s] 有效样本: %d (Responder: %d)\n", cohort, nrow(df), sum(df$y)))

  if (nrow(df) < 20 || length(unique(df$y)) < 2) {
    cat("样本不足或单类别, 跳过\n")
    return(NULL)
  }

  # --- AUC + 置换/精确 CI ---
  r <- pROC::roc(df$y, df$score, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(r))
  ci <- pROC::ci.auc(r, method = "delong")
  cat(sprintf("[%s] AUC = %.3f (95%% CI %.3f-%.3f)\n", cohort, auc_val, ci[1], ci[3]))

  # --- Wilcoxon: Responder vs NonResponder ---
  wt <- wilcox.test(score ~ y, data = df)
  cat(sprintf("[%s] Wilcoxon p = %.4g\n", cohort, wt$p.value))

  # --- 输出 ---
  res <- data.frame(
    cohort = cohort, n = nrow(df), responders = sum(df$y),
    auc = auc_val, auc_low = ci[1], auc_high = ci[3],
    wilcox_p = wt$p.value,
    median_score = median(df$score),
    train_median_cutoff = train_median,
    stringsAsFactors = FALSE
  )
  write.table(res, file.path(out_dir, paste0("validation_", cohort, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(df, file.path(out_dir, paste0("scores_", cohort, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(out_dir, paste0("boxplot_", cohort, ".pdf")), width = 5.5, height = 5)
  boxplot(score ~ y, data = df, col = c("#D95F02", "#1B9E77"),
          names = c("NonResponder", "Responder"),
          ylab = "ExoImmune Score", main = paste0(cohort, " (AUC=", round(auc_val, 3), ")"))
  stripchart(score ~ y, data = df, vertical = TRUE, method = "jitter",
             add = TRUE, pch = 16, col = rgb(0, 0, 0, 0.4))
  dev.off()

  invisible(list(result = res, scores = df, roc = r))
}


# ============================================================
# 路径 A: 通用模板 —— 使用你自己整理的外部队列
# ============================================================
# 示例 (取消注释并替换文件路径):
#
# expr_ext <- read.delim("你的队列_expr_symbol.tsv", sep = "\t", check.names = FALSE)
# rownames(expr_ext) <- expr_ext[[1]]; expr_ext <- expr_ext[, -1]
# expr_ext <- as.matrix(expr_ext); storage.mode(expr_ext) <- "numeric"
# expr_ext <- log2(expr_ext + 1)   # 若是原始 counts/TPM, 务必 log 转换
#
# resp_ext <- read.delim("你的队列_response.tsv", sep = "\t")
# # resp_ext 需含 sample_id 和 y (1/0) 两列
#
# validate_external(expr_ext, resp_ext, cohort = "MyCohort")


# ============================================================
# 路径 B: 自动下载 GSE91061 (Riaz 2017 黑色素瘤 anti-PD-1)
# ============================================================
run_gse91061 <- function() {

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    cat("需要 GEOquery: BiocManager::install('GEOquery')\n")
    return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(GEOquery))

  cat("下载 GSE91061 ...\n")
  gse <- getGEO("GSE91061", GSEMatrix = TRUE)
  eset <- gse[[1]]

  expr_ext <- exprs(eset)   # 已为 log2 转换的表达谱 (RNA-seq processed)
  pheno_ext <- pData(eset)

  cat("外部队列维度:", dim(expr_ext), "\n")
  cat("pheno 列名:\n"); print(colnames(pheno_ext))

  # --- 探针 -> gene symbol ---
  fdat <- fData(eset)
  sym_col <- colnames(fdat)[grepl("symbol|gene.?symbol", colnames(fdat), ignore.case = TRUE)][1]
  if (is.na(sym_col)) {
    cat("找不到 symbol 列, fData 列名:\n"); print(colnames(fdat))
    return(invisible(NULL))
  }
  syms <- as.character(fdat[[sym_col]])
  # 去除多探针注释 "A /// B" 只取第一个
  syms <- sapply(strsplit(syms, " /// "), `[`, 1)
  keep <- !is.na(syms) & syms != ""
  expr_ext <- expr_ext[keep, , drop = FALSE]
  syms <- syms[keep]

  # 按 symbol 合并 (取均值)
  idx <- split(seq_along(syms), syms)
  expr_sym <- t(vapply(idx, function(ii) {
    if (length(ii) == 1) as.numeric(expr_ext[ii, ]) else colMeans(expr_ext[ii, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(expr_ext))))
  rownames(expr_sym) <- names(idx)

  # --- 疗效标签 ---
  # Riaz 2017 的 pheno 里通常有 "response" 或 "clinical response" 相关列;
  # 列名可能随 GEO 版本变化, 先打印再映射
  resp_col <- colnames(pheno_ext)[grepl("response", colnames(pheno_ext), ignore.case = TRUE)]
  cat("候选疗效列:\n"); print(resp_col)
  if (length(resp_col) == 0) {
    cat("未找到疗效列, 请手动检查 pheno_ext\n")
    return(invisible(NULL))
  }
  resp_raw <- as.character(pheno_ext[[resp_col[1]]])
  resp_raw <- sub("^.*: ", "", resp_raw)  # 去掉 "field: value" 前缀
  print(table(resp_raw, useNA = "ifany"))

  # 注意: 该队列含 on-treatment 样本, 仅保留 baseline (pre-treatment)
  time_col <- colnames(pheno_ext)[grepl("time|biopsy|visit|pre|on", colnames(pheno_ext), ignore.case = TRUE)]
  cat("候选时间点列:\n"); print(time_col)
  if (length(time_col) > 0) {
    tp <- sub("^.*: ", "", as.character(pheno_ext[[time_col[1]]]))
    print(table(tp, useNA = "ifany"))
    pre_idx <- grepl("pre|baseline", tp, ignore.case = TRUE)
    cat("保留 baseline 样本:", sum(pre_idx), "/", length(tp), "\n")
  } else {
    pre_idx <- rep(TRUE, nrow(pheno_ext))
    cat("未识别时间点列, 使用全部样本 —— 请人工核对是否含 on-treatment!\n")
  }

  y_map <- c("CR" = 1, "PR" = 1, "CR/PR" = 1, "Responder" = 1, "R" = 1,
             "SD" = 0, "PD" = 0, "SD/PD" = 0, "NonResponder" = 0, "NR" = 0)
  y <- unname(y_map[resp_raw])

  resp_df <- data.frame(sample_id = rownames(pheno_ext)[pre_idx],
                        y = y[pre_idx], stringsAsFactors = FALSE)
  resp_df <- resp_df[!is.na(resp_df$y), , drop = FALSE]
  expr_pre <- expr_sym[, colnames(expr_sym) %in% resp_df$sample_id, drop = FALSE]

  validate_external(expr_pre, resp_df, cohort = "GSE91061_Riaz")
}

# 运行示例 (需要网络):
# run_gse91061()

cat("\n模块3 已加载. 使用方式:\n")
cat("  1) 通用: validate_external(expr_mat, resp_df, cohort='名字')\n")
cat("  2) 自动: run_gse91061()  (需要网络 + GEOquery)\n")
cat("注意: 自动下载路径的列名可能因 GEO 版本变化, 首次运行请检查打印的列名\n")
