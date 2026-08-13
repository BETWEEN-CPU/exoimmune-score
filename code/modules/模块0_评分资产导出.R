# ============================================================
# 模块0: ExoImmune Score 评分资产导出
# ============================================================
# 目的:
#   原 Phase 3 只保存了 LASSO 系数(lasso_coef_lambda_1se.tsv),
#   但没有保存训练集的基因均值/SD。外部验证时必须用【训练集】的
#   缩放参数对新数据做 z-score,否则评分不可比。
#   本模块补齐这一缺口,导出评分所需的全部资产。
#
# 依赖: Phase 2/3 的输出 (core_exosome_expr.rds, lasso_coef_lambda_1se.tsv)
# 产出: exoimmune_scoring_assets.rds / .tsv
# ============================================================

root_dir <- Sys.getenv("EXOIMMUNE_ROOT", unset = "D:/善意无穷倒退/李新月/EV")

in_dir_p2 <- file.path(root_dir, "IMvigor210_core_screen")
in_dir_p3 <- file.path(root_dir, "IMvigor210_lasso")
out_dir   <- file.path(root_dir, "ExoImmune_scoring_assets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1) 加载系数
# -------------------------
coef_file <- file.path(in_dir_p3, "lasso_coef_lambda_1se.tsv")
if (!file.exists(coef_file)) {
  coef_file <- file.path(in_dir_p3, "lasso_coef_lambda_min.tsv")
  cat("lambda.1se 系数不存在, 改用 lambda.min\n")
}
coef_df <- read.delim(coef_file, sep = "\t", header = TRUE, check.names = FALSE)
coef_df$gene <- as.character(coef_df$gene)

intercept <- 0
if ("(Intercept)" %in% coef_df$gene) {
  intercept <- coef_df$coef[coef_df$gene == "(Intercept)"]
  if (is.na(intercept)) intercept <- 0
}
coef_genes <- setdiff(coef_df$gene, "(Intercept)")
coef_vec   <- coef_df$coef[match(coef_genes, coef_df$gene)]

cat("评分基因数:", length(coef_genes), "\n")
cat("截距:", intercept, "\n")

# -------------------------
# 2) 从训练集表达矩阵计算均值/SD
#    注意: Phase 3 中 z-score 是在【响应分析样本子集】上按行(scale)做的,
#    这里保持一致: 使用 core_exosome_expr (已对齐响应样本) 的逐基因均值/SD
# -------------------------
core_expr <- readRDS(file.path(in_dir_p2, "core_exosome_expr.rds"))

present <- coef_genes[coef_genes %in% rownames(core_expr)]
missing <- setdiff(coef_genes, rownames(core_expr))
if (length(missing) > 0) {
  cat("警告: 以下评分基因不在 core_exosome_expr 中:", paste(missing, collapse = ", "), "\n")
}

scale_tbl <- data.frame(
  gene = present,
  mean = rowMeans(core_expr[present, , drop = FALSE], na.rm = TRUE),
  sd   = apply(core_expr[present, , drop = FALSE], 1, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
scale_tbl$sd[is.na(scale_tbl$sd) | scale_tbl$sd == 0] <- NA

if (any(is.na(scale_tbl$sd))) {
  cat("警告: 以下基因训练集 SD 为 0 或 NA, 外部验证时这些基因将贡献 0:\n")
  print(scale_tbl$gene[is.na(scale_tbl$sd)])
}

# -------------------------
# 3) 保存资产
# -------------------------
assets <- list(
  coef_df    = data.frame(gene = coef_genes, coef = coef_vec, stringsAsFactors = FALSE),
  intercept  = intercept,
  scale_tbl  = scale_tbl,
  trained_on = "IMvigor210 (core exosome genes, lambda.1se)",
  created    = as.character(Sys.Date())
)

saveRDS(assets, file.path(out_dir, "exoimmune_scoring_assets.rds"))

write.table(assets$coef_df, file.path(out_dir, "exoimmune_coef.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(scale_tbl, file.path(out_dir, "exoimmune_training_scale.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 4) 通用评分函数 (外部验证/模块1-3 共用)
# -------------------------
score_new_data <- function(expr_mat, assets) {
  # expr_mat: 基因 x 样本 的数值矩阵, rownames 为基因 symbol
  scale_tbl <- assets$scale_tbl
  overlap <- intersect(scale_tbl$gene, rownames(expr_mat))
  cat("评分基因覆盖: ", length(overlap), "/", nrow(scale_tbl), "\n")
  if (length(overlap) == 0) stop("表达矩阵与评分基因无重叠")

  expr_use <- expr_mat[overlap, , drop = FALSE]
  st <- scale_tbl[match(overlap, scale_tbl$gene), , drop = FALSE]

  z <- sweep(expr_use, 1, st$mean, "-")
  z <- sweep(z, 1, st$sd, "/")
  z[!is.finite(z)] <- 0   # 缺失基因/零方差贡献 0

  cf <- assets$coef_df
  cf <- cf[match(overlap, cf$gene), , drop = FALSE]
  if (any(is.na(cf$coef))) stop("部分重叠基因缺少系数")

  as.numeric(assets$intercept + t(z) %*% cf$coef)
}

saveRDS(score_new_data, file.path(out_dir, "score_new_data.rds"))

cat("\n模块0 完成. 资产保存于:", out_dir, "\n")
cat("下一步: 运行 模块1_生物标志物头对头比较.R\n")
