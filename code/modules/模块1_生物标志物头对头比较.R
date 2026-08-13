# ============================================================
# 模块1: ExoImmune Score vs 已知生物标志物 头对头比较
# ============================================================
# 目的:
#   回答审稿人必问的问题: 9基因评分相比 PD-L1、TMB、T cell-inflamed
#   GEP (Ayer 18基因) 是否有增量价值?
#
# 内容:
#   1) 从 phenoData 提取 PD-L1 IC Level / TMB
#   2) 从表达矩阵计算 Ayer TIG-GEP 评分
#   3) 单标志物 AUC 两两 DeLong 检验
#   4) 多变量 logistic: 评分是否在校正 PD-L1/TMB 后仍独立显著
#   5) 联合模型 AUC
#
# 依赖: Phase 1/3 输出 + 模块0 资产
# 产出: biomarker_headtohead_summary.tsv, AUC对比图, 多变量模型表
# ============================================================

root_dir <- Sys.getenv("EXOIMMUNE_ROOT", unset = "D:/善意无穷倒退/李新月/EV")

suppressPackageStartupMessages({
  library(pROC)
})

out_dir <- file.path(root_dir, "ExoImmune_biomarker_headtohead")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 0) 工具函数
# -------------------------
pick_col <- function(df, patterns) {
  cn <- colnames(df)
  hit <- cn[grepl(paste(patterns, collapse = "|"), cn, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

to_num <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "NE")] <- NA
  suppressWarnings(as.numeric(x))
}

# -------------------------
# 1) 加载数据
# -------------------------
score_df <- readRDS(file.path(root_dir, "IMvigor210_lasso", "exoimmune_score.rds"))
pheno    <- readRDS(file.path(root_dir, "IMvigor210_clean", "IMvigor210_pheno_aligned.rds"))
expr     <- readRDS(file.path(root_dir, "IMvigor210_clean", "IMvigor210_expr_symbol.rds"))

pheno$sample_id   <- as.character(pheno$sample_id)
score_df$sample_id <- as.character(score_df$sample_id)

df <- merge(score_df, pheno, by = "sample_id", all.x = TRUE, suffixes = c("", ".pheno"))

# 疗效标签
df$y <- ifelse(df$response_clean == "Responder", 1, 0)
df <- df[!is.na(df$y), , drop = FALSE]
cat("有效样本:", nrow(df), "  Responder:", sum(df$y), "\n")

# -------------------------
# 2) 提取 PD-L1 与 TMB
# -------------------------
cat("\nphenoData 列名 (供核对):\n")
print(colnames(pheno))

# IMvigor210CoreBiologies 中常见列:
#   "IC Level" (IC0/IC1/IC2+), "TC Level", "Immune Phenotype" (inflamed/excluded/desert)
#   "FMOne mutation burden per MB", "Neoantigen burden per MB"
ic_col  <- pick_col(pheno, c("^IC Level$", "IC.Level", "PD-L1.*IC", "ICLevel"))
tmb_col <- pick_col(pheno, c("FMOne mutation burden", "mutation burden per MB", "TMB"))
neo_col <- pick_col(pheno, c("Neoantigen burden"))

cat("检测到 IC 列:", ic_col, " TMB 列:", tmb_col, " Neo 列:", neo_col, "\n")

if (!is.na(ic_col)) {
  ic_raw <- as.character(pheno[[ic_col]])
  # 序数编码: IC0=0, IC1=1, IC2+=2; 兼容纯数字
  ic_ord <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ic_raw)))
  df$pdl1_ic <- ic_ord[match(df$sample_id, pheno$sample_id)]
} else {
  df$pdl1_ic <- NA
  cat("警告: 未找到 PD-L1 IC 列, 该标志物跳过\n")
}

if (!is.na(tmb_col)) {
  df$tmb <- to_num(pheno[[tmb_col]])[match(df$sample_id, pheno$sample_id)]
} else if (!is.na(neo_col)) {
  df$tmb <- to_num(pheno[[neo_col]])[match(df$sample_id, pheno$sample_id)]
  cat("用 neoantigen burden 替代 TMB\n")
} else {
  df$tmb <- NA
  cat("警告: 未找到 TMB 列, 该标志物跳过\n")
}

# -------------------------
# 3) Ayer T cell-inflamed GEP (18基因) z-score 均值
# -------------------------
tig_genes <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1",
               "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1",
               "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")

tig_present <- intersect(tig_genes, rownames(expr))
cat("TIG-GEP 基因覆盖:", length(tig_present), "/", length(tig_genes), "\n")

tig_mat <- expr[tig_present, , drop = FALSE]
tig_z <- t(scale(t(as.matrix(tig_mat))))
tig_z[!is.finite(tig_z)] <- NA
tig_score <- colMeans(tig_z, na.rm = TRUE)

df$tig_gep <- as.numeric(tig_score[match(df$sample_id, colnames(expr))])

# -------------------------
# 4) 单标志物 AUC + DeLong 检验
# -------------------------
biomarkers <- list(
  ExoImmune_Score = df$score,
  PDL1_IC         = df$pdl1_ic,
  TMB             = df$tmb,
  TIG_GEP         = df$tig_gep
)

auc_rows <- list()
roc_objs <- list()

for (nm in names(biomarkers)) {
  x <- biomarkers[[nm]]
  ok <- !is.na(x) & !is.na(df$y)
  if (sum(ok) < 30 || length(unique(df$y[ok])) < 2) next
  r <- pROC::roc(response = df$y[ok], predictor = x[ok],
                 levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_objs[[nm]] <- r
  ci <- pROC::ci.auc(r)
  auc_rows[[nm]] <- data.frame(
    biomarker = nm, n = sum(ok),
    auc = as.numeric(pROC::auc(r)),
    ci_low = ci[1], ci_high = ci[3],
    stringsAsFactors = FALSE
  )
}

auc_tbl <- do.call(rbind, auc_rows)
cat("\n=== 单标志物 AUC ===\n")
print(auc_tbl)

# DeLong: 每个标志物 vs ExoImmune Score
delong_rows <- list()
for (nm in setdiff(names(roc_objs), "ExoImmune_Score")) {
  # DeLong 要求相同样本集: 取两标志物都非缺失的样本重算
  ok <- !is.na(df$score) & !is.na(biomarkers[[nm]]) & !is.na(df$y)
  r1 <- pROC::roc(df$y[ok], df$score[ok], levels = c(0, 1), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(df$y[ok], biomarkers[[nm]][ok], levels = c(0, 1), direction = "<", quiet = TRUE)
  tst <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
  delong_rows[[nm]] <- data.frame(
    comparison = paste0("ExoImmune_Score vs ", nm),
    n_paired = sum(ok),
    auc_score = as.numeric(pROC::auc(r1)),
    auc_other = as.numeric(pROC::auc(r2)),
    delong_z = unname(tst$statistic),
    delong_p = tst$p.value,
    stringsAsFactors = FALSE
  )
}
delong_tbl <- do.call(rbind, delong_rows)
cat("\n=== DeLong 检验 (配对) ===\n")
print(delong_tbl)

# -------------------------
# 5) 多变量 logistic: 独立增量价值
# -------------------------
mv <- df[!is.na(df$score) & !is.na(df$pdl1_ic) & !is.na(df$tmb), , drop = FALSE]
if (nrow(mv) >= 50) {
  mv$score_z <- as.numeric(scale(mv$score))
  mv$tmb_z   <- as.numeric(scale(mv$tmb))

  fit_full <- glm(y ~ score_z + pdl1_ic + tmb_z, data = mv, family = binomial)
  fit_base <- glm(y ~ pdl1_ic + tmb_z, data = mv, family = binomial)

  cat("\n=== 多变量 logistic (y ~ score + PD-L1 IC + TMB) ===\n")
  print(summary(fit_full)$coefficients)

  lrt <- anova(fit_base, fit_full, test = "Chisq")
  cat("\n=== 似然比检验: 加入评分是否改善模型 ===\n")
  print(lrt)

  mv_tbl <- data.frame(
    term = rownames(summary(fit_full)$coefficients),
    OR = exp(coef(fit_full)),
    ci_low = exp(confint.default(fit_full)[, 1]),
    ci_high = exp(confint.default(fit_full)[, 2]),
    p = summary(fit_full)$coefficients[, 4],
    stringsAsFactors = FALSE
  )
  lrt_p <- lrt$`Pr(>Chi)`[2]

  # 联合模型 AUC
  mv$pred_full <- predict(fit_full, type = "response")
  mv$pred_base <- predict(fit_base, type = "response")
  r_full <- pROC::roc(mv$y, mv$pred_full, levels = c(0, 1), direction = "<", quiet = TRUE)
  r_base <- pROC::roc(mv$y, mv$pred_base, levels = c(0, 1), direction = "<", quiet = TRUE)
  tst_comb <- pROC::roc.test(r_full, r_base, method = "delong", paired = TRUE)

  comb_tbl <- data.frame(
    comparison = "Score+PDL1+TMB vs PDL1+TMB",
    auc_full = as.numeric(pROC::auc(r_full)),
    auc_base = as.numeric(pROC::auc(r_base)),
    delong_p = tst_comb$p.value,
    lrt_p = lrt_p,
    stringsAsFactors = FALSE
  )
  cat("\n=== 联合模型 AUC ===\n")
  print(comb_tbl)

  write.table(mv_tbl, file.path(out_dir, "multivariable_logistic.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(comb_tbl, file.path(out_dir, "combined_model_auc.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("\n多变量模型样本不足 (n=", nrow(mv), "), 跳过\n")
}

# -------------------------
# 6) 输出与图
# -------------------------
write.table(auc_tbl, file.path(out_dir, "single_biomarker_auc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(delong_tbl, file.path(out_dir, "delong_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(out_dir, "roc_comparison.pdf"), width = 6.5, height = 6)
cols <- c("#2C7FB8", "#D95F02", "#1B9E77", "#7570B3")
plot(roc_objs[[1]], col = cols[1], lwd = 2,
     main = "ExoImmune Score vs established biomarkers")
if (length(roc_objs) > 1) {
  for (i in 2:length(roc_objs)) plot(roc_objs[[i]], col = cols[i], lwd = 2, add = TRUE)
}
abline(a = 0, b = 1, lty = 2, col = "grey60")
legend("bottomright",
       legend = paste0(names(roc_objs), " (AUC=",
                       sapply(roc_objs, function(r) round(as.numeric(pROC::auc(r)), 3)), ")"),
       col = cols[seq_along(roc_objs)], lwd = 2, bty = "n", cex = 0.8)
dev.off()

cat("\n模块1 完成. 输出目录:", out_dir, "\n")
cat("解读提示: 若评分在 DeLong 中优于至少部分标志物, 或在多变量模型中独立显著,\n")
cat("即可主张增量价值; 若不显著, 可转向主张'外泌体来源的简便替代标志物'.\n")
