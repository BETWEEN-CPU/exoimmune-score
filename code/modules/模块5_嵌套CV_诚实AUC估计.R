# ============================================================
# 模块5: 嵌套交叉验证 —— 诚实 AUC 估计
# ============================================================
# 目的:
#   原管线在同一批样本上先后做 Wilcoxon 筛选 + LASSO 训练,
#   再用 CV/bootstrap 评估 —— 特征选择的信息已泄漏进评估,
#   AUC 系统性偏高 (双重 dipping)。
#
#   本模块用嵌套 CV 给出无偏估计:
#     外层 5 折: 留出测试折
#     内层: 仅在训练折内完成 Wilcoxon 筛选 + cv.glmnet 调参
#     用训练折的均值/SD 缩放测试折, 预测并计 AUC
#   重复 20 次取分布, 即"流程级"的诚实性能。
#
# 依赖: Phase 1/2 输出
# 产出: nested_cv_auc_summary.tsv, nested_cv_predictions.rds
# ============================================================

root_dir <- Sys.getenv("EXOIMMUNE_ROOT", unset = "D:/善意无穷倒退/李新月/EV")

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

out_dir <- file.path(root_dir, "ExoImmune_nested_cv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1) 加载数据 (与 Phase 3 相同的输入)
# -------------------------
expr <- readRDS(file.path(root_dir, "IMvigor210_core_screen", "core_exosome_expr.rds"))
resp <- readRDS(file.path(root_dir, "IMvigor210_clean", "IMvigor210_response_df.rds"))

resp <- resp[!is.na(resp$response_clean), , drop = FALSE]
resp$response_clean <- factor(resp$response_clean, levels = c("NonResponder", "Responder"))

common <- intersect(colnames(expr), resp$sample_id)
expr <- expr[, common, drop = FALSE]
resp <- resp[match(common, resp$sample_id), , drop = FALSE]

y <- ifelse(resp$response_clean == "Responder", 1, 0)
names(y) <- resp$sample_id

cat("样本:", ncol(expr), " 基因:", nrow(expr), " Responder:", sum(y), "\n")

# -------------------------
# 2) 单次嵌套 CV 流程
# -------------------------
nested_cv_once <- function(expr, y, n_outer = 5, seed = 1, top_n = 50) {

  set.seed(seed)
  # 分层外层折
  folds <- rep(NA_integer_, length(y))
  idx1 <- which(y == 1); idx0 <- which(y == 0)
  folds[idx1] <- sample(rep(1:n_outer, length.out = length(idx1)))
  folds[idx0] <- sample(rep(1:n_outer, length.out = length(idx0)))

  pred_all <- rep(NA_real_, length(y))

  for (k in 1:n_outer) {
    tr <- which(folds != k)
    te <- which(folds == k)

    expr_tr <- expr[, tr, drop = FALSE]
    y_tr <- y[tr]

    # --- 内层步骤1: Wilcoxon 筛选 (只用训练折) ---
    pvals <- apply(expr_tr, 1, function(x) {
      tryCatch(wilcox.test(x[y_tr == 1], x[y_tr == 0])$p.value,
               error = function(e) NA_real_)
    })
    padj <- p.adjust(pvals, method = "BH")
    ord <- order(padj, pvals, na.last = TRUE)
    sig <- which(padj < 0.1)
    cand <- unique(c(sig, head(ord, top_n)))
    cand <- cand[!is.na(pvals[cand])]
    if (length(cand) < 10) cand <- head(ord, min(100, length(ord)))
    genes_k <- rownames(expr_tr)[cand]

    # --- 内层步骤2: 训练折 z-score (记录缩放参数) ---
    x_tr <- expr_tr[genes_k, , drop = FALSE]
    mu <- rowMeans(x_tr, na.rm = TRUE)
    sds <- apply(x_tr, 1, sd, na.rm = TRUE)
    sds[is.na(sds) | sds == 0] <- 1   # 防除零
    xz_tr <- sweep(sweep(x_tr, 1, mu, "-"), 1, sds, "/")
    xz_tr[is.na(xz_tr)] <- 0
    keep <- apply(xz_tr, 1, sd) > 0
    xz_tr <- xz_tr[keep, , drop = FALSE]
    mu <- mu[keep]; sds <- sds[keep]; genes_use <- rownames(xz_tr)

    # --- 内层步骤3: cv.glmnet 调参 (只用训练折) ---
    set.seed(seed + k)
    cvfit <- cv.glmnet(t(xz_tr), y_tr, family = "binomial", alpha = 1,
                       nfolds = min(10, min(table(y_tr))), type.measure = "deviance",
                       standardize = FALSE)
    co <- as.matrix(coef(cvfit, s = "lambda.1se"))
    nz <- rownames(co)[co[, 1] != 0]
    if (length(setdiff(nz, "(Intercept)")) < 2) {
      co <- as.matrix(coef(cvfit, s = "lambda.min"))
    }

    # --- 测试折: 用训练折缩放参数 ---
    x_te <- expr[genes_use, te, drop = FALSE]
    xz_te <- sweep(sweep(x_te, 1, mu, "-"), 1, sds, "/")
    xz_te[is.na(xz_te)] <- 0

    pred_all[te] <- as.numeric(predict(cvfit, newx = t(xz_te), s = "lambda.1se", type = "response"))
  }

  if (any(is.na(pred_all)) || length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, pred_all, levels = c(0, 1), direction = "<", quiet = TRUE)))
}

# -------------------------
# 3) 重复运行
# -------------------------
N_REPEAT <- 20
cat("开始嵌套 CV (5折 x", N_REPEAT, "次重复) ...\n")
cat("提示: 每次重复含 5 折 x (Wilcoxon + cv.glmnet), 可能耗时较长\n")

aucs <- numeric(N_REPEAT)
t0 <- Sys.time()
for (i in 1:N_REPEAT) {
  aucs[i] <- nested_cv_once(expr, y, n_outer = 5, seed = 20260726 + i)
  cat(sprintf("  repeat %2d/%d: AUC = %.3f  (elapsed %.1f min)\n",
              i, N_REPEAT, aucs[i], as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

aucs <- aucs[!is.na(aucs)]
cat("\n=== 嵌套 CV 诚实 AUC ===\n")
cat(sprintf("mean = %.3f, sd = %.3f, 2.5%%-97.5%% = %.3f-%.3f\n",
            mean(aucs), sd(aucs), quantile(aucs, 0.025), quantile(aucs, 0.975)))

# -------------------------
# 4) 与表观 AUC 对比
# -------------------------
score_df <- readRDS(file.path(root_dir, "IMvigor210_lasso", "exoimmune_score.rds"))
app <- merge(data.frame(sample_id = names(y), y = y), score_df, by = "sample_id")
apparent_auc <- as.numeric(pROC::auc(pROC::roc(app$y, app$score, levels = c(0, 1),
                                               direction = "<", quiet = TRUE)))
cat(sprintf("\n表观 AUC (原管线) = %.3f\n", apparent_auc))
cat(sprintf("乐观偏差估计 = %.3f\n", apparent_auc - mean(aucs)))

summary_tbl <- data.frame(
  metric = c("apparent_auc", "nested_cv_auc_mean", "nested_cv_auc_sd",
             "nested_cv_auc_low", "nested_cv_auc_high", "optimism",
             "n_repeats", "n_samples", "n_genes_input"),
  value = c(apparent_auc, mean(aucs), sd(aucs),
            quantile(aucs, 0.025), quantile(aucs, 0.975),
            apparent_auc - mean(aucs),
            length(aucs), ncol(expr), nrow(expr)),
  stringsAsFactors = FALSE
)
write.table(summary_tbl, file.path(out_dir, "nested_cv_auc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(repeat_id = seq_along(aucs), auc = aucs),
            file.path(out_dir, "nested_cv_auc_per_repeat.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(out_dir, "nested_cv_auc_distribution.pdf"), width = 6, height = 5)
hist(aucs, breaks = 10, col = "#2C7FB8", border = "white",
     xlab = "Nested CV AUC", main = "Honest performance estimate")
abline(v = apparent_auc, col = "#D95F02", lwd = 2, lty = 2)
legend("topright", legend = c("Nested CV", "Apparent"),
       col = c("#2C7FB8", "#D95F02"), lwd = c(8, 2), lty = c(1, 2), bty = "n")
dev.off()

cat("\n模块5 完成. 输出目录:", out_dir, "\n")
cat("写作建议: 主文报告 nested CV AUC 作为内部性能, 表观 AUC 可放补充材料;\n")
cat("二者之差即乐观偏差, 主动报告反而增强可信度。\n")
