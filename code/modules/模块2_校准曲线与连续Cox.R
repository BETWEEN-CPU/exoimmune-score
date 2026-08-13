# ============================================================
# 模块2: 校准评估 + 连续变量 Cox 分析
# ============================================================
# 目的:
#   补齐原管线的两个统计缺口:
#   1) 校准: 预测概率是否与实际响应率一致 (校准曲线 + Brier + HL检验)
#      —— 原代码 import 了 brier_score_loss 却未使用
#   2) 生存: 主分析改用连续变量 Cox (HR per 1 SD), 中位数 KM 仅作展示;
#      输出 C-index 供与外部队列比较
#
# 依赖: Phase 3/4 输出
# 产出: calibration_plot.pdf, calibration_stats.tsv, cox_continuous.tsv
# ============================================================

root_dir <- Sys.getenv("EXOIMMUNE_ROOT", unset = "D:/善意无穷倒退/李新月/EV")

suppressPackageStartupMessages({
  library(survival)
})

out_dir <- file.path(root_dir, "ExoImmune_calibration_cox")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1) 加载与合并
# -------------------------
score_df <- readRDS(file.path(root_dir, "IMvigor210_lasso", "exoimmune_score.rds"))
pheno    <- readRDS(file.path(root_dir, "IMvigor210_clean", "IMvigor210_pheno_aligned.rds"))

pheno$sample_id    <- as.character(pheno$sample_id)
score_df$sample_id <- as.character(score_df$sample_id)

df <- merge(score_df, pheno, by = "sample_id", all.x = TRUE, suffixes = c("", ".pheno"))
df$y <- ifelse(df$response_clean == "Responder", 1, 0)

# ============================================================
# PART A: 校准分析 (响应预测)
# ============================================================

cal <- df[!is.na(df$y) & !is.na(df$score), , drop = FALSE]
cat("校准分析样本:", nrow(cal), "\n")

# 用评分拟合 logistic (表观拟合, 校准曲线的标准做法)
fit_cal <- glm(y ~ score, data = cal, family = binomial)
cal$pred <- predict(fit_cal, type = "response")

# --- 校准曲线: 分位数分箱 (每组约等量样本) ---
n_bins <- 10
cal$bin <- cut(cal$pred, breaks = quantile(cal$pred, probs = seq(0, 1, length.out = n_bins + 1)),
               include.lowest = TRUE, labels = FALSE)

cal_curve <- aggregate(cbind(pred, y) ~ bin, data = cal, FUN = mean)
cal_n     <- aggregate(y ~ bin, data = cal, FUN = length)
colnames(cal_n)[2] <- "n"
cal_curve <- merge(cal_curve, cal_n, by = "bin")

cat("\n=== 校准分箱表 ===\n")
print(cal_curve)

# --- Brier score ---
brier <- mean((cal$pred - cal$y)^2)
# 参考: 无信息模型 (患病率) 的 Brier
prev <- mean(cal$y)
brier_null <- mean((prev - cal$y)^2)
brier_skill <- 1 - brier / brier_null

cat(sprintf("\nBrier = %.4f  (null=%.4f, skill score=%.3f)\n", brier, brier_null, brier_skill))

# --- Hosmer-Lemeshow 检验 (手动实现, 避免额外依赖) ---
hl_group <- split(seq_len(nrow(cal)), cal$bin)
hl_stat <- 0
for (idx in hl_group) {
  obs1 <- sum(cal$y[idx]); exp1 <- sum(cal$pred[idx])
  n_g  <- length(idx); p_g <- mean(cal$pred[idx])
  if (p_g > 0 & p_g < 1) {
    hl_stat <- hl_stat + (obs1 - exp1)^2 / (n_g * p_g * (1 - p_g))
  }
}
hl_df <- length(hl_group) - 2
hl_p <- 1 - pchisq(hl_stat, df = hl_df)
cat(sprintf("Hosmer-Lemeshow: X2=%.2f, df=%d, p=%.4f\n", hl_stat, hl_df, hl_p))
cat("  (p>0.05 提示无显著失配; 注意 HL 检验功效有限, 以校准曲线目视为主)\n")

# --- 校准图 ---
pdf(file.path(out_dir, "calibration_plot.pdf"), width = 6, height = 6)
plot(cal_curve$pred, cal_curve$y, pch = 16, col = "#2C7FB8",
     xlim = 0:1, ylim = 0:1,
     xlab = "Predicted response probability",
     ylab = "Observed response proportion",
     main = sprintf("Calibration (Brier=%.3f)", brier))
abline(a = 0, b = 1, lty = 2, col = "grey60")
lines(cal_curve$pred, cal_curve$y, col = "#2C7FB8")
# 直方图背景: 预测概率分布
par(new = TRUE)
hist(cal$pred, breaks = 20, axes = FALSE, xlab = "", ylab = "", main = "",
     col = rgb(0.8, 0.8, 0.8, 0.4), border = NA, freq = FALSE)
dev.off()

cal_stats <- data.frame(
  metric = c("n", "brier", "brier_null", "brier_skill", "HL_chisq", "HL_df", "HL_p"),
  value = c(nrow(cal), brier, brier_null, brier_skill, hl_stat, hl_df, hl_p),
  stringsAsFactors = FALSE
)
write.table(cal_stats, file.path(out_dir, "calibration_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cal_curve, file.path(out_dir, "calibration_bins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# PART B: 连续变量 Cox 分析 (OS)
# ============================================================

df$os      <- as.numeric(df$os)
df$censOS  <- as.numeric(df$censOS)
df$event_os <- ifelse(df$censOS == 0, 1, 0)

sv <- df[!is.na(df$os) & !is.na(df$event_os) & !is.na(df$score), , drop = FALSE]
cat("\n生存分析样本:", nrow(sv), "  事件数:", sum(sv$event_os), "\n")

# --- 连续评分, HR per 1 SD ---
sv$score_z <- as.numeric(scale(sv$score))
cox_cont <- coxph(Surv(os, event_os) ~ score_z, data = sv)
s <- summary(cox_cont)

cox_cont_tbl <- data.frame(
  variable = "ExoImmune_Score (per 1 SD)",
  n = nrow(sv), events = sum(sv$event_os),
  HR = unname(s$coefficients[1, "exp(coef)"]),
  ci_low = unname(s$conf.int[1, "lower .95"]),
  ci_high = unname(s$conf.int[1, "upper .95"]),
  p = unname(s$coefficients[1, "Pr(>|z|)"]),
  c_index = unname(s$concordance[1]),
  stringsAsFactors = FALSE
)
cat("\n=== 连续变量 Cox (HR per 1 SD) ===\n")
print(cox_cont_tbl)

# --- 限制性立方样条检查线性假设 (若 rms 可用) ---
has_rms <- requireNamespace("rms", quietly = TRUE)
if (has_rms) {
  suppressPackageStartupMessages(library(rms))
  dd <- datadist(sv); options(datadist = "dd")
  fit_rcs <- cph(Surv(os, event_os) ~ rcs(score_z, 4), data = sv, x = TRUE, y = TRUE)
  pdf(file.path(out_dir, "cox_rcs_nonlinearity.pdf"), width = 6, height = 5)
  plot(Predict(fit_rcs), main = "HR vs score (restricted cubic spline)")
  dev.off()
  anova_rcs <- anova(fit_rcs)
  print(anova_rcs)
  write.table(as.data.frame(anova_rcs), file.path(out_dir, "cox_rcs_anova.tsv"),
              sep = "\t", quote = FALSE)
  cat("RCS 非线性检验已输出; 若 nonlinear p>0.05, 连续线性 Cox 假设可接受\n")
} else {
  cat("提示: 安装 rms 包可做非线性检查 (install.packages('rms'))\n")
}

# --- 中位数 KM 保留为展示图 (沿用原分析, 注明仅展示用) ---
med <- median(sv$score, na.rm = TRUE)
sv$score_group <- factor(ifelse(sv$score >= med, "High", "Low"), levels = c("Low", "High"))
fit_km <- survfit(Surv(os, event_os) ~ score_group, data = sv)
lr <- survdiff(Surv(os, event_os) ~ score_group, data = sv)
p_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

pdf(file.path(out_dir, "km_median_display_only.pdf"), width = 7, height = 6)
plot(fit_km, col = c("#D95F02", "#1B9E77"), lwd = 2,
     xlab = "Overall survival time", ylab = "Survival probability",
     main = paste0("KM by median split (display only)\nLog-rank p = ", signif(p_lr, 3)))
legend("bottomleft", legend = levels(sv$score_group),
       col = c("#D95F02", "#1B9E77"), lwd = 2, bty = "n")
dev.off()

write.table(cox_cont_tbl, file.path(out_dir, "cox_continuous.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n模块2 完成. 输出目录:", out_dir, "\n")
cat("写作建议: 主文报告 HR per 1 SD + C-index; 中位数 KM 放附图并注明展示用途\n")
