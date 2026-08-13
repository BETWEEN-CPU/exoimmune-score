# ExoImmune — Exosome-Associated Gene Signature for Immunotherapy Response

> A transcriptome-derived signature of exosome-associated genes (ExoImmune Score) for predicting immune checkpoint inhibitor (ICI) response, built on IMvigor210 (atezolizumab-treated urothelial carcinoma) and validated across single-cell, TCGA, and external ICI cohorts.

## 项目简介 (Overview)

本研究构建并验证了一个基于**外泌体相关基因**的转录组评分（ExoImmune Score，9 个 LASSO 基因），用于预测免疫检查点抑制剂（ICI）治疗应答，并系统评估其独立于 PD-L1 / TMB / TIG-GEP 的增量价值。

**分析流程全景：**

- **训练与内部验证**（R）：IMvigor210 队列 → ExoCarta/Vesiclepedia 外泌体基因交集 → Wilcoxon 筛选 → LASSO（10 折 CV, lambda.1se）→ 构建评分 → KM / Cox 生存分析
- **机制分析**（R）：评分在 TCGA-BLCA / TCGA-SKCM 中的预后关联与免疫模块相关性
- **单细胞来源与通讯**（Python）：GSE72056 / GSE115978 / GSE126044 单细胞数据 → 细胞注释 → 配体-受体通讯分析
- **统计严谨性优化模块**（code/modules/）：评分资产导出、头对头比较、校准 + 连续 Cox、外部 ICI 队列验证、通讯置换检验、嵌套交叉验证

## 仓库结构 (Repository Structure)

```
ExoImmune/
├── README.md                     # 本文件
├── data/
│   └── README.md                 # 数据获取说明（数据本身不入库，约 7.5 GB）
└── code/
    ├── python/
    │   └── ExoImmune_pipeline.py # Python 主管线：单细胞来源 + 细胞通讯 + TCGA 生物标志物 + 汇总
    ├── r/
    │   └── ExoImmune_pipeline.R  # R 主管线：IMvigor210 评分构建 + 生存分析 + TCGA 验证
    └── modules/                  # 统计严谨性优化模块包（0~5）
        ├── README.md             # 模块使用说明（执行顺序、注意事项、论文写法建议）
        ├── 模块0_评分资产导出.R
        ├── 模块1_生物标志物头对头比较.R
        ├── 模块2_校准曲线与连续Cox.R
        ├── 模块3_外部ICI队列验证.R
        ├── 模块4_通讯置换检验.py
        └── 模块5_嵌套CV_诚实AUC估计.R
```

## 快速开始 (Quick Start)

### 1. 准备数据

数据不入库。请参照 [`data/README.md`](data/README.md) 获取所需数据并设置数据根目录：

```bash
export EXOIMMUNE_ROOT=/path/to/data_root   # 未设置时使用默认路径
```

### 2. 运行主管线

```r
# RStudio 中打开 code/r/ExoImmune_pipeline.R 并 source()
```

```bash
# Python（单细胞 + TCGA + IMvigor210 内部验证）
python code/python/ExoImmune_pipeline.py
```

### 3. 运行优化模块（统计严谨性补强）

先阅读 [`code/modules/README.md`](code/modules/README.md)（含执行顺序与注意事项）：

```
模块0 → 模块1 → 模块2 → 模块5 → 模块3 → 模块4
```

## 环境依赖 (Dependencies)

**Python 3.8+**：`numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `scikit-learn`, `lifelines`（可选，缺失时自动跳过 Cox/对数秩检验）

**R ≥ 4.0**：`survival`, `glmnet`, `pROC`, `ggplot2`；优化模块额外需要 `GEOquery`（模块3）、`rms`（模块2，可选）

## 论文报告要点

- TCGA 部分统一用 **prognostic association**，IMvigor210/外部队列部分用 **predictive**，两者不混用。
- 主文报嵌套 CV AUC，表观 AUC 放补充材料，主动报告乐观偏差。
- 外部验证结论措辞建议用 "provides initial external support for..." 而非 "confirms"。

详细论文写法建议见 [`code/modules/README.md`](code/modules/README.md)。

## 许可证 (License)

本项目代码仅供学术研究使用，如有引用需求请联系作者。
