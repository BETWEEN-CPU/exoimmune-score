# 数据说明 (Data)

原始数据体积约 7.5 GB，**不入库**。请按下方说明获取数据，并保持与代码预期一致的目录结构。

代码通过环境变量 `EXOIMMUNE_ROOT` 定位数据根目录；未设置时默认 `D:/善意无穷倒退/李新月/EV`。

```bash
# Linux / macOS
export EXOIMMUNE_ROOT=/path/to/data_root

# Windows PowerShell
$env:EXOIMMUNE_ROOT = "D:\path\to\data_root"
```

数据根目录下需存在以下文件/文件夹：

## 训练与内部验证队列

| 数据 | 说明 | 获取方式 |
|---|---|---|
| `IMvigor210CoreBiologies.Rdata` | IMvigor210 尿路上皮癌 anti-PD-L1 (atezolizumab) 队列（训练集） | 通过 R 包 `IMvigor210CoreBiologies`（Bioconductor）获取：`BiocManager::install("IMvigor210CoreBiologies")`，包内自带 `IMvigor210CoreBiologies.Rdata`。原始出处：Mariathasan et al., *Nature* 2018, 554:544-548 |
| `ExoCarta_protein_mRNA_details_6.txt` | ExoCarta 数据库蛋白/mRNA 条目（外泌体基因来源） | [exocarta.org](http://www.exocarta.org) 下载 v6 全库 |
| `VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt` | Vesiclepedia 数据库蛋白/mRNA 条目（外泌体基因来源） | [microvesicles.org](http://www.microvesicles.org) 下载 v5.1 全库 |

## 单细胞数据（来源分析与细胞通讯）

| 数据 | 说明 | 获取方式 |
|---|---|---|
| `GSE72056_melanoma_single_cell_revised_v2.txt` | 黑色素瘤单细胞 RNA-seq（Tirosh et al., *Science* 2016, 352:189-196） | GEO: [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056) |
| `GSE115978_counts.csv.gz` | 单细胞计数矩阵 | GEO: [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978) |
| `GSE126044_counts.txt.gz` | 单细胞计数矩阵 | GEO: [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) |

## 外部验证候选队列（模块3 可选）

| 数据 | 说明 | 获取方式 |
|---|---|---|
| `GSE91061` | 黑色素瘤 anti-PD-1 队列（Riaz et al., *Cell* 2017） | GEO: [GSE91061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061)；模块3 提供自动下载脚本 |
| `GSE78220_series_matrix.txt` | 黑色素瘤 anti-PD-1 队列（Hugo et al., *Cell* 2016） | GEO: [GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220) |
| `GSE135222` | 黑色素瘤 anti-PD-1/anti-CTLA4 队列 | GEO: [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222) |
| Gide 2019 / IMmotion150/151 | 其他可选外部队列 | 需申请或从文献补充数据整理 |

## 机制分析（TCGA）

| 数据 | 说明 | 获取方式 |
|---|---|---|
| `EV_immunotherapy_project/TCGA_processed/` | TCGA-BLCA / TCGA-SKCM 处理后的表达与临床数据 | 原始数据来自 [GDC Data Portal](https://portal.gdc.cancer.gov)（GENCODE v36 注释的 STAR 定量），下载后按脚本步骤处理 |
| `msigdb_v2026.1.Hs_GMTs/` | MSigDB 人类基因集（免疫/非免疫模块） | [gsea-msigdb.org](https://www.gsea-msigdb.org/gsea/msigdb) 下载 v2026.1 |

## 建议的文件布局（与脚本默认一致）

```
EXOIMMUNE_ROOT/
├── IMvigor210CoreBiologies.Rdata         # 或 IMvigor210/5ff67-main/IMvigor210CoreBiologies/...
├── ExoCarta_protein_mRNA_details_6.txt
├── VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt
├── GSE72056_melanoma_single_cell_revised_v2.txt/
├── GSE115978_counts.csv.gz
├── GSE126044_counts.txt.gz
├── EV_immunotherapy_project/TCGA_processed/
└── msigdb_v2026.1.Hs_GMTs/
```

注：脚本运行过程中会**自动创建**各分析输出目录（`IMvigor210_clean/`、`IMvigor210_lasso/`、`ExoImmune_*/`、`singlecell_*/`、`TCGA_*/` 等），无需手动准备。
