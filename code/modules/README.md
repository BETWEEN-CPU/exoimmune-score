# ExoImmune 课题优化模块包 —— 使用说明

## 这套代码解决什么问题

针对原管线（`总.txt` / `总(1).txt`）的三个主要短板：

| 短板 | 对应模块 |
|---|---|
| 没有真正的外部验证 | 模块3 |
| 未与 PD-L1 / TMB / TIG-GEP 头对头比较 | 模块1 |
| 特征筛选与评估同源（双重 dipping），AUC 偏乐观 | 模块5 |
| 校准缺失、中位数切分丢信息 | 模块2 |
| 细胞通讯分析没有统计检验 | 模块4 |
| 外部验证缺训练集缩放参数 | 模块0 |

## 执行顺序

```
模块0 → 模块1 → 模块2 → 模块5 → 模块3 → 模块4
```

- **模块0 必须先跑**：它导出评分资产（系数 + 训练集均值/SD + 通用评分函数），模块1/3 依赖它。这也是补原管线的一个真实漏洞——Phase 3 没保存训练集缩放参数，此前任何"外推"其实在数学上都不可复现。
- R 模块直接在 RStudio 里 `source()` 即可；模块4 是 Python（`python 模块4_通讯置换检验.py`）。
- 所有模块的 `root_dir` 默认沿用 `D:/善意无穷倒退/李新月/EV`；如需更换数据根目录，设置环境变量 `EXOIMMUNE_ROOT` 即可（`Sys.setenv(EXOIMMUNE_ROOT="...")` 或 PowerShell `$env:EXOIMMUNE_ROOT="..."`），无需改代码。

## 各模块注意事项

**模块1（头对头比较）**
- PD-L1 IC / TMB 的列名自动匹配，但 IMvigor210CoreBiologies 的版本不同列名可能有差异。跑完先看打印的"检测到 IC 列 / TMB 列"是否正确，不对就手动改 `pick_col` 的 patterns。
- 三种结局都要诚实写进论文：若评分优于部分标志物或在多变量模型中独立显著 → 主张增量价值；若不显著 → 转向"外泌体来源的简便替代标志物"叙事，反而更有故事性。

**模块2（校准 + 连续 Cox）**
- HL 检验功效低，以校准曲线目视判断为主。
- 装了 `rms` 包会自动多做一步限制性立方样条非线性检查，建议装。

**模块5（嵌套 CV）**
- 计算量最大（20 次重复 × 5 折 × cv.glmnet），预计跑几十分钟到几小时，建议挂后台。
- 嵌套 CV AUC 大概率低于表观 AUC——这是**正常的**，主文报嵌套值、表观值放补充材料，主动报告乐观偏差是加分项不是减分项。

**模块3（外部验证）**
- 核心原则：**系数和缩放参数全部来自训练集**，外部数据上一个参数都不能重估。
- `run_gse91061()` 是自动下载示例，但 GEO 的 pheno 列名随版本变化，**第一次跑务必人工核对**打印出的列名、疗效映射和 baseline 筛选（Riaz 队列含 on-treatment 样本，必须剔除）。
- 如果有条件拿到 Gide 2019 (PRJEB23709) 或 IMmotion 系列数据，验证价值更高；哪怕只有一个外部队列也够用。

**模块4（通讯置换检验）**
- 默认 1000 次置换，GSE72056 约 4600 细胞，预计几分钟到十几分钟。
- 结果用法：只有 FDR<0.05 的边进入主文机制叙事；原来证据分高但置换不显著的边降级为描述性观察。

## 论文写法建议（与既有审慎风格一致）

- TCGA 部分统一用 prognostic association，IMvigor210/外部队列部分才用 predictive；两者不混用。
- 性能表述："The score showed an apparent AUC of X in the training cohort and a nested cross-validation AUC of Y, suggesting moderate optimism in the apparent estimate." 这类写法审稿人很吃。
- 外部验证若不显著，可用 "warrants validation in larger cohorts" 收尾；若显著，用 "provides initial external support for..." 而非 "confirms"。

## 已知遗留问题（下一步可做）

1. GSE115978 的 marker 注释较粗糙，可考虑换 scType/SingleR 重注释做稳健性。
2. 组织错配（膀胱癌评分用黑色素瘤单细胞解释）——TISCH 数据库查 BLCA 单细胞数据可补。
3. 若跑完模块1发现 TMB 数据可用，可加 score+TMB 联合模型作为二级卖点。
