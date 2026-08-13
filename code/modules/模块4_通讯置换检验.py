# -*- coding: utf-8 -*-
"""
模块4: 单细胞配体-受体通讯 置换检验
====================================
目的:
  原管线的 evidence_score (检出率乘积 x log表达量乘积) 没有统计检验。
  本模块对每条 LR 边做细胞标签置换检验 (label permutation):
    - 打乱细胞类型标签 N 次, 重算 evidence_score 得到零分布
    - 经验 p 值 = 置换中 score >= 观测值的比例
    - BH-FDR 校正
  输出显著边列表, 可直接替换/补强原 Part 4-5 的叙事。

依赖: 原 总.txt 中的数据文件 (GSE72056)
产出: singlecell_communication_analysis/lr_permutation_results.tsv

运行: python 模块4_通讯置换检验.py
"""

from pathlib import Path
import gzip
import math
import os
import numpy as np
import pandas as pd

# ============================================================
# 配置
# ============================================================
root = Path(os.environ.get("EXOIMMUNE_ROOT", r"D:/善意无穷倒退/李新月/EV"))
gse72056 = root / "GSE72056_melanoma_single_cell_revised_v2.txt" / "GSE72056_melanoma_single_cell_revised_v2.txt"
out_dir = root / "singlecell_communication_analysis"
out_dir.mkdir(exist_ok=True)

N_PERM = 1000          # 置换次数 (建议 >=1000; 跑不动可先设 100 试跑)
SEED = 20260726
MIN_DET_RATE = 0.05    # 配体/受体在该细胞群的最低检出率, 低于此值的边不检验

lr_pairs = [
    ("IFNG", "IFNGR1", "IFNG_response"),
    ("IFNG", "IFNGR2", "IFNG_response"),
    ("CXCL9", "CXCR3", "CXCL9_10_11_CXCR3"),
    ("CXCL10", "CXCR3", "CXCL9_10_11_CXCR3"),
    ("CXCL11", "CXCR3", "CXCL9_10_11_CXCR3"),
    ("CD274", "PDCD1", "PD1_PDL1"),
    ("PDCD1LG2", "PDCD1", "PD1_PDL1"),
    ("LGALS9", "HAVCR2", "Exhaustion_checkpoint"),
    ("NECTIN2", "TIGIT", "Exhaustion_checkpoint"),
    ("PVR", "TIGIT", "Exhaustion_checkpoint"),
    ("CD80", "CTLA4", "Exhaustion_checkpoint"),
    ("CD86", "CTLA4", "Exhaustion_checkpoint"),
    ("HLA-A", "CD8A", "MHC_I_CD8"),
    ("HLA-B", "CD8A", "MHC_I_CD8"),
    ("HLA-C", "CD8A", "MHC_I_CD8"),
    ("B2M", "CD8A", "MHC_I_CD8"),
]

genes_needed = sorted(set([g for p in lr_pairs for g in p[:2]]))


# ============================================================
# 数据读取 (沿用原 总.txt 的逻辑)
# ============================================================
def clean(x):
    return x.strip().strip('"')


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return math.nan


def gse72056_label(malignant_code, type_code):
    if malignant_code == "2":
        return "Malignant"
    if malignant_code == "1":
        mp = {"1": "T", "2": "B", "3": "Myeloid", "4": "Endothelial",
              "5": "CAF", "6": "NK", "0": "NonMalignant_Unknown", "": "NonMalignant_Unknown"}
        return mp.get(type_code, f"NonMalignant_{type_code}")
    return "Unresolved"


def read_data(path, genes):
    genes = set(genes)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        _tumor = f.readline()
        malignant_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        type_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        labels = np.array([gse72056_label(m, t) for m, t in zip(malignant_row[1:], type_row[1:])])
        gene_to_values = {}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            gene = clean(parts[0])
            if gene in genes:
                gene_to_values[gene] = np.array([safe_float(x) for x in parts[1:]], dtype=float)
    return labels, gene_to_values


# ============================================================
# 核心: 按标签聚合 + evidence score
# ============================================================
def aggregate_fast(gene_to_values, labels):
    """返回 {gene: {label: (mean_expr, det_rate)}} —— 向量化版本"""
    label_order = np.unique(labels)
    agg = {}
    for gene, vals in gene_to_values.items():
        agg[gene] = {}
        for lab in label_order:
            xs = vals[labels == lab]
            xs = xs[~np.isnan(xs)]
            if len(xs) == 0:
                continue
            agg[gene][lab] = (xs.mean(), (xs > 0).mean(), len(xs))
    return agg


def edge_score(ligand_stats, receptor_stats):
    lmean, ldet, _ = ligand_stats
    rmean, rdet, _ = receptor_stats
    return ldet * rdet * math.log1p(lmean) * math.log1p(rmean)


def compute_observed(agg, label_order):
    """计算所有候选边的观测分数"""
    rows = []
    for ligand, receptor, axis in lr_pairs:
        if ligand not in agg or receptor not in agg:
            continue
        for src in label_order:
            if src not in agg[ligand]:
                continue
            ls = agg[ligand][src]
            if ls[1] < MIN_DET_RATE:
                continue
            for tgt in label_order:
                if tgt not in agg[receptor]:
                    continue
                rs = agg[receptor][tgt]
                if rs[1] < MIN_DET_RATE:
                    continue
                rows.append({
                    "axis": axis, "ligand": ligand, "receptor": receptor,
                    "source_cell": src, "target_cell": tgt,
                    "obs_score": edge_score(ls, rs),
                })
    return pd.DataFrame(rows)


def permutation_null(gene_to_values, labels, edges_df, n_perm, rng):
    """对每条观测边, 置换标签算零分布"""
    label_order = np.unique(labels)
    n = len(labels)
    # 计数器: 每条边被超过的次数
    exceed = np.zeros(len(edges_df), dtype=int)
    valid_perm = 0

    # 预取需要的 (gene, label) 组合
    needed = set()
    for _, r in edges_df.iterrows():
        needed.add((r["ligand"], r["source_cell"]))
        needed.add((r["receptor"], r["target_cell"]))

    for b in range(n_perm):
        perm_labels = labels[rng.permutation(n)]
        # 重算所需统计量
        stats_cache = {}
        for gene, lab in needed:
            if gene not in gene_to_values:
                continue
            vals = gene_to_values[gene]
            xs = vals[perm_labels == lab]
            xs = xs[~np.isnan(xs)]
            if len(xs) == 0:
                continue
            stats_cache[(gene, lab)] = (xs.mean(), (xs > 0).mean(), len(xs))

        for i, (_, r) in enumerate(edges_df.iterrows()):
            ls = stats_cache.get((r["ligand"], r["source_cell"]))
            rs = stats_cache.get((r["receptor"], r["target_cell"]))
            if ls is None or rs is None:
                continue
            if edge_score(ls, rs) >= r["obs_score"]:
                exceed[i] += 1
        valid_perm += 1

        if (b + 1) % 200 == 0:
            print(f"  permutation {b+1}/{n_perm}")

    return exceed, valid_perm


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty(len(p))
    out[order] = q
    return out


# ============================================================
# 主流程
# ============================================================
def main():
    print("读取 GSE72056 ...")
    labels, gene_to_values = read_data(gse72056, genes_needed)
    print(f"细胞数: {len(labels)}, 基因覆盖: {len(gene_to_values)}/{len(genes_needed)}")

    label_order = np.unique(labels)
    agg = aggregate_fast(gene_to_values, labels)
    edges = compute_observed(agg, label_order)
    print(f"候选边数 (过检出率过滤): {len(edges)}")

    if len(edges) == 0:
        print("无边可检验, 退出")
        return

    rng = np.random.default_rng(SEED)
    print(f"开始置换检验 (n_perm={N_PERM}) ...")
    exceed, valid_perm = permutation_null(gene_to_values, labels, edges, N_PERM, rng)

    # 经验 p 值 (+1 校正, 置换检验标准做法)
    edges["perm_p"] = (exceed + 1) / (valid_perm + 1)
    edges["perm_fdr"] = bh_fdr(edges["perm_p"].values)
    edges = edges.sort_values(["perm_fdr", "perm_p", "obs_score"],
                              ascending=[True, True, False])

    out_path = out_dir / "lr_permutation_results.tsv"
    edges.to_csv(out_path, sep="\t", index=False)
    print(f"\n结果保存: {out_path}")

    sig = edges[edges["perm_fdr"] < 0.05]
    print(f"\nFDR<0.05 的显著边: {len(sig)}/{len(edges)}")
    if len(sig) > 0:
        print(sig[["axis", "ligand", "receptor", "source_cell", "target_cell",
                   "obs_score", "perm_p", "perm_fdr"]].head(30).to_string(index=False))

    print("\n解读: 只有置换显著的边才应进入机制叙事;")
    print("对原管线的高证据分但置换不显著的边, 建议在论文中降级为描述性观察。")


if __name__ == "__main__":
    main()
