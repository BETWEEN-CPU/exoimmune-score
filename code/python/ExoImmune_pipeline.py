# -*- coding: utf-8 -*-
"""
ExoImmune Project - Complete Integrated Analysis Pipeline
===========================================================
Merged from 13 step-by-step debug scripts. Later versions override earlier
conflicts. Run this script end-to-end in one pass.

Execution order:
  Part 1  Data exploration (from 1,2,3,4,9)
  Part 2  Single-cell source analysis (from 5)
  Part 3  Source summary & heatmap (from 6)
  Part 4  Cell-cell communication analysis (from 7)
  Part 5  High-confidence edges & narrative (from 8)
  Part 6  TCGA biomarker analysis (from 10)
  Part 7  TCGA Cox survival analysis (from 11, overrides Cox in 10)
  Part 8  IMvigor210 internal validation (from 13)
  Part 9  Project integrated summary (from 12)
"""

# ============================================================
#  IMPORTS
# ============================================================
from pathlib import Path
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")
import gzip
import math
import os

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

try:
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test
    HAS_LIFELINES = True
except Exception as e:
    HAS_LIFELINES = False
    print("lifelines not available; Cox/log-rank will be skipped:", repr(e))

# ============================================================
#  PATHS
# ============================================================
root = Path(os.environ.get("EXOIMMUNE_ROOT", r"D:/善意无穷倒退/李新月/EV"))

sc_source_dir = root / "singlecell_source_analysis"
sc_comm_dir = root / "singlecell_communication_analysis"
tcga_adj_dir = root / "TCGA_adjustment_biomarker_analysis"
imvigor_dir = root / "IMvigor210_internal_validation_biomarker"
summary_dir = root / "project_integrated_summary"

for d in [sc_source_dir, sc_comm_dir, tcga_adj_dir, imvigor_dir, summary_dir]:
    d.mkdir(exist_ok=True)

# Data file paths
gse72056 = root / "GSE72056_melanoma_single_cell_revised_v2.txt" / "GSE72056_melanoma_single_cell_revised_v2.txt"
gse115978 = root / "GSE115978_counts.csv.gz"
gse126044 = root / "GSE126044_counts.txt.gz"

# ============================================================
#  CONSTANTS  (later overrides earlier)
# ============================================================

# --- LASSO genes (universal) ---
lasso_coef = {
    "SH3TC1": -0.408018745460635,
    "CHST3": -0.386004323741037,
    "IFNG": 0.327810512159627,
    "NOMO3": 0.257418006143344,
    "KCND3": -0.228226990651115,
    "MMP10": -0.161545021847452,
    "LOC81691": 0.144321819042801,
    "UHRF1": 0.137177737980373,
    "MYBL2": 0.101269860842756,
}
lasso_genes = list(lasso_coef.keys())
gene_order = ["SH3TC1", "CHST3", "IFNG", "NOMO3", "KCND3", "MMP10", "LOC81691", "UHRF1", "MYBL2"]

# --- Marker sets for cell annotation (from 5.txt, used by 7.txt) ---
marker_sets = {
    "T_NK": ["CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY", "GZMB", "PRF1"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74"],
    "Myeloid": ["LYZ", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "CD68", "FCGR3A"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "ENG", "ESAM"],
    "CAF": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2", "PDPN"],
    "Malignant_Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "MLANA", "PMEL", "TYR", "MITF"],
    "Proliferation": ["MKI67", "TOP2A", "UBE2C", "HMGB2", "STMN1", "PCNA"],
    "Checkpoint": ["PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "Antigen_presentation": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5"],
    "IFNG_CXCL_axis": ["IFNG", "CXCL9", "CXCL10", "CXCL11", "CXCR3", "STAT1", "IRF1"],
}
all_needed_genes = sorted(set(lasso_genes + [g for gs in marker_sets.values() for g in gs]))

# --- Communication LR pairs (from 7.txt) ---
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
extra_genes = [
    "TAP1", "TAP2", "NLRC5", "STAT1", "IRF1",
    "GZMB", "PRF1", "NKG7", "CD3D", "CD3E",
    "LAG3", "TIGIT", "HAVCR2", "CTLA4", "PDCD1",
    "CD274", "PDCD1LG2",
]
comm_genes = sorted(set([x for p in lr_pairs for x in p[:2]] + extra_genes))

# --- Priority axes for filtering (from 8.txt) ---
priority_axes = [
    "IFNG_response", "CXCL9_10_11_CXCR3", "PD1_PDL1",
    "Exhaustion_checkpoint", "MHC_I_CD8",
]
priority_pairs = [
    ("IFNG_response", "IFNG", "IFNGR1"),
    ("IFNG_response", "IFNG", "IFNGR2"),
    ("CXCL9_10_11_CXCR3", "CXCL9", "CXCR3"),
    ("CXCL9_10_11_CXCR3", "CXCL10", "CXCR3"),
    ("CXCL9_10_11_CXCR3", "CXCL11", "CXCR3"),
    ("PD1_PDL1", "CD274", "PDCD1"),
    ("PD1_PDL1", "PDCD1LG2", "PDCD1"),
    ("Exhaustion_checkpoint", "LGALS9", "HAVCR2"),
    ("Exhaustion_checkpoint", "PVR", "TIGIT"),
    ("Exhaustion_checkpoint", "CD86", "CTLA4"),
    ("MHC_I_CD8", "HLA-A", "CD8A"),
    ("MHC_I_CD8", "HLA-B", "CD8A"),
    ("MHC_I_CD8", "HLA-C", "CD8A"),
    ("MHC_I_CD8", "B2M", "CD8A"),
]
axis_interpretation = {
    "IFNG_response": "T/NK-derived IFNG may activate IFNGR-positive myeloid and tumor/stromal cells, consistent with interferon-responsive inflamed TME.",
    "CXCL9_10_11_CXCR3": "Myeloid-derived CXCL9/CXCL10/CXCL11 may recruit or retain CXCR3-positive T/NK cells.",
    "PD1_PDL1": "Myeloid or stromal CD274/PDCD1LG2 signals toward PDCD1-positive T/NK cells, supporting checkpoint activation.",
    "Exhaustion_checkpoint": "LGALS9-HAVCR2, PVR-TIGIT and CD86-CTLA4 axes support T-cell exhaustion or inhibitory immune regulation.",
    "MHC_I_CD8": "HLA/B2M-CD8A signals indicate preserved antigen presentation or immune recognition activity.",
}

# --- TCGA gene sets (from 10.txt, renamed to avoid clash with 13.txt) ---
tcga_gene_sets = {
    "IFNG_signature": ["IFNG", "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "IDO1"],
    "CD8_score": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "cytolytic_score_expr": ["GZMA", "PRF1"],
    "checkpoint_expr_score": ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "T_cell_exhaustion_expr_score": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"],
    "MHC_I_expr_score": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5"],
    "myeloid_surrogate": ["LYZ", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "CD68", "FCGR3A"],
    "CAF_stromal_surrogate": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2", "PDPN"],
    "endothelial_surrogate": ["PECAM1", "VWF", "KDR", "ENG", "ESAM"],
    "proliferation_score": ["MKI67", "TOP2A", "UBE2C", "HMGB2", "STMN1", "PCNA"],
}
tcga_single_biomarkers = [
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
    "IFNG", "CXCL9", "CXCL10", "CXCL11", "CD8A", "CD8B", "GZMB", "PRF1",
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5",
    "COL1A1", "DCN", "PECAM1", "VWF", "LYZ", "AIF1", "MKI67", "TOP2A",
]

# --- IMvigor210 gene sets (from 13.txt, renamed) ---
imvigor_gene_sets = {
    "IFNG_signature": ["IFNG", "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "IDO1"],
    "CD8_score": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "cytolytic_score": ["GZMA", "PRF1"],
    "checkpoint_score": ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "T_cell_inflamed_score": ["CD8A", "CD8B", "CXCL9", "CXCL10", "IFNG", "GZMB", "PRF1", "IDO1", "STAT1"],
    "MHC_I_score": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5"],
    "exhaustion_score": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"],
    "myeloid_score": ["LYZ", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "CD68", "FCGR3A"],
    "CAF_stromal_score": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2", "PDPN"],
    "proliferation_score": ["MKI67", "TOP2A", "UBE2C", "HMGB2", "STMN1", "PCNA"],
}
imvigor_single_genes = [
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
    "IFNG", "CXCL9", "CXCL10", "CXCL11", "CD8A", "CD8B", "GZMB", "PRF1",
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5",
]

lasso_component_genes = set(lasso_genes)

# --- Biomarker flag mapping (from 12.txt) ---
biomarker_component_flag = {
    "IFNG": "score_component_single_gene",
    "IFNG_signature": "contains_score_component_IFNG",
    "checkpoint": "independent_module",
    "antigen_presentation": "independent_module",
    "cytolytic": "independent_module",
    "exhaustion": "independent_module",
    "t_cell_inflamed": "independent_module",
    "CD8_score": "independent_signature",
    "cytolytic_score_expr": "independent_signature",
    "checkpoint_expr_score": "independent_signature",
    "T_cell_exhaustion_expr_score": "independent_signature",
    "MHC_I_expr_score": "independent_signature",
    "CD274": "independent_single_gene",
    "PDCD1": "independent_single_gene",
    "CXCL9": "independent_single_gene",
    "CXCL10": "independent_single_gene",
    "CD8A": "independent_single_gene",
    "immune_surrogate": "adjustment_covariate",
    "stromal_surrogate": "adjustment_covariate",
    "proliferation_score": "adjustment_covariate",
}


# ============================================================
#  HELPER FUNCTIONS  (later overrides earlier conflicts)
# ============================================================

def clean(x):
    """Strip whitespace and quotes. Final version from 7.txt."""
    return x.strip().strip('"')


def open_text(path):
    """Open plain or gzipped text. Final version from 7.txt."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def safe_float(x):
    """Safe float conversion. Final version from 7.txt."""
    try:
        return float(x)
    except Exception:
        return math.nan


def gse72056_label(malignant_code, type_code):
    """Cell label mapping for GSE72056 annotation rows. Final version from 7.txt."""
    if malignant_code == "2":
        return "Malignant"
    if malignant_code == "1":
        mp = {
            "1": "T", "2": "B", "3": "Myeloid", "4": "Endothelial",
            "5": "CAF", "6": "NK", "0": "NonMalignant_Unknown", "": "NonMalignant_Unknown",
        }
        return mp.get(type_code, f"NonMalignant_{type_code}")
    return "Unresolved"


def zscore(s):
    """Z-score a pandas Series. From 10.txt (overrides earlier)."""
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean(skipna=True)) / sd


def zscore_series(s):
    """Z-score for Series, from 13.txt."""
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean(skipna=True)) / sd


def zscore_df_cols(df):
    """Z-score DataFrame columns, from 13.txt."""
    return df.apply(lambda x: (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) != 0 else x * np.nan, axis=0)


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR. Final version from 13.txt."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(len(pvals), np.nan)
    ok = ~np.isnan(pvals)
    p = pvals[ok]
    if len(p) == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    tmp = np.empty(len(p))
    tmp[order] = q
    out[ok] = tmp
    return out


def find_col(df, candidates):
    """Find column by candidate names (case-insensitive). From 13.txt."""
    cols_lower = {c.lower(): c for c in df.columns}
    for x in candidates:
        if x.lower() in cols_lower:
            return cols_lower[x.lower()]
    for c in df.columns:
        cl = c.lower()
        if any(x.lower() in cl for x in candidates):
            return c
    return None


def fmt_p(p):
    """Format p-value for display. From 12.txt."""
    try:
        p = float(p)
    except Exception:
        return "NA"
    if p == 0:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


# --- Single-cell source helpers (from 5.txt) ---

def read_selected_rows_matrix(path, sep, genes):
    """Read only selected gene rows from a delimited matrix."""
    genes = set(genes)
    with open_text(path) as f:
        header = f.readline().rstrip("\n").split(sep)
        header = [clean(x) for x in header]
        cell_names = header[1:]
        gene_to_values = {}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(sep)
            gene = clean(parts[0])
            if gene in genes:
                gene_to_values[gene] = [safe_float(x) for x in parts[1:]]
    return cell_names, gene_to_values


def summarize_gene_by_label(dataset, gene_to_values, labels, output_path):
    """Summarize gene expression per cell-type label. From 5.txt."""
    rows = []
    label_order = sorted(set(labels))
    for gene in lasso_genes:
        if gene not in gene_to_values:
            continue
        vals = gene_to_values[gene]
        for lab in label_order:
            idxs = [i for i, x in enumerate(labels) if x == lab]
            xs = [vals[i] for i in idxs if i < len(vals) and not math.isnan(vals[i])]
            if not xs:
                continue
            mean_expr = sum(xs) / len(xs)
            det_rate = sum(1 for x in xs if x > 0) / len(xs)
            rows.append({
                "dataset": dataset, "gene": gene, "coef": lasso_coef.get(gene, ""),
                "coef_direction": "positive" if lasso_coef.get(gene, 0) > 0 else "negative",
                "cell_type": lab, "n_cells": len(xs), "mean_expr": mean_expr,
                "det_rate": det_rate, "expressing_cells": sum(1 for x in xs if x > 0),
            })
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("dataset\tgene\tcoef\tcoef_direction\tcell_type\tn_cells\tmean_expr\tdet_rate\texpressing_cells\n")
        for r in rows:
            f.write(f"{r['dataset']}\t{r['gene']}\t{r['coef']}\t{r['coef_direction']}\t"
                    f"{r['cell_type']}\t{r['n_cells']}\t{r['mean_expr']:.6f}\t"
                    f"{r['det_rate']:.6f}\t{r['expressing_cells']}\n")
    return rows


def top_source_table(rows, output_path):
    """Pick top source per gene. From 5.txt."""
    best = {}
    for r in rows:
        key = (r["dataset"], r["gene"])
        score = (r["det_rate"], r["mean_expr"], r["n_cells"])
        if key not in best or score > best[key][0]:
            best[key] = (score, r)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("dataset\tgene\tcoef\tcoef_direction\ttop_cell_type\tn_cells\tmean_expr\tdet_rate\texpressing_cells\n")
        for key in sorted(best.keys()):
            r = best[key][1]
            f.write(f"{r['dataset']}\t{r['gene']}\t{r['coef']}\t{r['coef_direction']}\t"
                    f"{r['cell_type']}\t{r['n_cells']}\t{r['mean_expr']:.6f}\t"
                    f"{r['det_rate']:.6f}\t{r['expressing_cells']}\n")


# --- Source summary helpers (from 6.txt) ---

def coef_sign(x):
    return "positive" if x > 0 else "negative"


def broad_category(cell_type):
    ct = str(cell_type)
    if ct in ["T", "T_NK", "NK"]:
        return "Immune_T_NK"
    if ct in ["B"]:
        return "Immune_B"
    if ct in ["Myeloid", "Macro"]:
        return "Immune_Myeloid"
    if ct in ["Endothelial", "Endo"]:
        return "Stromal_Endothelial"
    if ct in ["CAF"]:
        return "Stromal_CAF"
    if ct in ["Malignant", "Malignant_Epithelial", "TumorOrStromal_CD45neg", "TumorOrStromal_CD45neg_CD90neg"]:
        return "Tumor_or_Unknown_CD45neg"
    if ct in ["NonMalignant_Unknown", "Unresolved", "Immune_CD45pos", "Stromal_CD45neg_CD90pos", "Unknown"]:
        return ct
    return ct


def pick_top_source(df, dataset):
    rows = []
    for gene in gene_order:
        sub = df[df["gene"] == gene].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(["det_rate", "mean_expr", "n_cells"], ascending=[False, False, False])
        top = sub.iloc[0]
        rows.append({
            "dataset": dataset, "gene": gene, "coef": float(top["coef"]),
            "coef_sign": coef_sign(float(top["coef"])),
            "top_cell_type": top["cell_type"],
            "top_broad_category": broad_category(top["cell_type"]),
            "top_n_cells": int(top["n_cells"]),
            "top_mean_expr": float(top["mean_expr"]),
            "top_det_rate": float(top["det_rate"]),
            "expressing_cells": int(top["expressing_cells"]),
        })
    return pd.DataFrame(rows)


# --- Communication helpers (from 7.txt) ---

def read_selected_gse72056(path, genes):
    genes = set(genes)
    with open_text(path) as f:
        header = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        tumor_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        malignant_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        type_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        cells = header[1:]
        labels = [gse72056_label(m, t) for m, t in zip(malignant_row[1:], type_row[1:])]
        gene_to_values = {}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            gene = clean(parts[0])
            if gene in genes:
                gene_to_values[gene] = [safe_float(x) for x in parts[1:]]
    return cells, labels, gene_to_values


def read_selected_delim(path, sep, genes):
    genes = set(genes)
    with open_text(path) as f:
        header = [clean(x) for x in f.readline().rstrip("\n").split(sep)]
        cells = header[1:]
        gene_to_values = {}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(sep)
            gene = clean(parts[0])
            if gene in genes:
                gene_to_values[gene] = [safe_float(x) for x in parts[1:]]
    return cells, gene_to_values


def aggregate_by_label(dataset, gene_to_values, labels):
    """Aggregate gene expression by cell-type label."""
    rows = []
    label_order = sorted(set(labels))
    for gene, vals in gene_to_values.items():
        for lab in label_order:
            idx = [i for i, x in enumerate(labels) if x == lab]
            xs = [vals[i] for i in idx if i < len(vals) and not math.isnan(vals[i])]
            if not xs:
                continue
            mean_expr = sum(xs) / len(xs)
            det_rate = sum(1 for x in xs if x > 0) / len(xs)
            rows.append({
                "dataset": dataset, "gene": gene, "cell_type": lab,
                "n_cells": len(xs), "mean_expr": mean_expr,
                "det_rate": det_rate, "expressing_cells": sum(1 for x in xs if x > 0),
            })
    return pd.DataFrame(rows)


def make_lr_table(dataset, agg):
    """Build ligand-receptor evidence table."""
    rows = []
    for ligand, receptor, axis in lr_pairs:
        ligand_df = agg[agg["gene"] == ligand]
        receptor_df = agg[agg["gene"] == receptor]
        if ligand_df.empty or receptor_df.empty:
            rows.append({
                "dataset": dataset, "axis": axis, "ligand": ligand, "receptor": receptor,
                "source_cell": "MISSING", "target_cell": "MISSING",
                "ligand_mean": np.nan, "ligand_det_rate": np.nan,
                "receptor_mean": np.nan, "receptor_det_rate": np.nan,
                "evidence_score": np.nan, "status": "missing_ligand_or_receptor",
            })
            continue
        for _, lrow in ligand_df.iterrows():
            for _, rrow in receptor_df.iterrows():
                score = (
                    float(lrow["det_rate"]) * float(rrow["det_rate"]) *
                    math.log1p(float(lrow["mean_expr"])) * math.log1p(float(rrow["mean_expr"]))
                )
                rows.append({
                    "dataset": dataset, "axis": axis, "ligand": ligand, "receptor": receptor,
                    "source_cell": lrow["cell_type"], "target_cell": rrow["cell_type"],
                    "ligand_mean": float(lrow["mean_expr"]), "ligand_det_rate": float(lrow["det_rate"]),
                    "receptor_mean": float(rrow["mean_expr"]), "receptor_det_rate": float(rrow["det_rate"]),
                    "evidence_score": score, "status": "ok",
                })
    return pd.DataFrame(rows)


# --- TCGA helpers (from 10.txt) ---

def residualize(y, covars):
    y = pd.to_numeric(y, errors="coerce")
    X = covars.copy().apply(pd.to_numeric, errors="coerce")
    dat = pd.concat([y.rename("y"), X], axis=1).dropna()
    if dat.shape[0] < 10 or X.shape[1] == 0:
        return pd.Series(index=y.index, data=np.nan)
    X2 = sm.add_constant(dat[X.columns], has_constant="add")
    fit = sm.OLS(dat["y"], X2).fit()
    resid = pd.Series(index=y.index, data=np.nan)
    resid.loc[dat.index] = fit.resid
    return resid


def partial_spearman(x, y, covars):
    dat = pd.concat([x.rename("x"), y.rename("y"), covars], axis=1).dropna()
    if dat.shape[0] < 20:
        return np.nan, np.nan, dat.shape[0]
    rx = residualize(dat["x"], dat[covars.columns])
    ry = residualize(dat["y"], dat[covars.columns])
    ok = rx.notna() & ry.notna()
    if ok.sum() < 20:
        return np.nan, np.nan, ok.sum()
    rho, p = stats.spearmanr(rx[ok], ry[ok])
    return rho, p, ok.sum()


def read_expression_matrix(path):
    expr = pd.read_csv(path, sep="\t")
    if "gene_symbol" not in expr.columns:
        raise ValueError(f"gene_symbol column not found in {path}")
    expr = expr.drop_duplicates("gene_symbol").set_index("gene_symbol")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr


def compute_tcga_signature_scores(expr):
    """Compute TCGA biomarker signature scores. From 10.txt."""
    samples = expr.columns
    out = pd.DataFrame(index=samples)
    present_report = []
    for name, genes in tcga_gene_sets.items():
        present = [g for g in genes if g in expr.index]
        missing = [g for g in genes if g not in expr.index]
        present_report.append({
            "signature": name, "n_genes": len(genes), "n_present": len(present),
            "present_genes": ",".join(present), "missing_genes": ",".join(missing),
        })
        if len(present) == 0:
            out[name] = np.nan
        else:
            mat = expr.loc[present].T
            out[name] = mat.apply(zscore, axis=0).mean(axis=1)
    for g in tcga_single_biomarkers:
        out[g] = expr.loc[g] if g in expr.index else np.nan
    out["immune_surrogate"] = out[["IFNG_signature", "CD8_score", "myeloid_surrogate"]].mean(axis=1)
    out["stromal_surrogate"] = out[["CAF_stromal_surrogate", "endothelial_surrogate"]].mean(axis=1)
    return out.reset_index().rename(columns={"index": "sample_id"}), pd.DataFrame(present_report)


# --- Cox survival helpers (from 11.txt, overrides 10.txt Cox) ---

# (uses zscore, HAS_LIFELINES from global)


# --- IMvigor210 helpers (from 13.txt) ---

def read_expr(path):
    expr = pd.read_csv(path, sep="\t")
    first = expr.columns[0]
    expr = expr.rename(columns={first: "gene"})
    expr["gene"] = expr["gene"].astype(str)
    expr = expr.drop_duplicates("gene").set_index("gene")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr


def read_response(path):
    df = pd.read_csv(path, sep="\t")
    if "sample_id" not in df.columns:
        raise ValueError(f"sample_id column not found in {path}. Columns={list(df.columns)}")
    if "response_clean" in df.columns:
        response_col = "response_clean"
    elif "response_binary" in df.columns:
        response_col = "response_binary"
    elif "response_best" in df.columns:
        response_col = "response_best"
    else:
        raise ValueError(f"Cannot infer response column from {path}. Columns={list(df.columns)}")
    out = df[["sample_id", response_col]].copy()
    out.columns = ["sample_id", "response_raw"]
    out["response_raw"] = out["response_raw"].astype(str).str.strip()
    response_map = {
        "Responder": 1, "NonResponder": 0, "CR/PR": 1, "SD/PD": 0,
        "CR": 1, "PR": 1, "SD": 0, "PD": 0,
    }
    out["y"] = out["response_raw"].map(response_map)
    out = out.dropna(subset=["y"]).copy()
    out["y"] = out["y"].astype(int)
    out["response"] = out["y"].map({1: "Responder", 0: "NonResponder"})
    return out[["sample_id", "response", "y"]].drop_duplicates("sample_id")


def read_score(path):
    df = pd.read_csv(path, sep="\t")
    if "sample_id" not in df.columns:
        raise ValueError(f"sample_id column not found in {path}. Columns={list(df.columns)}")
    if "score" not in df.columns:
        raise ValueError(f"score column not found in {path}. Columns={list(df.columns)}")
    out = df[["sample_id", "score"]].copy()
    out.columns = ["sample_id", "ExoImmune_Score"]
    out["ExoImmune_Score"] = pd.to_numeric(out["ExoImmune_Score"], errors="coerce")
    return out.drop_duplicates("sample_id")


def compute_imvigor_biomarkers(expr):
    """Compute IMvigor210 biomarker scores. From 13.txt."""
    samples = expr.columns
    biom = pd.DataFrame(index=samples)
    coverage = []
    for name, genes in imvigor_gene_sets.items():
        present = [g for g in genes if g in expr.index]
        missing = [g for g in genes if g not in expr.index]
        coverage.append({
            "biomarker": name, "type": "signature", "n_genes": len(genes),
            "n_present": len(present), "present_genes": ",".join(present),
            "missing_genes": ",".join(missing),
            "contains_lasso_component": any(g in lasso_component_genes for g in genes),
        })
        if len(present) == 0:
            biom[name] = np.nan
        else:
            mat = expr.loc[present].T
            biom[name] = zscore_df_cols(mat).mean(axis=1)
    for g in imvigor_single_genes:
        coverage.append({
            "biomarker": g, "type": "single_gene", "n_genes": 1,
            "n_present": int(g in expr.index), "present_genes": g if g in expr.index else "",
            "missing_genes": "" if g in expr.index else g,
            "contains_lasso_component": g in lasso_component_genes,
        })
        biom[g] = expr.loc[g] if g in expr.index else np.nan
    biom = biom.reset_index().rename(columns={"index": "sample_id"})
    return biom, pd.DataFrame(coverage)


def bootstrap_auc(y, score, n_boot=2000):
    y = np.asarray(y)
    score = np.asarray(score)
    aucs = []
    n = len(y)
    rng = np.random.default_rng(20260628)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], score[idx]))
    return np.mean(aucs), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), len(aucs)


def repeated_cv_auc(X, y, n_splits=5, n_repeats=50):
    aucs = []
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=20260628)
    for tr, te in rskf.split(X, y):
        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(solver="liblinear", penalty="l2", max_iter=1000))
        ])
        pipe.fit(X.iloc[tr], y.iloc[tr])
        pred = pipe.predict_proba(X.iloc[te])[:, 1]
        aucs.append(roc_auc_score(y.iloc[te], pred))
    return np.mean(aucs), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), len(aucs)


def fit_logistic_auc(df, features, y_col="y"):
    dat = df[[y_col] + features].replace([np.inf, -np.inf], np.nan).dropna()
    if dat.shape[0] < 30 or dat[y_col].nunique() < 2:
        return None
    X = dat[features]
    y = dat[y_col]
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(solver="liblinear", penalty="l2", max_iter=1000))
    ])
    pipe.fit(X, y)
    pred = pipe.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, pred)
    cv_mean, cv_low, cv_high, cv_n = repeated_cv_auc(X, y, n_splits=5, n_repeats=50)
    return {
        "n": dat.shape[0], "responders": int(y.sum()),
        "features": "+".join(features), "apparent_auc": auc,
        "cv_auc_mean": cv_mean, "cv_auc_low": cv_low, "cv_auc_high": cv_high, "cv_n": cv_n,
        "pred": pd.Series(pred, index=dat.index), "y": y,
    }


def decision_curve(y, pred, thresholds):
    y = np.asarray(y)
    pred = np.asarray(pred)
    n = len(y)
    prevalence = y.mean()
    rows = []
    for pt in thresholds:
        if pt <= 0 or pt >= 1:
            continue
        treat = pred >= pt
        tp = np.sum((treat == 1) & (y == 1))
        fp = np.sum((treat == 1) & (y == 0))
        nb = tp / n - fp / n * (pt / (1 - pt))
        nb_all = prevalence - (1 - prevalence) * (pt / (1 - pt))
        nb_none = 0.0
        rows.append({"threshold": pt, "net_benefit_model": nb, "net_benefit_all": nb_all, "net_benefit_none": nb_none})
    return pd.DataFrame(rows)


# ============================================================
#  PART 1: SINGLE-CELL SOURCE ANALYSIS  (from 5.txt)
# ============================================================

def run_part2_source_analysis():
    print("\n" + "=" * 80)
    print("Part 1. Single-cell source analysis of LASSO genes")
    print("=" * 80)

    # 1) GSE72056: annotated source analysis
    print("\n2.1 GSE72056 annotated source analysis")
    with open_text(gse72056) as f:
        header = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        tumor_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        malignant_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]
        type_row = [clean(x) for x in f.readline().rstrip("\n").split("\t")]

        cell_names_72056 = header[1:]
        malignant_codes = malignant_row[1:]
        type_codes = type_row[1:]
        labels_72056 = [gse72056_label(m, t) for m, t in zip(malignant_codes, type_codes)]

        gene_to_values_72056 = {}
        need = set(lasso_genes)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            gene = clean(parts[0])
            if gene in need:
                gene_to_values_72056[gene] = [safe_float(x) for x in parts[1:]]

    n_found_72056 = len(gene_to_values_72056)
    print(f"  GSE72056: {len(cell_names_72056)} single cells annotated by Tirosh et al. labels.")
    print(f"  Mean expression and detection rate computed per cell type for {n_found_72056}/9 LASSO genes.")

    rows_72056 = summarize_gene_by_label(
        "GSE72056", gene_to_values_72056, labels_72056,
        sc_source_dir / "GSE72056_lasso_gene_by_celltype.tsv"
    )
    top_source_table(rows_72056, sc_source_dir / "GSE72056_lasso_gene_top_source.tsv")

    # 2) GSE115978: marker-based rough annotation
    print("\n2.2 GSE115978 marker-based annotation")
    cell_names_115978, gene_to_values_115978 = read_selected_rows_matrix(gse115978, ",", all_needed_genes)
    n_cells = len(cell_names_115978)
    module_scores = {module: [0.0] * n_cells for module in marker_sets}
    module_gene_counts = {}
    for module, genes in marker_sets.items():
        present = [g for g in genes if g in gene_to_values_115978]
        module_gene_counts[module] = len(present)
        if not present:
            continue
        for i in range(n_cells):
            vals = []
            for g in present:
                x = gene_to_values_115978[g][i]
                if not math.isnan(x):
                    vals.append(math.log1p(x))
            module_scores[module][i] = sum(vals) / len(vals) if vals else 0.0

    # Conservative rule-based labels
    labels_115978 = []
    for i, cell in enumerate(cell_names_115978):
        c = cell.lower()
        explicit = None
        if "cd45pos" in c or "cd45_pos" in c or "cd45-pos" in c:
            explicit = "Immune_CD45pos"
        elif "cd45neg" in c or "cd45_neg" in c or "cd45-neg" in c:
            if "cd90pos" in c or "cd90_pos" in c or "cd90-pos" in c:
                explicit = "Stromal_CD45neg_CD90pos"
            elif "cd90neg" in c or "cd90_neg" in c or "cd90-neg" in c:
                explicit = "TumorOrStromal_CD45neg_CD90neg"
            else:
                explicit = "TumorOrStromal_CD45neg"

        identity_modules = ["T_NK", "B", "Myeloid", "Endothelial", "CAF", "Malignant_Epithelial"]
        scores = {m: module_scores[m][i] for m in identity_modules}
        winner = max(scores, key=scores.get)
        winner_score = scores[winner]

        if explicit == "Immune_CD45pos":
            label = winner if winner in ["T_NK", "B", "Myeloid"] else "Immune_CD45pos"
        elif explicit == "Stromal_CD45neg_CD90pos":
            label = winner if winner in ["CAF", "Endothelial"] else "Stromal_CD45neg_CD90pos"
        elif explicit in ["TumorOrStromal_CD45neg_CD90neg", "TumorOrStromal_CD45neg"]:
            label = winner if winner in ["Malignant_Epithelial", "CAF", "Endothelial"] else explicit
        else:
            label = winner if winner_score > 0 else "Unknown"
        labels_115978.append(label)

    print(f"  GSE115978: {n_cells} cells with marker-based rule annotation.")
    print(f"  Identity modules: { {k:v for k,v in Counter(labels_115978).items() if v>10} }")

    # Save cell annotation
    with open(sc_source_dir / "GSE115978_inferred_cell_annotation.tsv", "w", encoding="utf-8") as f:
        f.write("cell\tinferred_cell_type")
        for module in marker_sets:
            f.write(f"\t{module}_score")
        f.write("\n")
        for i, cell in enumerate(cell_names_115978):
            f.write(f"{cell}\t{labels_115978[i]}")
            for module in marker_sets:
                f.write(f"\t{module_scores[module][i]:.6f}")
            f.write("\n")

    rows_115978 = summarize_gene_by_label(
        "GSE115978_inferred",
        {g: gene_to_values_115978[g] for g in lasso_genes if g in gene_to_values_115978},
        labels_115978,
        sc_source_dir / "GSE115978_lasso_gene_by_inferred_celltype.tsv"
    )
    top_source_table(rows_115978, sc_source_dir / "GSE115978_lasso_gene_top_source.tsv")

    # 3) GSE126044: coverage only
    print("\n2.3 GSE126044 coverage check")
    cell_names_126044, gene_to_values_126044 = read_selected_rows_matrix(gse126044, "\t", lasso_genes)
    print(f"  GSE126044: {len(cell_names_126044)}-column matrix; sample-level bulk data, excluded from single-cell source evidence.")


# ============================================================
#  PART 3: SOURCE SUMMARY & HEATMAP  (from 6.txt)
# ============================================================

def run_part3_source_summary():
    print("\n" + "=" * 100)
    print("PART 3: SOURCE SUMMARY & HEATMAP")
    print("=" * 100)

    gse72056_df = pd.read_csv(sc_source_dir / "GSE72056_lasso_gene_by_celltype.tsv", sep="\t")
    gse115978_df = pd.read_csv(sc_source_dir / "GSE115978_lasso_gene_by_inferred_celltype.tsv", sep="\t")

    top72056 = pick_top_source(gse72056_df, "GSE72056")
    top115978 = pick_top_source(gse115978_df, "GSE115978_inferred")

    summary = pd.concat([top72056, top115978], ignore_index=True)
    summary["gene"] = pd.Categorical(summary["gene"], categories=gene_order, ordered=True)
    summary = summary.sort_values(["gene", "dataset"])
    summary_path = sc_source_dir / "singlecell_lasso_source_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(summary.to_string(index=False))

    # Category enrichment
    cat_summary = (
        summary.groupby(["coef_sign", "top_broad_category"])
        .agg(n=("gene", "count"), genes=("gene", lambda x: ",".join(list(x))),
             mean_det_rate=("top_det_rate", "mean"), mean_expr=("top_mean_expr", "mean"))
        .reset_index().sort_values(["coef_sign", "n"], ascending=[True, False])
    )
    cat_path = sc_source_dir / "singlecell_lasso_source_category_summary.tsv"
    cat_summary.to_csv(cat_path, sep="\t", index=False)
    print("\nCategory summary:")
    print(cat_summary.to_string(index=False))

    # Consensus
    consensus_rows = []
    for gene in gene_order:
        sub = summary[summary["gene"] == gene]
        if sub.empty:
            continue
        categories = list(sub["top_broad_category"])
        celltypes = list(sub["top_cell_type"])
        consensus_rows.append({
            "gene": gene, "coef": lasso_coef[gene], "coef_sign": coef_sign(lasso_coef[gene]),
            "dataset1_celltype": sub.iloc[0]["top_cell_type"] if len(sub) > 0 else "",
            "dataset1_category": sub.iloc[0]["top_broad_category"] if len(sub) > 0 else "",
            "dataset2_celltype": sub.iloc[1]["top_cell_type"] if len(sub) > 1 else "",
            "dataset2_category": sub.iloc[1]["top_broad_category"] if len(sub) > 1 else "",
            "category_overlap": len(set(categories)) == 1 if len(categories) > 1 else True,
            "celltype_set": " | ".join(celltypes),
        })
    consensus = pd.DataFrame(consensus_rows)
    consensus_path = sc_source_dir / "singlecell_lasso_gene_consensus.tsv"
    consensus.to_csv(consensus_path, sep="\t", index=False)
    print("\nGene consensus:")
    print(consensus.to_string(index=False))

    # Heatmap
    hm = summary.pivot(index="gene", columns="dataset", values="top_det_rate").reindex(gene_order)
    plt.figure(figsize=(7, 4.8))
    sns.heatmap(hm, annot=True, fmt=".2f", cmap="YlOrRd", cbar_kws={"label": "Top detection rate"})
    plt.title("LASSO genes: top source detection rate across single-cell datasets")
    plt.tight_layout()
    fig_path = sc_source_dir / "singlecell_lasso_top_source_heatmap.pdf"
    plt.savefig(fig_path)
    plt.close()

    # Quick readout
    print("\nQuick takeaway:")
    for gene in gene_order:
        sub = summary[summary["gene"] == gene]
        if sub.empty:
            continue
        txt = [f"{r['dataset']}->{r['top_cell_type']}({r['top_broad_category']})" for _, r in sub.iterrows()]
        print(f"{gene:10s} coef={lasso_coef[gene]:+.3f}  " + " ; ".join(txt))

    print(f"\n[Part 3 complete]  Saved: {summary_path}, {cat_path}, {consensus_path}, {fig_path}")


# ============================================================
#  PART 4: CELL-CELL COMMUNICATION  (from 7.txt)
# ============================================================

def run_part4_communication():
    print("\n" + "=" * 100)
    print("PART 4: CELL-CELL COMMUNICATION ANALYSIS")
    print("=" * 100)

    # GSE72056
    print("\n4.1 GSE72056")
    cells72056, labels72056, genes72056 = read_selected_gse72056(gse72056, comm_genes)
    print(f"Cells: {len(cells72056)}, Labels: {dict(Counter(labels72056))}")
    print(f"Genes found: {len(genes72056)}/{len(comm_genes)}")
    print(f"Missing: {sorted(set(comm_genes) - set(genes72056))}")

    agg72056 = aggregate_by_label("GSE72056", genes72056, labels72056)
    agg72056.to_csv(sc_comm_dir / "GSE72056_communication_gene_by_celltype.tsv", sep="\t", index=False)

    lr72056 = make_lr_table("GSE72056", agg72056)
    lr72056.to_csv(sc_comm_dir / "GSE72056_targeted_lr_evidence.tsv", sep="\t", index=False)

    top72056 = lr72056[lr72056["status"] == "ok"].sort_values("evidence_score", ascending=False).head(30)
    print("\nTop GSE72056 LR evidence:")
    print(top72056[["axis", "ligand", "receptor", "source_cell", "target_cell", "evidence_score"]].to_string(index=False))

    # GSE115978
    print("\n4.2 GSE115978 with inferred labels")
    ann115978 = pd.read_csv(sc_source_dir / "GSE115978_inferred_cell_annotation.tsv", sep="\t")
    labels115978 = ann115978["inferred_cell_type"].astype(str).tolist()
    cells115978, genes115978 = read_selected_delim(gse115978, ",", comm_genes)

    if len(cells115978) != len(labels115978):
        print("WARNING: cell count mismatch", len(cells115978), len(labels115978))
    else:
        mismatch = sum(c1 != c2 for c1, c2 in zip(cells115978, ann115978["cell"].astype(str)))
        print(f"Cell name mismatches vs annotation: {mismatch}")

    print(f"Cells: {len(cells115978)}, Labels: {dict(Counter(labels115978))}")
    print(f"Genes found: {len(genes115978)}/{len(comm_genes)}")

    agg115978 = aggregate_by_label("GSE115978_inferred", genes115978, labels115978)
    agg115978.to_csv(sc_comm_dir / "GSE115978_communication_gene_by_inferred_celltype.tsv", sep="\t", index=False)

    lr115978 = make_lr_table("GSE115978_inferred", agg115978)
    lr115978.to_csv(sc_comm_dir / "GSE115978_targeted_lr_evidence.tsv", sep="\t", index=False)

    top115978 = lr115978[lr115978["status"] == "ok"].sort_values("evidence_score", ascending=False).head(30)
    print("\nTop GSE115978 LR evidence:")
    print(top115978[["axis", "ligand", "receptor", "source_cell", "target_cell", "evidence_score"]].to_string(index=False))

    # Combined summary
    print("\n4.3 Combined targeted LR summary")
    combined = pd.concat([lr72056, lr115978], ignore_index=True)
    combined_ok = combined[combined["status"] == "ok"].copy()

    top_by_pair = (
        combined_ok.sort_values("evidence_score", ascending=False)
        .groupby(["dataset", "axis", "ligand", "receptor"], as_index=False).head(1)
        .sort_values(["dataset", "axis", "evidence_score"], ascending=[True, True, False])
    )
    top_by_pair.to_csv(sc_comm_dir / "targeted_lr_top_source_target_by_pair.tsv", sep="\t", index=False)

    axis_summary = (
        combined_ok.groupby(["dataset", "axis"])
        .agg(max_evidence=("evidence_score", "max"), median_evidence=("evidence_score", "median"),
             n_pairs=("evidence_score", "size"))
        .reset_index().sort_values(["dataset", "max_evidence"], ascending=[True, False])
    )
    axis_summary.to_csv(sc_comm_dir / "targeted_lr_axis_summary.tsv", sep="\t", index=False)
    print(axis_summary.to_string(index=False))

    # Heatmap
    hm = axis_summary.pivot(index="axis", columns="dataset", values="max_evidence").fillna(0)
    plt.figure(figsize=(7.5, 4.8))
    sns.heatmap(hm, annot=True, fmt=".3f", cmap="Reds", cbar_kws={"label": "Max LR evidence score"})
    plt.title("Targeted communication axes in single-cell datasets")
    plt.tight_layout()
    fig = sc_comm_dir / "targeted_communication_axis_heatmap.pdf"
    plt.savefig(fig)
    plt.close()

    print(f"\n[Part 4 complete]  Saved files: {list(sc_comm_dir.glob('*.tsv'))} + {fig}")


# ============================================================
#  PART 5: HIGH-CONFIDENCE EDGES & NARRATIVE  (from 8.txt)
# ============================================================

def run_part5_comm_edges():
    print("\n" + "=" * 100)
    print("PART 5: HIGH-CONFIDENCE EDGES & MECHANISM NARRATIVE")
    print("=" * 100)

    top_path = sc_comm_dir / "targeted_lr_top_source_target_by_pair.tsv"
    axis_path = sc_comm_dir / "targeted_lr_axis_summary.tsv"
    source_path = sc_source_dir / "singlecell_lasso_gene_consensus.tsv"

    top = pd.read_csv(top_path, sep="\t")
    axis = pd.read_csv(axis_path, sep="\t")
    source = pd.read_csv(source_path, sep="\t")

    top["pair_key"] = list(zip(top["axis"], top["ligand"], top["receptor"]))
    keep_keys = set(priority_pairs)
    high = top[top["pair_key"].isin(keep_keys)].copy()

    high["edge"] = high["source_cell"].astype(str) + " -> " + high["target_cell"].astype(str)
    high["lr_pair"] = high["ligand"].astype(str) + "-" + high["receptor"].astype(str)

    pair_dataset_counts = high.groupby(["axis", "ligand", "receptor"])["dataset"].nunique().reset_index(name="n_datasets")
    high = high.merge(pair_dataset_counts, on=["axis", "ligand", "receptor"], how="left")
    high["support_level"] = np.where(high["n_datasets"] >= 2, "cross_dataset", "single_dataset")

    high["axis_rank_in_dataset"] = high.groupby(["dataset", "axis"])["evidence_score"].rank(ascending=False, method="dense")
    high["evidence_tier"] = np.where(
        high["support_level"].eq("cross_dataset") & high["axis_rank_in_dataset"].le(3),
        "high",
        np.where(high["axis_rank_in_dataset"].le(3), "moderate", "supportive")
    )
    high = high.sort_values(["axis", "ligand", "receptor", "dataset"])

    high_path = sc_comm_dir / "targeted_lr_high_confidence_edges.tsv"
    high.drop(columns=["pair_key"]).to_csv(high_path, sep="\t", index=False)

    print("High confidence / priority edges:")
    print(high[["dataset", "axis", "lr_pair", "edge", "evidence_score", "support_level", "evidence_tier"]].to_string(index=False))

    # Mechanism narrative
    narrative_rows = []
    for ax in priority_axes:
        sub = high[high["axis"] == ax].copy()
        if sub.empty:
            continue
        best = sub.sort_values("evidence_score", ascending=False).iloc[0]
        datasets = ",".join(sorted(sub["dataset"].unique()))
        pairs = "; ".join(sub.sort_values("evidence_score", ascending=False).drop_duplicates(["ligand", "receptor"]).head(5)["lr_pair"].tolist())
        edges = "; ".join(sub.sort_values("evidence_score", ascending=False).drop_duplicates(["source_cell", "target_cell", "lr_pair"]).head(5)["edge"].tolist())
        narrative_rows.append({
            "axis": ax, "datasets_supporting": datasets, "best_pair": best["lr_pair"],
            "best_edge": best["edge"], "best_evidence_score": best["evidence_score"],
            "representative_pairs": pairs, "representative_edges": edges,
            "interpretation": axis_interpretation.get(ax, ""),
        })

    narrative = pd.DataFrame(narrative_rows)
    narrative_path = sc_comm_dir / "targeted_lr_mechanism_narrative.tsv"
    narrative.to_csv(narrative_path, sep="\t", index=False)
    print("\nMechanism narrative:")
    print(narrative.to_string(index=False))

    # Key results text
    notes = []
    notes.append("Single-cell source analysis:")
    notes.append("- IFNG was consistently localized to T/T_NK cells across GSE72056 and GSE115978, supporting a lymphocyte-derived positive ExoImmune Score component.")
    notes.append("- NOMO3 was consistently localized to malignant tumor cells, whereas UHRF1 and MYBL2 were enriched in tumor/stromal CD45-negative or proliferative compartments.")
    notes.append("- LOC81691 was enriched in endothelial/CAF compartments, suggesting a stromal contribution to the positive score.")
    notes.append("- Negative-score genes including SH3TC1, CHST3, KCND3 and MMP10 were mainly distributed across myeloid, endothelial, CAF or CD90-positive stromal compartments.")
    notes.append("")
    notes.append("Targeted cell-cell communication analysis:")
    for _, r in narrative.iterrows():
        notes.append(f"- {r['axis']}: {r['interpretation']} Best evidence: {r['best_pair']} on {r['best_edge']}.")
    notes.append("")
    notes.append("Interpretation:")
    notes.append("These findings support a model in which ExoImmune Score captures a multicellular exosome-associated program rather than a single-cell-lineage marker. The positive score component links T/NK-derived IFNG, malignant-cell programs, proliferative CD45-negative compartments and stromal/endothelial signals. Communication analysis further connects this program to IFNG response, CXCL9/10/11-CXCR3 inflammatory recruitment, MHC-I antigen presentation and PD-1/exhaustion checkpoint signaling.")

    notes_path = sc_comm_dir / "singlecell_mechanism_key_results.txt"
    with open(notes_path, "w", encoding="utf-8") as f:
        f.write("\n".join(notes))

    # Bubble plot
    plot_df = high[high["evidence_tier"].isin(["high", "moderate", "supportive"])].copy()
    non_mhc = plot_df[plot_df["axis"] != "MHC_I_CD8"].copy()
    mhc = plot_df[plot_df["axis"] == "MHC_I_CD8"].sort_values("evidence_score", ascending=False).head(8)
    plot_df = pd.concat([non_mhc, mhc], ignore_index=True).sort_values(["axis", "evidence_score"], ascending=[True, False])
    plot_df["y_label"] = plot_df["dataset"] + " | " + plot_df["edge"]
    plot_df["x_label"] = plot_df["axis"] + "\n" + plot_df["lr_pair"]

    plt.figure(figsize=(12, max(5, 0.25 * plot_df["y_label"].nunique())))
    sns.scatterplot(data=plot_df, x="x_label", y="y_label", size="evidence_score", hue="axis",
                    sizes=(40, 500), alpha=0.82, edgecolor="black", linewidth=0.3)
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("")
    plt.ylabel("")
    plt.title("Targeted single-cell communication evidence")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    plt.tight_layout()
    bubble_path = sc_comm_dir / "targeted_lr_priority_edges_bubbleplot.pdf"
    plt.savefig(bubble_path)
    plt.close()

    print(f"\n[Part 5 complete]  Saved: {high_path}, {narrative_path}, {notes_path}, {bubble_path}")


# ============================================================
#  PART 6: TCGA BIOMARKER ANALYSIS  (from 10.txt)
# ============================================================

def run_part6_tcga_biomarker():
    print("\n" + "=" * 100)
    print("PART 6: TCGA BIOMARKER ANALYSIS")
    print("=" * 100)

    tcga_dir = root / "EV_immunotherapy_project" / "TCGA_processed"
    mech_dir = root / "TCGA_mechanism_prognosis"
    cohorts = ["TCGA-BLCA", "TCGA-SKCM"]

    all_corr = []
    all_partial = []
    all_linear = []
    all_present = []

    for cohort in cohorts:
        print(f"\n--- {cohort} ---")
        expr_path = tcga_dir / f"{cohort}_log2TPM_symbol_tumor.tsv"
        mech_path = mech_dir / f"{cohort}_mechanism_scores.tsv"
        surv_path = tcga_dir / f"{cohort}_survival_matched.tsv"

        expr = read_expression_matrix(expr_path)
        biom, present = compute_tcga_signature_scores(expr)
        present["cohort"] = cohort
        all_present.append(present)

        mech = pd.read_csv(mech_path, sep="\t")
        surv = pd.read_csv(surv_path, sep="\t")

        df = mech.merge(biom, on="sample_id", how="left", suffixes=("", "_expr"))
        merge_on = "patient_id" if "patient_id" in surv.columns else "sample_id"
        surv2 = surv.drop_duplicates(merge_on)
        df = df.merge(surv2, on=merge_on, how="left", suffixes=("", "_clin"))

        df["cohort"] = cohort
        df["score_z"] = zscore(df["score"])

        # Age
        if "age_at_diagnosis" in df.columns:
            df["age_years"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce") / 365.25
        elif "age_at_index" in df.columns:
            df["age_years"] = pd.to_numeric(df["age_at_index"], errors="coerce")
        else:
            df["age_years"] = np.nan

        # Stage
        if "ajcc_pathologic_stage" in df.columns:
            df["stage_simple"] = df["ajcc_pathologic_stage"].astype(str).str.extract(r"(Stage [IVX]+)", expand=False)
        else:
            df["stage_simple"] = np.nan

        df.to_csv(tcga_adj_dir / f"{cohort}_biomarker_adjustment_merged.tsv", sep="\t", index=False)
        print(f"Merged shape: {df.shape}")

        # Correlations
        candidate_vars = [
            "checkpoint", "antigen_presentation", "cytolytic", "exhaustion", "t_cell_inflamed",
            "IFNG_signature", "CD8_score", "cytolytic_score_expr", "checkpoint_expr_score",
            "T_cell_exhaustion_expr_score", "MHC_I_expr_score",
            "immune_surrogate", "stromal_surrogate", "myeloid_surrogate",
            "CAF_stromal_surrogate", "endothelial_surrogate", "proliferation_score",
        ] + tcga_single_biomarkers

        for v in candidate_vars:
            if v not in df.columns:
                continue
            dat = df[["score", v]].dropna()
            if dat.shape[0] < 20:
                continue
            rho, p = stats.spearmanr(dat["score"], dat[v])
            all_corr.append({"cohort": cohort, "variable": v, "n": dat.shape[0], "spearman_rho": rho, "p": p})

        # Partial correlations
        covar_sets = {
            "immune_surrogate": ["immune_surrogate"],
            "stromal_surrogate": ["stromal_surrogate"],
            "immune_plus_stromal": ["immune_surrogate", "stromal_surrogate"],
            "immune_stromal_proliferation": ["immune_surrogate", "stromal_surrogate", "proliferation_score"],
        }
        for v in candidate_vars:
            if v not in df.columns or v in ["immune_surrogate", "stromal_surrogate", "proliferation_score"]:
                continue
            for covar_name, covars in covar_sets.items():
                use_covars = [c for c in covars if c in df.columns]
                if not use_covars:
                    continue
                rho, p, n = partial_spearman(df["score"], df[v], df[use_covars])
                all_partial.append({
                    "cohort": cohort, "variable": v, "adjustment": covar_name, "n": n,
                    "partial_spearman_rho": rho, "p": p,
                })

        # Linear models
        for v in [
            "checkpoint", "antigen_presentation", "cytolytic", "exhaustion", "t_cell_inflamed",
            "IFNG_signature", "CD8_score", "checkpoint_expr_score", "MHC_I_expr_score",
            "CD274", "PDCD1", "IFNG", "CXCL9", "CXCL10", "CD8A",
        ]:
            if v not in df.columns:
                continue
            covars = ["score_z", "immune_surrogate", "stromal_surrogate", "proliferation_score"]
            dat = df[[v] + covars].dropna()
            if dat.shape[0] < 30:
                continue
            X = sm.add_constant(dat[covars], has_constant="add")
            y = dat[v]
            fit = sm.OLS(y, X).fit()
            all_linear.append({
                "cohort": cohort, "outcome": v, "n": dat.shape[0],
                "model": "outcome ~ score_z + immune_surrogate + stromal_surrogate + proliferation_score",
                "score_beta": fit.params.get("score_z", np.nan),
                "score_p": fit.pvalues.get("score_z", np.nan),
                "score_ci_low": fit.conf_int().loc["score_z", 0] if "score_z" in fit.params.index else np.nan,
                "score_ci_high": fit.conf_int().loc["score_z", 1] if "score_z" in fit.params.index else np.nan,
                "r2": fit.rsquared,
            })

    # Save
    present_df = pd.concat(all_present, ignore_index=True)
    present_df.to_csv(tcga_adj_dir / "TCGA_signature_gene_coverage.tsv", sep="\t", index=False)

    corr = pd.DataFrame(all_corr)
    if not corr.empty:
        corr["fdr"] = bh_fdr(corr["p"])
        corr = corr.sort_values(["cohort", "p"])
        corr.to_csv(tcga_adj_dir / "TCGA_score_biomarker_spearman.tsv", sep="\t", index=False)

    partial = pd.DataFrame(all_partial)
    if not partial.empty:
        partial["fdr"] = bh_fdr(partial["p"])
        partial = partial.sort_values(["cohort", "adjustment", "p"])
        partial.to_csv(tcga_adj_dir / "TCGA_score_biomarker_partial_spearman.tsv", sep="\t", index=False)

    linear = pd.DataFrame(all_linear)
    if not linear.empty:
        linear["score_fdr"] = bh_fdr(linear["score_p"])
        linear = linear.sort_values(["cohort", "score_p"])
        linear.to_csv(tcga_adj_dir / "TCGA_adjusted_linear_models.tsv", sep="\t", index=False)

    # Compact summary
    key_vars = [
        "checkpoint", "antigen_presentation", "cytolytic", "exhaustion", "t_cell_inflamed",
        "IFNG_signature", "CD8_score", "checkpoint_expr_score", "MHC_I_expr_score",
        "CD274", "PDCD1", "IFNG", "CXCL9", "CXCL10", "CD8A",
        "immune_surrogate", "stromal_surrogate", "proliferation_score",
    ]
    summary_rows = []
    for cohort in cohorts:
        for v in key_vars:
            row = {"cohort": cohort, "variable": v}
            csub = corr[(corr["cohort"] == cohort) & (corr["variable"] == v)] if not corr.empty else pd.DataFrame()
            if not csub.empty:
                row["spearman_rho"] = csub.iloc[0]["spearman_rho"]
                row["spearman_p"] = csub.iloc[0]["p"]
                row["spearman_fdr"] = csub.iloc[0]["fdr"]
            psub = partial[(partial["cohort"] == cohort) & (partial["variable"] == v) & (partial["adjustment"] == "immune_stromal_proliferation")] if not partial.empty else pd.DataFrame()
            if not psub.empty:
                row["partial_rho_adj_immune_stromal_prolif"] = psub.iloc[0]["partial_spearman_rho"]
                row["partial_p_adj_immune_stromal_prolif"] = psub.iloc[0]["p"]
                row["partial_fdr_adj_immune_stromal_prolif"] = psub.iloc[0]["fdr"]
            lsub = linear[(linear["cohort"] == cohort) & (linear["outcome"] == v)] if not linear.empty else pd.DataFrame()
            if not lsub.empty:
                row["linear_score_beta_adj"] = lsub.iloc[0]["score_beta"]
                row["linear_score_p_adj"] = lsub.iloc[0]["score_p"]
                row["linear_score_fdr_adj"] = lsub.iloc[0]["score_fdr"]
            summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tcga_adj_dir / "TCGA_adjusted_association_summary.tsv", sep="\t", index=False)

    # Plots
    if not summary.empty:
        plot_df = summary[summary["variable"].isin(key_vars[:15])].copy()
        hm = plot_df.pivot(index="variable", columns="cohort", values="spearman_rho")
        plt.figure(figsize=(7, 7))
        sns.heatmap(hm, annot=True, fmt=".2f", cmap="vlag", center=0, cbar_kws={"label": "Spearman rho"})
        plt.title("ExoImmune Score vs biomarkers/modules")
        plt.tight_layout()
        plt.savefig(tcga_adj_dir / "TCGA_score_biomarker_spearman_heatmap.pdf")
        plt.close()

        hm2 = plot_df.pivot(index="variable", columns="cohort", values="partial_rho_adj_immune_stromal_prolif")
        plt.figure(figsize=(7, 7))
        sns.heatmap(hm2, annot=True, fmt=".2f", cmap="vlag", center=0, cbar_kws={"label": "Partial Spearman rho"})
        plt.title("Adjusted association with ExoImmune Score")
        plt.tight_layout()
        plt.savefig(tcga_adj_dir / "TCGA_score_biomarker_partial_heatmap.pdf")
        plt.close()

    print("\nSignature gene coverage:")
    print(present_df.to_string(index=False))
    print(f"\n[Part 6 complete]  Saved files: {list(tcga_adj_dir.glob('*.tsv'))}")


# ============================================================
#  PART 7: TCGA COX SURVIVAL ANALYSIS  (from 11.txt, overrides 10.txt Cox)
# ============================================================

def run_part7_tcga_cox():
    print("\n" + "=" * 100)
    print("PART 7: TCGA COX SURVIVAL ANALYSIS")
    print("=" * 100)

    if not HAS_LIFELINES:
        print("lifelines not available; skipping Cox analysis.")
        return

    cohorts = ["TCGA-BLCA", "TCGA-SKCM"]
    rows = []

    for cohort in cohorts:
        print(f"\n--- {cohort} ---")
        fp = tcga_adj_dir / f"{cohort}_biomarker_adjustment_merged.tsv"
        if not fp.exists():
            print(f"File not found: {fp}")
            continue
        df = pd.read_csv(fp, sep="\t")

        df["OS.time"] = pd.to_numeric(df["OS.time"], errors="coerce")
        df["OS.event"] = pd.to_numeric(df["OS.event"], errors="coerce")
        df["score_z"] = zscore(df["score"])

        base = df[["OS.time", "OS.event", "score", "score_z"]].replace([np.inf, -np.inf], np.nan).dropna()
        base = base[base["OS.time"] > 0]
        print(f"Cox usable n: {base.shape[0]}, Events: {int(base['OS.event'].sum())}")

        # Cox score-only
        for score_var in ["score", "score_z"]:
            dat = base[["OS.time", "OS.event", score_var]].dropna()
            if dat.shape[0] < 50 or dat["OS.event"].sum() < 10:
                print(f"Skip Cox: {score_var}")
                continue
            cph = CoxPHFitter()
            try:
                cph.fit(dat, duration_col="OS.time", event_col="OS.event")
                s = cph.summary.loc[score_var]
                rows.append({
                    "cohort": cohort, "model": f"{score_var}_only", "n": dat.shape[0],
                    "events": int(dat["OS.event"].sum()), "term": score_var,
                    "coef": s["coef"], "HR": s["exp(coef)"], "p": s["p"],
                    "ci_low": s["exp(coef) lower 95%"], "ci_high": s["exp(coef) upper 95%"],
                })
                print(f"Cox {score_var}: HR={s['exp(coef)']:.3f}, p={s['p']:.4f}")
            except Exception as e:
                print(f"Cox failed: {score_var}: {e}")

        # Adjusted Cox models
        covar_sets = {
            "score_z_plus_immune_stromal": ["score_z", "immune_surrogate", "stromal_surrogate"],
            "score_z_plus_immune_stromal_prolif": ["score_z", "immune_surrogate", "stromal_surrogate", "proliferation_score"],
            "score_z_plus_modules": ["score_z", "checkpoint", "antigen_presentation", "cytolytic", "exhaustion", "t_cell_inflamed"],
        }

        for model_name, covars in covar_sets.items():
            use = ["OS.time", "OS.event"] + covars
            dat = df[use].copy()
            dat["OS.time"] = pd.to_numeric(dat["OS.time"], errors="coerce")
            dat["OS.event"] = pd.to_numeric(dat["OS.event"], errors="coerce")
            for c in covars:
                dat[c] = zscore(dat[c])
            dat = dat.replace([np.inf, -np.inf], np.nan).dropna()
            dat = dat[dat["OS.time"] > 0]

            if dat.shape[0] < 50 or dat["OS.event"].sum() < 10:
                continue

            try:
                cph = CoxPHFitter(penalizer=0.01)
                cph.fit(dat, duration_col="OS.time", event_col="OS.event")
                for term in cph.summary.index:
                    s = cph.summary.loc[term]
                    rows.append({
                        "cohort": cohort, "model": model_name, "n": dat.shape[0],
                        "events": int(dat["OS.event"].sum()), "term": term,
                        "coef": s["coef"], "HR": s["exp(coef)"], "p": s["p"],
                        "ci_low": s["exp(coef) lower 95%"], "ci_high": s["exp(coef) upper 95%"],
                    })
                print(f"Cox {model_name}: n={dat.shape[0]}, events={int(dat['OS.event'].sum())}")
            except Exception as e:
                print(f"Adjusted Cox failed: {model_name}: {e}")

        # Log-rank high vs low
        lrdat = base.copy()
        med = lrdat["score"].median()
        low = lrdat[lrdat["score"] <= med]
        high = lrdat[lrdat["score"] > med]
        if low.shape[0] > 20 and high.shape[0] > 20:
            res = logrank_test(
                high["OS.time"], low["OS.time"],
                event_observed_A=high["OS.event"], event_observed_B=low["OS.event"],
            )
            print(f"\nLog-rank high vs low: median={med:.3f}, p={res.p_value:.4f}")
            rows.append({
                "cohort": cohort, "model": "median_split_logrank", "n": lrdat.shape[0],
                "events": int(lrdat["OS.event"].sum()), "term": "high_vs_low",
                "coef": np.nan, "HR": np.nan, "p": res.p_value, "ci_low": np.nan, "ci_high": np.nan,
            })

    res_df = pd.DataFrame(rows)
    res_path = tcga_adj_dir / "TCGA_cox_debug_score_models.tsv"
    res_df.to_csv(res_path, sep="\t", index=False)

    print("\nFinal Cox debug results:")
    print(res_df.to_string(index=False))
    print(f"\n[Part 7 complete]  Saved: {res_path}")


# ============================================================
#  PART 8: IMVIGOR210 INTERNAL VALIDATION  (from 13.txt)
# ============================================================

def run_part8_imvigor_validation():
    print("\n" + "=" * 100)
    print("PART 8: IMVIGOR210 INTERNAL VALIDATION")
    print("=" * 100)

    clean_dir = root / "IMvigor210_clean"
    lasso_dir = root / "IMvigor210_lasso"

    expr_path = clean_dir / "IMvigor210_response_expr.tsv"
    response_path = clean_dir / "IMvigor210_response_df.tsv"
    score_path = lasso_dir / "exoimmune_score.tsv"

    np.random.seed(20260628)

    # Load
    expr = read_expr(expr_path)
    response = read_response(response_path)
    score = read_score(score_path)
    biom, coverage = compute_imvigor_biomarkers(expr)

    df = response.merge(score, on="sample_id", how="inner").merge(biom, on="sample_id", how="left")
    df.to_csv(imvigor_dir / "IMvigor210_response_score_biomarker_merged.tsv", sep="\t", index=False)
    coverage.to_csv(imvigor_dir / "IMvigor210_biomarker_gene_coverage.tsv", sep="\t", index=False)

    print(f"Expression: {expr.shape}, Response: {response.shape}, Score: {score.shape}")
    print(f"Merged: {df.shape}, Response counts: {df['response'].value_counts().to_dict()}")
    print("\nBiomarker gene coverage:")
    print(coverage.to_string(index=False))

    # ExoImmune robustness
    base = df[["sample_id", "response", "y", "ExoImmune_Score"]].dropna()
    y = base["y"]
    score_vec = base["ExoImmune_Score"]

    auc = roc_auc_score(y, score_vec)
    boot_mean, boot_low, boot_high, boot_n = bootstrap_auc(y, score_vec, n_boot=2000)
    cv_mean, cv_low, cv_high, cv_n = repeated_cv_auc(base[["ExoImmune_Score"]], y, n_splits=5, n_repeats=100)

    fpr, tpr, thr = roc_curve(y, score_vec)
    youden_idx = np.argmax(tpr - fpr)
    youden_cutoff = thr[youden_idx]
    pred_youden = (score_vec >= youden_cutoff).astype(int)
    cm = confusion_matrix(y, pred_youden)

    robust = pd.DataFrame([{
        "model": "ExoImmune_Score", "n": base.shape[0], "responders": int(y.sum()),
        "nonresponders": int((1 - y).sum()), "apparent_auc": auc,
        "bootstrap_auc_mean": boot_mean, "bootstrap_auc_low95": boot_low,
        "bootstrap_auc_high95": boot_high, "bootstrap_n": boot_n,
        "repeated_5fold_cv_auc_mean": cv_mean, "repeated_5fold_cv_auc_low95": cv_low,
        "repeated_5fold_cv_auc_high95": cv_high, "cv_n": cv_n,
        "youden_cutoff": youden_cutoff, "youden_sensitivity": tpr[youden_idx],
        "youden_specificity": 1 - fpr[youden_idx],
        "tn": cm[0, 0], "fp": cm[0, 1], "fn": cm[1, 0], "tp": cm[1, 1],
    }])
    robust.to_csv(imvigor_dir / "IMvigor210_exoimmune_internal_robustness.tsv", sep="\t", index=False)
    print("\nExoImmune internal robustness:")
    print(robust.to_string(index=False))

    # Biomarker comparison
    biomarker_vars = list(imvigor_gene_sets.keys()) + imvigor_single_genes
    comparison_rows = []
    pred_store = {}

    exo_res = fit_logistic_auc(df, ["ExoImmune_Score"])
    if exo_res:
        comparison_rows.append({
            "model": "ExoImmune_Score", "model_type": "exoimmune_only", "biomarker": "",
            "n": exo_res["n"], "responders": exo_res["responders"],
            "apparent_auc": exo_res["apparent_auc"], "cv_auc_mean": exo_res["cv_auc_mean"],
            "cv_auc_low95": exo_res["cv_auc_low"], "cv_auc_high95": exo_res["cv_auc_high"],
        })
        pred_store["ExoImmune_Score"] = exo_res

    for v in biomarker_vars:
        if v not in df.columns:
            continue
        res1 = fit_logistic_auc(df, [v])
        if res1:
            comparison_rows.append({
                "model": v, "model_type": "biomarker_only", "biomarker": v,
                "n": res1["n"], "responders": res1["responders"],
                "apparent_auc": res1["apparent_auc"], "cv_auc_mean": res1["cv_auc_mean"],
                "cv_auc_low95": res1["cv_auc_low"], "cv_auc_high95": res1["cv_auc_high"],
            })
            pred_store[v] = res1

        res2 = fit_logistic_auc(df, ["ExoImmune_Score", v])
        if res2:
            comparison_rows.append({
                "model": f"ExoImmune_Score+{v}", "model_type": "exoimmune_plus_biomarker", "biomarker": v,
                "n": res2["n"], "responders": res2["responders"],
                "apparent_auc": res2["apparent_auc"], "cv_auc_mean": res2["cv_auc_mean"],
                "cv_auc_low95": res2["cv_auc_low"], "cv_auc_high95": res2["cv_auc_high"],
            })
            pred_store[f"ExoImmune_Score+{v}"] = res2

    comparison = pd.DataFrame(comparison_rows)
    comparison["contains_lasso_component"] = comparison["biomarker"].apply(
        lambda x: any(g in lasso_component_genes for g in imvigor_gene_sets.get(x, [x])) if x else False
    )
    comparison["note"] = ""
    comparison.loc[comparison["biomarker"].eq("IFNG"), "note"] = "IFNG is a LASSO score component; do not interpret as independent comparison."
    comparison.loc[comparison["biomarker"].eq("IFNG_signature"), "note"] = "Signature contains IFNG, a score component."
    comparison = comparison.sort_values(["model_type", "cv_auc_mean"], ascending=[True, False])
    comparison.to_csv(imvigor_dir / "IMvigor210_biomarker_model_comparison_auc.tsv", sep="\t", index=False)
    print("\nBiomarker model comparison:")
    print(comparison.to_string(index=False))

    # Univariate association
    assoc_rows = []
    for v in ["ExoImmune_Score"] + biomarker_vars:
        if v not in df.columns:
            continue
        dat = df[["y", v]].dropna().copy()
        if dat.shape[0] < 30 or dat["y"].nunique() < 2:
            continue
        dat[v + "_z"] = zscore_series(dat[v])
        X = dat[[v + "_z"]].values
        yv = dat["y"].values
        clf = LogisticRegression(solver="liblinear")
        clf.fit(X, yv)
        coef = clf.coef_[0][0]
        pred = clf.predict_proba(X)[:, 1]
        auc_v = roc_auc_score(yv, pred)
        try:
            X2 = sm.add_constant(dat[[v + "_z"]])
            fit = sm.Logit(dat["y"], X2).fit(disp=False)
            beta, se, p = fit.params[v + "_z"], fit.bse[v + "_z"], fit.pvalues[v + "_z"]
            ci_low, ci_high = beta - 1.96 * se, beta + 1.96 * se
        except Exception:
            beta, p, ci_low, ci_high = coef, np.nan, np.nan, np.nan

        responder_mean = dat.loc[dat["y"] == 1, v].mean()
        nonresponder_mean = dat.loc[dat["y"] == 0, v].mean()
        try:
            wt_p = stats.mannwhitneyu(dat.loc[dat["y"] == 1, v], dat.loc[dat["y"] == 0, v], alternative="two-sided").pvalue
        except Exception:
            wt_p = np.nan

        assoc_rows.append({
            "variable": v, "n": dat.shape[0], "responders": int(dat["y"].sum()),
            "OR_per_1SD": np.exp(beta), "OR_low95": np.exp(ci_low), "OR_high95": np.exp(ci_high),
            "wald_p": p, "auc": auc_v, "responder_mean": responder_mean,
            "nonresponder_mean": nonresponder_mean, "wilcoxon_p": wt_p,
            "contains_lasso_component": any(g in lasso_component_genes for g in imvigor_gene_sets.get(v, [v])) if v != "ExoImmune_Score" else True,
        })

    assoc = pd.DataFrame(assoc_rows)
    assoc["wald_fdr"] = bh_fdr(assoc["wald_p"])
    assoc["wilcoxon_fdr"] = bh_fdr(assoc["wilcoxon_p"])
    assoc = assoc.sort_values("auc", ascending=False)
    assoc.to_csv(imvigor_dir / "IMvigor210_response_univariate_biomarker_association.tsv", sep="\t", index=False)
    print("\nUnivariate response association:")
    print(assoc.to_string(index=False))

    # Calibration and DCA
    dat = df[["y", "ExoImmune_Score"]].dropna().copy()
    X = dat[["ExoImmune_Score"]]
    y = dat["y"]
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(solver="liblinear", penalty="l2", max_iter=1000))
    ])
    pred_cv = cross_val_predict(pipe, X, y, cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=20260628), method="predict_proba")[:, 1]

    cal_df = pd.DataFrame({"y": y.values, "pred_cv": pred_cv})
    cal_df["bin"] = pd.qcut(cal_df["pred_cv"], q=5, duplicates="drop")
    cal_tab = cal_df.groupby("bin").agg(n=("y", "size"), mean_pred=("pred_cv", "mean"), observed_response_rate=("y", "mean")).reset_index()
    cal_tab["bin"] = cal_tab["bin"].astype(str)
    cal_tab.to_csv(imvigor_dir / "IMvigor210_exoimmune_calibration_table.tsv", sep="\t", index=False)

    brier = brier_score_loss(y, pred_cv)
    print(f"\nCalibration Brier score: {brier:.4f}")
    print(cal_tab.to_string(index=False))

    thresholds = np.arange(0.05, 0.81, 0.01)
    dca = decision_curve(y, pred_cv, thresholds)
    dca.to_csv(imvigor_dir / "IMvigor210_exoimmune_decision_curve.tsv", sep="\t", index=False)

    # Plots
    plt.figure(figsize=(6, 5))
    for name in ["ExoImmune_Score", "IFNG_signature", "CD8_score", "checkpoint_score", "T_cell_inflamed_score"]:
        if name not in pred_store:
            continue
        yy = pred_store[name]["y"]
        pp = pred_store[name]["pred"]
        fpr, tpr, _ = roc_curve(yy, pp)
        plt.plot(fpr, tpr, lw=2, label=f"{name} AUC={roc_auc_score(yy, pp):.3f}")
    plt.plot([0, 1], [0, 1], "--", color="grey", lw=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("IMvigor210 response prediction")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(imvigor_dir / "IMvigor210_biomarker_roc_comparison.pdf")
    plt.close()

    plot_cmp = comparison[comparison["model_type"].isin(["exoimmune_only", "biomarker_only"])].copy()
    plot_cmp = plot_cmp.sort_values("cv_auc_mean", ascending=False).head(20)
    plt.figure(figsize=(8, max(4, 0.28 * len(plot_cmp))))
    sns.barplot(data=plot_cmp, y="model", x="cv_auc_mean", hue="model_type", dodge=False)
    plt.xlim(0.4, 1.0)
    plt.xlabel("Repeated 5-fold CV AUC")
    plt.ylabel("")
    plt.title("IMvigor210 biomarker comparison")
    plt.tight_layout()
    plt.savefig(imvigor_dir / "IMvigor210_biomarker_cv_auc_barplot.pdf")
    plt.close()

    plt.figure(figsize=(5, 5))
    plt.plot([0, 1], [0, 1], "--", color="grey")
    plt.plot(cal_tab["mean_pred"], cal_tab["observed_response_rate"], marker="o")
    plt.xlabel("Mean predicted probability")
    plt.ylabel("Observed response rate")
    plt.title(f"ExoImmune calibration, Brier={brier:.3f}")
    plt.tight_layout()
    plt.savefig(imvigor_dir / "IMvigor210_exoimmune_calibration_plot.pdf")
    plt.close()

    plt.figure(figsize=(6, 5))
    plt.plot(dca["threshold"], dca["net_benefit_model"], label="ExoImmune Score")
    plt.plot(dca["threshold"], dca["net_benefit_all"], "--", label="Treat all")
    plt.plot(dca["threshold"], dca["net_benefit_none"], "--", label="Treat none")
    plt.xlabel("Threshold probability")
    plt.ylabel("Net benefit")
    plt.title("Decision curve analysis")
    plt.legend()
    plt.tight_layout()
    plt.savefig(imvigor_dir / "IMvigor210_exoimmune_decision_curve.pdf")
    plt.close()

    plt.figure(figsize=(5, 4))
    sns.boxplot(data=df, x="response", y="ExoImmune_Score", order=["NonResponder", "Responder"])
    sns.stripplot(data=df, x="response", y="ExoImmune_Score", order=["NonResponder", "Responder"], color="black", alpha=0.35, size=2)
    plt.xlabel("")
    plt.ylabel("ExoImmune Score")
    plt.title("ExoImmune Score by response")
    plt.tight_layout()
    plt.savefig(imvigor_dir / "IMvigor210_exoimmune_response_boxplot.pdf")
    plt.close()

    print(f"\n[Part 8 complete]  Saved files: {list(imvigor_dir.glob('*'))}")


# ============================================================
#  PART 9: PROJECT INTEGRATED SUMMARY  (from 12.txt)
# ============================================================

def run_part9_integrated_summary():
    print("\n" + "=" * 100)
    print("PART 9: PROJECT INTEGRATED SUMMARY")
    print("=" * 100)

    # Load results
    source = pd.read_csv(sc_source_dir / "singlecell_lasso_gene_consensus.tsv", sep="\t")
    lr_narrative = pd.read_csv(sc_comm_dir / "targeted_lr_mechanism_narrative.tsv", sep="\t")
    lr_edges = pd.read_csv(sc_comm_dir / "targeted_lr_high_confidence_edges.tsv", sep="\t")
    assoc = pd.read_csv(tcga_adj_dir / "TCGA_adjusted_association_summary.tsv", sep="\t")
    cox = pd.read_csv(tcga_adj_dir / "TCGA_cox_debug_score_models.tsv", sep="\t")

    # 1) Biomarker table with leakage flag
    assoc["biomarker_flag"] = assoc["variable"].map(biomarker_component_flag).fillna("other")
    assoc["preferred_for_main_text"] = assoc["biomarker_flag"].isin([
        "independent_module", "independent_signature", "independent_single_gene", "contains_score_component_IFNG",
    ])
    assoc["interpretation_note"] = ""
    assoc.loc[assoc["variable"].eq("IFNG"), "interpretation_note"] = "Do not use as independent biomarker evidence because IFNG is part of the LASSO score."
    assoc.loc[assoc["variable"].eq("IFNG_signature"), "interpretation_note"] = "Useful immune-inflamed signature, but note it contains IFNG, a score component."
    assoc.loc[assoc["variable"].isin(["immune_surrogate", "stromal_surrogate", "proliferation_score"]), "interpretation_note"] = "Used as adjustment/surrogate covariate, not primary biomarker comparison."
    assoc.to_csv(summary_dir / "integrated_TCGA_biomarker_adjusted_summary.tsv", sep="\t", index=False)

    # 2) Compact key association table
    key_vars = [
        "checkpoint", "antigen_presentation", "cytolytic", "exhaustion", "t_cell_inflamed",
        "IFNG_signature", "CD8_score", "checkpoint_expr_score", "MHC_I_expr_score",
        "CD274", "PDCD1", "CD8A", "CXCL9", "CXCL10", "IFNG",
    ]
    key_assoc = assoc[assoc["variable"].isin(key_vars)].copy()
    key_assoc = key_assoc[[
        "cohort", "variable", "biomarker_flag", "spearman_rho", "spearman_p",
        "partial_rho_adj_immune_stromal_prolif", "partial_p_adj_immune_stromal_prolif",
        "linear_score_beta_adj", "linear_score_p_adj", "interpretation_note",
    ]]
    key_assoc = key_assoc.sort_values(["cohort", "variable"])
    key_assoc.to_csv(summary_dir / "main_text_TCGA_key_biomarker_table.tsv", sep="\t", index=False)

    # 3) Cox interpretation table
    score_cox = cox[cox["term"].isin(["score_z", "high_vs_low"])].copy()

    def cox_context(row):
        model = row["model"]
        if model == "score_z_only":
            return "univariate_continuous_score"
        if model == "median_split_logrank":
            return "median_split_logrank"
        if "immune_stromal" in model:
            return "adjusted_for_immune_stromal_surrogates"
        if "modules" in model:
            return "adjusted_for_immune_modules"
        return model

    def cox_note(row):
        if row["term"] == "high_vs_low":
            return "High vs low score survival separation is significant by log-rank." if row["p"] < 0.05 else "High vs low score survival separation is not significant by log-rank."
        if row["term"] == "score_z":
            if row["p"] < 0.05 and row["HR"] < 1:
                return "Higher score is associated with lower hazard in this model."
            if row["p"] < 0.05 and row["HR"] > 1:
                return "Higher score is associated with higher hazard in this model."
            return "Score is not independently significant in this model."
        return ""

    score_cox["context"] = score_cox.apply(cox_context, axis=1)
    score_cox["interpretation_note"] = score_cox.apply(cox_note, axis=1)
    score_cox.to_csv(summary_dir / "integrated_TCGA_score_survival_interpretation.tsv", sep="\t", index=False)

    # 4) Source + communication tables
    source.to_csv(summary_dir / "integrated_singlecell_lasso_source_summary.tsv", sep="\t", index=False)

    edge_keep = lr_edges[[
        "dataset", "axis", "ligand", "receptor", "source_cell", "target_cell",
        "ligand_det_rate", "receptor_det_rate", "evidence_score", "support_level", "evidence_tier",
    ]].copy().sort_values(["axis", "ligand", "receptor", "dataset"])
    edge_keep.to_csv(summary_dir / "integrated_singlecell_priority_communication_edges.tsv", sep="\t", index=False)
    lr_narrative.to_csv(summary_dir / "integrated_singlecell_communication_narrative.tsv", sep="\t", index=False)

    # 5) Human-readable Chinese summary
    def get_assoc(cohort, var):
        sub = assoc[(assoc["cohort"] == cohort) & (assoc["variable"] == var)]
        return sub.iloc[0] if not sub.empty else None

    lines = []
    lines.append("ExoImmune project integrated summary")
    lines.append("=" * 80)
    lines.append("")
    lines.append("1. Core interpretation")
    lines.append("- ExoImmune Score is best interpreted as an immunotherapy response-associated exosome-related immune state score, not as a robust independent OS biomarker.")
    lines.append("- The score tracks an inflamed, checkpoint-high, antigen-presentation/cytolytic and exhaustion-associated tumor microenvironment.")
    lines.append("- Single-cell evidence supports a multicellular origin: T/T_NK-derived IFNG, malignant-cell NOMO3, tumor/stromal proliferative components, and endothelial/CAF/stromal contributions.")
    lines.append("")
    lines.append("2. Single-cell source evidence")
    for _, r in source.iterrows():
        lines.append(f"- {r['gene']}: coef={float(r['coef']):+.3f}, sign={r['coef_sign']}, source={r['celltype_set']}, category_overlap={r['category_overlap']}")

    lines.append("")
    lines.append("3. Targeted communication evidence")
    for _, r in lr_narrative.iterrows():
        lines.append(f"- {r['axis']}: best={r['best_pair']} on {r['best_edge']}, datasets={r['datasets_supporting']}. {r['interpretation']}")

    lines.append("")
    lines.append("4. TCGA adjusted biomarker evidence")
    for cohort in ["TCGA-BLCA", "TCGA-SKCM"]:
        lines.append(f"- {cohort}:")
        for var in ["t_cell_inflamed", "IFNG_signature", "CD8_score", "checkpoint", "checkpoint_expr_score", "PDCD1", "CD274", "CD8A", "MHC_I_expr_score"]:
            r = get_assoc(cohort, var)
            if r is None:
                continue
            lines.append(f"  {var}: rho={r['spearman_rho']:.3f}, p={fmt_p(r['spearman_p'])}; partial_rho={r['partial_rho_adj_immune_stromal_prolif']:.3f}, partial_p={fmt_p(r['partial_p_adj_immune_stromal_prolif'])}; flag={r['biomarker_flag']}")
        r_ifng = get_assoc(cohort, "IFNG")
        if r_ifng is not None:
            lines.append(f"  IFNG single gene: rho={r_ifng['spearman_rho']:.3f}, but should be treated as a score component rather than an independent biomarker.")

    lines.append("")
    lines.append("5. TCGA survival interpretation")
    for cohort in ["TCGA-BLCA", "TCGA-SKCM"]:
        sub = score_cox[score_cox["cohort"] == cohort]
        lines.append(f"- {cohort}:")
        for _, r in sub.iterrows():
            if r["term"] == "score_z":
                lines.append(f"  {r['model']}: HR={r['HR']:.3f}, p={fmt_p(r['p'])}. {r['interpretation_note']}")
            elif r["term"] == "high_vs_low":
                lines.append(f"  median split log-rank: p={fmt_p(r['p'])}. {r['interpretation_note']}")

    lines.append("")
    lines.append("6. Suggested manuscript wording")
    lines.append("The ExoImmune Score showed strong positive associations with checkpoint, cytolytic, antigen-presentation, exhaustion and T cell-inflamed programs in both TCGA-BLCA and TCGA-SKCM. These associations remained partly preserved after adjustment for immune, stromal and proliferation surrogate scores, particularly for T cell-inflamed, CD8-related and checkpoint-related signals. Because IFNG is one of the score genes, single-gene IFNG correlations were considered component-driven and were not interpreted as independent validation.")
    lines.append("Single-cell analysis further indicated that the score reflects a multicellular exosome-associated program rather than a single lineage marker. IFNG localized consistently to T/T_NK cells, NOMO3 to malignant cells, and UHRF1/MYBL2 to tumor/stromal CD45-negative or proliferative compartments. Targeted ligand-receptor analysis supported IFNG-IFNGR, CXCL9/10/11-CXCR3, PD-1/PD-L1, exhaustion checkpoint and MHC-I-CD8 communication axes.")
    lines.append("In survival analyses, higher score was associated with lower hazard in univariate Cox models, especially in TCGA-SKCM. However, this association was largely attenuated after adjustment for immune/stromal surrogates or immune modules, supporting the interpretation that ExoImmune Score mainly captures immune microenvironment state and immunotherapy response potential rather than an independent prognostic biomarker.")

    summary_txt = summary_dir / "integrated_project_summary_notes.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # Print
    print("\nMain TCGA key biomarker table:")
    print(key_assoc.to_string(index=False))
    print("\nScore survival interpretation:")
    print(score_cox[["cohort", "model", "term", "n", "events", "HR", "p", "context", "interpretation_note"]].to_string(index=False))
    print(f"\n[Part 9 complete]  Saved files: {list(summary_dir.glob('*'))}")


# ============================================================
#  MAIN: Run all parts in order
# ============================================================

if __name__ == "__main__":
    print("=" * 100)
    print("ExoImmune Complete Analysis Pipeline")
    print("=" * 100)
    print(f"Root: {root}")
    print(f"Python: {os.sys.version}")
    print(f"Start time: {pd.Timestamp.now()}")

    run_part1_exploration()      # from 1,2,3,4,9
    run_part2_source_analysis()  # from 5
    run_part3_source_summary()   # from 6
    run_part4_communication()    # from 7
    run_part5_comm_edges()       # from 8
    run_part6_tcga_biomarker()   # from 10
    run_part7_tcga_cox()         # from 11
    run_part8_imvigor_validation()  # from 13
    run_part9_integrated_summary()  # from 12

    print("\n" + "=" * 100)
    print("ALL PARTS COMPLETE")
    print("=" * 100)
    print(f"End time: {pd.Timestamp.now()}")
    print("\nOutput directories:")
    for d in [sc_source_dir, sc_comm_dir, tcga_adj_dir, imvigor_dir, summary_dir]:
        files = list(d.glob("*"))
        print(f"  {d.name}: {len(files)} files")
