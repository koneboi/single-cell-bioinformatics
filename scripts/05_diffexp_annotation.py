"""
05_diffexp_annotation.py
========================

Annotate Leiden clusters to cell types using a curated marker-gene mapping and
run a differential-expression analysis between cell types (or clusters).

Steps:
    * Map each Leiden cluster to a cell type by scoring cluster marker genes
      against a curated panel (automated manual-style annotation).
    * Store the annotation in `obs['cell_type_pred']` and `obs['annotation']`.
    * Run differential expression (Wilcoxon) between e.g. T cells vs B cells,
      or in general the two largest clusters, to demonstrate `rank_genes_groups`.
    * Produce a volcano plot and export the DE table.

Outputs (output/):
    * annotated.h5ad            - AnnData with annotation + DE results stored
    * umap_celltype_labels.png  - UMAP colored by predicted cell type
    * volcano_DE.png            - volcano plot of the chosen contrast
    * de_table.png              - bar chart of top DE genes (log2FC)
    * DE_<A>_vs_<B>.csv         - full DE table
"""

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from scipy import stats

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
OUT_DIR = "output"
CLUSTERED_FILE = os.path.join(OUT_DIR, "leiden_clusters.h5ad")

# Curated marker panel used to annotate clusters (borrowed from 01).
MARKER_PANEL = {
    "T cells": ["CD3D", "CD3E", "CD3G", "CD2", "IL7R", "TRAC"],
    "B cells": ["CD19", "CD79A", "MS4A1", "CD22", "PAX5"],
    "Monocytes": ["CD14", "FCGR3A", "LYZ", "CSF1R", "S100A8"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "GZMB", "PRF1"],
    "Dendritic cells": ["CD1C", "FCER1A", "LILRA4", "CLEC9A", "ITGAX"],
    "Neutrophils": ["FCGR3B", "S100A12", "CSF3R", "CXCL8", "CEACAM8"],
}

# Differential-expression contrast: which pair to compare.
# If these cell types exist, use them; otherwise default to the two largest.
CONTRAST_A = "T cells"
CONTRAST_B = "B cells"

N_TOP_DE = 20
UPREG_COLOR = "#DD8452"
DOWNREG_COLOR = "#4C72B0"
NS_COLOR = "#BBBBBB"


def score_clusters(adata):
    """
    Annotate each Leiden cluster by scoring the mean expression of curated
    marker genes per cell type. For each cluster, the cell type whose marker
    genes are most highly expressed (relative to the other types) wins.
    """
    categories = adata.obs["leiden"].cat.categories
    # Only consider genes actually present in the (HVG-reduced) data.
    available = set(adata.var_names)
    type_genes = {
        ct: [m for m in markers if m in available]
        for ct, markers in MARKER_PANEL.items()
    }

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    var_index = {g: i for i, g in enumerate(adata.var_names)}

    # Baseline: global mean expression of each cell type's markers.
    global_scores = {}
    for ct, genes in type_genes.items():
        idxs = [var_index[g] for g in genes if g in var_index]
        global_scores[ct] = X[:, idxs].mean() if idxs else 0.0

    assigned = {}
    for cluster in categories:
        m = adata.obs["leiden"] == cluster
        scores = {}
        for ct, genes in type_genes.items():
            idxs = [var_index[g] for g in genes if g in var_index]
            # Mean expression in this cluster.
            scores[ct] = X[m, :][:, idxs].mean() if len(idxs) else 0.0
        # Competition score: how much this cluster over-expresses each type's
        # markers relative to the dataset-wide mean.
        comp = {ct: (scores[ct] - global_scores[ct]) for ct in scores}
        best_type = max(comp, key=comp.get)
        assigned[str(cluster)] = (best_type, round(comp[best_type], 2))
        print(f"  * cluster {cluster:<3} -> {best_type:<16} "
              f"(enrichment={assigned[str(cluster)][1]})")
    return assigned


def build_volcano(adata, groupA, groupB, n_top=20):
    """Run Wilcoxon DE between groupA and groupB; return DataFrame + plot-ready."""
    sub = adata[adata.obs["annotation"].isin([groupA, groupB])].copy()
    a_mask = sub.obs["annotation"] == groupA
    b_mask = sub.obs["annotation"] == groupB

    results = []
    for gene in sub.var_names:
        expA = np.asarray(sub.X[a_mask, :].toarray() if hasattr(sub.X, "toarray") else sub.X[a_mask, :])
        expB = np.asarray(sub.X[b_mask, :].toarray() if hasattr(sub.X, "toarray") else sub.X[b_mask, :])
        expA = expA[:, sub.var_names.get_loc(gene)]
        expB = expB[:, sub.var_names.get_loc(gene)]
        # Log2 fold change of means.
        meanA = expA.mean()
        meanB = expB.mean()
        fc = np.log2((meanA + 1e-6) / (meanB + 1e-6))
        # Mann-Whitney U (equivalent to Wilcoxon rank-sum).
        u, p = stats.mannwhitneyu(expA, expB, alternative="two-sided")
        results.append({"gene": gene, "log2FC": fc, "meanA": meanA, "meanB": meanB, "pval": p})

    df = pd.DataFrame(results)
    df["-log10p"] = -np.log10(df["pval"].clip(lower=1e-300))
    # Multiple-testing correction (Benjamini-Hochberg).
    n = len(df)
    df = df.sort_values("pval")
    df["qval"] = (df["pval"].cummax()) * n / np.arange(1, n + 1)
    df["qval"] = df["qval"].clip(upper=1.0)
    df = df.sort_index()

    sig_a = (df["qval"] < 0.05) & (df["log2FC"] > 0.5)
    sig_b = (df["qval"] < 0.05) & (df["log2FC"] < -0.5)
    df["sig_type"] = np.where(sig_a, "up", np.where(sig_b, "down", "ns"))

    # Rank by a combined score for the top-DE list.
    df["score"] = df["-log10p"] * np.abs(df["log2FC"])
    df = df.sort_values("score", ascending=False)
    top = df.head(n_top)
    return df, top, groupA, groupB, sub


def main():
    print("=" * 70)
    print(" 05_diffexp_annotation.py — Cell-type annotation & DE")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"\n[1/5] Loading {CLUSTERED_FILE} ...")
    adata = sc.read_h5ad(CLUSTERED_FILE)
    print(f"  * loaded: {adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs['leiden'].nunique()} clusters")

    print("\n[2/5] Annotating clusters using the marker panel ...")
    mapping = score_clusters(adata)
    adata.obs["cell_type_pred"] = [
        mapping[str(c)][0] for c in adata.obs["leiden"]
    ]
    adata.obs["cell_type_pred"] = adata.obs["cell_type_pred"].astype("category")
    # Alias annotation column for downstream use.
    adata.obs["annotation"] = adata.obs["cell_type_pred"].astype("category")

    # Sanity: cross-tab of predicted vs ground truth.
    ct = pd.crosstab(adata.obs["cell_type"], adata.obs["annotation"])
    print("\n  Cross-tab (rows=truth, cols=predicted):")
    print("  " + ct.to_string().replace("\n", "\n  "))

    print("\n[3/5] Choosing DE contrast ...")
    present = set(adata.obs["annotation"].cat.categories)
    if CONTRAST_A in present and CONTRAST_B in present:
        groupA, groupB = CONTRAST_A, CONTRAST_B
        print(f"  * using fixed contrast {groupA} vs {groupB}")
    else:
        sizes = adata.obs["annotation"].value_counts()
        groupA, groupB = str(sizes.index[0]), str(sizes.index[1])
        print(f"  * falling back to largest groups {groupA} vs {groupB}")

    print(f"\n[4/5] Running DE ({groupA} vs {groupB}, Wilcoxon/Mann-Whitney) ...")
    df_res, top, groupA, groupB, _sub = build_volcano(
        adata, groupA, groupB, n_top=N_TOP_DE
    )
    n_up = int((df_res["sig_type"] == "up").sum())
    n_down = int((df_res["sig_type"] == "down").sum())
    print(f"  * significant up (q<0.05, FC>0.5): {n_up}")
    print(f"  * significant down             : {n_down}")

    # Store raw DE stats back onto the AnnData (optional, for reuse).
    adata.uns["DE_contrast"] = [groupA, groupB]
    out_de = os.path.join(OUT_DIR, f"DE_{groupA.replace(' ', '_')}_vs_{groupB.replace(' ', '_')}.csv")
    df_res.to_csv(out_de, index=False)
    print(f"  * saved -> {out_de}")

    # Export the annotated AnnData.
    adata.write(os.path.join(OUT_DIR, "annotated.h5ad"))
    print("  * saved -> output/annotated.h5ad")

    print("\nPlotting ...")
    # ---- Plot 1: UMAP colored by predicted cell type ----------------------- #
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(
        adata,
        color="annotation",
        palette="tab20",
        ax=ax,
        frameon=False,
        show=False,
        title="UMAP — predicted cell types",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "umap_celltype_labels.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 2: volcano plot ---------------------------------------------- #
    fig, ax = plt.subplots(figsize=(7.5, 6))
    colors = df_res["sig_type"].map(
        {"up": UPREG_COLOR, "down": DOWNREG_COLOR, "ns": NS_COLOR}
    )
    ax.scatter(df_res["log2FC"], df_res["-log10p"], c=colors, s=6, alpha=0.6)
    ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=1)
    ax.axvline(0.5, color="grey", ls="--", lw=1)
    ax.axvline(-0.5, color="grey", ls="--", lw=1)
    # Label the top DE genes.
    for _, row in top.head(12).iterrows():
        if abs(row["log2FC"]) < 0.5:
            continue
        ax.annotate(
            row["gene"], (row["log2FC"], row["-log10p"]),
            fontsize=7, ha="center", va="bottom",
        )
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Volcano — {groupA} vs {groupB}")
    from matplotlib.lines import Line2D
    legend_el = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=UPREG_COLOR,
               markersize=8, label=f"up in {groupA}"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=DOWNREG_COLOR,
               markersize=8, label=f"up in {groupB}"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=NS_COLOR,
               markersize=8, label="not significant"),
    ]
    ax.legend(handles=legend_el, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "volcano_DE.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 3: top DE genes bar chart ------------------------------------ #
    top_pos = top[top["log2FC"] > 0].head(10)
    top_neg = top[top["log2FC"] < 0].head(10)
    top_vis = pd.concat([top_pos, top_neg])
    top_vis = top_vis.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 6))
    colors_bar = [UPREG_COLOR if v > 0 else DOWNREG_COLOR for v in top_vis["log2FC"]]
    ax.barh(top_vis["gene"], top_vis["log2FC"], color=colors_bar, edgecolor="white")
    ax.axvline(0, color="grey", lw=1)
    ax.set_xlabel("log2 fold change")
    ax.set_title(f"Top DE genes — {groupA} vs {groupB}")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "de_table.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 4: mean expression of the top DE genes across cell types ----- #
    top_genes = top["gene"].head(12).tolist()
    sc.pl.heatmap(
        adata,
        var_names=top_genes,
        groupby="annotation",
        use_raw=False,
        cmap="viridis",
        show=False,
        figsize=(10, 4.5),
    )
    fig = plt.gcf()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "de_heatmap.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("  * saved 4 figures to output/")
    print("\nDone. Run 06_trajectory_pseudotime.py next.\n")


if __name__ == "__main__":
    main()
