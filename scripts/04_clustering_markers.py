"""
04_clustering_markers.py
========================

Cluster cells using a kNN graph + Leiden algorithm, embed with UMAP, compare
clusters to the ground truth, and rank marker genes per cluster (Wilcoxon).

Steps:
    * Build the kNN graph on the PCA embedding.
    * Run Leiden clustering over a range of resolutions (0.1 - 1.0) and select
      the resolution that best matches the ground truth (highest ARI / NMI).
    * Compute a UMAP embedding.
    * Evaluate agreement between clusters and ground-truth cell types with
      Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI).
    * Rank marker genes per cluster with a Wilcoxon rank-sum test and export
      a top-marker table.

Outputs (output/):
    * leiden_clusters.h5ad - AnnData with clustering + UMAP
    * umap_clusters.png      - UMAP colored by Leiden cluster
    * umap_resolution_compare.png - UMAP for a few resolutions / ARI-NMI bar
    * marker_heatmap.png     - top marker expression heatmap by cluster
    * marker_dotplot.png     - dotplot of top markers
    * marker_genes.csv       - top markers per cluster
"""

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
OUT_DIR = "output"
NORMALIZED_FILE = os.path.join(OUT_DIR, "normalized.h5ad")

RESOLUTIONS = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
N_NEIGHBORS = 15
N_PCS = 30
N_TOP_MARKERS = 10    # markers per cluster (ranked)
SEED = 42


def choose_resolution(adata, resolutions):
    """Run Leiden at several resolutions, return best (ARI), data, and report."""
    true_labels = adata.obs["cell_type"].astype(str).values
    score_rows = []
    best_score, best_res, best_adata = -1.0, None, None

    for res in resolutions:
        tmp = adata.copy()
        res_key = f"leiden_{res}"
        sc.tl.leiden(tmp, resolution=res, key_added=res_key, random_state=SEED)
        pred = tmp.obs[res_key].astype(int).values
        ari = adjusted_rand_score(true_labels, pred)
        nmi = normalized_mutual_info_score(true_labels, pred)
        score_rows.append({"resolution": res, "n_clusters": pred.max() + 1, "ARI": ari, "NMI": nmi})
        print(f"  * res={res:<4} n_clusters={pred.max() + 1:<2} ARI={ari:.3f} NMI={nmi:.3f}")
        if ari > best_score:
            best_score, best_res, best_adata = ari, res, tmp
    report = pd.DataFrame(score_rows)
    return best_adata, best_res, report


def main():
    print("=" * 70)
    print(" 04_clustering_markers.py — Leiden clustering, UMAP & markers")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"\n[1/6] Loading {NORMALIZED_FILE} ...")
    adata = sc.read_h5ad(NORMALIZED_FILE)
    print(f"  * loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"  * using X_pca with {adata.obsm['X_pca'].shape[1]} PCs")

    print(f"\n[2/6] Building kNN graph (n_neighbors={N_NEIGHBORS}) ...")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)
    print("  * graph built.")

    print("\n[3/6] Leiden clustering across resolutions ...")
    adata, best_res, report = choose_resolution(adata, RESOLUTIONS)
    adata.obs["leiden"] = adata.obs[f"leiden_{best_res}"]
    n_clusters = adata.obs["leiden"].nunique()
    print(f"  * selected resolution={best_res} with {n_clusters} clusters")
    report.to_csv(os.path.join(OUT_DIR, "resolution_report.csv"), index=False)

    print("\n[4/6] Computing UMAP embedding ...")
    sc.tl.umap(adata, random_state=SEED)
    print("  * UMAP done.")

    print("\n[5/6] Comparing clusters to ground truth ...")
    true_labels = adata.obs["cell_type"].astype(str).values
    pred = adata.obs["leiden"].astype(int).values
    ari = adjusted_rand_score(true_labels, pred)
    nmi = normalized_mutual_info_score(true_labels, pred)
    print(f"  * ARI = {ari:.3f}, NMI = {nmi:.3f}")

    print(f"\n[6/6] Ranking marker genes (Wilcoxon, top {N_TOP_MARKERS}/cluster) ...")
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method="wilcoxon", use_raw=False
    )
    marker_rows = []
    for cluster in adata.obs["leiden"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df = df.head(N_TOP_MARKERS)
        df.insert(0, "cluster", cluster)
        marker_rows.append(df)
    markers_df = pd.concat(marker_rows, ignore_index=True)
    markers_df.to_csv(os.path.join(OUT_DIR, "marker_genes.csv"), index=False)
    print(f"  * saved -> output/marker_genes.csv ({len(markers_df)} rows)")

    # Export the clustered AnnData.
    adata.write(os.path.join(OUT_DIR, "leiden_clusters.h5ad"))
    print("  * saved -> output/leiden_clusters.h5ad")

    print("\nPlotting ...")
    # ---- Plot 1: UMAP colored by Leiden cluster ---------------------------- #
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(
        adata,
        color="leiden",
        palette="tab20",
        ax=ax,
        frameon=False,
        show=False,
        title=f"UMAP — Leiden clusters (res={best_res})",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "umap_clusters.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 2: ARI / NMI vs resolution + UMAP by ground truth ------------ #
    fig = plt.figure(figsize=(13, 5))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(report["resolution"], report["ARI"], marker="o", color="#4C72B0", label="ARI")
    ax1.plot(report["resolution"], report["NMI"], marker="s", color="#DD8452", label="NMI")
    ax1.axvline(best_res, color="grey", ls="--", lw=1)
    ax1.set_xlabel("Resolution")
    ax1.set_ylabel("Score")
    ax1.set_title(f"Cluster quality vs resolution (best res={best_res})")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(1, 2, 2)
    sc.pl.umap(
        adata,
        color="cell_type",
        palette="tab20",
        ax=ax2,
        frameon=False,
        show=False,
        title="UMAP — ground truth",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "umap_resolution_compare.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 3: top-marker heatmap ---------------------------------------- #
    tops = adata.obs["leiden"].cat.categories
    top_genes = markers_df.groupby("cluster")["names"].apply(list).to_dict()
    genes_for_heatmap = []
    for c in tops:
        genes_for_heatmap += [g for g in top_genes[c]][:8]
    genes_for_heatmap = list(dict.fromkeys(genes_for_heatmap))  # unique, ordered
    fig = plt.figure(figsize=(max(9, len(genes_for_heatmap) * 0.22), 6))
    sc.pl.heatmap(
        adata,
        var_names=genes_for_heatmap,
        groupby="leiden",
        use_raw=False,
        swap_axes=True,
        cmap="viridis",
        show=False,
        figsize=(max(9, len(genes_for_heatmap) * 0.22), 6),
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "marker_heatmap.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 4: dotplot of top markers ------------------------------------ #
    fig = plt.figure(figsize=(max(10, len(top_genes) * 1.6), max(6, len(genes_for_heatmap) * 0.3)))
    sc.pl.dotplot(
        adata,
        var_names=genes_for_heatmap,
        groupby="leiden",
        use_raw=False,
        show=False,
        figsize=(max(10, len(top_genes) * 1.6), max(6, len(genes_for_heatmap) * 0.3)),
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "marker_dotplot.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("  * saved 4 figures to output/")
    print("\nDone. Run 05_diffexp_annotation.py next.\n")


if __name__ == "__main__":
    main()
