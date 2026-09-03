"""
03_normalize_hvg_pca.py
=======================

Normalize the QC-filtered data, identify highly variable genes (HVGs), and run PCA.

Steps:
    * Normalize to a fixed library size (target_sum) with log1p transform.
    * Identify highly variable genes (HVGs) using the Seurat v3 flavor.
    * Restrict the matrix to HVGs and run PCA on the first ~30-40 PCs.
    * Save plots:
        - PCA colored by ground-truth cell type
        - Variance ratio (scree) plot
        - HVG dispersion/mean plot

Outputs (output/):
    * normalized.h5ad - AnnData with `X` log-normalized on HVGs, PCA stored
    * pca_groundtruth.png
    * pca_scree.png
    * hvg_plot.png
"""

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
OUT_DIR = "output"
FILTERED_FILE = os.path.join(OUT_DIR, "filtered.h5ad")

TARGET_SUM = 1e4      # library size normalization target
N_HVG = 1200          # number of highly variable genes
N_PCS = 30            # number of principal components
SEED = 42


def main():
    print("=" * 70)
    print(" 03_normalize_hvg_pca.py — Normalization, HVG & PCA")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"\n[1/4] Loading {FILTERED_FILE} ...")
    adata = sc.read_h5ad(FILTERED_FILE)
    print(f"  * loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    print("\n[2/4] Normalizing (target_sum + log1p) ...")
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    print("  * done.")

    print(f"\n[3/4] Identifying {N_HVG} highly variable genes ...")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat", n_bins=20)
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"  * flagged {n_hvg} HVGs")
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"  * reduced to {adata.shape[1]} genes")

    print(f"\n[4/4] Running PCA ({N_PCS} components) ...")
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
    var_ratio = adata.uns["pca"]["variance_ratio"]
    evr = np.cumsum(var_ratio)
    print(f"  * variance explained by PC1-{N_PCS}: {evr[-1]:.1%}")

    # ---- Export normalized data -------------------------------------------- #
    adata.write(os.path.join(OUT_DIR, "normalized.h5ad"))
    print("  * saved -> output/normalized.h5ad")

    print("\nPlotting ...")
    # ---- Plot 1: PCA colored by ground-truth cell type --------------------- #
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.pca(
        adata,
        color="cell_type",
        palette="tab20",
        ax=ax,
        frameon=False,
        show=False,
        title="PCA (ground-truth cell types)",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pca_groundtruth.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 2: Scre plot (variance ratio) -------------------------------- #
    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(1, N_PCS + 1)
    ax.bar(x, var_ratio, color="#4C72B0", alpha=0.7, label="per-PC variance ratio")
    ax.plot(x, evr, color="crimson", marker="o", label="cumulative ratio")
    ax.axhline(0.02, color="grey", ls="--", lw=1, alpha=0.6)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance ratio")
    ax.set_title("PCA scree plot")
    ax.legend()
    ax.set_xticks(x[::2])
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pca_scree.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 3: HVG dispersion plot --------------------------------------- #
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "hvg_plot.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 4: PCA with % variance from the parameters ------------------- #
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(
        np.arange(1, N_PCS + 1),
        var_ratio,
        color="#55A868",
        s=40,
        zorder=3,
    )
    ax.axvline(5, color="grey", ls=":", lw=1)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained")
    ax.set_title("Explained variance per PC (top PCs highlighted)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pca_variance.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("  * saved 4 figures to output/")
    print("\nDone. Run 04_clustering_markers.py next.\n")


if __name__ == "__main__":
    main()
