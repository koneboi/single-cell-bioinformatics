"""
01_simulate_scrnaseq.py
=======================

Simulate a realistic synthetic scRNA-seq count matrix using a simplified model of
gene expression and known cell types.

Generative model:
    * 6 cell types with distinct marker-gene expression profiles:
      T cells, B cells, Monocytes, NK cells, Dendritic cells, Neutrophils.
    * Each gene has a baseline (background) mean expression level.
    * Marker genes are boosted for their associated cell type.
    * Counts are sampled from a negative-binomial-like distribution
      (mean = expression, overdispersion = noise).
    * Dropout (zero-inflation): a fraction of entries are forced to zero,
      mimicking the low capture efficiency of scRNA-seq.

Outputs (written to output/):
    * counts_matrix.h5ad   - AnnData with raw counts + ground-truth cell type labels
    * qc_preplots_*        - preliminary QC overview figures
    * simulate_UMAP_groundtruth.png / simulate_marker_heatmap.png
"""

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
OUT_DIR = "output"
SEED = 42

N_CELLS_PER_TYPE = 300        # cells per cell type
N_BASE_GENES = 1200           # "housekeeping" / background genes
N_MARKERS_PER_TYPE = 60       # marker genes per cell type
DROPOUT_BASE = 0.35           # baseline dropout probability
NOISE_FACTOR = 0.4            # extra overdispersion

CELL_TYPES = [
    "T cells",
    "B cells",
    "Monocytes",
    "NK cells",
    "Dendritic cells",
    "Neutrophils",
]

# Canonical-ish marker genes used for annotation / DE sanity checks downstream.
# These don't need to be biologically exhaustive; the simulation just needs to
# reproducibly assign distinct expression to each type.
MARKERS = {
    "T cells": [
        "CD3D", "CD3E", "CD3G", "CD2", "CD28", "IL7R", "TRAC", "TRBC1", "TRBC2",
        "CD6", "LCK", "FYN", "ZAP70", "CD8A", "CD8B", "CD4", "CCR7", "SELL",
    ],
    "B cells": [
        "CD19", "CD79A", "CD79B", "MS4A1", "CD22", "BLK", "BANK1", "PAX5", "IGHM",
        "IGHD", "CD74", "HLA-DRA", "HLA-DPB1", "MZB1", "CD37", "FCRLA", "TCL1A",
    ],
    "Monocytes": [
        "CD14", "FCGR3A", "CSF1R", "LYZ", "S100A8", "S100A9", "ITGAM", "CD68",
        "MSR1", "MNDA", "FCN1", "VCAN", "TYROBP", "CTSS", "SERPINA1", "C1QA",
    ],
    "NK cells": [
        "NKG7", "GNLY", "KLRD1", "KLRF1", "NKG2A", "TYROBP", "GZMB", "GZMH",
        "PRF1", "FCGR3A", "KLRB1", "KLRC1", "CD160", "CX3CR1", "ZEB2", "SPON2",
    ],
    "Dendritic cells": [
        "CD1C", "FCER1A", "LILRA4", "CLEC9A", "BATF3", "IRF8", "FLT3", "CLEC10A",
        "GPR183", "CADM1", "ITGAX", "HLA-DQA1", "HLA-DQB1", "CD83", "CCR7", "LAMP3",
    ],
    "Neutrophils": [
        "FCGR3B", "CSF3R", "S100A12", "CXCL8", "CEACAM3", "CEACAM8", "ARG1",
        "MNDA", "FCGR2A", "G0S2", "S100P", "FPR1", "FPR2", "CD177", "PGLYRP1",
        "TCN1",
    ],
}

# Keep the simulation deterministic.
rng = np.random.default_rng(SEED)


def build_gene_names():
    """Return a list of gene names: markers first, then background genes."""
    genes = []
    for _type, markers in MARKERS.items():
        genes += markers
    # A few dedicated per-type markers that are not in the curated list
    # (guarantees each type has *some* unique signature even if markers overlap).
    extra = []
    for i in range(N_BASE_GENES):
        extra.append(f"BGENE_{i}")
    # Drop duplicates while preserving order.
    seen = set()
    uniq = []
    for g in genes + extra:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq


def gene_type_mask(gene_names, marker_genes):
    """Boolean mask marking genes that are markers of a given cell type."""
    return np.array([g in marker_genes for g in gene_names], dtype=bool)


def simulate_counts():
    """Generate the count matrix and cell-type labels."""
    gene_names = build_gene_names()
    n_genes = len(gene_names)

    n_cells = N_CELLS_PER_TYPE * len(CELL_TYPES)
    # Background mean expression per gene (log-normal distributed).
    background_mean = np.exp(rng.normal(0.0, 1.2, size=n_genes)) * 5.0

    # Build per-cell-type mean expression matrix (cells x genes).
    counts = np.zeros((n_cells, n_genes))
    labels = []

    idx = 0
    for ctype in CELL_TYPES:
        markers = set(MARKERS[ctype])
        is_marker = gene_type_mask(gene_names, markers)
        for _ in range(N_CELLS_PER_TYPE):
            mean = background_mean.copy()
            # Boost marker genes for this cell type.
            mean[is_marker] *= rng.uniform(8.0, 15.0)
            # Add per-gene overdispersion / noise.
            mean = mean * np.exp(rng.normal(0.0, NOISE_FACTOR, size=n_genes))
            # Sample Poisson counts with a random per-cell size factor (library size).
            size_factor = rng.uniform(0.6, 1.4)
            lam = np.clip(mean * size_factor, 0.0, 1e5)
            counts[idx] = rng.poisson(lam)
            labels.append(ctype)
            idx += 1

    # Dropout (zero-inflation): with probability depending on the mean,
    # small counts are set to zero. Higher expression is more likely to be kept.
    with np.errstate(divide="ignore"):
        keep_prob = 1.0 - DROPOUT_BASE * np.exp(-counts / 5.0)
    keep_prob = np.clip(keep_prob, 0.05, 1.0)
    mask = rng.random(counts.shape) < keep_prob
    counts = counts * mask
    counts = counts.astype(np.int32)

    # Cell barcodes & metadata.
    barcodes = [f"cell_{i:05d}" for i in range(n_cells)]
    labels = np.array(labels)

    adata = sc.AnnData(X=counts)
    adata.obs_names = barcodes
    adata.var_names = gene_names
    adata.obs["cell_type"] = pd.Categorical(labels)
    return adata


def main():
    print("=" * 70)
    print(" 01_simulate_scrnaseq.py — Synthetic scRNA-seq data")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    print("\n[1/3] Simulating counts ...")
    adata = simulate_counts()
    print(f"  * shape       : {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"  * sparsity    : {(adata.X == 0).mean():.1%} zeros")
    nz = adata.X.sum(axis=1)
    print(f"  * per-cell UMI: median {np.median(nz):.0f}, mean {np.mean(nz):.0f}")
    print("  * cell types  :")
    for ct in CELL_TYPES:
        print(f"      - {ct}: {(adata.obs['cell_type'] == ct).sum()} cells")

    # Store gene-level housekeeping metadata (marker status not needed downstream).
    adata.obs["n_counts"] = np.asarray(nz).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    print("\n[2/3] Saving raw AnnData ...")
    adata.write(os.path.join(OUT_DIR, "counts_matrix.h5ad"))
    print(f"  * saved -> {OUT_DIR}/counts_matrix.h5ad")

    print("\n[3/3] Making pre-QC plots ...")
    # A) Ground-truth structure: UMAP of library-size-normalized, log data.
    ad_tmp = adata.copy()
    sc.pp.normalize_total(ad_tmp, target_sum=1e4)
    sc.pp.log1p(ad_tmp)
    sc.pp.pca(ad_tmp, n_comps=20)
    sc.pp.neighbors(ad_tmp, n_neighbors=15)
    sc.tl.umap(ad_tmp, random_state=42)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc.pl.umap(
        ad_tmp,
        color="cell_type",
        ax=ax,
        palette="tab20",
        frameon=False,
        show=False,
        title="Simulated ground truth (UMAP of log-normalized counts)",
    )
    fig.savefig(os.path.join(OUT_DIR, "simulate_UMAP_groundtruth.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # B) UMI / gene-count histograms.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(np.asarray(adata.obs["n_counts"]), bins=60, color="#4C72B0", edgecolor="white")
    axes[0].set_title("Library size (n_counts)")
    axes[0].set_xlabel("Total counts")
    axes[1].hist(np.asarray(adata.obs["n_genes"]), bins=60, color="#55A868", edgecolor="white")
    axes[1].set_title("Genes detected (n_genes)")
    axes[1].set_xlabel("Num. genes")
    fig.suptitle("Synthetic scRNA-seq — QC preplots")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "simulate_histograms.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("  * saved figures to output/")
    print("\nDone. Run 02_qc_filtering.py next.\n")


if __name__ == "__main__":
    main()
