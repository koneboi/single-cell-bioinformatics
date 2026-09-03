"""
06_trajectory_pseudotime.py
===========================

Explore cell-state ordering with a diffusion map and diffusion pseudotime (DPT).

Steps:
    * Recompute/verify neighbors on the annotated data.
    * Compute a diffusion map (sc.tl.diffmap).
    * Compute diffusion pseudotime (sc.tl.dpt) with a data-driven root cell.
    * Project onto the existing UMAP.
    * Plot pseudotime and the expression of key lineage genes along pseudotime.

NOTE: This is an unsupervised *illustrative* trajectory. The synthetic data does
not embed a genuine developmental continuum, so the pseudotime should be read as
a demonstration of the scanpy DPT machinery rather than a real lineage inference.

Outputs (output/):
    * trajectory.h5ad          - AnnData with diffmap + dpt stored
    * umap_pseudotime.png      - UMAP colored by pseudotime
    * pseudotime_gene_trends.png - lineage gene expression vs pseudotime
    * diffmap_components.png   - diffusion-map scatter (components 1-2, 2-3)
    * pseudotime_distribution.png - histogram of pseudotime values
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
ANNOTATED_FILE = os.path.join(OUT_DIR, "annotated.h5ad")

# Genes whose expression we trace along pseudotime.
LINEAGE_GENES = [
    "CD3D", "CD2",      # T-cell lineage
    "CD79A", "MS4A1",   # B-cell lineage
    "CD14", "LYZ",      # monocyte lineage
    "NKG7", "GNLY",     # NK lineage
]
SEED = 42


def pick_root(adata, n_candidates=25):
    """
    Choose a DPT root cell using a robust heuristic.

    We try several candidate roots (distributed across the diffusion embedding /
    cell-type clusters) and keep the one that reaches the most cells with a
    *finite* pseudotime (i.e. the largest DPT-connected component). This makes
    DPT robust when the data contains disconnected clusters, which is common in
    real datasets and in our discrete simulated population.
    """
    rng = np.random.default_rng(SEED)

    # Candidate set: spread across the data so at least one root lands inside
    # the largest connected component of the nearest-neighbour graph.
    if "X_diffmap" in adata.obsm and adata.obsm["X_diffmap"].shape[1] >= 2:
        coord = adata.obsm["X_diffmap"][:, :2]
    else:
        coord = adata.obsm["X_pca"][:, :2]

    candidates = []
    # Include the extremes of the first two diffusion components (endpoints).
    for c in range(min(2, coord.shape[1])):
        candidates.append(int(np.argmin(coord[:, c])))
        candidates.append(int(np.argmax(coord[:, c])))
    candidates.append(int(np.argmax(np.linalg.norm(coord - coord.mean(0), axis=1))))
    # Add a random spread.
    candidates += list(rng.choice(adata.n_obs, size=n_candidates, replace=False))
    candidates = list(dict.fromkeys(int(c) for c in candidates))

    best_root, best_n = candidates[0], -1
    for root in candidates:
        adata.uns["iroot"] = root
        try:
            sc.tl.dpt(adata, n_dcs=10, n_branchings=0)
        except Exception:
            continue
        pt = np.asarray(adata.obs["dpt_pseudotime"].values.astype(float))
        n_fin = int(np.isfinite(pt).sum())
        if n_fin > best_n:
            best_n, best_root = n_fin, root
    print(f"  * best candidate root covers {best_n} cells with finite pseudotime")
    adata.uns["iroot"] = best_root
    return best_root


def main():
    print("=" * 70)
    print(" 06_trajectory_pseudotime.py — Diffusion map & pseudotime")
    print("=" * 70)
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"\n[1/4] Loading {ANNOTATED_FILE} ...")
    adata = sc.read_h5ad(ANNOTATED_FILE)
    print(f"  * loaded: {adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs['annotation'].nunique()} cell types")

    # Ensure graph is present; recompute neighbors on PCA if missing.
    if "neighbors" not in adata.uns:
        print("  * neighbors not found, recomputing (n_neighbors=15)...")
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)

    print("\n[2/4] Computing diffusion map ...")
    sc.tl.diffmap(adata, n_comps=15)
    print("  * diffusion map computed (15 components)")

    print("\n[3/4] Computing diffusion pseudotime (DPT) ...")
    root = pick_root(adata)
    adata.uns["iroot"] = root
    sc.tl.dpt(adata, n_dcs=10, n_branchings=0)
    print(f"  * root cell index: {root}")

    # Project DPT onto the existing UMAP (dpt is stored in obs['dpt_pseudotime']).
    if "X_umap" not in adata.obsm:
        print("  * UMAP missing, recomputing ...")
        sc.tl.umap(adata, random_state=SEED)

    adata.write(os.path.join(OUT_DIR, "trajectory.h5ad"))
    print("  * saved -> output/trajectory.h5ad")

    # DPT reports `inf` for cells disconnected from the chosen root; mark as NaN
    # and keep only the finite, ordered set for pseudotime-driven plots.
    pt_raw = adata.obs["dpt_pseudotime"].values.astype(float)
    n_finite = int(np.isfinite(pt_raw).sum())
    print(f"  * finite pseudotime cells: {n_finite}/{adata.n_obs}")
    if n_finite < adata.n_obs:
        print("  * (some cells are disconnected from the chosen root -> excluded from plots)")
    # Replace inf with NaN so scanpy/matplotlib handle disconnected cells as blank.
    adata.obs["dpt_pseudotime"] = np.where(np.isfinite(pt_raw), pt_raw, np.nan)
    pt = pt_raw[np.isfinite(pt_raw)]
    if len(pt):
        print(f"  * pseudotime range (finite cells): [{pt.min():.3f}, {pt.max():.3f}]")

    print("\nPlotting ...")
    # ---- Plot 1: UMAP colored by pseudotime -------------------------------- #
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(
        adata,
        color="dpt_pseudotime",
        cmap="viridis",
        ax=ax,
        frameon=False,
        show=False,
        title="UMAP — diffusion pseudotime",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "umap_pseudotime.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 2: lineage gene expression vs pseudotime --------------------- #
    genes = [g for g in LINEAGE_GENES if g in adata.var_names]
    if not genes:
        genes = list(adata.var_names[:4])
    pt_all = adata.obs["dpt_pseudotime"].values.astype(float)
    finite = np.isfinite(pt_all)
    order = np.argsort(pt_all[finite])
    pt_sorted = pt_all[finite][order]
    cell_idx = np.nonzero(finite)[0][order]

    n_genes = len(genes)
    ncols = 4
    nrows = int(np.ceil(n_genes / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 2.6))
    axes = np.atleast_1d(axes).ravel()
    X = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()
    for i, gene in enumerate(genes):
        gidx = adata.var_names.get_loc(gene)
        vec = np.asarray(X[cell_idx, gidx]).ravel()
        # Smooth with a rolling mean.
        window = max(5, len(vec) // 30)
        window = window + 1 if window % 2 == 0 else window
        kernel = np.ones(window) / window
        smooth = np.convolve(vec, kernel, mode="same")
        axes[i].plot(pt_sorted, smooth, color="#4C72B0", lw=1.6)
        axes[i].scatter(pt_sorted, vec, s=2, alpha=0.15, color="#4C72B0")
        axes[i].set_title(gene, fontsize=9)
        axes[i].set_xlabel("Pseudotime", fontsize=7)
        axes[i].tick_params(labelsize=6)
    # Hide empty axes.
    for j in range(n_genes, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle("Lineage gene expression along pseudotime", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pseudotime_gene_trends.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 3: diffusion-map scatter (components 1-2 and 2-3) ------------ #
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.diffmap(
        adata, color="annotation", palette="tab20", ax=axes[0],
        frameon=False, show=False, title="Diffusion map (DC1-DC2)",
    )
    sc.pl.diffmap(
        adata, color="dpt_pseudotime", cmap="viridis", ax=axes[1],
        frameon=False, show=False, title="Diffusion map (DC1-DC2) — pseudotime",
        components="1,2",
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "diffmap_components.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 4: pseudotime distribution ----------------------------------- #
    fig, ax = plt.subplots(figsize=(7, 4.5))
    finite_pt = np.asarray(adata.obs["dpt_pseudotime"].values.astype(float))
    finite_pt = finite_pt[np.isfinite(finite_pt)]
    if len(finite_pt) == 0:
        finite_pt = np.array([0.0])
    ax.hist(finite_pt, bins=60, color="#55A868", edgecolor="white")
    ax.axvline(np.median(finite_pt), color="crimson", ls="--", lw=1.2, label="median")
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Number of cells")
    ax.set_title("Distribution of diffusion pseudotime")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pseudotime_distribution.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("  * saved 4 figures to output/")
    print("\nDone. Pipeline complete.\n")


if __name__ == "__main__":
    main()
