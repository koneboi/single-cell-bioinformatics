# Single-Cell RNA-seq Analysis (scanpy)

**Comprehensive single-cell transcriptomics workflow — from raw counts to biological interpretation**

This repository provides a complete, end-to-end single-cell RNA-seq (scRNA-seq) analysis
pipeline built with [scanpy](https://scanpy.readthedocs.io/). It is designed to be run on a
**synthetic dataset generated automatically**, so the entire workflow can be reproduced
locally without downloading any real sequencing data. The same scripts can be pointed at a
real count matrix (e.g. Cell Ranger `.h5` / `.mtx` output) with minimal changes.

## Pipeline

```
QC → filtering → normalization → HVG → PCA → batch integration → clustering
→ marker genes → cell-type annotation → differential expression → trajectory
```

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `scripts/01_simulate_scrnaseq.py` | Generate a realistic synthetic count matrix + cell-type labels |
| 02 | `scripts/02_qc_filtering.py` | Compute QC metrics and filter low-quality cells |
| 03 | `scripts/03_normalize_hvg_pca.py` | Normalize, find HVGs, run PCA |
| 04 | `scripts/04_clustering_markers.py` | kNN graph, Leiden clustering, UMAP, marker genes |
| 05 | `scripts/05_diffexp_annotation.py` | Annotate cell types & run differential expression |
| 06 | `scripts/06_trajectory_pseudotime.py` | Diffusion map & pseudotime analysis |

## Skills

- **Python / scRNA-seq bioinformatics** (scanpy, anndata)
- **Clustering**: Leiden & Louvain (leidenalg)
- **Visualization**: UMAP, PCA, violin / scatter / volcano / dotplot / heatmap plots
- **Differential expression**: Wilcoxon rank-sum tests (`sc.tl.rank_genes_groups`)
- **Trajectory & pseudotime**: diffusion map (`sc.tl.diffmap`), diffusion pseudotime (`sc.tl.dpt`)
- Optional velocity scaffolds with **scvelo**

## Author

- **Boï Kone** (Master's-level bioinformatician)
- GitHub: [github.com/koneboi](https://github.com/koneboi)

---

## Dependencies

Core dependencies (see [`requirements.txt`](requirements.txt)):

| Package | Purpose |
|---------|---------|
| `scanpy` | Main analysis framework (>= 1.9) |
| `anndata` | In-memory single-cell data structure |
| `numpy` / `scipy` / `pandas` | Numerical & tabular data handling |
| `scikit-learn` | Silhouette score, metrics (ARI / NMI) |
| `matplotlib` / `seaborn` | Plotting |
| `leidenalg` | Leiden clustering |
| `umap-learn` | UMAP embedding |
| `scvelo` | Optional RNA velocity / dynamics |

Install with:

```bash
pip install -r requirements.txt
```

## How to run

Run the scripts **in order**, from the repository root. Each script only depends on the
output of the previous ones, and all expected input files are produced by stage 01.

```bash
# 1. Simulate the synthetic dataset (must run first)
python scripts/01_simulate_scrnaseq.py

# 2. QC & filtering
python scripts/02_qc_filtering.py

# 3. Normalization, HVG & PCA
python scripts/03_normalize_hvg_pca.py

# 4. Clustering & marker genes
python scripts/04_clustering_markers.py

# 5. Annotation & differential expression
python scripts/05_diffexp_annotation.py

# 6. Trajectory & pseudotime
python scripts/06_trajectory_pseudotime.py
```

Or run the whole pipeline in one shot:

```bash
for f in scripts/0*.py; do echo "=== Running $f ==="; python "$f"; done
```

## Output

All artifacts (`.h5ad` AnnData objects, CSV tables, and `*.png` figures) are written to the
[`output/`](output/) directory:

```
output/
├── counts_matrix.h5ad       # raw simulated data (from 01)
├── qc_metrics.csv           # per-cell QC metrics (from 02)
├── filtered.h5ad            # QC-filtered data (from 02)
├── normalized.h5ad          # normalized + HVG + PCA data (from 03)
├── leiden_clusters.h5ad     # clustered + UMAP data (from 04)
├── marker_genes.csv         # Wilcoxon top markers per cluster (from 04)
├── annotated.h5ad           # annotated + DE results (from 05)
├── DE_<celltypeA>_vs_<celltypeB>.csv  # DE table (from 05)
├── trajectory.h5ad          # diffmap + pseudotime data (from 06)
└── *.png                    # all figures
```

---

## Notes on the synthetic data

Stage 01 simulates a count matrix for 6 cell types (T cells, B cells, monocytes/clusters,
NK cells, dendritic cells, neutrophils) using a simplified generative model of gene
expression that includes:

- **Cell-type-specific marker genes** with elevated mean expression
- **Dropout** (zero-inflation) mimicking low capture efficiency
- **Technical noise** (overdispersion) via negative-binomial-like sampling

This produces realistic-looking data (shape, sparsity, gene-level mean/variance relation)
so that downstream QC, clustering, and DE all behave like real scRNA-seq data.

> **Warning:** The mitochondrial fraction is *simulated* (set in 02) since it is not
> explicitly modeled in 01 — this keeps the synthetic pipeline self-contained while still
> demonstrating the standard QC filtering workflow.
