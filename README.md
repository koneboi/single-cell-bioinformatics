# Single-Cell RNA-seq Analysis (scanpy)

**Comprehensive single-cell transcriptomics workflow — from raw counts to biological interpretation**

This repository provides a complete, end-to-end single-cell RNA-seq (scRNA-seq) analysis
pipeline built with [scanpy](https://scanpy.readthedocs.io/). It includes two complementary
workflows:

1. A **synthetic-data demonstration** (scripts `01`–`06`) that runs the entire pipeline on a
   simulated dataset, so the workflow can be exercised locally without any downloads.
2. A **real-data analysis** ([`scripts/real_10x_pbmc_pipeline.py`](scripts/real_10x_pbmc_pipeline.py))
   that runs the same standard scanpy workflow on the **10x Genomics PBMC 3k** dataset — ~2700
   real peripheral blood mononuclear cells from a healthy donor — and recovers the expected
   cell types (T cells, B cells, monocytes, NK cells, dendritic cells, platelets).

See [`RESULTS.md`](RESULTS.md) for the real-data results, dataset sourcing, and a guide to
running the pipeline on your own scRNA-seq data. The real-data workflow is also fully
reproducible: the raw dataset is cached under `data/real/`.

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
| Real | `scripts/real_10x_pbmc_pipeline.py` | Full scanpy workflow on real 10x PBMC 3k dataset |
| Real | `scripts/download_real_data.py` | Download/cache the PBMC 3k dataset (urllib) |

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

## Real-data workflow

The real-data pipeline downloads (on first run) and analyses the **10x PBMC 3k** dataset —
an actual public scRNA-seq dataset of peripheral blood mononuclear cells. This demonstrates
the same workflow on genuine data and produces figures + a biological interpretation.

```bash
# Optional: pre-download the dataset (otherwise the pipeline downloads it automatically)
python scripts/download_real_data.py

# Run the full real-data analysis (~2-4 minutes)
python scripts/real_10x_pbmc_pipeline.py
```

The real-data results (cell-type recovery, ARI, marker genes, figures) are documented in
[`RESULTS.md`](RESULTS.md), which also includes a guide for running the pipeline on your own
scRNA-seq data (10x `matrix.mtx` / `features.tsv` / `barcodes.tsv` or an `.h5ad` file).

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

### Real-data outputs (`output/real_*`)

The real-data pipeline writes all artifacts with a `real_` prefix:

```
output/
├── real_pbmc3k.h5ad         # full AnnData with annotations
├── real_qc_violin.png       # QC distributions (genes, counts, %mito)
├── real_qc_scatter.png      # counts vs genes, counts vs %mito
├── real_elbow.png           # PCA scree plot
├── real_umap_clusters.png   # UMAP by Leiden cluster
├── real_umap_annotated.png  # UMAP by annotated cell type
├── real_marker_dotplot.png  # canonical marker dotplot
├── real_marker_heatmap.png  # top-marker heatmap
├── real_marker_genes.csv    # top Wilcoxon markers per cluster
└── real_cluster_summary.csv # cluster -> cell-type mapping
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

## Real-data notes

The [`scripts/real_10x_pbmc_pipeline.py`](scripts/real_10x_pbmc_pipeline.py) workflow uses
the **10x Genomics PBMC 3k** dataset — real scRNA-seq from a healthy donor, obtained freely
from 10x Genomics (Chemistry 1.1). It is the same canonical dataset used throughout the
scanpy and Seurat tutorials. Full source/licensing/citation details, results, and a guide to
running the pipeline on your own data are in [`RESULTS.md`](RESULTS.md).
