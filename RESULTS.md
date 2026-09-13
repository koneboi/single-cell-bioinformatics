# Real-Data Results: 10x Genomics PBMC 3k

## Dataset

| Field | Value |
|-------|-------|
| **Name** | PBMC 3k (3k PBMCs from a Healthy Donor) |
| **Source** | [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k) |
| **Download URL** | `https://falexwolf.de/data/pbmc3k_raw.h5ad` (cached by scanpy) |
| **Original file** | `https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz` |
| **Chemistry** | 10x Chromium v2 (Chemistry 1.1) |
| **Sequencing** | Illumina, ~90k reads/cell |
| **Raw shape** | 2,700 cells x 32,738 genes |
| **License** | Freely available from 10x Genomics |
| **Citation** | Zheng et al. (2017), *PeerJ*; also used in the [scanpy PBMC3k tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) and [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) |

## Key Results

| Metric | Value |
|--------|-------|
| **Cells after QC** | 2,638 (62 removed) |
| **Genes after QC** | 13,656 (19,082 removed — lowly expressed) |
| **Highly variable genes** | 2,000 |
| **Leiden clusters** | 6 (resolution 0.5) |
| **ARI vs annotation** | 1.000 (expected — annotation is cluster-derived) |
| **NMI vs annotation** | 1.000 |

### Cell-type composition

| Cluster | Cell type | Cells | % of total | Canonical markers |
|---------|-----------|------:|----------:|-------------------|
| 0 | CD4+ T cells | 1,190 | 45.1% | CD3D, IL7R |
| 1 | B cells | 341 | 12.9% | MS4A1, CD79A |
| 2 | CD14+ Monocytes | 640 | 24.3% | CD14, LYZ, S100A8/A9 |
| 3 | Platelets | 12 | 0.5% | PPBP, PF4 |
| 4 | NK cells | 419 | 15.9% | GNLY, NKG7 |
| 5 | Dendritic cells | 36 | 1.4% | FCER1A, CST3 |

### Interpretation

The pipeline recovers the **major PBMC cell types** expected from a healthy donor blood sample:

- **CD4+ T cells** (45%) — the largest population, consistent with peripheral blood where naive/memory CD4+ T cells dominate. High IL7R expression confirms the memory/naive phenotype.
- **CD14+ Monocytes** (24%) — the second-largest population, expressing classical monocyte markers CD14 and LYZ. These are the main phagocytic population in PBMC.
- **NK cells** (16%) — high expression of cytotoxicity markers GNLY and NKG7. These cells are CD3-negative (or low) and functionally distinct from T cells.
- **B cells** (13%) — clearly separated, expressing B-cell receptor markers MS4A1 (CD20) and CD79A.
- **Dendritic cells** (1.4%) — rare antigen-presenting cells expressing FCER1A (high-affinity IgE receptor).
- **Platelets** (0.5%) — very rare, expressing platelet-specific markers PPBP and PF4. These are technically contaminants from blood processing.

The cell-type proportions are **biologically realistic** for human PBMC, where T cells typically comprise 50-70%, monocytes 10-30%, NK cells 5-15%, and B cells 5-15%.

## Generated Files

| File | Description |
|------|-------------|
| `output/real_qc_violin.png` | QC distributions: genes, counts, % mito |
| `output/real_qc_scatter.png` | Counts vs genes, counts vs % mito |
| `output/real_elbow.png` | PCA scree plot |
| `output/real_umap_clusters.png` | UMAP colored by Leiden cluster |
| `output/real_umap_annotated.png` | UMAP colored by annotated cell type |
| `output/real_marker_dotplot.png` | Dotplot of canonical marker genes |
| `output/real_marker_heatmap.png` | Heatmap of top markers per cluster |
| `output/real_marker_genes.csv` | Top 5 Wilcoxon markers per cluster |
| `output/real_cluster_summary.csv` | Cluster-to-cell-type mapping with counts |
| `output/real_pbmc3k.h5ad` | Full AnnData with all annotations |
| `data/real/pbmc3k_raw.h5ad` | Cached raw dataset (21 MB) |

## Reproduce

```bash
python scripts/real_10x_pbmc_pipeline.py
```

The script automatically downloads PBMC 3k from figshare on first run and caches it in `data/real/`.

## Run with your own scRNA-seq data

### Expected input

The pipeline accepts data in either of two formats:

**1. 10x Cell Ranger output** (recommended for raw data):

```
your_data/
  matrix.mtx      # sparse count matrix (Market Exchange format)
  features.tsv    # gene IDs and names (tab-separated)
  barcodes.tsv    # cell barcodes (one per line)
```

Load with:
```python
import scanpy as sc
adata = sc.read_10x_mtx("your_data/", var_names="gene_symbols", cache=True)
```

**2. AnnData .h5ad file** (recommended for processed data):

```python
import scanpy as sc
adata = sc.read_h5ad("your_data.h5ad")
```

### How to adapt the pipeline

1. Replace the data-loading step in `real_10x_pbmc_pipeline.py` with your file path
2. Adjust QC thresholds if needed (e.g., `MAX_PCT_MITO`, `MIN_GENES`, `MAX_GENES`)
3. Update `CELL_TYPE_MARKERS` with markers relevant to your tissue/cell types
4. Adjust `LEIDEN_RESOLUTION` to get the expected number of clusters
5. Re-interpret clusters based on marker expression and biological knowledge

### Quick checklist

- [ ] Are the gene names in HGNC symbol format? If using Ensembl IDs, set `var_names="gene_symbols"` in `read_10x_mtx`
- [ ] Is the matrix raw counts (integers)? If already log-normalised, skip `normalize_total` and `log1p`
- [ ] What is the expected number of cell types? Adjust Leiden resolution accordingly
- [ ] Do you have known marker genes for your tissue? Update `CELL_TYPE_MARKERS`
