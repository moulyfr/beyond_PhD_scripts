#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr

# ------------------------------------------------------------------
# 1) processing
# ------------------------------------------------------------------
adata = sc.read_csv('raw_counts/GSM5226574_C51ctr_raw_counts.csv').T
adata.var['mt'] = adata.var.index.str.startswith('MT-')
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.calculate_qc_metrics(adata, qc_vars=['pct_counts_mt'], inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20]

# ------------------------------------------------------------------
# 2) doublet removal - Scrublet
# ------------------------------------------------------------------
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_scores'] = doublet_scores
adata.obs['predicted_doublets'] = predicted_doublets
adata = adata[~adata.obs['predicted_doublets']]

# ------------------------------------------------------------------
# 3) norm'n
# ------------------------------------------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ------------------------------------------------------------------
# 4) clustering
# ------------------------------------------------------------------
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.pl.umap(adata, show=True)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'], show=True)

# ------------------------------------------------------------------
# 5) cell type annotation
# ------------------------------------------------------------------
sc.tl.leiden(adata, resolution=1)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=True)
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]
markers

# markers_scvi = model.differential_expression(groupby='leiden')
# markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > 0.5)]
# markers_scvi

# Cell type annotation
cell_type = {
    "0": "Macrophage",
    "1": "Fibroblast",
    "2": "CD4+ T-cell",
    "3": "AT2",
    "4": "AT1",
    "5": "CD8+ T-cell",
    "6": "Endothelial cell",
    "7": "Plasma cell",
    "8": "Macrophage",
    "9": "AT2",
    "10": "Fibroblast",
    "11": "Fibroblast",
    "12": "Macrophage",
    "13": "Macrophage",
    "14": "Airway epithelial",
    "15": "Airway epithelial",
    "16": "Monocyte",
    "17": "Airway epithelial",
    "18": "B-cell",
    "19": "Aerocyte",
    "20": "Airway epithelial",
    "21": "Smooth muscle cell",
    "22": "Cycling T/NK",
    "23": "Neuronal cell",
    "24": "Dendritic cell",
    "25": "Pericyte",
    "26": "Fibroblast",
    "27": "Erythroid-like",
    "28": "Macrophage"
}

adata.obs['cell_type'] = adata.obs.leiden.map(cell_type)
sc.pl.umap(adata, color=['cell_type'], frameon=False, show=True)
adata.write_h5ad('integrated.h5ad')
