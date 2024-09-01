#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Spyder Editor
scVelo script for python 3.8+ 

"""

import os
import scvelo as scv
import pandas as pd
import anndata as ad
import numpy as np

# change matplotlib settings to scv settings
scv.set_figure_params('scvelo')

# ---------------------------------------------------------------------------
# -----------------------------     DATA INPUT   ----------------------------
# ---------------------------------------------------------------------------
# read in data (loom, h5ad, csv. etc) w both unspliced+spliced counts

# ..........................................................................
# a) if working with looms files, convert to anndata:
# dirs
loom_dir = "/cloud-home/mfr/scrnaseq/ingested/ILD_Lung_GSE135893/looms" 
output_h5ad_file = "/cloud-home/mfr/scrnaseq/ingested/ILD_Lung_GSE135893/looms/combined_data.h5ad"  

# list to hold individual AnnData objects
adatas = []

# loop over each Loom file in the loom_dir
for loom_file in os.listdir(loom_dir):
    if loom_file.endswith(".loom"):
        loom_path = os.path.join(loom_dir, loom_file)
        sample_name = os.path.splitext(loom_file)[0] 
        # load Loom file into AnnData
        adata = scv.read(loom_path, cache=True)
        # add sample metadata
        adata.obs['sample'] = sample_name
        # append to list 
        adatas.append(adata)

# concat all AnnData objs
combined_adata = ad.concat(adatas, merge='same')

# save combined AnnData obj
combined_adata.write(output_h5ad_file)
print(f"Combined AnnData object saved to {output_h5ad_file}")

# bring in adata
adata = sc.read(output_h5ad_file, cache=True)

# ..........................................................................
# b) if working with scvelo dataset:
# or use built in dataset (4 cell fates; alpha, beta, delta, epsilon cells)
adata = scv.datasets.pancreas()
#adata = scv.datasets.dentategyrus()
adata


# ---------------------------------------------------------------------------
# -----------------------------     PROCESSING   ----------------------------
# ---------------------------------------------------------------------------
# check proportions of spliced vs unspliced counts; usually 10-25% is unspliced (containing introns)
scv.pl.proportions(adata)
# data should be processed for min umi/cell, min gene/cell, batch corr

# use same filters as was done for seurat obj for a particular ds
sc.pp.filter_cells(adata, min_counts=250)
sc.pp.filter_cells(adata, min_genes=500)
# filter genes expressed in fewer than 3 cells
sc.pp.filter_genes(adata, min_cells=3)
# mitochondrial gene related qc
adata.var['mito'] = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, adata.var['mito']].X, axis=1) / np.sum(adata.X, axis=1) * 100
sc.pp.filter_cells(adata, max_percent_mito=10) 
scv.pp.log1p(adata)

# additional scv filters
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# ---------------------------------------------------------------------------
# -----------------------------     ANALYSIS    -----------------------------
# ---------------------------------------------------------------------------
# estimate velocities, which will be stored in adata.layers
# velocities = vectors in GE spaces obtained by stochastic model of transcriptional dynamics
# other mode options = deterministic, dynamical
# to run dynamical, need to run the following first: scv.tl.recover_dynamics(adata, **params)
scv.tl.velocity(adata, mode='stochastic') # stochastic is default
# project velocities into lower-D UMAP embedding by translating to likely cell transitions 
# cell transition in accordance to direction of velocity vector inferred 
# i.e. transition probs of cell-cell transitions estimated
# transition prob calc'd w cosine correlation between cell-cell transitions and velocity vector and stored in velocity graph
scv.tl.velocity_graph(adata)
# access markov transition matrix
scv.utils.get_transition_matrix
scv.pl.velocity_embedding(adata, basis='umap')
scv.pl.velocity_embedding_grid(adata, basis='umap')
scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters')

# id genes that help explain inferred lienages
# to infer, test which genes have cluter-specific differentail velocity expression
# velocity t-test
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

# panel per gene
scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True)

#graph confidence of speed/rate of differentiation 
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

# velocity length/confidence for diff cell types 
df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

# velocity graph to portray infererred cell-cell transitions 
scv.pl.velocity_graph(adata, threshold=.1)

# best graph - paga trajectory inference 
# this is needed due to a current bug - bugfix is coming soon
#!pip install python-igraph --upgrade --quiet
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)
