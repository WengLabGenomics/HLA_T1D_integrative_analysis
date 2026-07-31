import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import doubletdetection
import os
import scanpy.external as sce

###___________________PBMC scRNA-seq___________________
date = 'PBMC'

file_path ='./processing_h5ad'

adata  = sc.read(file_path+'/PBMC_raw.h5ad')

clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)
for i in adata.obs['sample'].unique():
    temp = adata[adata.obs['sample'] == i]
    doublets = clf.fit(temp.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    temp.obs["doublet"] = doublets
    temp.obs["doublet_score"] = doublet_score
    if i == adata.obs['sample'].unique()[0]:
        adata3 = temp
    else:
        adata3 = ad.concat([adata3,temp],join = 'outer')

adata = adata3

adata = adata[adata.obs['doublet'] == 0]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts > 500, :]
adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.pct_counts_mt < 7, :]

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.write(file_path+'/{}/{}_merge_outer_QC_norm.h5ad'.format(date,date))
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sce.pp.harmony_integrate(adata, key=['batch'])
n_n = 35
n_pc = 30
sc.pp.neighbors(adata, n_neighbors=n_n, n_pcs=n_pc, use_rep='X_pca_harmony')
sc.tl.umap(adata)
res =1.3
sc.tl.leiden(adata, resolution=res)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
adata.write(file_path+'/{}/{}_harmony_{}_{}_{}_withcli.h5ad'.format(date,date,n_n,n_pc,res))

###___________________Pancreatic islet scRNA-seq___________________
date = 'Islet'

file_path ='./processing_h5ad'

adata  = sc.read(file_path+'/Islet_raw.h5ad')

clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)
for i in adata.obs['sample'].unique():
    temp = adata[adata.obs['sample'] == i]
    doublets = clf.fit(temp.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    temp.obs["doublet"] = doublets
    temp.obs["doublet_score"] = doublet_score
    if i == adata.obs['sample'].unique()[0]:
        adata3 = temp
    else:
        adata3 = ad.concat([adata3,temp],join = 'outer')

adata = adata3
adata = adata[adata.obs['doublet'] == 0]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts > 150, :]
adata = adata[adata.obs.n_genes_by_counts < 8500, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.write(file_path+'/{}/{}_merge_outer_QC_norm.h5ad'.format(date,date))
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sce.pp.harmony_integrate(adata, key=['batch'])
n_n = 35
n_pc = 30
sc.pp.neighbors(adata, n_neighbors=n_n, n_pcs=n_pc, use_rep='X_pca_harmony')
sc.tl.umap(adata)
res =1.3
sc.tl.leiden(adata, resolution=res)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
adata.write(file_path+'/{}/{}_harmony_{}_{}_{}_withcli.h5ad'.format(date,date,n_n,n_pc,res))