import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import doubletdetection
import os
import scanpy.external as sce
import scipy.stats as stats
import seaborn as sns
import statsmodels.api as sm
import warnings

import matplotlib.pyplot as plt

import matplotlib.ticker as mtick
import math
import scipy
import matplotlib as mpl

import statsmodels.stats.weightstats as sw

###___________________PBMC scRNA-seq___________________
raw = sc.read('./processing_h5ad/PBMC_raw.h5ad')
meta = pd.read_csv('./processing_h5ad/PBMC/PBMC_harmony_35_30_1.3_meta.csv',index_col=0)
raw = raw[meta.index]
raw.obs['Celltype1'] = meta['Celltype1']
raw.obs['Lineage'] = raw.obs['Celltype1']
raw.obs.loc[raw.obs['Celltype1'].isin(['CD4T_Mem','CD4T_Naive','CD4T_Reg','CD4T_ISG',
      'CD8T_Cytotoxic','CD8T_Mem','CD8T_Naive','MAIT']),'Lineage'] = 'T'
raw.obs.loc[raw.obs['Celltype1'].isin(['cDC','cMono','ncMono','pDC']),'Lineage'] = 'Myeloid'
raw.obs.loc[raw.obs['Celltype1'].isin(['B_Mem','B_Naive']),'Lineage'] = 'B'
sc.pp.normalize_total(raw, target_sum=1e4)
raw.obs['sample'] = raw.obs['sample'].astype('category')

for ct in raw.obs['Celltype1'].unique():
    temp = raw[raw.obs['Celltype1'] == ct]
    res = pd.DataFrame(columns=temp.var_names, index=temp.obs['sample'].cat.categories)
    for clust in temp.obs['sample'].cat.categories: 
        res.loc[clust] = temp[temp.obs['sample'].isin([clust]),:].X.mean(0)
    res.to_csv('./eQTL-input/PBMC-pseudobulk/{}.csv'.format(ct))
    
for ct in ['T','Myeloid','B']:
    temp = raw[raw.obs['Lineage'] == ct]
    res = pd.DataFrame(columns=temp.var_names, index=temp.obs['sample'].cat.categories)
    for clust in temp.obs['sample'].cat.categories: 
        res.loc[clust] = temp[temp.obs['sample'].isin([clust]),:].X.mean(0)
    res.to_csv('./eQTL-input/PBMC/PBMC-pseudobulk/{}.csv'.format(ct))

    
###___________________Pancreatic islet scRNA-seq___________________
raw = sc.read('./processing_h5ad/Islet_raw.h5ad')
meta = pd.read_csv('./processing_h5ad/Islet/Islet_harmony_35_30_1.3_meta.csv',index_col=0)
raw = raw[meta.index]
raw.obs['Celltype1'] = meta['Celltype1']
sc.pp.normalize_total(raw, target_sum=1e4)
raw.obs['sample'] = raw.obs['sample'].astype('category')

for ct in raw.obs['Celltype1'].unique():
    temp = raw[raw.obs['Celltype1'] == ct]
    res = pd.DataFrame(columns=temp.var_names, index=temp.obs['sample'].cat.categories)
    for clust in temp.obs['sample'].cat.categories: 
        res.loc[clust] = temp[temp.obs['sample'].isin([clust]),:].X.mean(0)
    res.to_csv('./eQTL-input/Islet/Islet-pseudobulk/{}.csv'.format(ct))
