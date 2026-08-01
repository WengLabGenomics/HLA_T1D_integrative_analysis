
##---------------------step1-------------------------------------
magma_dir=~/software/magma_dir
${magma_dir}/magma \
    --annotate window=1000,1000 \ # 1000kb
    --snp-loc ~/software/MHC-ref/SNP2HLA_package_v1.0.3/SNP2HLA/HM_CEU_REF.bim \
    --gene-loc ${magma_dir}/NCBI36.3.gene.loc \
    --out out/step1
    
##---------------------step2-------------------------------------
${magma_dir}/magma \
    --bfile ~/software/MHC-ref/SNP2HLA_package_v1.0.3/SNP2HLA/HM_CEU_REF \
    --pval T1D_HLA_summary.pval use='SNP,P' ncol='N' \
    --gene-annot out/step1.genes.annot \
    --out out/step2/T1D_HLA

##---------------------step3-------------------------------------   
scdrs munge-gs \
    --out-file ./out/step3/T1D_HLA_gene.gs \
    --zscore-file ./out/step2/T1D_HLA.z_file.tsv \
    --weight zscore \
    --n-max 1000

##---------------------step4------------------------------------- 
# PBMC
scdrs compute-score \
    --h5ad-file PBMC_raw_QC.h5ad \
    --h5ad-species human \
    --gs-file ./out/step3/T1D_HLA_gene.gs \
    --gs-species human \
    --flag-filter-data True \
    --flag-raw-count True \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder out/step4-PBMC

# Pancreatic islet
scdrs compute-score \
    --h5ad-file Islet_raw_QC.h5ad \
    --h5ad-species human \
    --gs-file ./out/step3/T1D_HLA_gene.gs \
    --gs-species human \
    --flag-filter-data True \
    --flag-raw-count True \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder out/step4-Islet
    
##---------------------step5------------------------------------- 
# PBMC
scdrs perform-downstream \
        --h5ad-file PBMC_raw_QC.h5ad  \
        --score-file ./out/step4-PBMC/T1D.full_score.gz \
        --out-folder ./out/step5-PBMC/ \
        --group_analysis celltype \
        --flag-filter-data True \
        --flag-raw-count True

# Pancreatic islet        
scdrs perform-downstream \
        --h5ad-file Islet_raw_QC.h5ad  \
        --score-file ./out/step4-Islet/T1D.full_score.gz \
        --out-folder ./out/step5-Islet/ \
        --group_analysis celltype \
        --flag-filter-data True \
        --flag-raw-count True