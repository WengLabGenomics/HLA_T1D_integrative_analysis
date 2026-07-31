arg=commandArgs(T)
date = arg[1] # PBMC or Islet
ct = arg[2] # celltypes in PBMC or pancreatic islet

library(MatrixEQTL)
base.dir = paste0('./eQTL-input/',date,'/',date,'-pseudobulk/')
save.dir = paste0('./eQTL-res/',date,'/')

cov = read.table(paste0(base.dir,'/',ct, "/Covariates_cli.txt"))
re_cov = c('sex','age','disease','race')
cov = cov[re_cov,]
write.table(cov,paste0(base.dir,'/',ct, "/Covariates.txt"),quote = FALSE)
        
useModel = modelLINEAR; 
        
# Genotype file name
SNP_file_name = paste0(base.dir,'/',ct, "/SNP.txt");
snps_location_file_name = paste0(base.dir,'/',ct, "/snpsloc.txt");

# Gene expression file name
expression_file_name = paste0(base.dir,'/',ct, "/GE.txt");
gene_location_file_name = paste0(base.dir,'/',ct, "/geneloc.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir,'/',ct,"/Covariates.txt");
        
# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
        
## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";    
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name,header = TRUE, stringsAsFactors = FALSE);
snpspos$snpid = rownames(snpspos)
snpspos = snpspos[c('snpid','chr','pos')]
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos$geneid = rownames(genepos)
genepos = genepos[c('geneid','chr','left','right')]


me = Matrix_eQTL_main(
        snps = snps, 
        gene = gene, 
        cvrt = cvrt,
        output_file_name  = output_file_name_tra,
        pvOutputThreshold = pvOutputThreshold_tra,
        useModel = useModel, 
        errorCovariance = errorCovariance, 
        verbose = TRUE, 
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos, 
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);
cis = me$cis$eqtls
cis = cis[cis$FDR < 0.05,]
write.csv(cis,paste0(save.dir,'/',ct,'_cis_fdr005_all.csv'))
