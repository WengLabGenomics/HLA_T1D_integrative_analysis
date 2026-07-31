library(dplyr)
library(data.table)
library(stringr)

####______________________PBMC________________________
cts = c( 'NK','B_Mem','B_Naive','CD4T_Mem','CD4T_Naive','CD4T_ISG','CD4T_Reg',
      'CD8T_Cytotoxic','CD8T_Mem','CD8T_Naive','MAIT','cDC','cMono','ncMono','pDC')
#cts = c('B','T','Myeloid')

geno = data.frame(t(read.csv('./eQTL-input/PBMC/HLA_dosage_for_eQTL.csv',row.names = 1)))
colnames(geno) <-  gsub("\\.", "-", colnames(geno) )
cli = read.csv('./eQTL-input/PBMC/cli_for_QTL.csv',row.names = 1)
ref_clean = read.csv('./eQTL-input/gene_ref_clean.csv',row.names = 1)
loci = c('A','B','C','DRB1','DQA1','DQB1','DPB1','DPA1')
pos = c(29944214,31355516,31270415,32584287,32640529,32663284,33081591,33069309)
df = data.frame('loci'=loci,pos = pos)
x= data.frame(str_split_fixed(rownames(geno),'_',3))
colnames(x) <- c('X1','loci','X3')
rownames(x) = rownames(geno)
library(dplyr)
x = left_join(x,df,by = 'loci')
rownames(x) = rownames(geno)
x$chr = 'chr6'
snploc = x[c('chr','pos')]
for (ct in cts){
data = data.frame(fread(paste0('./eQTL-input/PBMC/PBMC-pseudobulk/',ct,'.csv'),header = FALSE))
var_names = data.frame(t(data[1,]))$X[2:length(data.frame(t(data[1,]))$X)]
data = data.frame(fread(paste0('./eQTL-input/PBMC/PBMC-pseudobulk/',ct,'.csv')))
rownames(data) <- data$V1
data = data[,2:ncol(data)]
colnames(data) = var_names
keep = colSums(data > 0 ) >= nrow(data)*0.1
print(paste0(ct,':',dim(data[,keep])[2]))
data = data[,keep]
data  = log2(data+1)

re = intersect(rownames(cli),rownames(data))
re = intersect(re,colnames(geno))
geno_re = geno[,re]
data_re = data.frame(t(data[re,]))
colnames(data_re) <- re
gene_loc = ref_clean[colnames(data),]
cli_re = cli[re,]
cli_re = data.frame(t(cli_re))
colnames(cli_re) = re 

folder_path <- paste0("./eQTL-input/PBMC/PBMC-pseudobulk/",ct)
if (!file.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created successfully."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

folder_path <- paste0("./eQTL-input/PBMC/PBMC-pseudobulk/",ct)
write.table(data_re,paste0(folder_path,'/GE.txt'),quote = FALSE)
write.table(geno_re,paste0(folder_path,'/SNP.txt'),quote = FALSE)
write.table(snploc,paste0(folder_path,'/snpsloc.txt'),quote = FALSE)
write.table(cli_re,paste0(folder_path,'/Covariates_cli.txt'),quote = FALSE)
write.table(gene_loc,paste0(folder_path,'/geneloc.txt'),quote = FALSE) 
    }

####______________________Pancreatic islet________________________
cts = c('Acinar','Activated stellate','Alpha','Beta','CFTR+ ductal','Delta','Endothelial','Gamma',
       'Immune cell','MUC5B+ ductal','Proliferative alpha','Quiescent stellate','Schwann')
geno = data.frame(t(read.csv('./eQTL-input/Islet/HLA_dosage_for_eQTL.csv',row.names = 1)))
colnames(geno) <-  gsub("\\.", "-", colnames(geno) )
cli = read.csv('./eQTL-input/Islet/cli_for_QTL.csv',row.names = 1)
ref_clean = read.csv('./eQTL-input/gene_ref_clean.csv',row.names = 1)
for (ct in cts){
data = data.frame(fread(paste0('./eQTL-input/Islet/Islet-pseudobulk/',ct,'.csv'),header = FALSE))
var_names = data.frame(t(data[1,]))$X[2:length(data.frame(t(data[1,]))$X)]
data = data.frame(fread(paste0('./eQTL-input/Islet/Islet-pseudobulk/',ct,'.csv')))
rownames(data) <- data$V1
data = data[,2:ncol(data)]
colnames(data) = var_names
keep = colSums(data > 0 ) >= nrow(data)*0.1
print(paste0(ct,':',dim(data[,keep])[2]))
data = data[,keep]
data  = log2(data+1)

re = intersect(rownames(cli),rownames(data))
re = intersect(re,colnames(geno))
geno_re = geno[,re]
data_re = data.frame(t(data[re,]))
colnames(data_re) <- re
gene_loc = ref_clean[colnames(data),]
cli_re = cli[re,]
cli_re = data.frame(t(cli_re))
colnames(cli_re) = re 

folder_path <- paste0("./eQTL-input/Islet/Islet-pseudobulk/",ct)
if (!file.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created successfully."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

folder_path <- paste0("./eQTL-input/Islet/Islet-pseudobulk/",ct)
write.table(data_re,paste0(folder_path,'/GE.txt'),quote = FALSE)
write.table(geno_re,paste0(folder_path,'/SNP.txt'),quote = FALSE)
write.table(snploc,paste0(folder_path,'/snpsloc.txt'),quote = FALSE)
write.table(cli_re,paste0(folder_path,'/Covariates_cli.txt'),quote = FALSE)
write.table(gene_loc,paste0(folder_path,'/geneloc.txt'),quote = FALSE) 
    }
