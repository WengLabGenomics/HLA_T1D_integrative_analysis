library(data.table)
# Chinese cohort
hap1 = fread('Chinese-DR-DQ-hap.csv')
r_col  = colnames(hap1)
hap1 = data.frame(hap1)
colnames(hap1) = r_col 
hap2 =  fread('Chinese-DP-hap.csv')
r_col  = colnames(hap2)
hap2 = data.frame(hap2)
colnames(hap2) = r_col 

rownames(hap1) <- hap1$sample
rownames(hap2) <- hap2$sample

re = intersect(rownames(hap1),rownames(hap2))

hap = cbind(hap1[re,],hap2[re,])

res_hap = read.csv('./res/T1D-haplotype-association-meta-res.csv',row.names = 1)
hap_df = hap[,res_hap$name]

HLA <- read.csv('./Chinese_HLA_dosage.csv',row.names = 1)
res = read.csv('./res/T1D-HLA-association-meta-res.csv',row.names = 1)
res = res[res$HLA_gene %in% c('A','B','C'),]
HLA = HLA[c(res$name)]
df = cbind(HLA[rownames(hap),],hap)
cli = read.csv('./Chinese_HLA_association_covariates.csv',row.names = 1)
df = cbind(df,cli[rownames(df),])

bt <- function(x) paste0("`", x, "`")
name = c()
p <- c()
coe <- c()
OR <- c()
CI_L = c()
CI_H = c()
SE = c()
for(i in res$name){
#for( i in seq(1:4)){
    temp = df[c(i,'sex','disease_status','PC1','PC2',colnames(hap_df))]
    colnames(temp) <- c('loci','sex','disease_status','PC1','PC2',colnames(hap_df))
    ff <- as.formula(paste( "disease_status ~", "loci", "+ sex + PC1 + PC2 + ", paste(bt(colnames(hap_df)), collapse = " + ")))
    nf <- as.formula(paste( "disease_status ~", "loci", "+ sex + PC1 + PC2"))
    model <- glm(ff, family = binomial(link = "logit"), data = temp)
      
  
  model_sum <- summary(model)$coefficients
  
  p_value <- model_sum["loci", "Pr(>|z|)"]
  beta <- model_sum["loci", "Estimate"]
  se <- model_sum["loci", "Std. Error"]
  
  or <- exp(beta)
  ci_low <- exp(beta - 1.96 * se)
  ci_high <- exp(beta + 1.96 * se)
  
  name <- append(name, i)
  p <- append(p, p_value)
  OR <- append(OR, or)
  coe <- append(coe, beta)
  CI_L <- append(CI_L, ci_low)
  CI_H <- append(CI_H, ci_high)
  SE <- append(SE, se)
    }
HLA_res = data.frame(name,p,
                    coe,OR,CI_L,CI_H,SE)

c_res = HLA_res 

# European cohort
hap1 = fread('European-DR-DQ-hap.csv')
r_col  = colnames(hap1)
hap1 = data.frame(hap1)
colnames(hap1) = r_col 
hap2 =  fread('European-DP-hap.csv')
r_col  = colnames(hap2)
hap2 = data.frame(hap2)
colnames(hap2) = r_col 

rownames(hap1) <- hap1$sample
rownames(hap2) <- hap2$sample

re = intersect(rownames(hap1),rownames(hap2))

hap = cbind(hap1[re,],hap2[re,])

res_hap = read.csv('./res/T1D-haplotype-association-meta-res.csv',row.names = 1)
hap_df = hap[,res_hap$name]

HLA <- read.csv('./UKB_HLA_dosage.csv',row.names = 1)
res = read.csv('./res/T1D-HLA-association-meta-res.csv',row.names = 1)
res = res[res$HLA_gene %in% c('A','B','C'),]
HLA = HLA[c(res$name)]
df = cbind(HLA[rownames(hap),],hap)
cli = read.csv('./UKB_HLA_association_covariates.csv',row.names = 1)
df = cbind(df,cli[rownames(df),])

bt <- function(x) paste0("`", x, "`")
name = c()
p <- c()
coe <- c()
OR <- c()
CI_L = c()
CI_H = c()
SE = c()
for(i in res$name){
#for( i in seq(1:4)){
    temp = df[c(i,'sex','disease_status','PC1','PC2',colnames(hap_df))]
    colnames(temp) <- c('loci','sex','disease_status','PC1','PC2',colnames(hap_df))
    ff <- as.formula(paste( "disease_status ~", "loci", "+ sex + PC1 + PC2 + ", paste(bt(colnames(hap_df)), collapse = " + ")))
    nf <- as.formula(paste( "disease_status ~", "loci", "+ sex + PC1 + PC2"))
    model <- glm(ff, family = binomial(link = "logit"), data = temp)
      
  
  model_sum <- summary(model)$coefficients
  
  p_value <- model_sum["loci", "Pr(>|z|)"]
  beta <- model_sum["loci", "Estimate"]
  se <- model_sum["loci", "Std. Error"]
  
  or <- exp(beta)
  ci_low <- exp(beta - 1.96 * se)
  ci_high <- exp(beta + 1.96 * se)
  
  name <- append(name, i)
  p <- append(p, p_value)
  OR <- append(OR, or)
  coe <- append(coe, beta)
  CI_L <- append(CI_L, ci_low)
  CI_H <- append(CI_H, ci_high)
  SE <- append(SE, se)
    }
HLA_res = data.frame(name,p,
                    coe,OR,CI_L,CI_H,SE)

e_res = HLA_res

# meta-analysis
library(metafor)
meta_data <- merge(c_res, e_res, by = "name", suffixes = c("_CHI", "_EUR"))

meta_data$Beta_META <- NA
meta_data$SE_META <- NA
meta_data$P_META <- NA
meta_data$CI_Lower <- NA
meta_data$CI_Upper <- NA
meta_data$Q_META <- NA
meta_data$Q_pval <- NA

for (i in 1:nrow(meta_data)) {
  betas <- c(meta_data$coe_CHI[i], meta_data$coe_EUR[i])
  ses <- c(meta_data$SE_CHI[i], meta_data$SE_EUR[i])

  tryCatch({
    res <- rma(yi = betas, sei = ses, method = "FE")
    
    meta_data$Beta_META[i] <- res$b
    meta_data$SE_META[i] <- res$se
    meta_data$P_META[i] <- res$pval
    meta_data$CI_Lower[i] <- res$ci.lb
    meta_data$CI_Upper[i] <- res$ci.ub
    meta_data$Q_META[i] <- res$QE
    meta_data$Q_pval[i] <- res$QEp
  }, error = function(e) {
    meta_data$Beta_META[i] <- NA
  })
}

meta_data = meta_data[meta_data$P_META < 0.05,]