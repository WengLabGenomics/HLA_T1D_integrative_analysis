library(stringr)
library(metafor)

# European cohort
HLA <- read.csv('UKB_HLA_dosage.csv',row.names = 1)
test = HLA
test = test[,colSums(test)>= (nrow(test)*2*.005)]
pca = prcomp(scale(test))
test[c('PC1','PC2')] = data.frame(pca$x)[c('PC1','PC2')]
cli <-  read.csv('UKB_cli.csv',row.names = 1)
rownames(cli) <- cli$eid
cli = cli[rownames(test),]
write.csv(test[c('sex','disease_status','PC1','PC2')],'UKB_HLA_association_covariates.csv')

name = c()
p <- c()
OR <- c()
coe <- c()
CI_L = c()
CI_H = c()
SE = c()
for( i in seq(1:(length(test)-length(c('sex','disease_status','PC1','PC2'))))){
#for( i in seq(1:4)){
    temp = test[c(colnames(test)[i],'sex','disease_status','PC1','PC2')]
    colnames(temp) <- c('loci','sex','disease_status','PC1','PC2')
    model <- glm(disease_status ~loci+sex+PC1+PC2, family = binomial(link = "logit"), data = temp)
    null_model <- glm(disease_status ~ sex+PC1+PC2, family = binomial(link = "logit"), data = temp)
    likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
    p_value <- likelihood_ratio_test$Pr[2]
    or <- exp(coef(model))['loci']
    se <- sqrt(diag(vcov(model)))['loci']
    ci_low =exp(confint(model, "loci"))[1]
    ci_high =exp(confint(model, "loci"))[2]
    name = append(name,colnames(test)[i])
    p <- append(p,p_value)
    OR <- append(OR,or)
    coe <- append(coe,coef(model)["loci"])
    CI_L = append(CI_L,ci_low)
    CI_H = append(CI_H,ci_high)
    SE = append(SE,se)
    }
HLA_res = data.frame(name,p,
                    coe,OR,CI_L,CI_H,SE)
e_res = HLA_res

# Chinese cohort
HLA = read.csv('Chinese_HLA_dosage.csv',row.names = 1)
test = HLA
test = test[,colSums(test)>= (nrow(test)*2*.005)]
pca = prcomp(scale(test))
test[c('PC1','PC2')] = data.frame(pca$x)[c('PC1','PC2')]
cli = read.csv('Chinese_cli.csv',row.names = 1)
cli = cli[rownames(test),]
test$disease_status = as.numeric(cli$disease_status)
test$sex = cli$sex
write.csv(test[c('sex','disease_status','PC1','PC2')],'Chinese_HLA_association_covariates.csv')
name = c()
p <- c()
OR <- c()
coe <- c()
CI_L = c()
CI_H = c()
SE = c()
for( i in seq(1:(length(test)-length(c('sex','disease_status','PC1','PC2'))))){
#for( i in seq(1:4)){
    temp = test[c(colnames(test)[i],'sex','disease_status','PC1','PC2')]
    colnames(temp) <- c('loci','sex','disease_status','PC1','PC2')
    model <- glm(disease_status ~loci+sex+PC1+PC2, family = binomial(link = "logit"), data = temp)
    null_model <- glm(disease_status ~ sex+PC1+PC2, family = binomial(link = "logit"), data = temp)
    likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
    p_value <- likelihood_ratio_test$Pr[2]
    or <- exp(coef(model))['loci']
    se <- sqrt(diag(vcov(model)))['loci']
    ci_low =exp(confint(model, "loci"))[1]
    ci_high =exp(confint(model, "loci"))[2]
    name = append(name,colnames(test)[i])
    p <- append(p,p_value)
    OR <- append(OR,or)
    coe <- append(coe,coef(model)["loci"])
    CI_L = append(CI_L,ci_low)
    CI_H = append(CI_H,ci_high)
    SE = append(SE,se)
    }
HLA_res = data.frame(name,p,
                    coe,OR,CI_L,CI_H,SE)
c_res =  HLA_res

#meta-analysis
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
meta_data$HLA_gene = str_split_fixed(meta_data$name,'_',3)[,2]
meta_data2 = data.frame()
for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')){
    temp = meta_data[meta_data$HLA_gene == i,]
    temp$Padj_META = p.adjust(temp$P_META, method = "bonferroni")
    meta_data2  = rbind(meta_data2,temp)
    }
meta_data2 = meta_data2[meta_data2$Padj_META < 0.05,]
write.csv(meta_data2 ,'./res/T1D-HLA-association-meta-res.csv')

# T1D-associated HLA allele LD patterns
library(pheatmap)
calc_r2 <- function(x,y){
  
  idx <- complete.cases(x,y)
  
  r <- cor(
    x[idx],
    y[idx],
    method="pearson"
  )
  
  return(r^2)
}

HLA <- read.csv('UKB_HLA_dosage.csv',row.names = 1)
res = read.csv('./res/T1D-HLA-association-meta-res.csv',row.names = 1)
hla_dosage_sig = HLA[,res$name]

alleles <- colnames(hla_dosage_sig)
ld_result <- data.frame()


for(i in 1:(length(alleles)-1)){
  
  for(j in (i+1):length(alleles)){
    
    r2 <- calc_r2(
      hla_dosage_sig[,i],
      hla_dosage_sig[,j]
    )
    
    ld_result <- rbind(
      ld_result,
      data.frame(
        allele1=alleles[i],
        allele2=alleles[j],
        r2=r2
      )
    )
  }
}


alleles <- colnames(hla_dosage_sig)

ld_matrix <- matrix(
  0,
  nrow = length(alleles),
  ncol = length(alleles),
  dimnames = list(alleles, alleles)
)

diag(ld_matrix) <- NA

for(i in 1:nrow(ld_result)){
  
  a1 <- ld_result$allele1[i]
  a2 <- ld_result$allele2[i]
  r2 <- ld_result$r2[i]
  
  ld_matrix[a1, a2] <- r2
  ld_matrix[a2, a1] <- r2
}
bk <- c(0, 0.2, 0.5, 0.8, 1)
col <- c("white", "#CADAE8", "#AAC4E2",'#8AA9D6', "#7883BA")
pheatmap(ld_matrix,cluster_rows =FALSE,cluster_cols   = FALSE,border_color = 'lightgrey',clustering_method =  'average',na_col = "#F5F5F5", color =col,breaks = bk ,filename = 'European_allele_correlation.pdf')

HLA <- read.csv('Chinese_HLA_dosage.csv',row.names = 1)
res = read.csv('./res/T1D-HLA-association-meta-res.csv',row.names = 1)
hla_dosage_sig = HLA[,res$name]

alleles <- colnames(hla_dosage_sig)
ld_result <- data.frame()


for(i in 1:(length(alleles)-1)){
  
  for(j in (i+1):length(alleles)){
    
    r2 <- calc_r2(
      hla_dosage_sig[,i],
      hla_dosage_sig[,j]
    )
    
    ld_result <- rbind(
      ld_result,
      data.frame(
        allele1=alleles[i],
        allele2=alleles[j],
        r2=r2
      )
    )
  }
}


alleles <- colnames(hla_dosage_sig)

ld_matrix <- matrix(
  0,
  nrow = length(alleles),
  ncol = length(alleles),
  dimnames = list(alleles, alleles)
)

diag(ld_matrix) <- NA

for(i in 1:nrow(ld_result)){
  
  a1 <- ld_result$allele1[i]
  a2 <- ld_result$allele2[i]
  r2 <- ld_result$r2[i]
  
  ld_matrix[a1, a2] <- r2
  ld_matrix[a2, a1] <- r2
}

bk <- c(0,0.2, 0.5, 0.8, 1)
col <- c("white", "#CC85B1", "#B280B2",'#9270A5', "#755A91")
pheatmap(ld_matrix,cluster_rows =FALSE,cluster_cols   = FALSE,border_color = 'lightgrey',clustering_method =  'average',na_col = "#F5F5F5", color =col,breaks = bk ,filename = 'Chinese_allele_correlation.pdf')