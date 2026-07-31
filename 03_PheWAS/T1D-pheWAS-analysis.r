library(data.table)


# allele level
## Chinese cohort
HLA <- read.csv('Chinese_HLA_dosage.csv',row.names = 1)
cli <- read.csv('Chinese-T1D-phenotype.csv',row.names = 1)

re = intersect(rownames(HLA),rownames(cli))

HLA = HLA[re,]
cli = cli[re,]

test1 = HLA
test1 = test1[,colSums(test1)>= nrow(test1)*2*.005]
pca = prcomp(test1)
test1[c('PC1','PC2','PC3')] = data.frame(pca$x)[c('PC1','PC2','PC3')]
merge_df2 = cbind(test1,cli)

line = c('PLT','MON','MCHC','BAS','EOS','MCV',
                  'LYM','HCT','RBC','HGB','NEU','HDL','LDL','TG','CHOL','SCr',
                  'UA','UREA','IBIL','TBIL','DBIL','TP','ALB',
                   'ALT','GLOB','GCT','AST','AG','GLU','HbA1c','TPOAb',
                   'ACR','ALb','Cr','onset_age')

library(stringr)
m = data.frame(str_split_fixed(colnames(test1),'_',3))
m$loci = colnames(test1)

library(RNOmni)
HLA_res_all_phe <- data.frame()

for (y in line){
HLA_res_all <- data.frame()
    for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')){
        test = test1[m[m$X2 == i,'loci']]
        name = c()
        p <- c()
        coe <- c()
        CI_L = c()
        CI_H = c()  
        SE = c()
        ff = as.formula(paste0(y,' ~ loci+sex+PC1+PC2'))
        nf = as.formula(paste0(y,' ~ sex+PC1+PC2'))
    
        for (j in seq(1,length(colnames(test)))){
            temp = merge_df2[c(colnames(test)[j],c('PC1','PC2'),colnames(cli))]
            colnames(temp) <- c('loci',c('PC1','PC2'),colnames(cli))
            temp = temp[!is.na(temp[y]),]
            mean_value <- mean(temp[,y], na.rm = TRUE)
            sd_value <- sd(temp[,y], na.rm = TRUE)
            lower_bound <- mean_value - 5 * sd_value
            upper_bound <- mean_value + 5 * sd_value
            temp = temp[temp[y] > lower_bound &  temp[y] <upper_bound, ]
            temp[y] =  RankNorm(temp[,y])
            model <- lm(ff, data = temp)
            null_model <- lm(nf, data = temp)
            likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
            p_value <- likelihood_ratio_test$Pr[2]
            ci_low =confint(model, "loci")[1]
            ci_high =confint(model, "loci")[2]
            se <- sqrt(diag(vcov(model)))['loci']
            name = append(name,colnames(test)[j])
            p <- append(p,p_value)
            coe <- append(coe,coef(model)["loci"])
            CI_L = append(CI_L,ci_low)
            CI_H = append(CI_H,ci_high)
            SE = append(SE,se)
    }
        HLA_res = data.frame(name,p,coe,CI_L,CI_H,SE)
        HLA_res$HLA_gene = i
        HLA_res_all = rbind(HLA_res_all,HLA_res)
    }
    print(paste0(y,':',length(temp[,y])))
    HLA_res_all$trial = y
    HLA_res_all_phe = rbind(HLA_res_all_phe,HLA_res_all)
    }

c_res = HLA_res_all_phe

## European cohort
HLA = read.csv('./UKB_HLA_dosage.csv',row.names = 1)
cli <- read.csv('UKB-T1D-phenotype.csv',row.names = 1)

re = intersect(rownames(HLA),rownames(cli))

HLA = HLA[re,]
cli = cli[re,]

test1 = HLA
test1 = test1[,colSums(test1)>= nrow(test1)*2*.005]
pca = prcomp(test1)
test1[c('PC1','PC2')] = data.frame(pca$x)[c('PC1','PC2')]
merge_df2 = cbind(test1,cli)

line <- c('Erythrocyte_Count','Haemoglobin_concentration','Mean_corpuscular_volume',
        'Mean_corpuscular_haemoglobin_concentration','Platelet_count','Neutrophill_count','Eosinophill_count',
         'Lymphocyte_percentage','Monocyte_percentage','Basophill_percentage','Albumin','Alanine_aminotransferase','Aspartate_aminotransferase',
             'Direct_bilirubin','Urea','Cholesterol','Serum_creatinine','Gamma_glutamyltransferase','Glucose',
         'HbA1c','HDLC','LDLC','Total_bilirubin','Total_protein','Triglycerides','Age_diabetes_diagnosed','Urinary_albumin_to_creatinine_ratio'
        )

library(stringr)
m = data.frame(str_split_fixed(colnames(test1),'_',3))
m$loci = colnames(test1)

library(RNOmni)
HLA_res_all_phe <- data.frame()

for (y in line){
HLA_res_all <- data.frame()
    for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')){
        test = test1[m[m$X2 == i,'loci']]
        name = c()
        p <- c()
        coe <- c()
        CI_L = c()
        CI_H = c()  
        SE = c()
        ff = as.formula(paste0(y,' ~ loci+sex+PC1+PC2'))
        nf = as.formula(paste0(y,' ~ sex+PC1+PC2'))
    
        for (j in seq(1,length(colnames(test)))){
            temp = merge_df2[c(colnames(test)[j],c('PC1','PC2'),colnames(cli))]
            colnames(temp) <- c('loci',c('PC1','PC2'),colnames(cli))
            temp = temp[!is.na(temp[y]),]
            mean_value <- mean(temp[,y], na.rm = TRUE)
            sd_value <- sd(temp[,y], na.rm = TRUE)
            lower_bound <- mean_value - 5 * sd_value
            upper_bound <- mean_value + 5 * sd_value
            temp = temp[temp[y] > lower_bound &  temp[y] <upper_bound, ]
            temp[y] =  RankNorm(temp[,y])
            model <- lm(ff, data = temp)
            null_model <- lm(nf, data = temp)
            likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
            p_value <- likelihood_ratio_test$Pr[2]
            ci_low =confint(model, "loci")[1]
            ci_high =confint(model, "loci")[2]
            se <- sqrt(diag(vcov(model)))['loci']
            name = append(name,colnames(test)[j])
            p <- append(p,p_value)
            coe <- append(coe,coef(model)["loci"])
            CI_L = append(CI_L,ci_low)
            CI_H = append(CI_H,ci_high)
            SE = append(SE,se)
    }
        HLA_res = data.frame(name,p,coe,CI_L,CI_H,SE)
        HLA_res$HLA_gene = i
        HLA_res_all = rbind(HLA_res_all,HLA_res)
    }
    print(paste0(y,':',length(temp[,y])))
    HLA_res_all$trial = y
    HLA_res_all_phe = rbind(HLA_res_all_phe,HLA_res_all)
    }

e_res = HLA_res_all_phe

## meta-analysis
w_re <- c('Erythrocyte_Count','Haemoglobin_concentration','Mean_corpuscular_volume',
        'Mean_corpuscular_haemoglobin_concentration','Platelet_count','Neutrophill_count','Eosinophill_count',
         'Lymphocyte_percentage','Monocyte_percentage','Basophill_percentage','Albumin','Alanine_aminotransferase','Aspartate_aminotransferase',
             'Direct_bilirubin','Urea','Cholesterol','Serum_creatinine','Gamma_glutamyltransferase','Glucose',
         'HbA1c','HDLC','LDLC','Total_bilirubin','Total_protein','Triglycerides','Age_diabetes_diagnosed','Urinary_albumin_to_creatinine_ratio'
        )
c_re <- c('RBC','HGB','MCV',
  'MCHC','PLT','NEU','EOS',
  'LYM','MON','BAS','ALB','ALT','AST',
  'DBIL','UREA','CHOL','SCr','GCT','GLU',
  'HbA1c','HDL','LDL','TBIL','TP',
  'TG','onset_age','ACR'
 )
uni <- c('Erythrocyte_Count','Haemoglobin_concentration','Mean_corpuscular_volume',
        'Mean_corpuscular_haemoglobin_concentration','Platelet_count','Neutrophill_count','Eosinophill_count',
         'Lymphocyte_percentage','Monocyte_percentage','Basophill_percentage','Albumin','Alanine_aminotransferase','Aspartate_aminotransferase',
             'Direct_bilirubin','Urea','Cholesterol','Serum_creatinine','Gamma_glutamyltransferase','Glucose',
         'HbA1c','HDL-C','LDL-C','Total_bilirubin','Total_protein','Triglycerides','Age_at_diagnosis','Urinary_albumin_to_creatinine_ratio'
        )

e_res <- e_res[e_res$trial %in% w_re,]
c_res <- c_res[c_res$trial %in% c_re,]

library(dplyr)
temp1 <- data.frame(trial = w_re,trial2 = uni)
temp2 <- data.frame(trial = c_re,trial2 = uni)
e_res2 = left_join(e_res,temp1,by = 'trial')
c_res2 = left_join(c_res,temp2,by = 'trial')

res_df <- read.csv('./res/T1D-HLA-association-meta-res.csv',row.names = 1)
e_res2 = e_res2[e_res2$name %in% res_df$name,]
c_res2 = c_res2[c_res2$name %in% res_df$name,]

library(metafor)

meta_data_all = data.frame()
for (y in uni){
#y = 'Neutrophill_count' 
e_res3 = e_res2[e_res2$trial2 == y,]
c_res3 = c_res2[c_res2$trial2 == y,]
meta_data <- merge(c_res3, e_res3, by = "name", suffixes = c("_CHI", "_EUR"))

# 初始化输出列
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
  
  # 有时数据质量不好，要加 tryCatch 避免报错中断
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
    # 如果报错就留空
    meta_data$Beta_META[i] <- NA
  })
}
    meta_data$gene = str_split_fixed(meta_data$name,'_',3)[,2]
    meta_data2 = data.frame()
for (i in c('A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1')){
    temp = meta_data[meta_data$gene == i,]
    temp$Padj_META = p.adjust(temp$P_META, method = "bonferroni")
    meta_data2  = rbind(meta_data2,temp)
    }
    meta_data_all = rbind(meta_data_all,meta_data2)
    }

write.csv(meta_data_all[meta_data_all$Padj_META < 0.05,],'./res/T1D-pheWAS-meta-res.csv')

# haplotype level
## Chinese cohort
line = c('PLT','MON','MCHC','BAS','EOS','MCV',
                  'LYM','HCT','RBC','HGB','NEU','HDL','LDL','TG','CHOL','SCr',
                  'UA','UREA','IBIL','TBIL','DBIL','TP','ALB',
                   'ALT','GLOB','GCT','AST','AG','GLU','HbA1c','TPOAb',
                   'ACR','ALb','Cr','onset_age')
library(RNOmni)
HLA_res_all_phe <- data.frame()
cli <- read.csv('Chinese-T1D-phenotype.csv',row.names = 1)
covariates = read.csv('./Chinese_HLA_association_covariates.csv',row.names = 1)
for (y in line){
HLA_res_all <- data.frame()
    for (hap in c('DR-DQ','DP')){
        df = fread(paste0('Chinese-',hap,'-hap.csv'))
        hp = colnames(df)[3:ncol(df)]
        df = data.frame(df)
        rownames(df) = df$sample
        df = df[,3:ncol(df)]
        colnames(df) = hp
        re = intersect(rownames(df),rownames(cli))
        df = df[re,]
        cli = cli[re,]
        covariates = covariates[re,]
        merge_df2 = cbind(df,cbind(covariates[c('PC1','PC2')],cli))
        test = df
        name = c();p <- c();coe <- c();CI_L = c();CI_H = c() ;SE = c()
        ff = as.formula(paste0(y,' ~ loci+sex+PC1+PC2'))
        nf = as.formula(paste0(y,' ~ sex+PC1+PC2'))
        for (j in seq(1,length(colnames(test)))){
            temp = merge_df2[c(colnames(test)[j],c('PC1','PC2'),colnames(cli))]
            colnames(temp) <- c('loci',c('PC1','PC2'),colnames(cli))
            temp = temp[!is.na(temp[y]),]
            mean_value <- mean(temp[,y], na.rm = TRUE)
            sd_value <- sd(temp[,y], na.rm = TRUE)
            lower_bound <- mean_value - 5 * sd_value
            upper_bound <- mean_value + 5 * sd_value
            temp = temp[temp[y] > lower_bound &  temp[y] <upper_bound, ]
            temp[y] =  RankNorm(temp[,y])
            model <- lm(ff, data = temp)
            null_model <- lm(nf, data = temp)
            likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
            p_value <- likelihood_ratio_test$Pr[2]
             ci_low =confint(model, "loci")[1]
            ci_high =confint(model, "loci")[2]
             se <- sqrt(diag(vcov(model)))['loci']
            name = append(name,colnames(test)[j])
            p <- append(p,p_value)
            coe <- append(coe,coef(model)["loci"])
            CI_L = append(CI_L,ci_low)
            CI_H = append(CI_H,ci_high)
             SE = append(SE,se)
            }
        
        Padj = p.adjust(p, method = "bonferroni")
        HLA_res = data.frame(name,p,Padj,coe,CI_L,CI_H,SE)
        HLA_res$loci = hap
        HLA_res_all = rbind(HLA_res_all,HLA_res)
        }
    print(paste0(y,':',length(temp[,y])))
    HLA_res_all$trial = y
    HLA_res_all_phe = rbind(HLA_res_all_phe,HLA_res_all)
    }
c_res = HLA_res_all_phe
## European cohort
line <- c('Erythrocyte_Count','Haemoglobin_concentration','Mean_corpuscular_volume','Mean_corpuscular_haemoglobin_concentration','Platelet_count','Neutrophill_count','Eosinophill_count','Lymphocyte_percentage','Monocyte_percentage','Basophill_percentage','Albumin','Alanine_aminotransferase','Aspartate_aminotransferase','Direct_bilirubin','Urea','Cholesterol','Serum_creatinine','Gamma_glutamyltransferase','Glucose','HbA1c','HDLC','LDLC','Total_bilirubin','Total_protein','Triglycerides','Age_diabetes_diagnosed','Urinary_albumin_to_creatinine_ratio')
HLA_res_all_phe <- data.frame()
cli <- read.csv('UKB-T1D-phenotype.csv',row.names = 1)
covariates = read.csv('UKB_HLA_association_covariates.csv',row.names = 1)
for (y in line){
HLA_res_all <- data.frame()
    for (hap in c('DR-DQ','DP')){
        df = fread(paste0('European-',hap,'-hap.csv'))
        hp = colnames(df)[3:ncol(df)]
        df = data.frame(df)
        rownames(df) = df$sample
        df = df[,3:ncol(df)]
        colnames(df) = hp
        re = intersect(rownames(df),rownames(cli))
        df = df[re,]
        cli = cli[re,]
        covariates = covariates[re,]
        merge_df2 = cbind(df,cbind(covariates[c('PC1','PC2')],cli))
        test = df
        name = c();p <- c();coe <- c();CI_L = c();CI_H = c() ;SE = c()
        ff = as.formula(paste0(y,' ~ loci+sex+PC1+PC2'))
        nf = as.formula(paste0(y,' ~ sex+PC1+PC2'))
        for (j in seq(1,length(colnames(test)))){
            temp = merge_df2[c(colnames(test)[j],c('PC1','PC2'),colnames(cli))]
            colnames(temp) <- c('loci',c('PC1','PC2'),colnames(cli))
            temp = temp[!is.na(temp[y]),]
            mean_value <- mean(temp[,y], na.rm = TRUE)
            sd_value <- sd(temp[,y], na.rm = TRUE)
            lower_bound <- mean_value - 5 * sd_value
            upper_bound <- mean_value + 5 * sd_value
            temp = temp[temp[y] > lower_bound &  temp[y] <upper_bound, ]
            temp[y] =  RankNorm(temp[,y])
            model <- lm(ff, data = temp)
            null_model <- lm(nf, data = temp)
            likelihood_ratio_test <- anova(null_model,model, test = "Chisq")
            p_value <- likelihood_ratio_test$Pr[2]
             ci_low =confint(model, "loci")[1]
            ci_high =confint(model, "loci")[2]
             se <- sqrt(diag(vcov(model)))['loci']
            name = append(name,colnames(test)[j])
            p <- append(p,p_value)
            coe <- append(coe,coef(model)["loci"])
            CI_L = append(CI_L,ci_low)
            CI_H = append(CI_H,ci_high)
             SE = append(SE,se)
            }
        
        Padj = p.adjust(p, method = "bonferroni")
        HLA_res = data.frame(name,p,Padj,coe,CI_L,CI_H,SE)
        HLA_res$loci = hap
        HLA_res_all = rbind(HLA_res_all,HLA_res)
        }
    print(paste0(y,':',length(temp[,y])))
    HLA_res_all$trial = y
    HLA_res_all_phe = rbind(HLA_res_all_phe,HLA_res_all)
    }
e_res =  HLA_res_all_phe

library(stringr)

w_re <- c('Erythrocyte_Count','Haemoglobin_concentration','Mean_corpuscular_volume','Mean_corpuscular_haemoglobin_concentration','Platelet_count','Neutrophill_count','Eosinophill_count','Lymphocyte_percentage','Monocyte_percentage','Basophill_percentage','Albumin','Alanine_aminotransferase','Aspartate_aminotransferase', 'Direct_bilirubin','Urea','Cholesterol','Serum_creatinine','Gamma_glutamyltransferase','Glucose', 'HbA1c','HDLC','LDLC','Total_bilirubin','Total_protein','Triglycerides','Age_diabetes_diagnosed','Urinary_albumin_to_creatinine_ratio')
c_re <- c('RBC','HGB','MCV',
  'MCHC','PLT','NEU','EOS',
  'LYM','MON','BAS','ALB','ALT','AST',
  'DBIL','UREA','CHOL','SCr','GCT','GLU',
  'HbA1c','HDL','LDL','TBIL','TP',
  'TG','onset_age','ACR')
uni <- c('Erythrocyte_Count','Haemoglobin_concentration','Mean_corpuscular_volume','Mean_corpuscular_haemoglobin_concentration','Platelet_count','Neutrophill_count','Eosinophill_count','Lymphocyte_percentage','Monocyte_percentage','Basophill_percentage','Albumin','Alanine_aminotransferase','Aspartate_aminotransferase','Direct_bilirubin','Urea','Cholesterol','Serum_creatinine','Gamma_glutamyltransferase','Glucose','HbA1c','HDL-C','LDL-C','Total_bilirubin','Total_protein','Triglycerides','Age_at_diagnosis','Urinary_albumin_to_creatinine_ratio')

e_res <- e_res[e_res$trial %in% w_re,]
c_res <- c_res[c_res$trial %in% c_re,]
library(dplyr)
temp1 <- data.frame(trial = w_re,trial2 = uni)
temp2 <- data.frame(trial = c_re,trial2 = uni)
e_res2 = left_join(e_res,temp1,by = 'trial')
c_res2 = left_join(c_res,temp2,by = 'trial')
c_res2 = c_res2[!is.na(c_res2$p),]
e_res2 = e_res2[!is.na(e_res2$p),]

meta_data_all = data.frame()
for (y in uni){
e_res3 = e_res2[e_res2$trial2 == y,]
c_res3 = c_res2[c_res2$trial2 == y,]
meta_data <- merge(c_res3, e_res3 , by = "name", suffixes = c("_CHI", "_EUR"))

# 初始化输出列
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
  
  # 有时数据质量不好，要加 tryCatch 避免报错中断
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
    # 如果报错就留空
    meta_data$Beta_META[i] <- NA
  })
}
    meta_data$gene = str_split_fixed(meta_data$name,'_',3)[,2]
    meta_data2 = data.frame()
for (i in c('DR-DQ','DP')){
    temp = meta_data[meta_data$loci_EUR == i,]
    temp$Padj_META = p.adjust(temp$P_META, method = "bonferroni")
    meta_data2  = rbind(meta_data2,temp)
    }
    meta_data_all = rbind(meta_data_all,meta_data2)
    }

write.csv(meta_data_all[meta_data_all$Padj_META < 0.05,],'./res/T1D-pheWAS-haplotype-meta-res.csv')


