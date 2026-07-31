library(stringr)
library(haplo.stats)
library(dplyr)
library(tidyr)
library(data.table)


# Obtaining haplotypes using the EM algorithm
run_hla_em <- function(
    df,
    loci,
    id = NULL,
    min_freq = 0.01,
    seed = 123
) {
  

  if (is.null(id)) {
    if (is.null(rownames(df))) {
      sample_id <- as.character(seq_len(nrow(df)))
    } else {
      sample_id <- rownames(df)
    }
  } else {
    sample_id <- as.character(df[[id]])
  }
  
    
  allele_cols <- unlist(
    lapply(
      loci,
      function(x) {
        paste0(x, c(".Allele1", ".Allele2"))
      }
    )
  )
  
  missing_cols <- setdiff(allele_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "missing data:",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  geno <- df[, allele_cols, drop = FALSE]

  geno[] <- lapply(
    geno,
    function(x) {
      x <- trimws(as.character(x))
      
      x[x %in% c("", ".", "NA", "N/A", "0")] <- NA
      
      sub("^.+\\*", "", x)
    }
  )
  
  geno <- as.matrix(geno)
  

  keep <- rowSums(!is.na(geno)) > 0
  
  geno <- geno[keep, , drop = FALSE]
  sample_id <- sample_id[keep]
  

  set.seed(seed)
  
  em_control <- haplo.em.control(
    insert.batch.size = length(loci),
    min.posterior = 0,
    tol = 1e-8,
    max.iter = 5000,
    random.start = 1,
    n.try = 20,
    iseed = seed,
    verbose = 0
  )
  
  em_fit <- haplo.em(
    geno = geno,
    locus.label = loci,
    miss.val = c(0, NA),
    control = em_control
  )
  
  if (em_fit$converge != 1) {
    warning("The EM algorithm has not converged. Please increase max.iter.")
  }
  

  hap_matrix <- as.data.frame(
    em_fit$haplotype,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  colnames(hap_matrix) <- loci
  

  hap_name <- apply(
    hap_matrix,
    1,
    function(x) {
      paste0(
        loci,
        "*",
        x,
        collapse = "-"
      )
    }
  )
  
  hap_frequency <- data.frame(
    hap_code = seq_len(nrow(hap_matrix)),
    haplotype = hap_name,
    hap_matrix,
    frequency = em_fit$hap.prob,
    expected_count = 2 * nrow(geno) * em_fit$hap.prob,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  hap_frequency$frequency_percent <-
    100 * hap_frequency$frequency
  
  hap_frequency$common <-
    hap_frequency$frequency >= min_freq
  
  hap_frequency <- hap_frequency[
    order(hap_frequency$frequency, decreasing = TRUE),
  ]
  
  rownames(hap_frequency) <- NULL
  

  diplotype_posterior <- data.frame(
    ID = sample_id[em_fit$indx.subj],
    subject_index = em_fit$indx.subj,
    hap1_code = em_fit$hap1code,
    hap2_code = em_fit$hap2code,
    haplotype_1 = hap_name[em_fit$hap1code],
    haplotype_2 = hap_name[em_fit$hap2code],
    posterior = em_fit$post,
    stringsAsFactors = FALSE
  )
  
  diplotype_posterior <- diplotype_posterior[
    order(
      diplotype_posterior$subject_index,
      -diplotype_posterior$posterior
    ),
  ]
  
  rownames(diplotype_posterior) <- NULL
  

  best_diplotype <- diplotype_posterior[
    !duplicated(diplotype_posterior$subject_index),
  ]
  
  best_diplotype$high_confidence <-
    best_diplotype$posterior >= 0.80
  
  rownames(best_diplotype) <- NULL
  

  dosage_matrix <- matrix(
    0,
    nrow = nrow(geno),
    ncol = nrow(hap_matrix)
  )
  
  for (i in seq_along(em_fit$post)) {
    
    subject_i <- em_fit$indx.subj[i]
    hap1_i <- em_fit$hap1code[i]
    hap2_i <- em_fit$hap2code[i]
    post_i <- em_fit$post[i]
    
    dosage_matrix[subject_i, hap1_i] <-
      dosage_matrix[subject_i, hap1_i] + post_i
    
    dosage_matrix[subject_i, hap2_i] <-
      dosage_matrix[subject_i, hap2_i] + post_i
  }
  
  colnames(dosage_matrix) <- hap_name
  
  haplotype_dosage <- data.frame(
    ID = sample_id,
    dosage_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  common_haplotypes <- hap_name[
    em_fit$hap.prob >= min_freq
  ]
  
  common_haplotype_dosage <- haplotype_dosage[
    ,
    c("ID", common_haplotypes),
    drop = FALSE
  ]
  

  list(
    em_fit = em_fit,
    convergence = em_fit$converge,
    log_likelihood = em_fit$lnlike,
    haplotype_frequency = hap_frequency,
    diplotype_posterior = diplotype_posterior,
    best_diplotype = best_diplotype,
    haplotype_dosage = haplotype_dosage,
    common_haplotype_dosage = common_haplotype_dosage,
    genotype_matrix = geno
  )
}
# chinese cohort
df = read.csv('./Chinese_HLA_data.csv',row.names = 1)
result_drdq <- run_hla_em(
  df = df,
  loci = c("DRB1", "DQA1", "DQB1"),
  min_freq = 0.005
)
best_drdq <- subset(
  result_drdq$best_diplotype,
  posterior >= 0.80
)

hap_long <- best_drdq %>%
  dplyr::select(sample = ID, haplotype_1, haplotype_2) %>%
  pivot_longer(cols = c(haplotype_1, haplotype_2),
               names_to = "hap_no",
               values_to = "DQ_haplotype")

hap_count <- hap_long %>%
  count(sample, DQ_haplotype) %>%
  tidyr::pivot_wider(names_from = DQ_haplotype, values_from = n, values_fill = 0)
write.csv(hap_count,'Chinese-DR-DQ-hap.csv')
result_dp <- run_hla_em(
  df = df,
  loci = c("DPA1", "DPB1"),
  min_freq = 0.005
)

best_dp <- subset(
  result_dp$best_diplotype,
  posterior >= 0.80
)

hap_long <- best_dp %>%
  dplyr::select(sample = ID, haplotype_1, haplotype_2) %>%
  pivot_longer(cols = c(haplotype_1, haplotype_2),
               names_to = "hap_no",
               values_to = "DP_haplotype")
hap_count <- hap_long %>%
  count(sample, DP_haplotype) %>%
  tidyr::pivot_wider(names_from = DP_haplotype, values_from = n, values_fill = 0)
write.csv(hap_count,'Chinese-DP-hap.csv')

# European cohort
df = read.csv('./UKB_HLA_data.csv',row.names = 1)

result_drdq <- run_hla_em(
  df = df,
  loci = c("DRB1", "DQA1", "DQB1"),
  min_freq = 0.005
)
best_drdq <- subset(
  result_drdq$best_diplotype,
  posterior >= 0.80
)

hap_long <- best_drdq %>%
  dplyr::select(sample = ID, haplotype_1, haplotype_2) %>%
  pivot_longer(cols = c(haplotype_1, haplotype_2),
               names_to = "hap_no",
               values_to = "DQ_haplotype")

hap_count <- hap_long %>%
  count(sample, DQ_haplotype) %>%
  tidyr::pivot_wider(names_from = DQ_haplotype, values_from = n, values_fill = 0)
write.csv(hap_count,'European-DR-DQ-hap.csv')
result_dp <- run_hla_em(
  df = df,
  loci = c("DPA1", "DPB1"),
  min_freq = 0.005
)

best_dp <- subset(
  result_dp$best_diplotype,
  posterior >= 0.80
)

hap_long <- best_dp %>%
  dplyr::select(sample = ID, haplotype_1, haplotype_2) %>%
  pivot_longer(cols = c(haplotype_1, haplotype_2),
               names_to = "hap_no",
               values_to = "DP_haplotype")
hap_count <- hap_long %>%
  count(sample, DP_haplotype) %>%
  tidyr::pivot_wider(names_from = DP_haplotype, values_from = n, values_fill = 0)
write.csv(hap_count,'European-DP-hap.csv')

# haplotype association analysis
cli = read.csv('./Chinese_HLA_association_covariates.csv',row.names = 1)
HLA_res_all = data.frame()
for (hap in c('DR-DQ','DP')){
    df = fread(paste0('Chinese-',hap,'-hap.csv'))
    hp = colnames(df)[3:ncol(df)]
    df = data.frame(df)
    rownames(df) = df$sample
    df = df[,3:ncol(df)]
    colnames(df) = hp
    df2 = df[,colSums(df)>= (9126*2*.005)]
    test = cbind(df2,cli[rownames(df2),c('sex','disease_status','PC1','PC2')])
    name = c();p <- c();coe <- c();OR <- c();CI_L = c();CI_H = c();SE = c()
    for( i in seq(1:(length(test)-length(c('sex','disease_status','PC1','PC2'))))){
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
    HLA_res$HLA_haplotype = hap
    HLA_res_all = rbind(HLA_res_all,HLA_res)
    }
c_res = HLA_res_all

cli = read.csv('./UKB_HLA_association_covariates.csv',row.names = 1)
HLA_res_all = data.frame()
for (hap in c('DR-DQ','DP')){
    df = fread(paste0('European-',hap,'-hap.csv'))
    hp = colnames(df)[3:ncol(df)]
    df = data.frame(df)
    rownames(df) = df$sample
    df = df[,3:ncol(df)]
    colnames(df) = hp
    df2 = df[,colSums(df)>= (11827*2*.005)]
    test = cbind(df2,cli[rownames(df2),c('sex','disease_status','PC1','PC2')])
    name = c();p <- c();coe <- c();OR <- c();CI_L = c();CI_H = c();SE = c()
    for( i in seq(1:(length(test)-length(c('sex','disease_status','PC1','PC2'))))){
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
    HLA_res$HLA_haplotype = hap
    HLA_res_all = rbind(HLA_res_all,HLA_res)
    }
e_res = HLA_res_all
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
meta_data2 = data.frame()
for (i in c('DR-DQ','DP')){
    temp = meta_data[meta_data$HLA_haplotype_EUR == i,]
    temp$Padj_META = p.adjust(temp$P_META, method = "bonferroni")
    meta_data2  = rbind(meta_data2,temp)
    }

write.csv(meta_data2[meta_data2$Padj_META < 0.05,],'./res/T1D-haplotype-association-meta-res.csv')