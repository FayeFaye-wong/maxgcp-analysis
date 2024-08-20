#   INPUT VARIBLES: 
#	ell, Mx1 vector of LD scores; 
#	z.1, Mx1 vector of estimated marginal per-normalized-genotype effects on trait 1
#   	(or Z scores; invariant to scaling);  
#	z.2, Mx2 vector of effects on trait2;  
# no.blocks, number of jackknife blocks;
#	crosstrait.intercept, 0 if crosstrait LDSC intercept should be fixed and 1 otherwise;
# ldsc.intercept, 0 if LDSC intercept should be fixed and 1 otherwise (1 is recommended);
# weights, Mx1 vector of regression weights; 
#	sig.threshold, threshold above which to discard chisq statistics for the purpose of estimating
#   	the LDSC intercept; large-effect SNPs discarded above sig_threshold*mean(z.x^2);
# n1, sample size for trait 1, only needed if ldsc.intercept=1; 
# n2, sample size for trait 2, only needed if ldsc.intercept=1; 
#	intercept.12, covariance between sampling errors for Z1 and Z2, only needed if
# crosstrait_intercept=0. If the 2 GWAS are disjoint, this can be set to zero.
#
#   OUTPUT VARIABLES: 
#	lcv.output, a list with named entries:
#   "zscore", Z score for partial genetic causality. zscore>>0 implies gcp>0.
#	pval.gcpzero.2tailed, 2-tailed p-value for null that gcp=0.
#   "gcp.pm", posterior mean gcp (gcp=1: trait 1 -> trait 2; gcp=-1: trait 2-> trait 1); 
#   "gcp.pse", posterior standard error for gcp; 
#   "rho.est", estimated genetic correlation; 
#   "rho.err", standard error of rho estimate;
#   "pval.fullycausal" [2 entries], p-values for null that gcp=1 or that gcp=-1, respectively; 
#   "h2.zscore" [2 entries], z scores for trait 1 and trait 2 being heritable, respectively;
#   	we recommend reporting results for h2.zscore > 7 (a very stringent threshold).


library(dplyr)


# github: https://github.com/lukejoconnor/LCV/blob/master/R/MomentFunctions.R
source("LCV/R/RunLCV.R")

# Load traits data and calculate Zs
ldsc_path <- "../Plink/sample_data/ldsc"
results_path <- "../Plink/sample_data/phenotypes_causal"
phenotypes <- c("pheno_0", "pheno_1", "pheno_2", "pheno_3", "pheno_4", "pheno_5", "pheno_6", "pheno_7", "pheno_8", "pheno_9")

results <- data.frame(trait1=character(), trait2=character(), gcp_pm=numeric(), gcp_pse=numeric(), log10_p=numeric(), stringsAsFactors=FALSE)

for (i in 1:(length(phenotypes)-1)) {
  for (j in (i+1):length(phenotypes)) {
    trait1 <- phenotypes[i]
    trait2 <- phenotypes[j]
    
    trait1File <- file.path(ldsc_path, paste0("naive.", trait1, "_munged.sumstats.gz"))
    trait2File <- file.path(ldsc_path, paste0("naive.", trait2, "_munged.sumstats.gz"))
    
    d0 <- na.omit(read.table(gzfile(trait1File), header=TRUE, sep="\t", stringsAsFactors=FALSE))
    d1 <- na.omit(read.table(gzfile(trait2File), header=TRUE, sep="\t", stringsAsFactors=FALSE))
    
    # Load LD scores
    ldscoresFile <- file.path(ldsc_path, "kpg.l2.ldscore")
    ldscoresFile <- read.table(ldscoresFile, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    
    # Merge datasets
    m <- merge(ldscoresFile, d0, by="SNP")
    data <- merge(m, d1, by="SNP")
    
    # Sort by position 
    data$CHR[data$CHR == "X"] <- 23
    data$CHR[data$CHR == "Y"] <- 24
    data$CHR <- as.numeric(data$CHR)
    data$BP <- as.numeric(data$BP)
    data <- data[order(data$CHR, data$BP), ]
    
    mhc_start <- 25500000
    mhc_end <- 34000000
    data <- data %>%
      filter(!(CHR == 6 & BP >= mhc_start & BP <= mhc_end))
    
    N1 = data$N.x[1]
    N2 = data$N.y[1]
    # Run LCV
    LCV <- RunLCV(data$L2, data$Z.x, data$Z.y, n.1=N1, n.2=N2, ldsc.intercept=0)
    # Store results
    result <- data.frame(trait1=trait1, trait2=trait2, 
                         gcp_pm=LCV$gcp.pm, 
                         gcp_pse=LCV$gcp.pse, 
                         log10_p=log(LCV$pval.gcpzero.2tailed)/log(10),
                         stringsAsFactors=FALSE)
    
    results <- bind_rows(results, result)
  }
}


write.csv(results, file.path(results_path, "lcv_naive_result.csv"), row.names=FALSE)


results <- data.frame(trait1=character(), trait2=character(), gcp_pm=numeric(), gcp_pse=numeric(), log10_p=numeric(), stringsAsFactors=FALSE)

for (i in 1:(length(phenotypes)-1)) {
  for (j in (i+1):length(phenotypes)) {
    trait1 <- phenotypes[i]
    trait2 <- phenotypes[j]
    
    trait1File <- file.path(ldsc_path, paste0("maxgcp.", trait1, "_munged.sumstats.gz"))
    trait2File <- file.path(ldsc_path, paste0("maxgcp.", trait2, "_munged.sumstats.gz"))
    
    d0 <- na.omit(read.table(gzfile(trait1File), header=TRUE, sep="\t", stringsAsFactors=FALSE))
    d1 <- na.omit(read.table(gzfile(trait2File), header=TRUE, sep="\t", stringsAsFactors=FALSE))
    
    # Load LD scores
    ldscoresFile <- file.path(ldsc_path, "kpg.l2.ldscore")
    ldscoresFile <- read.table(ldscoresFile, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    
    # Merge datasets
    m <- merge(ldscoresFile, d0, by="SNP")
    data <- merge(m, d1, by="SNP")
    
    # Sort by position 
    data$CHR[data$CHR == "X"] <- 23
    data$CHR[data$CHR == "Y"] <- 24
    data$CHR <- as.numeric(data$CHR)
    data$BP <- as.numeric(data$BP)
    data <- data[order(data$CHR, data$BP), ]
    
    mhc_start <- 25500000
    mhc_end <- 34000000
    data <- data %>%
      filter(!(CHR == 6 & BP >= mhc_start & BP <= mhc_end))
    
    N1 = data$N.x[1]
    N2 = data$N.y[1]
    # Run LCV
    LCV <- RunLCV(data$L2, data$Z.x, data$Z.y, n.1=N1, n.2=N2, ldsc.intercept=0)
    # Store results
    result <- data.frame(trait1=trait1, trait2=trait2, 
                         gcp_pm=LCV$gcp.pm, 
                         gcp_pse=LCV$gcp.pse, 
                         log10_p=log(LCV$pval.gcpzero.2tailed)/log(10),
                         stringsAsFactors=FALSE)
    
    results <- bind_rows(results, result)
  }
}


write.csv(results, file.path(results_path, "lcv_maxgcp_result.csv"), row.names=FALSE)

