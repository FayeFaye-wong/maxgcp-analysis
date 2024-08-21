library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
input_pca_path <- args[1]
output_pca_path <- args[2]

pc_results <- list.files(path = input_pca_path, pattern = "*.formatted.txt", full.names = TRUE)
gwas_list <- lapply(pc_results, read.table, header = TRUE)

K <- 2  # number of top PCs included in the first group
N <- length(gwas_list)  # Total number of PCs

combine_pvalues_modified <- function(pvalues, K, N) {
  pvalues_top_K <- pvalues[1:K]
  pvalues_remaining <- pvalues[(K+1):N]
  
  chi_square_top_K <- -2 * sum(log(pvalues_top_K))
  chi_square_remaining <- -2 * sum(log(pvalues_remaining))

  combined_statistic <- -2 * (log(1 - pchisq(chi_square_top_K, df = 2 * K)) +
                                log(1 - pchisq(chi_square_remaining, df = 2 * (N - K))))
  
  combined_pvalue <- pchisq(combined_statistic, df = 4, lower.tail = FALSE)
  
  return(combined_pvalue)
}

combined_results <- gwas_list[[1]] %>%
  select(CHR,POS,SNP) %>%
  mutate(combined_pvalue = apply(
    sapply(gwas_list, function(x) x$P), 
    1, 
    function(p) combine_pvalues_modified(p, K, N)
  ))


write.table(combined_results, file = output_pca_path, sep = "\t", row.names = FALSE, quote = FALSE)

