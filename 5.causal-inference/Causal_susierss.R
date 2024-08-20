library(susieR)
library(dplyr)
# z: z score
ldsc_path <- "../Plink/sample_data/ldsc"
gwas_path <- "../Plink/sample_data/gwas"
maxgcp_path <- "../Plink/sample_data/maxgcp_gwas"
results_path <- "../Plink/sample_data/causal_snps"
naive <- file.path(gwas_path, "naive.pheno_0.glm.linear.ldak")
naive <- na.omit(read.table(naive,header=TRUE,sep="",stringsAsFactors = FALSE))
colnames(naive)[colnames(naive) == "Predictor"] <- "SNP"
snps_of_interest <- unique(naive$SNP)

ldFile <- file.path(ldsc_path, "filtered_kpg.ld")
ld_all <- read.table(ldFile,header=TRUE,sep='',stringsAsFactors=FALSE)
ld_naive <- ld_all[ld_all$SNP_A %in% snps_of_interest & ld_all$SNP_B %in% snps_of_interest, ]

ld_naive <- ld_naive %>% 
  arrange(CHR_A, BP_A) %>% 
  mutate(group = NA)
bp_threshold <- 1000000
current_group <- 1
ld_naive$group[1] <- current_group
for (i in 2:nrow(ld_naive)) {
  if (ld_naive$CHR_A[i] != ld_naive$CHR_A[i-1] || 
      abs(ld_naive$BP_A[i] - ld_naive$BP_A[i-1]) > bp_threshold) {
    current_group <- current_group + 1
  }
  ld_naive$group[i] <- current_group
}


ld_naive_matrices <- list()
group_ids <- unique(ld_naive$group)
for (group_id in group_ids) {
  ld_group <- ld_naive[ld_naive$group == group_id, ]
  snps_group <- unique(c(ld_group$SNP_A, ld_group$SNP_B))
  ld_matrix_group <- matrix(0, nrow=length(snps_group), ncol=length(snps_group))
  diag(ld_matrix_group) <- 1
  
  dimnames(ld_matrix_group) <- list(as.character(snps_group), as.character(snps_group))
  row_indices <- match(ld_group$SNP_A, snps_group)
  col_indices <- match(ld_group$SNP_B, snps_group)
  ld_matrix_group[cbind(row_indices, col_indices)] <- ld_group$R2
  ld_matrix_group[cbind(col_indices, row_indices)] <- ld_group$R2 
  ld_naive_matrices[[paste0("chr", ld_group$CHR_A[1], "_group", group_id)]] <- ld_matrix_group
}

pip_naive_results <- list()

for (chr in names(ld_naive_matrices)) {
  ld_matrix <- ld_naive_matrices[[chr]]
  snps_chr <- rownames(ld_matrix)
  z_chr <- naive$Z[match(snps_chr, naive$SNP)]
  fit <- susie_rss(z_chr, ld_matrix, n = 2405)
  pip_naive_results[[chr]] <- fit$pip
}

# SNPs with higher PIPs are more likely to be causal. 
# Researchers often use a threshold (e.g., PIP > 0.95) to select potential causal variants.

pip_naive_table <- data.frame(CHR = character(), BP = numeric(), SNP = character(), PIP = numeric(), stringsAsFactors = FALSE)

for (chr in names(pip_naive_results)) {
  pip_chr <- pip_naive_results[[chr]]
  snps_chr <- names(pip_chr)
  bp_chr_A <- ld_naive$BP_A[match(snps_chr, ld_naive$SNP_A)]
  bp_chr_B <- ld_naive$BP_B[match(snps_chr, ld_naive$SNP_B)]
  bp_chr <- ifelse(!is.na(bp_chr_A), bp_chr_A, bp_chr_B)
  cleaned_chr <- gsub("_group.*", "", chr)
  temp_df <- data.frame(CHR = cleaned_chr, BP = bp_chr, SNP = snps_chr, PIP = pip_chr, stringsAsFactors = FALSE)
  pip_naive_table <- rbind(pip_naive_table, temp_df)
}

write.csv(pip_naive_table, file.path(results_path,"naive_pip_results.csv"), row.names = FALSE)

maxgcp_path <- "/Users/huangfeifei/Documents/TLAB/Plink/sample_data/maxgcp_gwas"
maxgcp_file <- file.path(maxgcp_path, "maxgcp.pheno_0.glm.linear")

maxgcp_data <- na.omit(read.table(maxgcp_file, header = FALSE, sep = "", stringsAsFactors = FALSE))
colnames(maxgcp_data) <- c("#CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", 
                           "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", 
                           "BETA", "SE", "T_STAT", "P", "ERRCODE")


maxgcp_data <- data.frame(
  SNP = maxgcp_data$ID,
  A1 = maxgcp_data$A1,
  A2 = maxgcp_data$ALT,
  n = maxgcp_data$OBS_CT,
  Z = maxgcp_data$BETA / maxgcp_data$SE
)
maxgcp_data <- maxgcp_data[abs(maxgcp_data$Z)>1.96, ]

snps_of_interest <- unique(maxgcp_data$SNP)

ld_maxgcp <- ld_all[ld_all$SNP_A %in% snps_of_interest & ld_all$SNP_B %in% snps_of_interest, ]

ld_maxgcp <- ld_maxgcp %>% 
  arrange(CHR_A, BP_A) %>% 
  mutate(group = NA)
current_group <- 1
ld_maxgcp$group[1] <- current_group
for (i in 2:nrow(ld_maxgcp)) {
  if (ld_maxgcp$CHR_A[i] != ld_maxgcp$CHR_A[i-1] || 
      abs(ld_maxgcp$BP_A[i] - ld_maxgcp$BP_A[i-1]) > bp_threshold) {
    current_group <- current_group + 1
  }
  ld_maxgcp$group[i] <- current_group
}


ld_maxgcp_matrices <- list()
group_ids <- unique(ld_maxgcp$group)
for (group_id in group_ids) {
  ld_group <- ld_maxgcp[ld_maxgcp$group == group_id, ]
  snps_group <- unique(c(ld_group$SNP_A, ld_group$SNP_B))
  ld_matrix_group <- matrix(0, nrow=length(snps_group), ncol=length(snps_group))
  diag(ld_matrix_group) <- 1
  
  dimnames(ld_matrix_group) <- list(as.character(snps_group), as.character(snps_group))
  row_indices <- match(ld_group$SNP_A, snps_group)
  col_indices <- match(ld_group$SNP_B, snps_group)
  ld_matrix_group[cbind(row_indices, col_indices)] <- ld_group$R2
  ld_matrix_group[cbind(col_indices, row_indices)] <- ld_group$R2 
  ld_maxgcp_matrices[[paste0("chr", ld_group$CHR_A[1], "_group", group_id)]] <- ld_matrix_group
}

pip_maxgcp_results <- list()

for (chr in names(ld_maxgcp_matrices)) {
  ld_matrix <- ld_maxgcp_matrices[[chr]]
  snps_chr <- rownames(ld_matrix)
  z_chr <- maxgcp_data$Z[match(snps_chr, maxgcp_data$SNP)]
  fit <- susie_rss(z_chr, ld_matrix, n = 2405)
  pip_maxgcp_results[[chr]] <- fit$pip
}

pip_maxgcp_table <- data.frame(CHR = character(), BP = numeric(), SNP = character(), PIP = numeric(), stringsAsFactors = FALSE)

for (chr in names(pip_maxgcp_results)) {
  pip_chr <- pip_maxgcp_results[[chr]]
  snps_chr <- names(pip_chr)
  bp_chr_A <- ld_maxgcp$BP_A[match(snps_chr, ld_maxgcp$SNP_A)]
  bp_chr_B <- ld_maxgcp$BP_B[match(snps_chr, ld_maxgcp$SNP_B)]
  bp_chr <- ifelse(!is.na(bp_chr_A), bp_chr_A, bp_chr_B)
  cleaned_chr <- gsub("_group.*", "", chr)
  temp_df <- data.frame(CHR = cleaned_chr, BP = bp_chr, SNP = snps_chr, PIP = pip_chr, stringsAsFactors = FALSE)
  pip_maxgcp_table <- rbind(pip_maxgcp_table, temp_df)
}

write.csv(pip_maxgcp_table, file.path(results_path,"maxgcp_pip_results.csv"), row.names = FALSE)

