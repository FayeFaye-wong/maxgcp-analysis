#!/bin/bash

# Need to be run in ldsc_env 

source activate ldsc_env

features=("pheno_0" "pheno_1" "pheno_2" "pheno_3" "pheno_4" "pheno_5" "pheno_6" "pheno_7" "pheno_8" "pheno_9")

base_dir="/Users/huangfeifei/Documents/TLAB/Plink/sample_data"
ldsc_dir="${base_dir}/ldsc/ldsc-master"
for pheno in "${features[@]}"; do
  python ${ldsc_dir}/munge_sumstats.py \
  --sumstats ${base_dir}/ldsc/naive.${pheno}.formatted.txt \
  --out ${base_dir}/ldsc/naive.${pheno}_munged \
  --N 2504
done

for pheno in "${features[@]}"; do
  python ${ldsc_dir}/munge_sumstats.py \
  --sumstats ${base_dir}/ldsc/maxgcp.${pheno}.formatted.txt \
  --out ${base_dir}/ldsc/maxgcp.${pheno}_munged \
  --N 2504
done

python ${ldsc_dir}/ldsc.py \
--bfile ${base_dir}/qc/filtered_kpg \
--l2 \
--ld-wind-cm 1 \
--out ${base_dir}/ldsc/kpg \
--yes-really

