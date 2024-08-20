Causal Inference Analysis Following GWAS

After performing GWAS, we implement various methods to evaluate the potential power of MaxGCP in understanding 
the causal mechanisms of diseases. The methods compared include Naive, PCA, MaxH, and MaxGCP.

Analysis Workflow

Causality Assessment with LCV:

We utilize the Latent Causal Variable (LCV) model to determine whether the correlated MaxGCP phenotypes exhibit 
causality. This helps in identifying potential causal relationships between phenotypes and the disease.

SNP Prioritization with Susie_RSS:

Using the Sum of Single Effects (SuSiE) with Summary Statistics (RSS), we calculate the Posterior Inclusion 
Probabilities (PIPs) for Single Nucleotide Polymorphisms (SNPs) within linkage disequilibrium regions. SNPs with 
high PIPs are strong candidates for being causal variants, indicating their potential role in disease etiology.

Functional Mechanism Testing with coloc:

The coloc method tests whether the GWAS interest association is driven by an expression Quantitative Trait Locus 
(eQTL). This could suggest a functional mechanism linking genetic variants to gene expression and, ultimately, 
to the disease.

Due to the lack of specific GWAS results for each phenotype in PCA and MaxH, these methods could not perform 
causal inference. 
