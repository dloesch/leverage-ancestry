# Data Description
_Data is from 1000 genomes Project_

1. **Local ancestry** file has been generated from chromosome 21 using RFMix-2.

2. **Global ancestry** files have been generated using ADMIXTURE. There are two files: **K3.results.txt** and **K5.results.txt**, the first includes only AFR, AMR, and EUR super-populations and the second containing all samples. 

3. **Simulated phenotype** files were generated using the liability threshold model as implemented by GCTA. SNPs were selected randomly; snp effect sizes and the age variable were randomly generated from a normal distribution. Sex and population labels are provided by 1000 Genomes. 
   1. The file **sims.pheno.csv** is a basic CSV file.
   2. The files **sims.pheno**, **sims.covar**, and **sims.sex.txt** are provided to use with PLINK.
   3. The file **summary_stats.txt** is a simulated GWAS summary stats file for PRS calculation.  
