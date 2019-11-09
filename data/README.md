# Data Description
_Data is from 1000 genomes Project_
## Genotype data can be accessed [here](https://drive.google.com/open?id=1O8QLlMd9gqd0voGAy-HiKZq4MO0igniV)

1. **Local ancestry** file has been generated from chromosome 21 using RFMix (Maples et al. 2013).

2. **Global ancestry** files have been generated using ADMIXTURE (Alexander et al. 2009). There are two files: **K3.results.txt** and **K5.results.txt**, the first includes only AFR, AMR, and EUR super-populations and the second containing all samples. 

3. **Simulated phenotype** files were generated using the liability threshold model as implemented by GCTA (Yang et al. 2011). SNPs were selected randomly; snp effect sizes and the age variable were randomly generated from a normal distribution. Sex and population labels are provided by 1000 Genomes. 
   1. The file **sims.pheno.csv** is a basic CSV file.
   2. The files **sims.pheno**, **sims.covar**, and **sims.sex.txt** are provided to use with PLINK (Chang et al. 2015).
   3. The file **summary_stats.txt** is a simulated GWAS summary stats file for PRS calculation.  
4. **MISC files** include script for generated phenotypes, summary stats from GWAS using simulated ata, a population information file, and an allele frequency file for SNPS used in simulation. 
