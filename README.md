# Leveraging Ancestry in Genetic Studies
## Knowing who is in your dataset
_Principal component analysis (PCA) and relationship inference using IBD or kinsihp coefficients are common practice. As an alternative, estimating ancestry proportions, a quantitative trait, can be more informative._
1. The script **run_admix.sh** runs ADMIXTURE; it does 10 replicates for each provided K. 
2. The script **admix_summary.sh** selects best replicate by log-likelihood and generates summary files such as CV error. 
3. The script **ave_admixture.R** calculates average admixture for each K per population. Population labels need to be provided as a text document with 2 columns: ID and Population. 
4. The script **admixture_workflow.sh** simply outlines the above procedure. 

## Adxmiture Mapping
_If there is evidence for differential disease risk by population, AM can be used to associate local ancestry segments with the disease._
1. The script **admix_mapping.sh** performs entire admix mapping procedure per chromosome. Provide global ancestry files (ADMIXTURE or RFMix output), local ancestry files (RFMix), KING unrelated output, and phenotype file.
2. The script **parse.local_ancestry.sh** parses the output of RFMix into separate files per reference ancestry. 
3. The script **admix_mapping.R** performs a likelihood ratio test (LRT) to determine if the inclusion of local ancestry segements improve the model. 

## Polygenic Risk Scores
_Polygenic Risk Scores (PRS) are the linear summation of GWAS summary statistics. The PRS is then used a variable in a predictive model (typically logistic regression)._
1. The script **PRS.sh** performs the entire PRS procedure. Provide a PLINK format file of target genotypes and a file continaing GWAS summary statistics. The script parses the genotype file into an R-readable format.
2. The script **calculate_PRS.R** calculates the PRS.
3. The script **test_PRS.R** tests the PRS in a logistic regression framework using 10-fold CV, Pseudo R2, and area under the reciever-operator curve.

## Data
Sample data is from phase 3 of the 1000 Genomes Project. Phenotype data is simulated using the liability threshold model as implemented by GCTA. 
Citation: “A global reference for human genetic variation” Nature 526 68-74 2015

## Software
1. PLINK 1.9: Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4. URL: https://www.cog-genomics.org/plink2/
2. ADMIXTURE: D.H. Alexander, J. Novembre, and K. Lange. Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664, 2009. URL: http://software.genetics.ucla.edu/admixture/
3. RFMix: Maples BK, Gravel S, Kenny EE, Bustamante CD. RFMix: a discriminative modeling approach for rapid and robust local-ancestry inference. Am J Hum Genet. 2013;93(2):278–288. URL https://github.com/slowkoni/rfmix
4. R:  R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
5. GCTA: Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011;88(1):76–82. URL https://cnsgenomics.com/software/gcta/#Overview
6. KING: Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) Robust relationship inference in genome-wide association studies. Bioinformatics 26(22):2867-2873. URL http://people.virginia.edu/~wc9c/KING/

