#set prefix
prefix=test

#ADMIXTURE can be run in supervised and unpservied manner
#recomendation is to merge with references (i.e. 1000 Genomes) and run unsupervised
#the inclusion of references makes results more interpretable

##run admixture analysis
#specifcy bed file
data=$prefix.bed

#specify K
K=3
#runs 10 replicates 
#bash run_admix.sh $data $K


#summarize results. Very fast, just run when admix jobs are finished
fam=$prefix.fam
#summarizes by population/site 
#identifies replicate with great log-likelihood
bash admix_summary.sh $K $fam