
for chr in {1..22}; do 
	plink --vcf chr$chr.1kg.phase3.v5a.vcf.gz --maf 0.05 --threads 2 --out chr$chr --indep-pairwise 50 5 0.8
	plink --vcf chr$chr.1kg.phase3.v5a.vcf.gz --maf 0.05 --threads 2 --out chr$chr --exclude chr$chr.prune.out --make-bed --snps-only 'just-agct'
	
done

touch merge_files
for chr in {1..22}; do echo chr$chr >> merge_files; done

plink --merge-list merge_files --make-bed --out 1KG.sim_data

#above leaves 2 million SNPS, create arbitrarily thinned file and pruned file

#thin with probability of 0.1
plink --bfile 1KG.sim_data --thin 0.1 --make-bed --out 1KG.sim_data.thinned


#prune with r2 of 0.2
plink --bfile 1KG.sim_data --indep-pairwise 50 5 0.2 --out exclude 
plink --bfile 1KG.sim_data --exclude exclude.prune.out --out 1KG.sim_data.pruned --make-bed 

#clean up 
for chr in {1..22}; do rm chr$chr.bed chr$chr.log chr$chr.bim chr$chr.fam; done
rm *nosex


#simulate pheno


#randomly select 20 SNPs to be causal
cut -f2 1KG.sim_data.thinned.bim | shuf -n 20 > snp_list

#get pop-specific frequency of randomly-selected SNPS
KG=integrated_call_samples_v3.20130502.ALL.panel
for pop in $pops;do
	cat $KG | grep $pop | cut -f1 > $pop.txt
	plink --bfile 1KG.sim_data.thinned --freq --keep-fam $pop.txt --extract snp_list --out $pop
done

paste AFR.frq AMR.frq EAS.frq EUR.frq SAS.frq > temp

echo -e "SNP\tAFR\tAMR\tEAS\tEUR\tSAS" > snp.freq_by_pop.txt
cat temp | awk -v OFS='\t' 'NR>1 {print $2, $5, $11, $17, $23, $29}' >> snp.freq_by_pop.txt

#clean up

rm AFR.* EUR.* EAS.* AMR.* SAS.* 


#simulate
N=2504
cases=$((N/2))
controls=$((N-cases))

gcta64  --bfile 1KG.sim_data  --simu-cc $cases $controls  --simu-causal-loci snp_list  --simu-hsq 0.3  --simu-k 0.5  --simu-rep 3  --out sims


#generate PCAs

plink --bfile 1KG.sim_data.pruned --pca 10 header tabs --threads 2 --out 1KG


#create pheno files

Rscript prep_sim_pheno.R 

#run association on traits
pheno=sims.pheno
covar=sims.covar
trait=T1
plink --bfile 1KG.sim_data.thinned --update-sex sims.sex.txt --logistic --pheno $pheno --pheno-name $trait --covar $covar --covar-name 'AGE-PC5'  --maf 0.01 --out sims.$trait --threads 2  

head -1 sims.$trait.assoc.logistic > sims.$trait.results.txt
cat sims.$trait.assoc.logistic | grep ADD | grep -v NA >> sims.$trait.results.txt

gzip sims.$trait.results.txt
gzip sims.$trait.assoc.logistic


#plot GWAS


Rscript plot_GWAS.R sims.$trait.results.txt.gz $trait TRUE #arguments are filename traitname thin==TRUE/FALSE



#make file with three populations for local ancestry
KG=integrated_call_samples_v3.20130502.ALL.panel

cat $KG | grep EAS | cut -f1 > EAS.txt
cat $KG | grep SAS | cut -f1 > SAS.txt

cat EAS.txt SAS.txt > EAS_SAS.txt

plink --bfile 1KG.sim_data.thinned --remove-fam EAS_SAS.txt --make-bed --out 1KG.sim_data.thinned.AFR_AMR_EUR
