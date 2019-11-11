####calculate PRS#######
###step 1: parse files
###step 2: calculate PRS
###step 3: test PRS


#get prefix
prefix=$1

#file containing summary stats
summary_stats=summary_stats.txt

#plink format file containing target genotype data
geno=1KG.sim_data.thinned

##PARSE FILES###

##hyperparamter selection###
if [ $2 = "TRUE" ];
then

	#hyper paramters
	thresh=$3 # select pvalue threshold
	r2=$4    #select valaue of r2 for clumping

	##create directory
	mkdir -p LD-${r2}
	cd LD-${r2}
	
	#reference data for LD clumping. Recomended to use LD structure of population...
	# most similar to that of GWAS summary stats
	ref=REF_DATA
	plink --bfile $ref --clump $gwas --clump-p1 $thresh --clump-r2 $r2 --out thresh_$thresh.r2_$r2 

	cat thresh_${thresh}.r2_${r2}.clumped | awk '{print $3}' > clumped_variants_${thresh}_$r2.txt
	
	#remove unneeded files
	rm thresh_${thresh}.r2_${r2}.clumped*
	
	#extract clumped variants from target genotype file
	plink --bfile $geno --extract clumped_variants_${thresh}_$r2.txt --make-bed --out $prefix.clumped_${thresh}_$r2

	new_geno=$prefix.clumped_${thresh}_$r2
	
	#specify effect allele
	a1=summary_stats.txt
	
	#convert to text format for R
	plink --bfile $new_geno --a1-allele $a1 2 1 'Q' --recode A --allow-no-sex --out $prefix.clumped_${thresh}_$r2
	
	target=$prefix.clumped_${thresh}_$r2.raw

else
	#if not training hyperparamters, then using GWAS signif results	
	snps=$(awk 'NR>1 {print $1}' summary_stats.txt)
	thresh=GWAS_signif
	plink --bfile $geno --snps $snps --make-bed --out $prefix.$thresh

	new_geno=$prefix.$thresh
	
	#specify effect allele
	a1=summary_stats.txt

	#convert to text format for R
	plink --bfile $new_geno --a1-allele $a1 2 1 'Q' --recode A --allow-no-sex --out $new_geno

	target=$new_geno.raw
fi



#####Calculate PRS and perform associations######


#calculates PRS
#SNP name column
SNP=QTL

#effect size column
effect=Effect


Rscript calculate_PRS.R $prefix $thresh $target $summary_stats $SNP $effect "beta" #beta or OR

#test PRS

PRS=$prefix.PRS.$thresh.txt

pheno=sims.pheno.csv

trait=T1

Rscript test_PRS.R $prefix $thresh $pheno $PRS $trait


if [ $2 = "TRUE" ]
then

	## move back up to PRS directory
	cd ..
else
	exit 0

fi
