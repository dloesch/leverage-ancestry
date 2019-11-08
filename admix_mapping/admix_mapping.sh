#!/bin/bash

#specify chromosome
chr=$1

#specify number of reference groups
groups=$2

#specify prefix
prefix=$3

#create subdirectories
mkdir -p admix_mapping

#go to admix_mapping sudirectory
cd admix_mapping

#create rest of subdirectories
mkdir -p admix_files
mkdir -p ./admix_files/${groups}_groups
mkdir -p ${groups}_groups
mkdir -p results
mkdir -p results/${groups}_groups


#RFmix results
rfmix=../local_ancestry/1KG.local_ancestry.${groups}_groups.chr$chr.msp.tsv.gz 

#get list of subjects used in RFmix analysis
if [ ! -f $prefix.subjects.${groups}_groups.txt ]
then
    gunzip -c $rfmix | head -2 | tail -1 | colrm 1 48 | sed 's/\.0//g' | sed 's/\.1//g' | tr '\t' '\n' | uniq > $prefix.subjects.${groups}_groups.txt

fi

#parse filles
echo "parsing local ancestry file"
subjects=$prefix.subjects.${groups}_groups.txt
Rscript ../parse.local_ancestry.R $chr $subjects $groups $prefix $rfmix


for i in AMR EAS AFR SAS; 
do 
    file=$prefix.${i}.chr${chr}.${groups}_groups.txt 
    if [ -f $file ]
    then
        gzip $file
        mv $file.gz ./admix_files/${groups}_groups
    fi
done

#admixture mapping
#pheno file
pheno=../sims.pheno.csv


#global admixture file
#can use average output of RFmix
#in this case, using ADMIXTURE output
admix=../K${groups}.results.txt

#related subjects to be removed, generated via KING's unrelated algorithim 

related=../simsunrelated_toberemoved.txt


###perform LRT######
trait="T1"
Rscript ../admix_mapping.R $chr $pheno $admix $related $groups $prefix $trait
results=$prefix.AM_results.chr${chr}.${groups}_groups.txt
gzip $results
mv $results.gz ./results/${groups}_groups

#go back to working directory and exit
cd ..
exit