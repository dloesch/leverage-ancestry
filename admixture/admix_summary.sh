
K=$1

#fam file
fam=$2

#go to directory
cd K${K}

echo "CV_error" > K${K}.CV_error.txt
echo "loglikelihood" > K${K}.loglikelihood.txt

for run in {1..10};
do
#CV error
cat ./run${run}/run${run}.K${K}.o* | grep CV >> K${K}.CV_error.txt

#loglikelihood
cat ./run${run}/run${run}.K${K}.o* | grep Log | tail -1 >> K${K}.loglikelihood.txt

done

#best loglikelihood
best=$(cat K${K}.loglikelihood.txt | grep Log | cut -d" " -f2 | sort -n | tail -1 | sed 's/-//g')

#summary file
paste K${K}.random.txt K${K}.CV_error.txt K${K}.loglikelihood.txt > K${K}.summary.txt

#best run
cat K${K}.summary.txt | grep "$best" > K${K}.best_run.txt


#create outputfile from best run with column ids
#column id files were made from IID column of PED file

cat $fam | cut -f1 -d" " > K${K}.subjects.txt
 
best=$(cat K${K}.best_run.txt | cut -f1 -d" ")
results=$(ls ./$best | grep Q)

paste K${K}.subjects.txt ./$best/$results > K${K}.results.txt


##now calculate average ancestry per population group/site
pops=populations.txt #path to file containing group/site information. 
#file format: two tab-sepearted columns - ID and POPULATION
Rscript ../ave_ancestry.R K${K}.results.txt $pops $K

cd ..


exit
