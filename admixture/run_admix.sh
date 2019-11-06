#specify bed file
data=$1

#specify number of Ks
K=$2

#program path
admix=/usr/local/bin/admixture

#organize files
mkdir -p K${K}
cd K${K}

echo "run   seed" > K${K}.random.txt

#runs 10 replicates
for run in {1..10};
do
    mkdir -p run${run}
    cd run${run}

    #set random seed
    seed=$RANDOM

    #ideally, submit as batch jobs. We use SGE
    $admix --cv=10 $data $K -j4 -s $seed 

    cd ..
    echo "run${run} $seed" >> K${K}.random.txt
done

cd ..
exit