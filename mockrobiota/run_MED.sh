source variables.env

declare -a mocks=("mock4" "mock5" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23")
declare -a Runs=("Merged" "Single")



for run in ${Runs[@]}; do 
workspace=${trunk}/${run}

input=${workspace}/input_fastqs

mkdir -p ${workspace}/MED

cp $input/*.fastq $workspace/MED

#convert files into fasta
for file in $workspace/MED/*.fastq; do 
$mothur "#fastq.info(fastq=${file},qfile=F)" 
rm $file
done

for mock in ${mocks[@]}; do 
mkdir -p $workspace/MED/${mock}
mv $workspace/MED/${mock}_*.fasta $workspace/MED/${mock}

for file in $workspace/MED/${mock}/*.fasta; do
sample=$(basename $file | cut -d'.' -f1 | cut -d'_' -f2)
sed -i -- "s/\(\-\)/:/g" $file
sed -i -- "s/\(\_\)/:/g" $file
sed -i -- "s/\(^>\)/\1Sample-${sample}_/g" $file
done
cat $workspace/MED/${mock}/*.fasta >$workspace/MED/${mock}.fasta > tmp && mv tmp $workspace/MED/${mock}/${mock}.fasta
done
done

####-----------------------------------------------------------------------------------------
##becareful to change environment according to working machine 
source activate py27

for run in ${Runs[@]}; do 
workspace=${trunk}/${run}

for mock in ${mocks[@]}; do 
cd $workspace/MED/${mock}
echo ${mock}
if [ ${mock} == 'mock4' ] || [ ${mock} == 'mock5' ] 
then 
decompose $workspace/MED/${mock}/${mock}.fasta -V 8 -M 1 --skip-removing-outliers --relocate-outliers -o $workspace/MED/${mock}
else
decompose $workspace/MED/${mock}/${mock}.fasta -V 8 -M 1 --skip-check-input-file --skip-removing-outliers --relocate-outliers -o $workspace/MED/${mock}
fi

done
done

cd ~



#parse output into mothur formate (.fasta and .count files)
for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
run_anlys=${anlystrunk}/${run}

for mock in ${mocks[@]}; do 
echo $mock
mkdir -p ${run_anlys}/MED/${mock}
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $workspace/MED/${mock}/MATRIX-COUNT.txt | sed -r -e 's/samples/Representative_Sequence\ttotal/' | sed -E 's/[[:space:]]Sample/\tSample/g' | awk '{if(NR==1) print; else {for(i=2;i<=NF;i++) sum+=$i; print "N"$1"\t"sum"\t"$2"\t"$3"\t"$4"\t"$5; sum=0} }' > ${run_anlys}/MED/${mock}/${mock}_output.count

sed -r -e 's/^>/>N/' $workspace/MED/${mock}/NODE-REPRESENTATIVES.fasta | sed -r -e 's/\|size\:[0-9]*//' >${run_anlys}/MED/${mock}/${mock}.fasta

done



done
