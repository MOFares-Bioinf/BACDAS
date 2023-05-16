source variables.env

declare -a mocks=("mockP" )
declare -a Runs=("Single" "Merged")


for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
dada2_output=${workspace}/DADA2
run_anlys=${anlystrunk}/${run}

for mock in ${mocks[@]}; do
mkdir -p ${workspace}/DADA2/${mock}
cp ${workspace}/input_fastqs/${mock}* ${workspace}/DADA2/${mock}

done

if [ $run = "Single" ]; then 
Rscript dada2_single.r
else
Rscript dada2_merge.r
fi

#parse output into .fasta & .count files
for mock in ${mocks[@]}; do
mkdir -p ${run_anlys}/DADA2/${mock}
echo $mock
sed -r -e 's/\"\"\,/Representative_Sequence\ttotal\t/' $workspace/DADA2/output_${mock}/transposed.csv | sed 's/\"//g' | sed 's/,/\t/g' | awk '{if(NR==1) print; else {for(i=1;i<=NF;i++) sum+=$i; print "ASV_"NR-1"\t"sum"\t"$2"\t"$3"\t"$4"\t"$5; sum=0} }' > ${run_anlys}/DADA2/${mock}/${mock}_output.count
awk -F"," 'NR!=1{print ">ASV_"NR-1"\n"$1}' $workspace/DADA2/output_${mock}/transposed.csv | sed 's/\"//g' > ${run_anlys}/DADA2/${mock}/${mock}.fasta
done

done
