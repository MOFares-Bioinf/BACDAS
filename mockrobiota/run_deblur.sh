source variables.env

declare -a mocks=("mock4" "mock5" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23")
declare -a Runs=("Single" "Merged")

mkdir -p $workspace/deblur

source activate deblurenv


for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
deblur_output=${workspace}/deblur
run_anlys=${anlystrunk}/${run}

#run workflow
for mock in ${mocks[@]};
do
mkdir -p $workspace/deblur/${mock}
cp ${workspace}/input_fastqs/${mock}_*.fastq $workspace/deblur/${mock}
deblur workflow --seqs-fp $workspace/deblur/${mock} --output-dir $workspace/deblur/${mock}/output -t -1 -w --keep-tmp-files --min-size 1 --min-reads 0 
done



#re-build biom tables using denoised reads before chimera removal -- skip chimera removal step--
for mock in ${mocks[@]}; do
echo $mock
mkdir -p $workspace/deblur/${mock}/compined_biome
for file in $workspace/deblur/${mock}/output/deblur_working_dir/*.fastq.trim.derep.no_artifacts.msa.deblur; do
sample=$(basename $file | cut -d'.' -f1);
cp ${file} $workspace/deblur/${mock}/compined_biome/${sample}.fasta;
done;
deblur build-biom-table $workspace/deblur/${mock}/compined_biome $workspace/deblur/${mock}/compined_biome --min-reads 1 --file_type '.fasta'
done;


#convert biom files into tsv
mkdir -p $workspace/deblur/biom_tables
for mock in ${mocks[@]}; do
echo $mock
biom convert -i /$workspace/deblur/${mock}/compined_biome/all.biom -o /$workspace/deblur/${mock}/compined_biome/table.from_biom.txt --to-tsv

cp $workspace/deblur/${mock}/compined_biome/table.from_biom.txt $workspace/deblur/biom_tables/${mock}_table.from.biome.txt
done



#Parse output into mothur format (.fasta and .count files)
for mock in ${mocks[@]}; do
mkdir -p ${run_anlys}/deblur/${mock}
echo $mock
awk 'NR!=1' $workspace/deblur/biom_tables/${mock}_table.from.biome.txt | sed -r -e 's/\#OTU ID/Representative_Sequence\ttotal/' | awk '{if(NR==1) print; else {for(i=1;i<=NF;i++) sum+=$i; print "ASV_"NR-1"\t"sum"\t"$2"\t"$3"\t"$4"\t"$5; sum=0} }' > ${run_anlys}/deblur/${mock}/${mock}_output.count
awk 'NR>2' $workspace/deblur/biom_tables/${mock}_table.from.biome.txt | awk '{print ">ASV_"NR"\n"$1}' > ${run_anlys}/deblur/${mock}/${mock}.fasta

done



done



conda deactivate

