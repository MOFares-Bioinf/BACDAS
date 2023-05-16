source variables.env

declare -a mocks=("mockP")

declare -a Runs=("Merged" "Single")


#==================================================================
#=================== Mothur pipeline ==============================
#==================================================================

for run in ${Runs[@]}; do 
workspace=${trunk}/${run}

input=${workspace}/input_fastqs

mkdir -p ${workspace}/mothur
cp $input/*.fastq $workspace/mothur
#convert files into fasta
for file in $workspace/mothur/*.fastq; do 
$mothur "#fastq.info(fastq=${file}, qfile=F)" 
rm $file
done


#make group files
for mock in ${mocks[@]}; do 
echo $mock
mkdir -p $workspace/mothur/${mock}
mv $workspace/mothur/${mock}_*.fasta $workspace/mothur/${mock}
fasta=""
groups=""
count=0
for file in ${workspace}/mothur/${mock}/${mock}_*.fasta; do 
if [ $count == 0 ] 
then
fasta=$file 
groups=$(basename $file | cut -d'.' -f1)
else
fasta=${fasta}-$file
groups=${groups}-$(basename $file | cut -d'.' -f1)
fi
echo $count
((count++))
done;
echo $fasta
echo $groups
$mothur "#make.group(fasta=${fasta}, groups=${groups}); merge.files(input=${fasta}, output=${workspace}/mothur/${mock}/${mock}.merge.fasta)"
cat $workspace/mothur/${mock}/*group* >$workspace/mothur/${mock}/${mock}.grp
rm $workspace/mothur/${mock}/*group*
done

done




# start pipeline
for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
for mock in ${mocks[@]}; do 
echo ${mock}

#preclustering
$mothur "#unique.seqs(fasta=${workspace}/mothur/${mock}/${mock}.merge.fasta); \
summary.seqs(name=current); \
count.seqs(name=current, group=$workspace/mothur/${mock}/${mock}.grp); \
pre.cluster(fasta=current, count=current,align=gotoh, processors=4); \
dist.seqs(fasta=current, cutoff=0.1)" 


#clustering step
echo opti 
$mothur "#set.dir(output=$workspace/mothur/Opti); cluster(column=$workspace/mothur/${mock}/${mock}.merge.unique.precluster.dist, count=$workspace/mothur/${mock}/${mock}.merge.unique.precluster.count_table, method=opti);\
make.shared(count=current, list=current, label=0.03)"
echo average 
$mothur "#set.dir(output=$workspace/mothur/AN); cluster(column=$workspace/mothur/${mock}/${mock}.merge.unique.precluster.dist, count=$workspace/mothur/${mock}/${mock}.merge.unique.precluster.count_table, method=average, cutoff=0.15);\
make.shared(count=current, list=current, label=0.03)" 
echo vsearch 
$mothur "#set.dir(output=$workspace/mothur/DGC); cluster(fasta=$workspace/mothur/${mock}/${mock}.merge.unique.precluster.fasta, count=$workspace/mothur/${mock}/${mock}.merge.unique.precluster.count_table, method=dgc, cutoff=0.03);\
make.shared(count=current, list=current, label=0.03)" 

done
done



#####extract otu representative seuquences and calculate error rates####
for run in ${Runs[@]}; do 
workspace=${trunk}/${run}

for mock in ${mocks[@]}; do   
$mothur "#set.dir(output=$workspace/mothur/AN); \
get.oturep(count=${workspace}/mothur/${mock}/${mock}.merge.unique.precluster.count_table,\
list=${workspace}/mothur/AN/${mock}.merge.unique.precluster.an.list,\
fasta=${workspace}/mothur/${mock}/${mock}.merge.unique.precluster.fasta, \
method=abundance, rename=True, cutoff=0.03)"


$mothur "#set.dir(output=$workspace/mothur/DGC); \
get.oturep(count=${workspace}/mothur/${mock}/${mock}.merge.unique.precluster.count_table,\
list=${workspace}/mothur/DGC/${mock}.merge.unique.precluster.dgc.list,\
fasta=${workspace}/mothur/${mock}/${mock}.merge.unique.precluster.fasta, \
method=abundance, rename=True, cutoff=0.03)"


$mothur "#set.dir(output=$workspace/mothur/Opti); 
get.oturep(count=${workspace}/mothur/${mock}/${mock}.merge.unique.precluster.count_table,\
list=${workspace}/mothur/Opti/${mock}.merge.unique.precluster.opti_mcc.list,\
fasta=${workspace}/mothur/${mock}/${mock}.merge.unique.precluster.fasta, \
method=abundance, rename=True, cutoff=0.03)"
done
done






for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
run_anlys=${anlystrunk}/${run}


#generate .fasta and .count files
for mock in ${mocks[@]}; do 
mkdir -p ${run_anlys}/AN/${mock}
cp $workspace/mothur/AN/${mock}.merge.unique.precluster.an.0.03.rep.fasta ${run_anlys}/AN/${mock}/${mock}.fasta
cp $workspace/mothur/AN/${mock}.merge.unique.precluster.an.0.03.rep.count_table ${run_anlys}/AN/${mock}/${mock}_output.count
#cp $workspace/mothur/AN/${mock}.merge.unique.precluster.an.0.03.rep.error.chimera ${run_anlys}/AN/${mock}/${mock}.error.chimera


mkdir -p ${run_anlys}/DGC/${mock}
cp $workspace/mothur/DGC/${mock}.merge.unique.precluster.dgc.0.03.rep.fasta ${run_anlys}/DGC/${mock}/${mock}.fasta
cp $workspace/mothur/DGC/${mock}.merge.unique.precluster.dgc.0.03.rep.count_table ${run_anlys}/DGC/${mock}/${mock}_output.count
#cp $workspace/mothur/DGC/${mock}.merge.unique.precluster.dgc.0.03.rep.error.chimera ${run_anlys}/DGC/${mock}/${mock}.error.chimera


mkdir -p ${run_anlys}/Opti/${mock}
cp $workspace/mothur/Opti/${mock}.merge.unique.precluster.opti_mcc.0.03.rep.fasta ${run_anlys}/Opti/${mock}/${mock}.fasta
cp $workspace/mothur/Opti/${mock}.merge.unique.precluster.opti_mcc.0.03.rep.count_table ${run_anlys}/Opti/${mock}/${mock}_output.count
#cp $workspace/mothur/Opti/${mock}.merge.unique.precluster.opti_mcc.0.03.rep.error.chimera ${run_anlys}/Opti/${mock}/${mock}.error.chimera
done
done
