source variables.env

declare -a mocks=("mock4" "mock5" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23")
declare -a Runs=("Single" "Merged")


for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
usearch_output=${workspace}/usearch

mkdir -p $workspace/usearch/{otus,zotus,uniques,uparse,tabbedout,chimeras,input_fastqs}
cp -r ${workspace}/input_fastqs ${usearch_output}/

#add sample names to read headers
for i in {1..4}; do 
sed -i -r -- "s/(^@M141.*)/\1;sample=mock5${i};/g" ${usearch_output}/input_fastqs/mock5_${i}.trim.filt.fastq
done

for i in {1..4}; do 
sed -i -r -- "s/(^@M00112.*)/\1;sample=mock4${i};/g" ${usearch_output}/input_fastqs/mock4_${i}.trim.filt.fastq
done


sed -i -r -- "s/(\slength=2..)/\1;sample=mock12;/g" ${usearch_output}/input_fastqs/mock12_1.trim.filt.fastq
sed -i -r -- "s/(\s1\:N\:0\:1)/\1;sample=mock13;/g" ${usearch_output}/input_fastqs/mock13_1.trim.filt.fastq
sed -i -r -- "s/(\s1\:N\:0\:2)/\1;sample=mock14;/g" ${usearch_output}/input_fastqs/mock14_1.trim.filt.fastq
sed -i -r -- "s/(\s1\:N\:0\:3)/\1;sample=mock15;/g" ${usearch_output}/input_fastqs/mock15_1.trim.filt.fastq
sed -i -r -- "s/(\s1\:N\:0\:................)/\1;sample=mock16;/g" ${usearch_output}/input_fastqs/mock16_1.trim.filt.fastq
sed -i -r -- "s/(\slength=2..)/\1;sample=mock18;/g" ${usearch_output}/input_fastqs/mock18_1.trim.filt.fastq
sed -i -r -- "s/(\slength=2..)/\1;sample=mock19;/g" ${usearch_output}/input_fastqs/mock19_1.trim.filt.fastq
sed -i -r -- "s/(\slength=3..)/\1;sample=mock20;/g" ${usearch_output}/input_fastqs/mock20_1.trim.filt.fastq
sed -i -r -- "s/(\slength=3..)/\1;sample=mock21;/g" ${usearch_output}/input_fastqs/mock21_1.trim.filt.fastq
sed -i -r -- "s/(\slength=3..)/\1;sample=mock22;/g" ${usearch_output}/input_fastqs/mock22_1.trim.filt.fastq
sed -i -r -- "s/(\slength=3..)/\1;sample=mock23;/g" ${usearch_output}/input_fastqs/mock23_1.trim.filt.fastq
sed -i -r -- "s/(\slength=2..)/\1;sample=mock24;/g" ${usearch_output}/input_fastqs/mock24_1.trim.filt.fastq
sed -i -r -- "s/(\slength=2..)/\1;sample=mock25;/g" ${usearch_output}/input_fastqs/mock25_1.trim.filt.fastq
sed -i -r -- "s/(\slength=2..)/\1;sample=mock27;/g" ${usearch_output}/input_fastqs/mock27_1.trim.filt.fastq
sed -i -r -- "s/(\slength=2..)/\1;sample=mock28;/g" ${usearch_output}/input_fastqs/mock28_1.trim.filt.fastq

#remove illegal charachters in read headers (replace underscore with hyphen)
sed -i -r -- "s/(\@M00967\_47\_000000000)\_/\1\-/g" ${usearch_output}/input_fastqs/mock13_1.trim.filt.fastq
sed -i -r -- "s/(\@M00967\_47\_000000000)_/\1\-/g" ${usearch_output}/input_fastqs/mock14_1.trim.filt.fastq
sed -i -r -- "s/(\@M00967\_47\_000000000)_/\1\-/g" ${usearch_output}/input_fastqs/mock15_1.trim.filt.fastq
sed -i -r -- "s/(\@MISEQ\_17\_000000000)_/\1\-/g" ${usearch_output}/input_fastqs/mock16_1.trim.filt.fastq


for file in ${usearch_output}/input_fastqs/*.fastq; do
sample=$(basename $file | cut -d'.' -f1 | cut -d'_' -f2)
$usearch -fastq_filter $file -sample ${sample} -fastqout ${file}_labeld
rm $file
mv ${file}_labeld ${file}
done

for mock in ${mocks[@]}; do

cat ${usearch_output}/input_fastqs/${mock}_*  > ${usearch_output}/input_fastqs/${mock}merged.fastq;
rm ${usearch_output}/input_fastqs/${mock}_*
mv ${usearch_output}/input_fastqs/${mock}merged.fastq ${usearch_output}/input_fastqs/${mock}_merged.fastq

# Find unique read sequences and abundances
$usearch \
 -fastx_uniques ${usearch_output}/input_fastqs/${mock}_merged.fastq \
 -sizeout \
 -relabel Uniq \
 -fastqout ${usearch_output}/uniques/${mock}.uniques.fastq
done


####-----------------------------------------###########
####-----------   Uparse workflow -----------###########
####-----------------------------------------###########

for mock in ${mocks[@]}; do
$usearch\
 -cluster_otus ${usearch_output}/uniques/${mock}.uniques.fastq \
 -minsize 1 \
 -otus ${usearch_output}/otus/${mock}.otus.fa \
 -uparseout ${usearch_output}/uparse/${mock}.uparse.txt
done


# Retain chimeric sequences
for mock in ${mocks[@]}; do
echo $mock
count=0
cat ${usearch_output}/otus/${mock}.otus.fa  > ${usearch_output}/otus/${mock}.addchim.otu.fa
while [ -f ${usearch_output}/uparse/${mock}.uparse.txt ]
do
echo 'uparse file exists'
((count++))
chimeras=$(grep '_chimera' ${usearch_output}/uparse/${mock}.uparse.txt)
if [[ $chimeras = *[!\ ]* ]]; then
echo 'chimeras detected'
grep '_chimera' ${usearch_output}/uparse/${mock}.uparse.txt | awk '{print $1}' > ${usearch_output}/chimeras/${mock}_chim.accons
$mothur "#get.seqs(accnos=${usearch_output}/chimeras/${mock}_chim.accons, fastq=${usearch_output}/uniques/${mock}.uniques.fastq)"
mv ${usearch_output}/uniques/${mock}.uniques.pick.fastq ${usearch_output}/uniques/${mock}.chime.uniques.fastq
mv ${usearch_output}/uparse/${mock}.uparse.txt ${usearch_output}/uparse/${mock}.uparse.txt.$count
$usearch\
 -cluster_otus ${usearch_output}/uniques/${mock}.chime.uniques.fastq \
 -minsize 1 \
 -otus ${usearch_output}/otus/${mock}_chime.otus.fa \
 -uparseout ${usearch_output}/uparse/${mock}.uparse.txt 
cat ${usearch_output}/otus/${mock}_chime.otus.fa | sed "s/>Uniq/>Uniq0/" >> ${usearch_output}/otus/${mock}.addchim.otu.fa
else 
echo 'chimeras finished'
mv ${usearch_output}/uparse/${mock}.uparse.txt ${usearch_output}/uparse/${mock}.uparse.txt.$count
fi
done
done


##generate the otu table 
for mock in ${mocks[@]}; do
echo ${mock}
$usearch \
 -otutab ${usearch_output}/input_fastqs/${mock}_merged.fastq \
 -otus ${usearch_output}/otus/${mock}.addchim.otu.fa \
 -otutabout ${usearch_output}/otus/${mock}.otutab.txt \
 -mapout ${usearch_output}/otus/${mock}.map.txt 
done

done # end of Run loop
 

 
####-----------------------------------------###########
####-----------   Unoise workflow -----------###########
####-----------------------------------------###########

for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
usearch_output=${workspace}/usearch


 
# Make zotus
for mock in ${mocks[@]}; do
$usearch\
 -unoise3 ${usearch_output}/uniques/${mock}.uniques.fastq \
 -minsize 1 \
 -tabbedout ${usearch_output}/tabbedout/${mock}.tabbedout \
 -ampout ${usearch_output}/tabbedout/${mock}.ampout \
 -zotus ${usearch_output}/zotus/${mock}.zotus.fa
done


for mock in ${mocks[@]}; do
echo $mock
count=0
cat ${usearch_output}/zotus/${mock}.zotus.fa  > ${usearch_output}/zotus/${mock}.addchim.zotu.fa
while [ -f ${usearch_output}/tabbedout/${mock}.tabbedout ]
do
echo 'tabbedout file exists'
((count++))
chimeras=$(grep 'chimera' ${usearch_output}/tabbedout/${mock}.tabbedout)
if [[ $chimeras = *[!\ ]* ]]; then
echo 'chimeras detected'
grep 'chimera' ${usearch_output}/tabbedout/${mock}.tabbedout | awk '{print $1}' > ${usearch_output}/chimeras/${mock}_zchim.accons
$mothur "#get.seqs(accnos=${usearch_output}/chimeras/${mock}_zchim.accons, fastq=${usearch_output}/uniques/${mock}.uniques.fastq)"
mv ${usearch_output}/uniques/${mock}.uniques.pick.fastq ${usearch_output}/uniques/${mock}.zchime.uniques.fastq
mv ${usearch_output}/tabbedout/${mock}.tabbedout ${usearch_output}/tabbedout/${mock}.tabbedout.$count
$usearch\
 -unoise3 ${usearch_output}/uniques/${mock}.zchime.uniques.fastq \
 -minsize 1 \
 -tabbedout ${usearch_output}/tabbedout/${mock}.tabbedout \
 -zotus ${usearch_output}/zotus/${mock}_chime.zotus.fa
cat ${usearch_output}/zotus/${mock}_chime.zotus.fa | sed "s/>Zotu/>Zotu0${count}00/" >> ${usearch_output}/zotus/${mock}.addchim.zotu.fa
else 
echo 'chimeras finished'
mv ${usearch_output}/tabbedout/${mock}.tabbedout ${usearch_output}/tabbedout/${mock}.tabbedout.$count
fi
done
done



for mock in ${mocks[@]}; do
##generate the otu table 
$usearch \
 -otutab ${usearch_output}/input_fastqs/${mock}_merged.fastq \
 -otus ${usearch_output}/tabbedout/${mock}.ampout \
 -otutabout ${usearch_output}/zotus/${mock}.zotutab.txt \
 -mapout ${usearch_output}/chimeras/${mock}.map.txt 

done
done 
#end of Run loop ---------------------------------------------------------------------------------------


## Parse output into mothur format (.fasta and .count files)
for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
run_anlys=${anlystrunk}/${run}
usearch_output=${workspace}/usearch


for mock in ${mocks[@]}; do
mkdir -p ${run_anlys}/Usearch/${mock}
cp ${usearch_output}/otus/${mock}.addchim.otu.fa ${run_anlys}/Usearch/${mock}/${mock}.fasta
cp ${usearch_output}/otus/${mock}.otutab.txt ${run_anlys}/Usearch/${mock}/${mock}_output.count
sed -i -r -e 's/\;size.+;//' ${run_anlys}/Usearch/${mock}/${mock}.fasta

sed -r -e 's/\#OTU ID/Representative_Sequence/' ${run_anlys}/Usearch/${mock}/${mock}_output.count | awk 'NR==1 {print $1,"total",$2,$3,$4,$5} NR!=1 {print $1,$2+$3+$4+$5,$2,$3,$4,$5}' | tr -s " " '\t' > tmp && mv tmp ${run_anlys}/Usearch/${mock}/${mock}_output.count


mkdir -p ${run_anlys}/Unoise3/${mock}

sed -r -e 's/^(>Uniq[0-9]+);.+/\1/' ${usearch_output}/tabbedout/${mock}.ampout >${run_anlys}/Unoise3/${mock}/${mock}.fasta
cp ${usearch_output}/zotus/${mock}.zotutab.txt ${run_anlys}/Unoise3/${mock}/${mock}_output.count

sed -r -e 's/\#OTU ID/Representative_Sequence/' ${run_anlys}/Unoise3/${mock}/${mock}_output.count | awk 'NR==1 {print $1,"total",$2,$3,$4,$5} NR!=1 {print $1,$2+$3+$4+$5,$2,$3,$4,$5}' | tr -s " " '\t' > tmp && mv tmp ${run_anlys}/Unoise3/${mock}/${mock}_output.count


$mothur "#set.dir(output=${run_anlys}/Usearch/${mock});
set.logfile(name=${run_anlys}/Usearch/usearch_${run}_log, append=T);
seq.error(fasta=${run_anlys}/Usearch/${mock}/${mock}.fasta, count=${run_anlys}/Usearch/${mock}/${mock}_output.count, reference=~/expected_seqs/${mock}.fasta, aligned=F, ignorechimeras=T);
seq.error(fasta=${run_anlys}/Usearch/${mock}/${mock}.fasta, count=${run_anlys}/Usearch/${mock}/${mock}_output.count, reference=~/expected_seqs/${mock}.fasta, aligned=F, ignorechimeras=F)"


$mothur "#set.dir(output=${run_anlys}/Unoise3/${mock}); 
set.logfile(name=${run_anlys}/Unoise3/Unoise3_${run}_log, append=T);
seq.error(fasta=${run_anlys}/Unoise3/${mock}/${mock}.fasta, count=${run_anlys}/Unoise3/${mock}/${mock}_output.count, reference=~/expected_seqs/${mock}.fasta, aligned=F, ignorechimeras=T);
seq.error(fasta=${run_anlys}/Unoise3/${mock}/${mock}.fasta, count=${run_anlys}/Unoise3/${mock}/${mock}_output.count, reference=~/expected_seqs/${mock}.fasta, aligned=F, ignorechimeras=F)"
done

done

#end of Run loop -------------------------------------------------------------------------------------


#checking correctness of count files and that they contain all generated otus
###############################################
#check sequence correctness 
################################################

for mock in ${mocks[@]}; do 
echo $mock
awk 'NR != 1 {print $3}' ${run_anlys}/Unoise3/${mock}/${mock}.error.summary | paste -sd+ | bc
awk 'NR != 1 {print $2}' ${run_anlys}/Unoise3/${mock}/${mock}_output.count | paste -sd+ | bc
done
done


for mock in ${mocks[@]}; do 
echo $mock
comm -3 <(grep '>' ${run_anlys}/Unoise3/${mock}/${mock}.fasta | sed 's/>//' | sort -u) <(awk 'NR != 1 {print $1}' ${run_anlys}/Unoise3/${mock}/${mock}_output.count | sort -u) #> Unoise3_${mock}.missing

done
