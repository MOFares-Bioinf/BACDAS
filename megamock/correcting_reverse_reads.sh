source variables.env



declare -a samples=("mockP_1" "mockP_2" "mockP_3")


for sample in ${samples[@]}; do 
echo ${sample}
$mothur "#set.logfile(name=${workspace}/${sample}.log);
fastq.info(fastq=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.fastq);
align.seqs(candidate=current, template=~/installs/silva.bacteria.fasta);
summary.seqs(fasta=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.align)"
done

for sample in ${samples[@]}; do 
echo ${sample}
#rm ${workspace}/${run}/input_fastqs/*
$mothur "#screen.seqs(fasta=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.align, start=6428, end=16351)"
rm ${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.good.align
done


for sample in ${samples[@]}; do 
echo ${sample}
$mothur "#remove.seqs(accnos=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.bad.accnos, fastq=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.fastq); 
remove.seqs(accnos=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.bad.accnos, fastq=${workspace}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq); 
remove.seqs(accnos=${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.bad.accnos, fastq=${workspace}/input_fastqs/${sample}.trim.filt.fastq);"
done

for sample in ${samples[@]}; do 
rm ${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.fastq
mv ${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.pick.fastq ${workspace}/input_fastqs_unm/${sample}.forward.trim.filt.fastq

rm ${workspace}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq
mv ${workspace}/input_fastqs_unm/${sample}.reverse.trim.filt.pick.fastq ${workspace}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq

rm ${workspace}/input_fastqs/${sample}.trim.filt.fastq
mv ${workspace}/input_fastqs/${sample}.trim.filt.pick.fastq ${workspace}/input_fastqs/${sample}.trim.filt.fastq

done

