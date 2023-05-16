source variables.env

declare -a samples=("mock12_1" "mock13_1" "mock14_1" "mock15_1" "mock16_1" "mock18_1" "mock19_1" "mock20_1" "mock21_1" "mock22_1" "mock23_1" "mock24_1" "mock25_1" "mock4_1" "mock4_2" "mock4_3" "mock4_4" "mock5_1" "mock5_2" "mock5_3" "mock5_4")



for sample in ${samples[@]}; do 
echo ${sample}
$mothur "#set.logfile(name=${trunk}/${sample}.log);
fastq.info(fastq=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.fastq);
align.seqs(candidate=current, template=~/installs/silva.bacteria.fasta);
summary.seqs(fasta=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.align)"
done



$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock12_1.forward.trim.filt.align, start=13862, end=22581)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock18_1.forward.trim.filt.align, start=13862, end=22535)" 
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock19_1.forward.trim.filt.align, start=13862, end=22535)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock20_1.forward.trim.filt.align, start=13862, end=26147)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock21_1.forward.trim.filt.align, start=13862, end=26147)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock22_1.forward.trim.filt.align, start=13862, end=25434)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock23_1.forward.trim.filt.align, start=13862, end=25434)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock4_1.forward.trim.filt.align, start=13876, end=21768)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock4_2.forward.trim.filt.align, start=13876, end=21768)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock4_3.forward.trim.filt.align, start=13876, end=21768)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock4_4.forward.trim.filt.align, start=13876, end=21768)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock5_1.forward.trim.filt.align, start=13876, end=22585)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock5_2.forward.trim.filt.align, start=13876, end=22585)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock5_3.forward.trim.filt.align, start=13876, end=22585)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align
$mothur "#screen.seqs(fasta=${trunk}/input_fastqs_unm/mock5_4.forward.trim.filt.align, start=13876, end=22585)"
rm ${trunk}/input_fastqs_unm/mock*.trim.filt.good.align


for sample in ${samples[@]}; do 
echo ${sample}
$mothur "#remove.seqs(accnos=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.bad.accnos, fastq=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.fastq); 
remove.seqs(accnos=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.bad.accnos, fastq=${trunk}/input_fastqs/${sample}.trim.filt.fastq);"

cp ${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.bad.accnos ${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.bad.accnos
sed -E -i 's/(SRR[0-9]+\.[0-9]+\.)1/\12/' ${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.bad.accnos
$mothur "#remove.seqs(accnos=${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.bad.accnos, fastq=${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq)"


rm ${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.fastq
mv ${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.pick.fastq ${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.fastq

rm ${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq
mv ${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.pick.fastq ${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq

rm ${trunk}/input_fastqs/${sample}.trim.filt.fastq
mv ${trunk}/input_fastqs/${sample}.trim.filt.pick.fastq ${trunk}/input_fastqs/${sample}.trim.filt.fastq


$mothur "#fastq.info(fastq=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.fastq, qfile=F); 
sub.sample(fasta=current, size=25000); list.seqs(fasta=current); 
get.seqs(accnos=current, fastq=${trunk}/input_fastqs_unm/${sample}.forward.trim.filt.fastq); 
get.seqs(accnos=current, fastq=${trunk}/input_fastqs_unm/${sample}.reverse.trim.filt.fastq); 
get.seqs(accnos=current, fastq=${trunk}/input_fastqs/${sample}.trim.filt.fastq)"
done