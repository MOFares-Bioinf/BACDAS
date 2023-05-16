source variables.env

export raw_path="/home/mohamed/16s_4/Peter_mock"

mkdir -p ${trunk}/clean
mkdir -p ${trunk}/sorted
mkdir -p ${trunk}/trimmed
mkdir -p ${trunk}/dada_paired


declare -a samples=("A" "B" "C")


for sample in ${samples[@]}; do

#unzip and sort
zcat ${raw_path}/${sample}_R1_001.fastq.gz | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/trimmed/${sample}_S1_L001_R1_001.srt.fastq
zcat ${raw_path}/${sample}_R2_001.fastq.gz | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/trimmed/${sample}_S1_L001_R2_001.srt.fastq

#find best trimming params
figaro -i ${trunk}/sorted -o ${trunk}/trimmed -a 252 -f 17 -r 21 -F illumina -n ${sample}_trimparams.json

#check for forward and reverse primers
for file in ${raw_path}/*.fastq.gz; do 
#echo $file; 
zcat $file | head -n 40 | grep -E "CCTACGGG.GGC.GCAG" #forward original 
zcat $file | head -n 40 | grep -E "GACTAC..GGGTATCTAATCC" #reverse orig
done;


echo "stripping primer $sample"
python $cutprim --readsFile_r1  ${trunk}/trimmed/${sample}_S1_L001_R1_001.srt.fastq  --primersFileR1_5 ${trunk}/forward_primers  -tr1 ${trunk}/trimmed/${sample}-forward.cptrm.fastq -utr1 ${trunk}/trimmed/${sample}.fwd.utr.fastq
python $cutprim --readsFile_r1  ${trunk}/trimmed/${sample}_S1_L001_R2_001.srt.fastq  --primersFileR1_5 ${trunk}/reverse_primers  -tr1 ${trunk}/trimmed/${sample}-reverse.cptrm.fastq -utr1 ${trunk}/trimmed/${sample}.rvs.utr.fastq
head ${trunk}/trimmed/${sample}-forward.cptrm.fastq -n 40 | awk 'NR%4==2 {print;}' | wc -c

rm ${trunk}/trimmed/${sample}_S1_L001_R1_001.srt.fastq ${trunk}/trimmed/${sample}_S1_L001_R2_001.srt.fastq


echo "filtering by read length, sample ${sample},,,,,,,,,,,, "
$prinseq -fastq ${trunk}/trimmed/${sample}-forward.cptrm.fastq -out_good ${trunk}/trimmed/${sample}-forward.prfilt -min_len 284 -max_len 284
$prinseq -fastq ${trunk}/trimmed/${sample}-reverse.cptrm.fastq -out_good ${trunk}/trimmed/${sample}-reverse.prfilt -min_len 280 -max_len 280

rm ${trunk}/trimmed/${sample}-forward.cptrm.fastq  ${trunk}/trimmed/${sample}-reverse.cptrm.fastq

# Filter out sequences with bad quality and ambiguity characters
$usearch \
 -fastq_filter  ${trunk}/trimmed/${sample}-forward.prfilt.fastq \
 -fastq_maxns 0 \
 -fastq_maxee_rate 0.01 \
 -fastqout  ${trunk}/trimmed/${sample}-forward.usfilt.fastq
 
$usearch \
 -fastq_filter  ${trunk}/trimmed/${sample}-reverse.prfilt.fastq \
 -fastq_maxns 0 \
 -fastq_maxee_rate 0.01 \
 -fastqout  ${trunk}/trimmed/${sample}-reverse.usfilt.fastq
 
rm ${trunk}/trimmed/${sample}-forward.prfilt.fastq ${trunk}/trimmed/${sample}-reverse.prfilt.fastq


done

for sample in ${samples[@]}; do

$mothur "#set.logfile(name=${trunk}/${sample}.log);
fastq.info(fastq=${trunk}/trimmed/${sample}-forward.usfilt.fastq, qfile=F);
align.seqs(candidate=current, template=~/installs/silva.bacteria.fasta);
summary.seqs(fasta=current);
screen.seqs(fasta=current, start=6428, end=16351);
remove.seqs(accnos=${trunk}/trimmed/${sample}-forward.usfilt.bad.accnos, fastq=${trunk}/trimmed/${sample}-forward.usfilt.fastq); 
remove.seqs(accnos=${trunk}/trimmed/${sample}-forward.usfilt.bad.accnos, fastq=${trunk}/trimmed/${sample}-reverse.usfilt.fastq)"

rm ${trunk}/trimmed/${sample}-forward.usfilt.good.align ${trunk}/trimmed/${sample}-forward.usfilt.align
done

for sample in ${samples[@]}; do
comm -12 <(awk 'NR%4==1 {print $0}' ${trunk}/trimmed/${sample}-forward.usfilt.pick.fastq | sed -E 's/([ \t])1/\12/' | sort -u) <(awk 'NR%4==1 {print $0}' ${trunk}/trimmed/${sample}-reverse.usfilt.pick.fastq | sort -u) | awk '{print $1}' | sed -E 's/^@//' > ${trunk}/trimmed/${sample}-common.accnos
$mothur "#get.seqs(accnos=${trunk}/trimmed/${sample}-common.accnos, fastq=${trunk}/trimmed/${sample}-forward.usfilt.pick.fastq); 
get.seqs(accnos=${trunk}/trimmed/${sample}-common.accnos, fastq=${trunk}/trimmed/${sample}-reverse.usfilt.pick.fastq)"

rm ${trunk}/trimmed/${sample}-forward.usfilt.pick.fastq ${trunk}/trimmed/${sample}-reverse.usfilt.pick.fastq ${trunk}/trimmed/${sample}-common.accnos

done


for sample in ${samples[@]}; do
cat ${trunk}/trimmed/${sample}-forward.usfilt.pick.pick.fastq | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/trimmed/${sample}-forward.ornt.srt.fastq
cat ${trunk}/trimmed/${sample}-reverse.usfilt.pick.pick.fastq | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/trimmed/${sample}-reverse.ornt.srt.fastq

$usearch -fastq_mergepairs ${trunk}/trimmed/${sample}-forward.ornt.srt.fastq -reverse ${trunk}/trimmed/${sample}-reverse.ornt.srt.fastq -fastqout  ${trunk}/trimmed/${sample}.merged.fastq

$usearch -fastx_truncate ${trunk}/trimmed/${sample}.merged.fastq  -trunclen 400 -fastqout ${trunk}/trimmed/${sample}.merged.ustrm.fastq


rm ${trunk}/trimmed/${sample}-forward.usfilt.pick.pick.fastq ${trunk}/trimmed/${sample}-reverse.usfilt.pick.pick.fastq
done

for sample in ${samples[@]}; do
$mothur "#fastq.info(fastq=${trunk}/trimmed/${sample}.merged.ustrm.fastq, qfile=F); list.seqs(fasta=current); get.seqs(accnos=current, fastq=${trunk}/trimmed/${sample}-forward.ornt.srt.fastq); get.seqs(accnos=current, fastq=${trunk}/trimmed/${sample}-reverse.ornt.srt.fastq)"


rm ${trunk}/trimmed/${sample}-forward.ornt.srt.fastq ${trunk}/trimmed/${sample}-reverse.ornt.srt.fastq ${trunk}/trimmed/${sample}.merged.ustrm.fasta


cp ${trunk}/trimmed/${sample}.merged.ustrm.fastq ${trunk}/clean/${sample}.merged.prcss.fastq
cp ${trunk}/trimmed/${sample}-forward.ornt.srt.pick.fastq ${trunk}/clean/${sample}-forward.prcss.fastq
cp ${trunk}/trimmed/${sample}-reverse.ornt.srt.pick.fastq ${trunk}/clean/${sample}-reverse.prcss.fastq

done

export usearch=~/installs/usearch
declare -a samples=("mockP_1" "mockP_2" "mockP_3")
for sample in ${samples[@]}; do
#$usearch -fastx_truncate ${sample}.reverse.trim.filt.fastq  -stripleft 160 -fastqout ${sample}.reverse.trim.filt.trim.fastq
#$usearch -fastx_truncate ${sample}.trim.filt.fastq  -trunclen 284 -fastqout ${sample}.trim.filt.trim.fastq
$usearch -fastq_mergepairs ${sample}.forward.trim.filt.fastq -reverse ${sample}.reverse.trim.filt.trim.fastq -fastqout  ${sample}.trim.filt.fastq
#head -n 400 ${sample}.trim.filt.trim.fastq | awk 'NR%4==2' | wc -c
done