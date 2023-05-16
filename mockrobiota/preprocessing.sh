source ~/cleaned_scripts/mockrobiota/Merged/variables.env

export raw_path="/media/mohamed/Prog/Mohamed/USB_Drive/raw_fastq_files"

mkdir -p ${trunk}/clean
mkdir -p ${trunk}/sorted
mkdir -p ${trunk}/trimmed
mkdir -p ${trunk}/dada_paired

declare -a mocks=("mock4_1" "mock4_2" "mock4_3" "mock4_4" "mock5_1" "mock5_2" "mock5_3" "mock5_4" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23")



for mock in ${mocks[@]}; do

#unzip and sort
zcat ${raw_path}/${mock}-forward-read.fastq.gz | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq
zcat ${raw_path}/${mock}-reverse-read.fastq.gz | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq

#find best trimming params
figaro -i ${trunk}/sorted -o ${trunk}/trimmed -a 252 -f 1 -r 1 -F illumina -n ${mock}_trimparams.json


#check for retained forward and reverse primers
if [ ${mock} == 'mock16' ] || [ ${mock} == 'mock18' ] || [ ${mock} == 'mock19' ] || [ ${mock} == 'mock22' ] || [ ${mock} == 'mock23' ] 
then 
echo "forward"
#zcat ${raw_path}/${mock}-forward-read.fastq.gz | head -n 40 | grep -E "GTGCCAGC.GCCGCGGTAA" #forward original 
echo "reverse_reads"
#zcat ${raw_path}/${mock}-reverse-read.fastq.gz | head -n 40 | grep -E "ACTAC..GGGT.TCTAAT" #reverse orig
fi




if [ ${mock} == 'mock16' ] || [ ${mock} == 'mock18' ] || [ ${mock} == 'mock19' ] || [ ${mock} == 'mock22' ] || [ ${mock} == 'mock23' ] 
then 
#echo "stripping primer $mock"
python $cutprim --readsFile_r1  ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq  --primersFileR1_5 ${trunk}/forward_primers  -tr1 ${trunk}/trimmed/${mock}-forward-read.fastq -utr1 ${trunk}/trimmed/${mock}.fwd.utr.fastq
python $cutprim --readsFile_r1  ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq  --primersFileR1_5 ${trunk}/reverse_primers  -tr1 ${trunk}/trimmed/${mock}-reverse-read.fastq -utr1 ${trunk}/trimmed/${mock}.rvs.utr.fastq

head ${trunk}/trimmed/${mock}-forward-read.fastq -n 40 | awk 'NR%4==2 {print;}' | wc -c

rm ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq
rm ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq

mv ${trunk}/trimmed/${mock}-forward-read.fastq ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq
mv ${trunk}/trimmed/${mock}-reverse-read.fastq ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq

fi

#rm ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq
#rm ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq


#echo "filtering by read length, sample ${mock},,,,,,,,,,,, "
$prinseq -fastq ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq -out_good ${trunk}/clean/${mock}_S1_L001_R1_001 -min_len 145 -max_len 258
$prinseq -fastq ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq -out_good ${trunk}/clean/${mock}_S1_L001_R2_001 -min_len 145 -max_len 258

echo "trimming read length, sample ${mock},,,,,,,,,,,, "
# global trimming in merged paired end reads is not necessary as properly merged reads should end at exact locus in the gene 
$usearch\
 -fastx_truncate ${trunk}/${trunk}/clean/${mock}_S1_L001_R1_001.fastq \
 -trunclen 297 \
 -fastqout ${trunk}/clean/${mock}-forward-read.trim.fastq 

rm ${trunk}/sorted/${mock}_S1_L001_R1_001.fastq
rm ${trunk}/clean/${mock}_S1_L001_R1_001.fastq

$usearch\
 -fastx_truncate ${trunk}/clean/${mock}_S1_L001_R2_001.fastq \
 -trunclen 130 \
 -fastqout ${trunk}/clean/${mock}-reverse-read.trim.fastq 

rm ${trunk}/sorted/${mock}_S1_L001_R2_001.fastq
rm ${trunk}/clean/${mock}_S1_L001_R2_001.fastq


# Filter out sequences with bad quality and ambiguity characters
#forward
$usearch \
 -fastq_filter  ${trunk}/clean/${mock}-forward-read.trim.fastq \
 -fastq_maxns 0 \
 -fastq_maxee_rate 0.01 \
 -fastqout  ${trunk}/clean/${mock}-forward-read.trim.filt.fastq
 
rm ${trunk}/clean/${mock}-forward-read.trim.fastq

#reverse
$usearch \
 -fastq_filter  ${trunk}/clean/${mock}-reverse-read.trim.fastq \
 -fastq_maxns 0 \
 -fastq_maxee_rate 0.01 \
 -fastqout  ${trunk}/clean/${mock}-reverse-read.trim.filt.fastq
 
rm ${trunk}/clean/${mock}-reverse-read.trim.fastq

done



if [ ${mock} == 'mock12' ]
then
comm -12 <(awk 'NR%4==1 {print $0}' ${trunk}/clean/${mock}-forward-read.trim.filt.fastq | sed -E 's/(SRR[0-9]+\.[0-9]+\.)1/\12/'|sort -u) <(awk 'NR%4==1 {print $0}' ${trunk}/clean/${mock}-reverse-read.trim.filt.fastq | sort -u) | awk '{print $1}' | sed -E 's/^@//' > ${trunk}/clean/${mock}-common.accnos
$mothur "#get.seqs(accnos=${trunk}/clean/${mock}-common.accnos, fastq=${trunk}/clean/${mock}-reverse-read.trim.filt.fastq)"
cp ${trunk}/clean/${mock}-common.accnos ${trunk}/clean/${mock}-common_f.accnos
sed -E -i 's/(SRR[0-9]+\.[0-9]+\.)2/\11/' ${trunk}/clean/${mock}-common_f.accnos
$mothur "#get.seqs(accnos=${trunk}/clean/${mock}-common_f.accnos, fastq=${trunk}/clean/${mock}-forward-read.trim.filt.fastq);"
else
comm -12 <(awk 'NR%4==1 {print $0}' ${trunk}/clean/${mock}-forward-read.trim.filt.fastq | sed -E 's/([ \t])1/\12/' | sort -u) <(awk 'NR%4==1 {print $0}' ${trunk}/clean/${mock}-reverse-read.trim.filt.fastq | sort -u) | awk '{print $1}' | sed -E 's/^@//' > ${trunk}/clean/${mock}-common.accnos
$mothur "#get.seqs(accnos=${trunk}/clean/${mock}-common.accnos, fastq=${trunk}/clean/${mock}-forward-read.trim.filt.fastq); get.seqs(accnos=${trunk}/clean/${mock}-common.accnos, fastq=${trunk}/clean/${mock}-reverse-read.trim.filt.fastq)"

fi

rm ${trunk}/clean/${mock}-forward-read.trim.filt.fastq
rm ${trunk}/clean/${mock}-reverse-read.trim.filt.fastq
rm ${trunk}/clean/${mock}-common.accnos



cat ${trunk}/clean/${mock}-forward-read.trim.filt.pick.fastq | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/clean/${mock}-forward-sorted.fastq
cat ${trunk}/clean/${mock}-reverse-read.trim.filt.pick.fastq | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/clean/${mock}-reverse-sorted.fastq
$usearch -fastq_mergepairs ${trunk}/clean/${mock}-forward-sorted.fastq -reverse ${trunk}/clean/${mock}-reverse-sorted.fastq -fastqout  ${trunk}/clean/${mock}_merged.fastq

rm ${trunk}/clean/${mock}-forward-read.trim.filt.pick.fastq
rm ${trunk}/clean/${mock}-reverse-read.trim.filt.pick.fastq



if [ ${mock} == 'mock12' ]
then
$mothur "#fastq.info(fastq=${trunk}/clean/${mock}_merged.fastq, qfile=F); 
sub.sample(fasta=current, size=30000); 
list.seqs(fasta=current); 
get.seqs(accnos=current, fastq=${trunk}/clean/${mock}_merged.fastq); 
get.seqs(accnos=current, fastq=${trunk}/clean/${mock}-forward-sorted.fastq)"

cp ${trunk}/clean/${mock}_merged.subsample.accnos ${trunk}/clean/${mock}_merged.subsample_r.accnos
sed -E -i 's/(SRR[0-9]+\.[0-9]+\.)1/\12/' ${trunk}/clean/${mock}_merged.subsample_r.accnos
$mothur "#get.seqs(accnos=${trunk}/clean/${mock}_merged.subsample_r.accnos, fastq=${trunk}/clean/${mock}-reverse-sorted.fastq)"
else

$mothur "#fastq.info(fastq=${trunk}/clean/${mock}_merged.fastq, qfile=F); 
sub.sample(fasta=current, size=30000); 
list.seqs(fasta=current); 
get.seqs(accnos=current, fastq=${trunk}/clean/${mock}_merged.fastq); 
get.seqs(accnos=current, fastq=${trunk}/clean/${mock}-forward-sorted.fastq); 
get.seqs(accnos=current, fastq=${trunk}/clean/${mock}-reverse-sorted.fastq)"
fi

rm ${trunk}/clean/${mock}-forward-sorted.fastq
rm ${trunk}/clean/${mock}-reverse-sorted.fastq
rm ${trunk}/clean/${mock}_merged.fastq
rm ${trunk}/clean/${mock}_merged.fasta
rm ${trunk}/clean/${mock}_merged.subsample.fasta
rm ${trunk}/clean/${mock}_merged.subsample.accnos

mv ${trunk}/clean/${mock}_merged.pick.fastq ${trunk}/clean/${mock}_merged.fastq
mv ${trunk}/clean/${mock}-forward-sorted.pick.fastq ${trunk}/clean/${mock}-forward-reads.fastq
mv ${trunk}/clean/${mock}-reverse-sorted.pick.fastq ${trunk}/clean/${mock}-reverse-reads.fastq


$usearch -fastq_mergepairs ${trunk}/clean/${mock}-forward-reads.fastq -reverse ${trunk}/clean/${mock}-reverse-reads.fastq -fastqout  ${trunk}/clean/${mock}_merged_test.fastq


done

mkdir -p ${trunk}/{input_fastqs,input_fastqs_unm}

declare -a mocks=("mock4_1" "mock4_2" "mock4_3" "mock4_4" "mock5_1" "mock5_2" "mock5_3" "mock5_4")
for mock in ${mocks[@]}; do
cp ${trunk}/clean/${mock}_merged.fastq ${trunk}/input_fastqs/${mock}.trim.filt.fastq 
cp ${trunk}/clean/${mock}-forward-reads.fastq ${trunk}/input_fastqs_unm/${mock}.forward.trim.filt.fastq
cp ${trunk}/clean/${mock}-reverse-reads.fastq ${trunk}/input_fastqs_unm/${mock}.reverse.trim.filt.fastq 
done

declare -a mocks=("mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23" "mock24" "mock25" "mock27" "mock28")
for mock in ${mocks[@]}; do
cp ${trunk}/clean/${mock}_merged.fastq ${trunk}/input_fastqs/${mock}_1.trim.filt.fastq
cp ${trunk}/clean/${mock}-forward-reads.fastq ${trunk}/input_fastqs_unm/${mock}_1.forward.trim.filt.fastq
cp ${trunk}/clean/${mock}-reverse-reads.fastq ${trunk}/input_fastqs_unm/${mock}_1.reverse.trim.filt.fastq 
done


#######
