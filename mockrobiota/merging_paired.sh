
declare -a mocks=("mock4_1" "mock4_2" "mock4_3" "mock4_4" "mock5_1" "mock5_2" "mock5_3" "mock5_4" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23" "mock24" "mock25" "mock27" "mock28" )
declare -a mocks=("mock16")
export usearch=~/installs/usearch
export mothur=~/installs/mothur
export trunk=/media/mohamed/8416E13916E12CBC/raw_fastq_files

mkdir -p ${trunk}/sorted
mkdir -p ${trunk}/merged
for mock in ${mocks[@]}; do

zcat ${trunk}/${mock}-forward-read.fastq.gz | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/sorted/${mock}-forward-read.fastq
zcat ${trunk}/${mock}-reverse-read.fastq.gz | paste - - - - | sort -V -k1,1 -t " " | tr '\t' '\n' > ${trunk}/sorted/${mock}-reverse-read.fastq
$usearch -fastq_mergepairs ${trunk}/sorted/${mock}-forward-read.fastq -reverse ${trunk}/sorted/${mock}-reverse-read.fastq -fastqout  ${trunk}/merged/${mock}.fastq 2>>merging.log
rm ${trunk}/sorted/${mock}-forward-read.fastq
rm ${trunk}/sorted/${mock}-reverse-read.fastq
$mothur "#fastq.info(fastq=${trunk}/merged/${mock}.fastq, qfile=F); summary.seqs(fasta=current);"
$mothur "#set.logfile(name=sequence_summary.log, append=T); summary.seqs(fasta=${trunk}/merged/${mock}.fasta);"
done
