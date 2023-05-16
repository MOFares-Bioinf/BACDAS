source variables.env

declare -a mocks=("mock4" "mock5" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23")
declare -a Runs=("Merged" "Single")



clustSpace=${anlystrunk}/clustered_refs
blastn=/home/mohamed/installs/ncbi-blast-2.10.0+/bin/blastn
makeblastdb=/home/mohamed/installs/ncbi-blast-2.10.0+/bin/makeblastdb



#### creating the blast reference databases for mock referrence--excuted once only 
mkdir -p ${trunk}/blast_dbs
for mock in ${mocks[@]}; do
cp ${trunk}/expected_seqs/${mock}.fasta ${trunk}/blast_dbs/ref_${mock}.fasta

#edit reference fasta sequence ids to add sequencial ids-ensure uniqueness-
awk 'NR % 2 {print ">Reference__",(NR+1)/2} !(NR % 2) {print}' ${trunk}/blast_dbs/ref_${mock}.fasta > tmp && mv tmp ${trunk}/blast_dbs/ref_${mock}.fasta
sed -i -- 's/ //g' ${trunk}/blast_dbs/ref_${mock}.fasta
#creat blast database from fasta file
$makeblastdb -in ${trunk}/blast_dbs/ref_${mock}.fasta -parse_seqids -title ${mock}_db -dbtype nucl
done
$makeblastdb -in ${trunk}/blast_dbs/silva.fasta -parse_seqids -title silva_db -dbtype nucl



#### creating clustered mock reference --excuted once only 

for mock in ${mocks[@]}; do 
mkdir -p ${clustSpace}/${mock}
cp ${trunk}/blast_dbs/ref_${mock}.fasta ${clustSpace}/${mock}
done


for mock in ${mocks[@]}; do 
echo ${mock}
$mothur "#set.logfile(name=${clustSpace}/clust_log, append=T);\
pcr.seqs(fasta=${clustSpace}/${mock}/ref_${mock}.fasta, oligos=${trunk}/expected_seqs/16sprimers.oligos);\
align.seqs(fasta=current, reference=${trunk}/expected_seqs/silva_seed/silva.seed_v138_1.align );\
filter.seqs(fasta=current, trump=., vertical=F);\
unique.seqs(fasta=current); \
summary.seqs(name=current); \
count.seqs(name=current); \
pre.cluster(fasta=current, count=current,align=gotoh); \
dist.seqs(fasta=current, cutoff=0.15); \
cluster(column=current, count=current, method=average, cutoff=0.15);\
make.shared(count=current, list=current, label=0.03);" 
done

for mock in ${mocks[@]}; do 
echo ${mock}
$mothur "#set.logfile(name=${clustSpace}/log.moth, append=T);\
get.oturep(count=${clustSpace}/${mock}/ref_${mock}.pcr.filter.unique.precluster.count_table, list=${clustSpace}/${mock}/ref_${mock}.pcr.filter.unique.precluster.an.list, fasta=${clustSpace}/${mock}/ref_${mock}.pcr.filter.unique.precluster.fasta, method=abundance, cutoff=0.03);" 
done

mkdir -p ${clustSpace}/reps/
for mock in ${mocks[@]}; do 
echo ${mock}
cp ${clustSpace}/${mock}/ref_${mock}.pcr.filter.unique.precluster.an.0.03.rep.fasta ${clustSpace}/reps/ref_${mock}.fasta
done


#### creating the blast reference databases clustered mock reference --excuted once only 
makeblastdb=/home/mohamed/installs/ncbi-blast-2.10.0+/bin/makeblastdb

for mock in ${mocks[@]}; do
sed -i -E 's/-//g' ${clustSpace}/reps/ref_${mock}.fasta 
sed -i -E 's/(>Reference__[0-9]+)..*/\1/g' ${clustSpace}/reps/ref_${mock}.fasta 
#creat blast database from fasta file
$makeblastdb -in ${clustSpace}/reps/ref_${mock}.fasta -parse_seqids -title ${mock}_db -dbtype nucl
done


###############  start of the analysis   ###############
declare -a tools=("DADA2" "deblur" "MED" "Unoise3" "Usearch" "DGC"  "AN" "Opti" "input")

for Run in ${Runs[@]}; do 
for tool in ${tools[@]}; do
rm -r ${anlystrunk}/$Run/${tool}/blast
done 
done

###### Composition  breakdown
for Run in ${Runs[@]}; do 
for tool in ${tools[@]}; do
mkdir -p ${anlystrunk}/$Run/${tool}/blast

for mock in ${mocks[@]}; do
mkdir -p ${anlystrunk}/$Run/${tool}/blast/${mock}

echo ${tool}__${mock}

#remove chimera
awk '$8 > 1 && NR != 1{print $1}' ${anlystrunk}/${Run}/${tool}/${mock}/${mock}.error.chimera > ${anlystrunk}/${Run}/${tool}/blast/${mock}/${mock}_chimeras.accnos
$mothur "#remove.seqs(accnos=${anlystrunk}/${Run}/${tool}/blast/${mock}/${mock}_chimeras.accnos, fasta=${anlystrunk}/$Run/$tool/${mock}/${mock}.fasta);"
mv ${anlystrunk}/$Run/$tool/${mock}/${mock}.pick.fasta ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.fasta

#align to mock ref exact match only
$blastn -db ${trunk}/blast_dbs/ref_${mock}.fasta \
-query ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.fasta \
-out ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_exact_mockref.blast \
-perc_identity 100 -num_alignments 1 -qcov_hsp_perc 100 \
-outfmt '6 mismatch pident length qcovhsp qseqid sseqid'
cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_exact_mockref.blast | cut -f5 -d$'\t' > ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exactmatch_only.accnos

#align to mock ref 97% similarity
$blastn -db ${trunk}/blast_dbs/ref_${mock}.fasta \
-query ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.fasta \
-out ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref.blast \
-perc_identity 97 -num_alignments 1 -qcov_hsp_perc 100 \
-outfmt '6 mismatch pident length qcovhsp qseqid sseqid'
cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref.blast | cut -f5 -d$'\t' > ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_all_97_mockref.accnos

#exclude exact match from 97%
comm -23 <(sort -u ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_all_97_mockref.accnos) <(sort -u ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exactmatch_only.accnos) | awk '{print $1}'> ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref_exactexcluded.accnos
$mothur "#remove.seqs(accnos=${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_all_97_mockref.accnos, fasta=${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.fasta);"
mv ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.pick.fasta ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.norefmatch.fasta

$mothur "#get.seqs(accnos=${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_all_97_mockref.accnos, fasta=${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.fasta)"
mv ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.pick.fasta ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exact_97_refmatch.fasta

done
done
done




### Identify non mock reference hits (Align to silva database)
for Run in ${Runs[@]}; do 
for tool in ${tools[@]}; do
for mock in ${mocks[@]}; do
#align and remove silva hits 
$blastn -db ${trunk}/expected_seqs/silva/SILVA_132_SSURef_Nr99_tax_silva.fasta \
-query ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.norefmatch.fasta \
-out ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.silva.blast \
-perc_identity 97 -num_alignments 1 -qcov_hsp_perc 100 -num_threads 4 \
-outfmt '6 mismatch pident length qcovhsp qseqid sseqid'
cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.silva.blast | cut -f5 -d$'\t' > ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_silva.accnos

#get other sequences (remaining non identified)
$mothur "#remove.seqs(accnos=${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_silva.accnos, fasta=${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.norefmatch.fasta); \
list.seqs(fasta=current)"

mv ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.norefmatch.pick.fasta ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.others.fasta
mv ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.nochime.norefmatch.pick.accnos ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.others.accnos

done
done
done



for Run in ${Runs[@]}; do 
for tool in ${tools[@]}; do
printf "%b" "mock\texactmatch\t97_mockref\t97_silva\tchimeras\tothers\ttotal_#_of_OTUs\n" > ${anlystrunk}/$Run/$tool/blast/${tool}_stats
for mock in ${mocks[@]}; do

#print per tool stats report 
printf "%b" "${mock}\t" "$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exactmatch_only.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref_exactexcluded.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_silva.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/${Run}/${tool}/blast/${mock}/${mock}_chimeras.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.others.accnos | wc -l)\t" \
"$(grep '>' ${anlystrunk}/$Run/$tool/${mock}/${mock}.fasta | wc -l)" \
"\n" >> ${anlystrunk}/$Run/$tool/blast/${tool}_stats

done
done
done


#print per tool stats report 
for Run in ${Runs[@]}; do 
printf "%b" "mock\ttool\texactmatch\t97_mockref\t97_silva\tchimeras\tothers\ttotal_#_of_OTUs\n" > ${anlystrunk}/$Run/compositionality_mocrobiota_otunumber
for mock in ${mocks[@]}; do
printf "%b" "${mock}\t" "\n" >> ${anlystrunk}/$Run/compositionality_mocrobiota_otunumber
echo ${mock}
for tool in ${tools[@]}; do

printf "%b" "${mock}\t" "${tool}\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exactmatch_only.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref_exactexcluded.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_silva.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/${Run}/${tool}/blast/${mock}/${mock}_chimeras.accnos | wc -l)\t" \
"$(cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.others.accnos | wc -l)\t" \
"$(grep '>' ${anlystrunk}/$Run/$tool/${mock}/${mock}.fasta | wc -l)" "\n" >> ${anlystrunk}/$Run/compositionality_mocrobiota_otunumber

done 
done
done



#####  Outpput reads abundance composition ######
for Run in ${Runs[@]}; do 
printf "%b" "mock\ttool\texactmatch\t97_mockref\t97_silva\tchimeras\tothers\t#_of_Reads\n" > ${anlystrunk}/$Run/compositionality_mocrobiota_readscount
for mock in ${mocks[@]}; do
printf "%b" "${mock}\t" "\n" >> ${anlystrunk}/$Run/compositionality_mocrobiota_readscount
echo ${mock}
for tool in ${tools[@]}; do



exact=$(awk 'NR==FNR{a[$1]=$2;next} {if($1 in a){print a[$1]} else {print 0}}' ${anlystrunk}/$Run/$tool/${mock}/${mock}_output.count ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exactmatch_only.accnos | awk '{s+=$1} END {if(s > 0){print s} else {print 0}}')
ref_97=$(awk 'NR==FNR{a[$1]=$2;next} {if($1 in a){print a[$1]} else {print 0}}' ${anlystrunk}/$Run/$tool/${mock}/${mock}_output.count ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref_exactexcluded.accnos | awk '{s+=$1} END {if(s > 0){print s} else {print 0}}')
silva_97=$(awk 'NR==FNR{a[$1]=$2;next} {if($1 in a){print a[$1]} else {print 0}}' ${anlystrunk}/$Run/$tool/${mock}/${mock}_output.count ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_silva.accnos | awk '{s+=$1} END {if(s > 0){print s} else {print 0}}')
chimera=$(awk 'NR==FNR{a[$1]=$2;next} {if($1 in a){print a[$1]} else {print 0}}' ${anlystrunk}/$Run/$tool/${mock}/${mock}_output.count ${anlystrunk}/${Run}/${tool}/blast/${mock}/${mock}_chimeras.accnos | awk '{s+=$1} END {if(s > 0){print s} else {print 0}}')
others=$(awk 'NR==FNR{a[$1]=$2;next} {if($1 in a){print a[$1]} else {print 0}}' ${anlystrunk}/$Run/$tool/${mock}/${mock}_output.count ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.others.accnos | awk '{s+=$1} END {if(s > 0){print s} else {print 0}}')
total=$(( $exact + $ref_97 + $silva_97 + $chimera + $others )) #$(grep '>' ${anlystrunk}/$Run/$tool/${mock}/${mock}.fasta | wc -l)

printf "%b" "${mock}\t" "${tool}\t" "${exact}\t" "${ref_97}\t" "${silva_97}\t" "${chimera}\t" "${others}\t" "${total}" "\n" >> ${anlystrunk}/$Run/compositionality_mocrobiota_readscount

done 
done
done



#######################################
# extract ref hits frequencies for specificity and sensitivity calculations
#
######################################


for Run in ${Runs[@]}; do 
for mock in ${mocks[@]}; do 
echo $mock
grep '^>' ${trunk}/blast_dbs/ref_${mock}.fasta | sed 's/>//' > ${anlystrunk}/$Run/${mock}_frequency
#| awk -v mock=$mock '{print mock,$0}'
for tool in ${tools[@]}; do

cat ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}_97_mockref.blast | cut -f6 -d$'\t' > tmp
awk 'NR==FNR {h[$1]++; next} h[$1] > 0 {print $0,h[$1]} ! h[$1] {print $0, 0}' tmp ${anlystrunk}/$Run/${mock}_frequency > tmp2 && mv tmp2 ${anlystrunk}/$Run/${mock}_frequency
done
awk -v mock=$mock '
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
        str=mock" "str
        print str
    }
}' ${anlystrunk}/$Run/${mock}_frequency > tmp && mv tmp ${anlystrunk}/$Run/${mock}_frequency
done
paste ${anlystrunk}/$Run/mock4_frequency \
${anlystrunk}/$Run/mock5_frequency \
${anlystrunk}/$Run/mock12_frequency \
${anlystrunk}/$Run/mock13_frequency \
${anlystrunk}/$Run/mock14_frequency \
${anlystrunk}/$Run/mock15_frequency \
${anlystrunk}/$Run/mock16_frequency \
${anlystrunk}/$Run/mock18_frequency \
${anlystrunk}/$Run/mock19_frequency \
${anlystrunk}/$Run/mock20_frequency \
${anlystrunk}/$Run/mock21_frequency \
${anlystrunk}/$Run/mock22_frequency \
${anlystrunk}/$Run/mock23_frequency > ${anlystrunk}/$Run/collective_frequency_$Run


for mock in ${mocks[@]}; do 
rm ${anlystrunk}/$Run/${mock}_frequency
done
done


############  Hits frequency to clustered mock reference analysis    #################



## in case of previous run clear the files 
for Run in ${Runs[@]}; do 
for tool in ${tools[@]}; do
rm -r ${anlystrunk}/$Run/${tool}/blast2
done 
done


for Run in ${Runs[@]}; do 
for tool in ${tools[@]}; do
mkdir -p ${anlystrunk}/$Run/${tool}/blast2
printf "%b" "mock\texactmatch\t97_mockref\t97_silva\tchimeras\tothers\ttotal\n" > ${anlystrunk}/$Run/$tool/blast2/${tool}_stats

for mock in ${mocks[@]}; do
mkdir -p ${anlystrunk}/$Run/${tool}/blast2/${mock}

echo ${tool}__${mock}

#align to clustered mock ref
$blastn -db ${clustSpace}/reps/ref_${mock}.fasta \
-query ${anlystrunk}/$Run/$tool/blast/${mock}/${mock}.exact_97_refmatch.fasta \
-out ${anlystrunk}/$Run/$tool/blast2/${mock}/${mock}_97_mockref.blast \
-num_alignments 1 \
-outfmt '6 mismatch pident length qcovhsp qseqid sseqid'


done
done
done



for Run in ${Runs[@]}; do 
for mock in ${mocks[@]}; do 
echo $mock
grep '^>' ${clustSpace}/reps/ref_${mock}.fasta | sed 's/>//' | sed -e 's/\tOtu.*//' > ${anlystrunk}/$Run/${mock}_clust_freq
#| awk -v mock=$mock '{print mock,$0}'
for tool in ${tools[@]}; do

cat ${anlystrunk}/$Run/$tool/blast2/${mock}/${mock}_97_mockref.blast | cut -f6 -d$'\t' > tmp
awk 'NR==FNR {h[$1]++; next} h[$1] > 0 {print $0,h[$1]} ! h[$1] {print $0, 0}' tmp ${anlystrunk}/$Run/${mock}_clust_freq > tmp2 && mv tmp2 ${anlystrunk}/$Run/${mock}_clust_freq
done
awk -v mock=$mock '
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
        str=mock" "str
        print str
    }
}' ${anlystrunk}/$Run/${mock}_clust_freq > tmp && mv tmp ${anlystrunk}/$Run/${mock}_clust_freq
done
paste ${anlystrunk}/$Run/mock4_frequency2 \
${anlystrunk}/$Run/mock5_frequency2 \
${anlystrunk}/$Run/mock12_frequency2 \
${anlystrunk}/$Run/mock13_frequency2 \
${anlystrunk}/$Run/mock14_frequency2 \
${anlystrunk}/$Run/mock15_frequency2 \
${anlystrunk}/$Run/mock16_frequency2 \
${anlystrunk}/$Run/mock18_frequency2 \
${anlystrunk}/$Run/mock19_frequency2 \
${anlystrunk}/$Run/mock20_frequency2 \
${anlystrunk}/$Run/mock21_frequency2 \
${anlystrunk}/$Run/mock22_frequency2 \
${anlystrunk}/$Run/mock23_frequency2 > ${anlystrunk}/$Run/collective_clust_freq_$Run

for mock in ${mocks[@]}; do 
rm ${anlystrunk}/$Run/${mock}_clust_freq
done
done


