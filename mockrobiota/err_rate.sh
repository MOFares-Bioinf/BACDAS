source variables.env

declare -a mocks=("mock4" "mock5" "mock12" "mock13" "mock14" "mock15" "mock16" "mock18" "mock19" "mock20" "mock21" "mock22" "mock23")
declare -a tools=("DADA2" "deblur" "MED" "Unoise3" "Usearch" "DGC"  "AN" "Opti" )
declare -a Runs=("Merged" "Single")



for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
run_anlys=${anlystrunk}/${run}

for tool in ${tools[@]}; do

for mock in ${mocks[@]}; do
$mothur "#set.dir(output=${run_anlys}/${tool}/${mock}); \
set.logfile(name=${run_anlys}/${tool}/${tool}_${run}_log, append=T);\
seq.error(fasta=${run_anlys}/${tool}/${mock}/${mock}.fasta, count=${run_anlys}/${tool}/${mock}/${mock}_output.count, reference=~/expected_seqs/${mock}.fasta, aligned=F, ignorechimeras=T)\;
seq.error(fasta=${run_anlys}/${tool}/${mock}/${mock}.fasta, count=${run_anlys}/${tool}/${mock}/${mock}_output.count, reference=~/expected_seqs/${mock}.fasta, aligned=F, ignorechimeras=F)"
done
done
done


declare -a tools=("DADA2" "deblur" "MED" "Unoise3" "Usearch" "DGC"  "AN" "Opti" )
for run in ${Runs[@]}; do 
workspace=${trunk}/${run}
run_anlys=${anlystrunk}/${run}

for i in {1..26}; do
for tool in ${tools[@]}; do
cat ${run_anlys}/${tool}/${tool}_${run}_log | grep 'Overall error rate:' | sed 's/Overall error rate\:\t//' | awk -v i=$i 'NR==i' >> ${run_anlys}/combined_err
done
done

awk '
  BEGIN {
    # Numbers of rows to print
    n=16;
  }
  {
    # Add to array with key = 0, 1, 2, 3, 0, 1, 2, ..
    l[(NR-1)%n] = l[(NR-1)%n] " " $0
  };
  END {
    # print the array
    for (i = 0; i < length(l); i++) {
      print l[i];
    }
  }
' ${run_anlys}/combined_err > ${run_anlys}/combined_err_sep
done