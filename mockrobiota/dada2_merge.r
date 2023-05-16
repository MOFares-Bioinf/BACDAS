#!/usr/bin/Rscript

library(dada2)
packageVersion("dada2")


path <- "~/16s_4/Merged/DADA2" #directory containing the fastq files after unzipping.
list.files(path)

mocks <- c("mock4", "mock5", "mock12", "mock13", "mock14", "mock15", "mock16", "mock18", "mock19", "mock20", "mock21", "mock22", "mock23")


for (mock in mocks) {
print(mock)


# Forward and reverse fastq filenames have format: SAMPLENAME.forward.trim.filt.fastq and SAMPLENAME.reverse.trim.filt.fastq
fnFs <- sort(list.files(paste(path, mock, sep="/"), pattern=".forward.trim.filt.fastq", full.names = TRUE))
fnRs <- sort(list.files(paste(path, mock, sep="/"), pattern=".reverse.trim.filt.fastq", full.names = TRUE))
  
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "\\.forward"), `[`, 1)
  
  
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, mock, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, mock, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))



names(fnFs) <- sample.names
names(fnRs) <- sample.names
names(filtFs) <- sample.names
names(filtRs) <- sample.names

  
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, truncQ=0, compress=TRUE, multithread=TRUE) 
head(out)
  
# Learn error rates
  
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
 
setDadaOpt(OMEGA_C = 0)
	
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
  
  
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
  
  
# Construct sequence table and write to disk
subdir <- paste("output",mock, sep="_")
outputfolder <- paste(path, subdir, sep="/")
dir.create(file.path(path, subdir))

#rds_fp <- paste(outputfolder, "seqtab.rds", sep="/")
#saveRDS(seqtab, rds_fp) 

tcsv_fp <-paste(outputfolder, "transposed.csv", sep="/")
write.csv(t(seqtab), file = tcsv_fp)
}
