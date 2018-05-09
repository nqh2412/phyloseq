#Bioconductor 3.7 R 3.5.0
version
#install bioconductor and independencies 
source("http://bioconductor.org/biocLite.R")
biocLite()
#install "phyloseq" package
biocLite("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")

#Examples 1
data(GlobalPatterns)
data(esophagus)
data(enterotype)
data(soilrep)
?GlobalPatterns
example("enterotype", ask=FALSE)
#### DADA2 pipeline tutorial 
#install dada2 package required Rcpp package
biocLite("Rcpp")
biocLite("dada2")
library("dada2")
#import data / + tab it will show the directories
path <- "/Users/lananh/Documents/nqhuynh/MiSeq_SOP/"
list.files(path)

## Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
list.files(path)
## Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
##quality profiles of reads
##forwards (R1= forward)
plotQualityProfile(fnFs[1:2])
##reverse (R2= reverse)
plotQualityProfile(fnRs[1:2])
###filter and trim the sq 
## place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), 
#truncQ=2, rm.phix=TRUE and maxEE=2. 
#The maxEE parameter sets the maximum number of “expected errors” allowed in a read, 
#which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
head(out)
#Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ = TRUE)
# De-replication combine all identical sequencing read to unique squences 
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Nam the derep-class objects by sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
#Sample interference
dadaFs <- dada(derepFs, err=errF, multithread = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread = TRUE)
#inspecting the return data class
dadaFs[[1]]
#merge paired reads 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[[1]])
##construct sq table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#inspect distribution of sq lenghts
table(nchar(getSequences(seqtab)))
#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "mergerd", "nochim")
rownames(track) <- sample.names
head(track)
