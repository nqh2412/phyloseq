source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
biocLite()
biocLite('phyloseq')
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

path <- "/Users/claudiaximenarestrepo-ortiz/Desktop/Lizards"
list.files(path)


##Filter and Trim
# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1.fastq and SAMPLENAME_L001_R2.fastq
fnFs <- sort(list.files(path, pattern="_L001_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

##Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

##Perform filtering and trimming
#Assign the filenames for the filtered fastq.gz files.
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,195),
                     maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(25,20) , rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)



##Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Dereplication
#Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Sample Inference
#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Inspecting the dada-class object returned by dada
dadaFs[[1]]
dadaRs[[1]]

##Merge paired reads
#Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab) # 1= 0% seq from chimeras removed, 0,9=1% seq from chimeras removed...etc

##Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

##Assign taxonomy
taxaRC <- assignTaxonomy(seqtab.nochim, "/Users/claudiaximenarestrepo-ortiz/Desktop/silva_nr_v132_train_set.fa.gz", tryRC=TRUE)
taxaSp <- addSpecies(taxaRC, "/Users/claudiaximenarestrepo-ortiz/Desktop/silva_species_assignment_v132.fa.gz")


#taxonomic assignments
taxa.print <- taxaSp # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#fichier Physeq
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxaSp))

#change nom seq/OTU
taxa_names(physeq)
taxa_names(physeq) <- paste0("OTU", seq(ntaxa(physeq)))


# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

# Extract abundance matrix from the phyloseq object
TAX = as(tax_table(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){TAX <- t(TAX)}
# Coerce to data.frame
TAXdf = as.data.frame(TAX)

write.table(OTU1, "/Users/claudiaximenarestrepo-ortiz/Desktop/Lizards/OTU1.txt", sep="\t")
write.table(TAX, "/Users/claudiaximenarestrepo-ortiz/Desktop/Lizards/TAX.txt", sep="\t")
