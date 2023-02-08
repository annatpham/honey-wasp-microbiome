# Sources: 
## https://benjjneb.github.io/dada2/tutorial_1_6.html
## http://web.stanford.edu/class/bios221/MicrobiomeWorkflowII.html

# Install dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

getwd()
setwd("/Users/lexiemartin/Desktop")

library(dada2); packageVersion("dada2")

# Getting to the Fastq Files ----------------------------------------------
## They should all be in a folder, set the "path" to be the path get to the files
mypath <- "/Users/lexiemartin/Desktop/5 April 2022 Files"
list.files(mypath) # Check to see all the files are there

# Filter and Trim ---------------------------------------------------------
# Forward and reverse fastq filenames are in the format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
## We want to  perform string manipulations to get matched lists of the forward and reverse fastq files
fnFs <- sort(list.files(mypath, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(mypath, pattern = "_R2_001.fastq", full.names = TRUE))
fnFs[1:3]
fnRs[1:3]

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,1)

# Visualize the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])
# Visualize the quality profiles of the reverse reads
plotQualityProfile(fnRs[1:2])

# Assign the filenames for the filtered fastq.gz files
myfilt_path <- file.path(mypath, "filtered") # Place filtered files in filtered/subdirectory
if(!file_test("-d", myfilt_path)) dir.create(myfilt_path)
filtFs <- file.path(mypath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(mypath, "filtered", paste0(sample.names, "_R_filt_fastq.gz"))

# Use standard filtering parameters
# Truncate forward reads at position 230 and reverse reads at position 240 based on quality profiles
# Remove 515F and 806R primers using truncLen
# Adapters already removed by GSAF
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,240), trimLeft= c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
# On Windows set multithread=FALSE 
head(out)

# Infer Sequence Variants -------------------------------------------------
# Dereplication 
## Combines all identical sequencing reads into "unique sequences" with a corresponding "abundance": the # of reads with that unique sequence
## Reduces computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
## Name the derep-class objects by the sample names 
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Error Rates
## Estimate error rates in forward and reverse reads
errF <-learnErrors(filtFs, multithread=TRUE)
errR <-learnErrors(filtRs, multithread=TRUE)
## Plot fitted error rate and observed error rate for each DNA substitution
plotErrors(errF)
plotErrors(errR) 

# Sample Inference
## Applying the core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
## Inspecting the returned dada-class object
dadaFs[[1]]

# Construct Sequence Table and Remove Chimeras ----------------------------
# Merge paired reads  
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# Create ASV table
seqtabAll <- makeSequenceTable(mergers) 
dim(seqtabAll)
table(nchar(getSequences(seqtabAll))) # Inspect the distribution of sequence lengths

#Remove Chimeras
seqtabNoC <- removeBimeraDenovo(seqtabAll)  

# Track Read Counts Through Pipeline --------------------------------------
getN <- function(x) sum(getUniques(x)) # Function to get unique reads
track <- 
  cbind(out, 
        sapply(dadaFs, getN), 
        sapply(dadaRs, getN), 
        sapply(mergers, getN), 
        rowSums(seqtabNoC))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Download as a data frame
dada2df <- data.frame(track) 
write.csv(dada2df, "./dada2df.csv") 

# Assign Taxonomy ---------------------------------------------------------
fastaRef <- "./silva_nr99_v138.1_train_set.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

# Install DECIPHER and phyloseq
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
library(DECIPHER)
library(phyloseq)

# Install packages
install.packages("phangorn")
library(phangorn)
install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library("gridExtra")

# Construct Phylogenetic Tree ---------------------------------------------
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phangAlign)  
treeNJ <- NJ(dm) #Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma = TRUE, rearrangment = "stochastic", control= pml.control(trace=0))
detach("package:phangorn",  unload=TRUE)

# Combine Data into a Phyloseq Object -------------------------------------
# Load metadata
honey <- read.csv(file = "Metadata.csv")
rownames(honey) <- honey$New_ID # Change row names to the sample labels used in metadata
all(rownames(seqtabNoC) %in% honey$New_ID) # Check that sample labels match

# Create phyloseq object
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(honey), 
               tax_table(taxTab),
               phy_tree(fitGTR$tree))
ps

# Make short names for ASVs while storing sequence in a new column
tax <- data.frame(tax_table(ps))
tax$Sequence = rownames(tax) # Make a Sequence column
tax_table(ps) <- as.matrix(tax) # Replace with new taxonomy table
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # Make short names for ASVs
