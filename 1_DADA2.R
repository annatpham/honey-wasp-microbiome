if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

getwd()
setwd("C:/Users/apham/Documents") #setwd("/Users/lexiemartin/Desktop")

library(dada2); packageVersion("dada2")
#Getting to the files
#They should all be in a folder on the computer, set the "path" to be the path get to  the files
mypath <- "./5 April 2022 Files"
list.files(mypath) #check to see all the files are there

#Forward and reverse fastq filenames are in the format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#May have to change this code of the filename format is different
#We want to  perform string manipulations to get matched lists of the forward and reverse fastq files
fnFs <- sort(list.files(mypath, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(mypath, pattern = "_R2_001.fastq", full.names = TRUE))
fnFs[1:3]
fnRs[1:3]

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <-sapply(strsplit(basename(fnFs), "_"), `[`,1)

#Visualize the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])
#Visualize the quality profiles of the reverse reads
plotQualityProfile(fnRs[8:18])

#Filter and Trim
#Assign the filenmaes for the filtered fastq.gz files
#Place filtered files in filtered/ subdirectory
myfilt_path <-file.path(mypath, "filtered") #place filtered files in filtered/subdirectory
if(!file_test("-d", myfilt_path)) dir.create(myfilt_path)
filtFs <- file.path(mypath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(mypath, "filtered",paste0(sample.names, "_R_filt_fastq.gz"))

#Adapters already removed by GSAF
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,240), trimLeft= c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
#On Mac set multithread=TRUE  
head(out)

##Infer Sequence Variants
#Dereplication
#depreplication combines all identical sequencing reads into "unique sequences" with a corresponding "abundance": the # of reads with that unique sequence
#reduces computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
#name the derep-class objects by the sample names 
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Error Rates
errF <-learnErrors(filtFs, multithread=TRUE)
errR <-learnErrors(filtRs, multithread=TRUE)
plotErrors(errF)
plotErrors(errR) 

#Sample Inference
#Applying the core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#Inspecting the returned dada-class object
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

#Create ASV table
seqtabAll <- makeSequenceTable(mergers[!grepl("Blank", "qblnk1", "qblnk2", names(mergers))])
dim(seqtabAll)
table(nchar(getSequences(seqtabAll)))

#Remove Chimeras
seqtabNoC <- removeBimeraDenovo(seqtabAll)  

# Track read counts through pipeline
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

#Download as a data frame
dada2df <- data.frame(track) 
write.csv(dada2df, "./dada2df.csv") 

# Assign taxonomy
fastaRef <- "./silva_nr99_v138.1_train_set.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

## to install DECIPHER and phyloseq
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")

#installing packages
install.packages("phangorn")
library(phangorn)
install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library("gridExtra")
library(DECIPHER)
library(phyloseq)

#Construct Phylogenetic Tree
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs #This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phangAlign)  
treeNJ <- NJ(dm) #Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma = TRUE, rearrangment = "stochastic", control= pml.control(trace=0))
detach("package:phangorn",  unload=TRUE)
fitGTR

#Combine data into a  phyloseq object
##duplicates were manually removed in excel because the code was not working
###honey$sample.id <- paste0(gsub("", "", honey$sample.id), "D", honey$Gut_portion)
###honey <- honey[!duplicated(honey$sample.id),] # Remove dupicate entries for reverse reads
honey <- read.csv(file = "Metadata.csv")
rownames(honey) <- honey$New_ID
all(rownames(seqtabNoC) %in% honey$New_ID) # TRUE
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(honey), 
               tax_table(taxTab),
               phy_tree(fitGTR$tree))
ps

#Make short names for ASVs while storing sequence in a new column
tax <- data.frame(tax_table(ps))
tax$Sequence = rownames(tax) # Make a Sequence column
tax_table(ps) <- as.matrix(tax) # Replace with new taxonomy table
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # Make short names for ASVs

#Remove samples with few reads
toRemove <- c("Pomp1", "Pomp2", "Sphec3", "Crab3", "Scol1", "Sphec6", "Sphec5") 
ps.remove <- prune_samples(!sample_names(ps) %in% toRemove, ps) # New ps object

#Using Decontam Package
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")

##Inspect library size
df <- as.data.frame(sample_data(ps.remove)) #Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.remove)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

##Identify contaminants using prevalence method
sample_data(ps.remove)$is.neg <- sample_data(ps.remove)$Sample_or_Control == "Control"
###Using the conservative threshold of 0.5
contamdf.prev05 <- isContaminant(ps,remove, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
###Number of times taxa were observed in neg controls and samples
#Make phyloseq object of presence-absence in neg controls and samples
ps.pa <- transform_sample_counts(ps.remove, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
#Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Samples)")
write.csv(contamdf.prev05, "contamdf.csv") # True/false contaminant classification data frame

#Remove chloroplast, mitochondria, and unassigned phyla
ps1 <- subset_taxa(ps.remove, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
filterorder = c("Chloroplast")
filterfamily = c("Mitochondria")
ps1 = subset_taxa(ps1, !Order %in% filterorder) 
ps1 = subset_taxa(ps1, !Family %in% filterfamily) # New ps object

#Remove ASVs classified as contaminants
contamTaxa = c("ASV9", "ASV17", "ASV21", "ASV39", "ASV55", "ASV104", "ASV125", 
               "ASV127", "ASV272", "ASV408", "ASV590", "ASV1324", "ASV1532")
allTaxa = taxa_names(ps1)
allTaxa <- allTaxa[!(allTaxa %in% contamTaxa)] # Remove contamTaxa from allTaxa
ps.clean = prune_taxa(allTaxa, ps1) # New ps object

#Make ASV, taxonomy, and read table of contaminants
ps.contam = prune_taxa(contamTaxa, ps1)
taxtab.contam <- data.frame(tax_table(ps.contam))
seqtab.contam <- data.frame(t(otu_table(ps.contam))) 
merge.contam <- merge(seqtab.contam, taxtab.contam, by="row.names", all = FALSE)
write.csv(merge.contam, "./merge.contam.csv")

#Remove blank samples
blank <- c("Blank", "qblnk1", "qblnk2") 
ps.clean <- prune_samples(!sample_names(ps.clean) %in% blank, ps.clean)

#Presence/Absence Filtering
BiocManager::install("genefilter")
library(genefilter)
library(vegan)
microbe <- filter_taxa(ps.clean, filterfun(kOverA(2, 2)), TRUE) # New ps object

ps.clean.sum <- sample_sums(ps.clean)
microbe.sum <- sample_sums(microbe)
filter1 <- microbe.sum/ps.clean.sum
filterdf <- data.frame(filter1)

#Merge taxonomy and ASV table
taxtab.microbe <- data.frame(tax_table(microbe))
seqtab.microbe <- data.frame(t(otu_table(microbe))) 
merge.microbe <- merge(seqtab.microbe, taxtab.microbe, by="row.names", all = FALSE)
write.csv(merge.microbe, "./merge.microbe.csv")



############################# Old Code ######################################
##New
###k=2, A=2
microbe1 <- filter_taxa(ps2, filterfun(kOverA(2, 2)), TRUE)
pssum<- sample_sums(ps2)
microbesum1 <-sample_sums(microbe1)
filter1 <-microbesum1/pssum
filter1.df <- as.data.frame(filter1) #14 less than 0.9 with Silva
View(filter1.df)
write.csv(filter1, "/Users/lexiemartin/Desktop/prevfiilter1.csv", row.names = FALSE)
rarecurve(otu_table(microbe1), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(ps), step=50, cex=0.5, label=FALSE)
##k=3, A = 2
microbe2 <- filter_taxa(ps2, filterfun(kOverA(3,2)), TRUE)
pssum<- sample_sums(ps2)
microbesum2 <-sample_sums(microbe2)
filter2 <-microbesum2/pssum
filter2.df <- as.data.frame(filter2) #17 less than 0.9 with Silva
View(filter2.df)
write.csv(filter2, "/Users/lexiemartin/Desktop/prevfiilter2.csv", row.names = FALSE)
rarecurve(otu_table(microbe2), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)
##k=2, A = 3
microbe3 <- filter_taxa(ps2, filterfun(kOverA(2, 3)), TRUE)
pssum<- sample_sums(ps2)
microbesum3 <-sample_sums(microbe3)
filter3 <-microbesum3/pssum
filter3.df <- as.data.frame(filter3) #14 less than 0.9
View(filter3.df)
write.csv(filter3, "/Users/lexiemartin/Desktop/prevfiilter3.csv", row.names = FALSE)
filter3/filter1
rarecurve(otu_table(microbe3), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)


############################## DON'T USE ###################################
#library(genefilter)
##k=2, A=0.001
#microbe1 <- filter_taxa(ps1, filterfun(kOverA(2, 0.001)), TRUE)
#pssum<- sample_sums(ps1)
#microbesum1 <-sample_sums(microbe1)
#filter1 <-microbesum1/pssum
#view(filter1)
#filter1
#write.csv(filter1, "/Users/lexiemartin/Desktop/prevfiilter1.csv", row.names = FALSE)
#rarecurve(otu_table(microbe1), step=50, cex=0.5, label=FALSE)
#rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)
##k=5, A = 0.01
#microbe2 <- filter_taxa(ps1, filterfun(kOverA(5, 0.01)), TRUE)
#pssum<- sample_sums(ps1)
#microbesum2 <-sample_sums(microbe2)
#filter2 <-microbesum2/pssum
#view(filter2)
#write.csv(filter2, "/Users/lexiemartin/Desktop/prevfiilter2.csv", row.names = FALSE)
#rarecurve(otu_table(microbe2), step=50, cex=0.5, label=FALSE)
#rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)
##k=2, A = 0.01
#microbe3 <- filter_taxa(ps1, filterfun(kOverA(2, 0.01)), TRUE)
#pssum<- sample_sums(ps1)
#microbesum3 <-sample_sums(microbe3)
#filter3 <-microbesum3/pssum
#view(filter3)
#write.csv(filter3, "/Users/lexiemartin/Desktop/prevfiilter3.csv", row.names = FALSE)
#filter3/filter1
#rarecurve(otu_table(microbe3), step=50, cex=0.5, label=FALSE)
#rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)
##k=3, A = 0.001
#microbe4 <- filter_taxa(ps1, filterfun(kOverA(3, 0.001)), TRUE)
#pssum<- sample_sums(ps1)
#microbesum4 <-sample_sums(microbe4)
#filter4 <-microbesum4/pssum
#view(filter4)
#write.csv(filter4, "/Users/lexiemartin/Desktop/prevfiilter4.csv", row.names = FALSE)
#rarecurve(otu_table(microbe4), step=50, cex=0.5, label=FALSE)
#rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)


#Rarecurves
rarecurve(otu_table(microbe3), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(microbe3), step=50, cex=0.5, label=FALSE,xlim=c(0, 10000))
rarecurve(otu_table(microbe2), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(microbe2), step=50, cex=0.5, label=FALSE,xlim=c(0, 20000))
rarecurve(otu_table(microbe4), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(microbe4), step=50, cex=0.5, label=FALSE,xlim=c(0, 20000))
rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE)
rarecurve(otu_table(ps1), step=50, cex=0.5, label=FALSE,xlim=c(0, 20000))

#Rarefy the data without replacement and create abundance plots (set.seed(1) used)
library("phyloseq")
library("ggplot2")
library("vegan")
install.packages("cowplot")
library(cowplot)
theme_set(theme_cowplot())
##prevfilter1 -USED THIS ONE
ps.rarefied1 = rarefy_even_depth(microbe1, rngseed=1, sample.size=0.9*min(sample_sums(microbe1)), replace=F) #From the tutorial: "the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset"
sampledepth.ps1 <- sample_sums(ps2)
sampledepth.microbe <- sample_sums(microbe1)
sampledepth.rare <- sample_sums(ps.rarefied1)
filter.reads <- data.frame(sampledepth.ps1, sampledepth.microbe, sampledepth.rare) # Convert values to dataframe

install.packages("dplyr")
library("dplyr")
#Removing unassigned reads
seqtabAll.df <- as.data.frame(t(ps.rarefied1@otu_table))
taxTab.df <- as.data.frame((ps.rarefied1@tax_table))
mergedf = merge(seqtabAll.df, taxTab.df, 
                by="row.names", all=FALSE)
NAGenusdf = mergedf %>% filter(is.na(Genus))
write.csv(NAGenusdf, "/Users/lexiemartin/Desktop/NAGenusdf.csv", row.names = FALSE)
write.csv(ps.rarefied1@tax_table, "/Users/lexiemartin/Desktop/rareTaxTab.csv", row.names = TRUE)
write.csv(mergedf, "/Users/lexiemartin/Desktop/mergedf.csv", row.names = FALSE)


##Replacing the tax table
rownames(NewTaxTab1) <- NewTaxTab1$X
NewTaxTab1$X <- NULL
tax_table(ps.rarefied1) <- as.matrix(NewTaxTab1)



microbe1_relabundance <- transform_sample_counts(ps.rarefied1, function(OTU) OTU/sum(OTU))
# Top 20
top25 <- names(sort(taxa_sums(ps.rarefied1), decreasing=TRUE))[1:25]
ps.top25 <- transform_sample_counts(ps.rarefied1, function(OTU) OTU/sum(OTU))
ps.top25 <- prune_taxa(top25, ps.top25)
ps.top25.rel <- microbiome::transform(ps.top25, "compositional")
mypalette=c(brewer.pal(name="Spectral", n=11), brewer.pal(name="PuOr", n=9))
mypalette = c("seagreen", "maroon4", "lightskyblue", "aquamarine4", "darkolivegreen2", "slateblue4", "pink1", "royalblue3", "springgreen1", "orchid1", "turquoise4", "aquamarine1", "chocolate1", "coral4", "antiquewhite", "dodgerblue1", "gold")
mypalette1 = c("lightcyan1", "turquoise4", "skyblue1", "orchid1", "tomato", "darkolivegreen","darkolivegreen1", "coral4", "gold", "maroon", "darkorange1", "palegreen", "lightpink", "green4", "mediumpurple1", "slateblue4", "bisque3")
mypalette2 = c("#001219", "#005F73","#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226")
mypalette3 = c("#7D5BA6", "#8D6A9F", "#005F73","#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03", "#BC4749")
mypalette4 = c("#f02a81","#fdc082", "#beaed5", "#fffc94", "#386cb2","#1a9e75", "#bf5b0e", "#7fc97b", "#666666") #Family
mypalette5 = c("#fffc94","#fdc082", "#beaed5", "#756fb5", "#386cb2","#a6cee4","#bf5b0e", "#ffc2f8", "#f02a81","#FF9233","#66a600", "#7fc97b","#1a9e75", "#e6ab00", "#a67612", "#333333", "#666666") #genus
mypalette6 = c("#beaed5", "#1a9e75", "#7fc97b") #Phylum

#USE THIS TO MAKE ABUNDANCE PLOTS
plot_bar(ps.top25.rel, x="New_ID", fill="Genus") + xlab("Sample Name") +ylab("Relative Abundance")  + geom_bar(stat="identity") +scale_y_continuous(expand = c(0, 0))+  theme(panel.grid = element_blank(),
        panel.border = element_blank(), axis.text.x = element_text(size = 8,angle = 75, hjust = 1)) #+ scale_fill_manual(values= mypalette, na.value ="black")
library(RColorBrewer)
n <- 17
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

mypalette2 <- c("#d73027", "#a50026", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1","#4575b4", "#313695", "#8e0152", "#de77ae", "#fde0ef", "#e6f5d0", "#7fbc41", "#276419")

#+facet_wrap(~Nest_name, scales="free_x", nrow = 1 ) to facet wrap with everything in 1 row

plot_bar(microbe1_relabundance, x="New_ID", fill="Genus") + theme(legend.position = "none") + xlab("Sample Name") +ylab("Relative Abundance") + geom_bar(stat="identity") +scale_y_continuous(expand = c(0, 0))+  theme(panel.grid = element_blank(),panel.border = element_blank()) #+scale_fill_manual(values = c("grey", "grey1", "grey2", "grey3","grey4","grey5","grey6","grey7","grey8","grey9","grey10","grey11","grey12","grey13","grey14","indianred3","grey15",  "grey16","grey17","grey18","grey19","grey20","grey21","grey22","grey23","grey24","grey25","grey26","grey27","grey28","grey29","grey30","grey31","grey32","grey33","grey34","grey35","grey36","grey37","grey38","grey39","grey40","grey41","grey42","grey43","grey44","lightseagreen", "grey45", "grey46","grey47","grey48","grey49","grey50","grey51","grey52","grey53","grey54","grey55","grey56","grey57","grey58","darkmagenta","grey59", "grey60","grey61","grey62","grey63","grey64","grey65","grey66","grey67","grey68","grey69","grey70","grey71","grey72","grey73","grey74","grey75","grey76","grey77","grey78","grey79","grey80","grey81","grey82","grey83","grey84","grey85","grey86","grey87","grey88","grey89","grey90","grey91","grey92", "grey93"), na.value ="mediumseagreen")#to get the plot
plot_bar(microbe1_relabundance, x="New_ID", fill="Genus") + theme(legend.title = element_text(size = 10),legend.text = element_text(size = 6)) + geom_bar(stat="identity") +scale_y_continuous(expand = c(0, 0))+  theme(panel.grid = element_blank(),panel.border = element_blank())#+scale_fill_manual(values = c("grey", "grey1", "grey2", "grey3","grey4","grey5","grey6","grey7","grey8","grey9","grey10","grey11","grey12","grey13","grey14","indianred3","grey15",  "grey16","grey17","grey18","grey19","grey20","grey21","grey22","grey23","grey24","grey25","grey26","grey27","grey28","grey29","grey30","grey31","grey32","grey33","grey34","grey35","grey36","grey37","grey38","grey39","grey40","grey41","grey42","grey43","grey44","lightseagreen", "grey45", "grey46","grey47","grey48","grey49","grey50","grey51","grey52","grey53","grey54","grey55","grey56","grey57","grey58","darkmagenta","grey59", "grey60","grey61","grey62","grey63","grey64","grey65","grey66","grey67","grey68","grey69","grey70","grey71","grey72","grey73","grey74","grey75","grey76","grey77","grey78","grey79","grey80","grey81","grey82","grey83","grey84","grey85","grey86","grey87","grey88","grey89","grey90","grey91","grey92", "grey93"), na.value ="mediumseagreen")#to get the plot
#plot_bar(microbe1_relabundance, x="New_ID", fill="Genus") + scale_fill_grey()#to get the legend

#c("grey", "grey1", "grey2", "grey3","grey4","grey5","grey6","grey7","grey8","grey9","grey10","grey11","grey12","grey13","grey14","indianred3","grey15",  "grey16","grey17","grey18","grey19","grey20","grey21","grey22","grey23","grey24","grey25","grey26","grey27","grey28","grey29","grey30","grey31","grey32","grey33","grey34","grey35","grey36","grey37","grey38","grey39","grey40","grey41","grey42","grey43","grey44","lightseagreen", "grey45", "grey46","grey47","grey48","grey49","grey50","grey51","grey52","grey53","grey54","grey55","grey56","grey57","grey58","darkmagenta","grey59", "grey60","grey61","grey62","grey63","grey64","grey65","grey66","grey67","grey68","grey69","grey70","grey71","grey72","grey73","grey74","grey75","grey76","grey77","grey78","grey79","grey80","grey81","grey82","grey83","grey84","grey85","grey86","grey87","grey88","grey89","grey90","grey91","grey92","grey93","grey94","grey95","black"))#to get the plot


#prevfilter2 - colors are different
ps.rarefied2 = rarefy_even_depth(microbe2, rngseed=1, sample.size=0.9*min(sample_sums(microbe2)), replace=F) #From the tutorial: "the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset"
microbe2_relabundance <- transform_sample_counts(ps.rarefied2, function(OTU) OTU/sum(OTU))
##prevfilter2 plot
top20 <- names(sort(taxa_sums(ps.rarefied2), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps.rarefied2, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
ps.top20.rel <- microbiome::transform(ps.top20, "compositional")
plot_bar(ps.top20.rel, x="New_ID", fill="Genus")

#plot_bar(microbe2_relabundance, x="New_ID", fill="Genus") + theme(legend.position = "none") + xlab("Sample Name") +ylab("Relative Abundance") #to get the plot
#plot_bar(microbe2_relabundance, x="New_ID", fill="Genus") #to get the legend

##prevfilter4
ps.rarefied4 = rarefy_even_depth(microbe4, rngseed=1, sample.size=0.9*min(sample_sums(microbe4)), replace=F) #From the tutorial: "the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset"
microbe4_relabundance <- transform_sample_counts(ps.rarefied4, function(OTU) OTU/sum(OTU))
plot_bar(microbe4_relabundance, x="New_ID", fill="Genus") + theme(legend.position = "none") + xlab("Sample Name") +ylab("Relative Abundance") 


palette = c("grey", "grey1", "grey2", "grey3","grey4","grey5","grey6","grey7","grey8","grey9","grey10","grey11","grey12","grey13","grey14","grey15", "indianred3", "grey16","grey17","grey18","grey19","grey20","grey21","grey22","grey23","grey24","grey25","grey26","grey27","grey28","grey29","grey30","grey31","grey32","grey33","grey34","grey35","grey36","grey37","grey38","grey39","grey40","grey41","grey42","grey43","grey44","grey45","lightseagreen", "grey46","grey47","grey48","grey49","grey50","grey51","grey52","grey53","grey54","grey55","grey56","grey57","grey58","grey59","darkmagenta", "grey60","grey61","grey62","grey63","grey64","grey65","grey66","grey67","grey68","grey69","grey70","grey71","grey72","grey73","grey74","grey75","grey76","grey77","grey78","grey79","grey80","grey81","grey82","grey83","grey84","grey85","grey86","grey87","grey88","grey89","grey90","grey91","grey92","grey93","grey94","grey95","grey0")

palette = c("grey")

##ps1
ps.rarefied = rarefy_even_depth(ps1, rngseed=1, sample.size=0.9*min(sample_sums(ps1)), replace=F) #From the tutorial: "the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset"
ps_relabundance <- transform_sample_counts(ps.rarefied, function(OTU) OTU/sum(OTU))
plot_bar(ps_relabundance, x="New_ID", fill="Genus") + theme(legend.position = "none") + xlab("Sample Name") +ylab("Relative Abundance") #to get the plot
plot_bar(ps_relabundance, x="New_ID", fill="Genus") #to get the legend

#ps
og.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F) #From the tutorial: "the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset"
og_relabundance <- transform_sample_counts(og.rarefied, function(OTU) OTU/sum(OTU))
plot_bar(og_relabundance, x="New_ID", fill="Genus") + theme(legend.position = "none") + xlab("Sample Name") +ylab("Relative Abundance") #to get the plot
plot_bar(og_relabundance, x="New_ID", fill="Genus") #to get the legend

#NMDS Plot
ps.prop1 <- transform_sample_counts(ps.rarefied1, function(otu) otu/sum(otu))
ord.nmds.bray1 <- ordinate(ps.prop1, method="NMDS", distance="bray")
plot_ordination(ps.prop1, ord.nmds.bray1, color="Genus",title="NMDS") +geom_point(size=5) +theme(axis.title = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15)) #+scale_color_manual(values =c("turquoise4", "skyblue", "indianred3","darkorchid", "royalblue", "palegreen", "cornflowerblue", "slateblue1", "darkolivegreen", "mediumorchid","springgreen3", "olivedrab", "forestgreen", "cyan3"))

#NEW NMDS Plot
# Used bquote() to italicize
tiff("NMDSFull.tif", units="in", width=12, height=10, res=300)
ord.nmds.bray1 <- ordinate(ps.rarefied1, method="NMDS", distance="bray")
plot_ordination(ps.rarefied1, ord.nmds.bray1, color="Genus", shape="Wasp.Type") +
  ggtitle("Bray-Curtis") +
  theme_bw() +
  geom_point(size=3, stroke=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(values=c(19, 3), name="Type") +scale_color_manual(values = mypalette7, labels = c(bquote(italic("Anoplius")), bquote(italic("Apis (Bee)")),bquote(bolditalic("Brachygastra")), bquote(italic("Cercis")), bquote(italic("Chalybion")), bquote(italic("Dasymutilla")), bquote(italic("Euodynerus")), bquote(italic("Isodontia")), bquote(italic("Myzinium")),  bquote(italic("Parancistrocerus")),  bquote(italic("Polistes")),  bquote(italic("Pseudomethoca")),  bquote(italic("Sceliphron")),  bquote(italic("Trielis")),  bquote(italic("Trogomorpha")),  bquote(italic("Trypoxylon"))), name="Insect Genus")+
  theme(axis.title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) 
dev.off()

mypalette7 = c("#666666","#fdc082", "#beaed5", "#756fb5", "#386cb2","#a6cee4","#bf5b0e", "#ffc2f8", "#f02a81","#FF9233","#337300", "#7fc97b","#1a9e75", "#e6ab00", "#a67612", "#333333") #genus

#NEW NMDS - Honey Wasp Only
ord.nmds.bray1 <- ordinate(honeywasponlyps, method="NMDS", distance="bray")
plot_ordination(honeywasponlyps, ord.nmds.bray1, color="Nest_name") +
  ggtitle("Bray-Curtis") +
  theme_bw() +
  geom_point(size=3, stroke=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = mypalette7, name = "Nest Name") +
  theme(axis.title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) 

#Heat Map
tiff("top50heatmap2.tif", units="in", width=10, height=10, res=300)
heat.genus <-tax_glom(ps.top50, "Genus")
p <-plot_heatmap(heat.genus,
                 sample.order="New_ID", sample.label="New_ID", taxa.label = "Genus", low="#66CCFF", high="#FF3300", na.value="black")
p$scales$scales[[1]]$name <- "Sample ID"
p$scales$scales[[2]]$name <- "Bacterial Genus"
p
dev.off()



#Creating a phyloseq object with only the honeywasps 
SamplesToRemove <- c("Mel1", "Mel2", "Vesp1", "Vesp2", "Crab1", "Thynn1", "Sphec2", "Mut1", "Vesp3", "Sphec3", "Ich1", "Vesp17", "Pomp1", "Vesp4", "Vesp5", "Vesp6", "Vesp7", "Crab2", "Vesp8", "Vesp10", "Vesp11", "Sphec4", "Vesp12", "Vesp13", "Crab3", "Vesp15", "Scol1", "Vesp16", "Thynn2", "Crab5", "Crab7", "Mut2", "Pomp2", "Sphec5", "Sphec6")
honeywasponlyps <- subset_samples(ps.rarefied1, !(New_ID %in% SamplesToRemove))

top25 <- names(sort(taxa_sums(honeywasponlyps), decreasing=TRUE))[1:25]
ps.top25 <- transform_sample_counts(honeywasponlyps, function(OTU) OTU/sum(OTU))
ps.top25 <- prune_taxa(top25, ps.top25)
ps.top25.rel <- microbiome::transform(ps.top25, "compositional")
plot_bar(ps.top25.rel, x="New_ID", fill="Phylum") + xlab("Sample Name") +ylab("Relative Abundance")  + geom_bar(stat="identity") +facet_wrap(~Wasp.Type, scales="free_x") +scale_y_continuous(expand = c(0, 0))+  theme(panel.grid = element_blank(),
                                                                                                                                                                                                                        panel.border = element_blank(), axis.text.x = element_text(size = 8,angle = 75, hjust = 1)) + scale_fill_manual(values= col_vector, na.value ="black")
plot_bar(ps.top25.rel, x="New_ID", fill="Genus") + xlab("Sample Name") +ylab("Relative Abundance")  + geom_bar(stat="identity") +facet_wrap(~Nest_name, scales="free_x", nrow = 1) +scale_y_continuous(expand = c(0, 0))+  theme(panel.grid = element_blank(),
      panel.border = element_blank(), axis.text.x = element_text(size = 8,angle = 75, hjust = 1)) + scale_fill_manual(values= col_vector, na.value ="black")










##Everything below this line is EXTRA
#EXTRA
ps.top20 <- transform_sample_counts(ps.rarefied, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="New_ID", fill="Genus") + theme(legend.position = "none") + xlab("Sample Name") +ylab("Relative Abundance")
plot_bar(ps.top20, x="New_ID", fill="Genus")


rare1 = rarefy_even_depth(microbe1)
rareog= rarefy_even_depth(ps2)

rare1sum = sample_sums(rare1)
rareotu=otu_table(rare1)
view(rareotu)
rare1divide = rareotu/rare1sum
view(rare1divide)
view(taxtab(ps1))
view(taxtab(ps))

sample_names(ps1)
rank_names(ps1)
sample_variables(ps1)
total = median(sample_sums(ps1))
standf = function(x, t=total) round(t * (x/sum(x)))
ps2 = transform_sample_counts(ps1, standf)
ps2

ps2_abund <- filter_taxa(ps2, function(x) sum(x > total*0.05) > 0, TRUE)
ps2_abund



##Rarefy
library("phyloseq")
library("ggplot2")
library("vegan")
#Rarefaction Curve
rarecurve(otu_table(ps), step=50, cex=0.5, label=FALSE)
#Rarefy the data without replacement
ps.rarefied = rarefy_even_depth(ps2_abund, rngseed=1, sample.size=0.9*min(sample_sums(ps2_abund)), replace=F) #From the tutorial: "the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset"
plot_bar(ps.rarefied, fill="Phylum")

#EXTRA
## Prevalence Filtering
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
##There are additional steps where you have to choose the filtering criteria, but I didn't understand how you are supposed to choose, so I didn't complete this filtering step


