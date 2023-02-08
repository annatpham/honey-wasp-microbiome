# After Processing Sequences via DADA2 ------------------------------------
# Remove samples with significant read cuts after the DADA2 pipeline
# Remove samples with >5000 unassigned ASVs
toRemove <- c("JCH644", "JCH713", "JCH750", "JCH623", "JCH625", "JCH642",
              "JCH648", "JCH653", "JCH655", "JCH664", "JCH769","JCH817",
              "Crab6", "Crab9", "Crab8", "Mut3", "Pomp3", "Pomp4", "Pomp5",
              "Crab4", "Vesp14", "Vesp9", "Sphec1", "Mut4") 
ps.remove <- prune_samples(!sample_names(ps) %in% toRemove, ps) # New phyloseq object

# Use Decontam Package ----------------------------------------------------
# Source: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")

# Inspect library size
df <- as.data.frame(sample_data(ps.remove)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.remove)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

# Identify contaminants using the prevalence method
sample_data(ps.remove)$is.neg <- sample_data(ps.remove)$Sample_or_Control == "Control"
## Use the conservative threshold of 0.5
contamdf.prev05 <- isContaminant(ps,remove, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
## Download the true/false contaminant classification data
write.csv(contamdf.prev05, "contamdf.csv")

# Number of times taxa were observed in neg controls and samples
## Make phyloseq object of presence-absence in neg controls and samples
ps.pa <- transform_sample_counts(ps.remove, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
## Make dataframe of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Samples)")

# Remove Chloroplast, Mitochondria, and Unassigned Phyla ------------------
ps1 <- subset_taxa(ps.remove, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
filterorder = c("Chloroplast")
filterfamily = c("Mitochondria")
ps1 = subset_taxa(ps1, !Order %in% filterorder) 
ps1 = subset_taxa(ps1, !Family %in% filterfamily) # New phyloseq object

# Remove ASVs Classified as Contaminants by Decontam ----------------------
contamTaxa = c("ASV9", "ASV17", "ASV21", "ASV39", "ASV55", "ASV104", "ASV125", 
               "ASV127", "ASV272", "ASV408", "ASV590", "ASV1324", "ASV1532")
allTaxa = taxa_names(ps1)
allTaxa <- allTaxa[!(allTaxa %in% contamTaxa)] # Remove contamTaxa from allTaxa
ps.clean = prune_taxa(allTaxa, ps1) # New phyloseq object
View(ps.clean@tax_table)

# Make ASV, taxonomy, and read table of contaminants
ps.contam = prune_taxa(contamTaxa, ps1)
taxtab.contam <- data.frame(tax_table(ps.contam))
seqtab.contam <- data.frame(t(otu_table(ps.contam))) 
merge.contam <- merge(seqtab.contam, taxtab.contam, by="row.names", all = FALSE)
write.csv(merge.contam, "./merge.contam.csv")

# Remove Blank Samples from Dataset ---------------------------------------
blank <- c("Blank", "qblnk1", "qblnk2") 
ps.clean <- prune_samples(!sample_names(ps.clean) %in% blank, ps.clean)

# Presence/Absence Filtering 1 ---------------------------------------------
# Remove ASVs with fewer than 2 reads in at least 2 samples. kOverA(X samples, X reads)
BiocManager::install("genefilter")
library(genefilter)
library(vegan)
microbe <- filter_taxa(ps.clean, filterfun(kOverA(2, 2)), TRUE) # New phyloseq object

# Check number of reads
# Check lines 208-213 for longer comment
ps.clean.sum <- sample_sums(ps.clean)
microbe.sum <- sample_sums(microbe) # microbe.sum <- data.frame(sample_sums(microbe))
filter1 <- microbe.sum/ps.clean.sum
filterdf <- data.frame(filter1)

# Remove Samples with <1000 Reads after P/A Filtering ---------------------
remove <- c("Pomp1", "Pomp2", "Sphec3", "Crab3", "Scol1", "Sphec6", "Sphec5") 
ps.clean <- prune_samples(!sample_names(ps.clean) %in% remove, ps.clean) 

# Presence-Absence Filtering 2 --------------------------------------------
# Redo p/a filtering
# This time using the final dataset, which does not include the samples from line 82
microbe <- filter_taxa(ps.clean, filterfun(kOverA(2, 2)), TRUE)

# Check number of reads
# Check lines 208-213 for longer comment
ps.clean.sum <- sample_sums(ps.clean)
microbe.sum <- sample_sums(microbe) # microbe.sum <- data.frame(sample_sums(microbe))
filter1 <- microbe.sum/ps.clean.sum
filterdf <- data.frame(filter1)

# Verify Taxonomy Assignments via BLAST -----------------------------------
#Merge taxonomy and ASV table
taxtab.microbe <- data.frame(tax_table(microbe))
seqtab.microbe <- data.frame(t(otu_table(microbe))) 
merge.microbe1 <- merge(seqtab.microbe, taxtab.microbe, by="row.names", all = FALSE)
write.csv(merge.microbe1, "./merge.microbe1.csv") # Download and BLAST sequences with >100 reads 

# Remove ASVs identified as non-bacteria via BLAST
badTaxa = c("ASV1077", "ASV1138", "ASV1290", "ASV1498", "ASV1799", "ASV233", 
            "ASV623", "ASV637", "ASV665", "ASV688", "ASV810", "ASV878", "ASV961", 
            "ASV396","ASV643","ASV879", "ASV966") 
          # Vector of ASVs identified as non-bacteria
allTaxa = taxa_names(microbe)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)] # Remove badTaxa from allTaxa
ps.taxa = prune_taxa(allTaxa, microbe) # New phyloseq object
View(ps.taxa@tax_table)

# Reassign taxon name for sequences >97% similar to a genus 
# Elevate to the family level sequences with high similarity to multiple genera 
tax <- as.data.frame(tax_table(ps.taxa)) 
tax$ASV <- rownames(tax)
tax$Genus[tax$ASV == "ASV1"] <- "Bifidobacterium"
tax$Genus[tax$ASV == "ASV10"] <- "Gluconobacter"
tax$Genus[tax$ASV == "ASV102"] <- "Pantoea"
tax$Genus[tax$ASV == "ASV108"] <- "Lactobacillaceae"
tax$Genus[tax$ASV == "ASV1094"] <- "Rhizobiaceae"
tax$Genus[tax$ASV == "ASV11"] <- "Ca. Stammerula"
tax$Genus[tax$ASV == "ASV1207"] <- "Pirellulaceae"
tax$Genus[tax$ASV == "ASV131"] <- "Lactobacillaceae"
tax$Genus[tax$ASV == "ASV138"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV14"] <- "Ca. Schmidhempelia"
tax$Genus[tax$ASV == "ASV145"] <- "Uncultured bacterium"
tax$Family[tax$ASV == "ASV145"] <- "Uncultured bacterium"
tax$Genus[tax$ASV == "ASV15"] <- "Acinetobacter"
tax$Genus[tax$ASV == "ASV16"] <- "Convivina"
tax$Genus[tax$ASV == "ASV200"] <- "Methylorubrum"
tax$Genus[tax$ASV == "ASV210"] <- "Aureimonas"
tax$Genus[tax$ASV == "ASV235"] <- "Convivina"
tax$Genus[tax$ASV == "ASV239"] <- "Gluconobacter"
tax$Genus[tax$ASV == "ASV245"] <- "Luteimonas"
tax$Genus[tax$ASV == "ASV246"] <- "Mammiliicoccus"
tax$Genus[tax$ASV == "ASV25"] <- "Convivina"
tax$Genus[tax$ASV == "ASV250"] <- "Yersiniaceae"
tax$Genus[tax$ASV == "ASV267"] <- "Acetobacteraceae"
tax$Genus[tax$ASV == "ASV27"] <- "Neisseriaceae"
tax$Genus[tax$ASV == "ASV28"] <- "Ca. Schmidhempelia"
tax$Genus[tax$ASV == "ASV281"] <- "Solirubrobacterales"
tax$Genus[tax$ASV == "ASV3"] <- "Bifidobacterium"
tax$Genus[tax$ASV == "ASV30"] <- "Lactobacillaceae"
tax$Genus[tax$ASV == "ASV308"] <- "Sphingomonas"
tax$Genus[tax$ASV == "ASV309"] <- "Microlunatus"
tax$Genus[tax$ASV == "ASV31"] <- "Acinetobacter"
tax$Genus[tax$ASV == "ASV34"] <- "Asaia"
tax$Genus[tax$ASV == "ASV36"] <- "Pectobacteriaceae"
tax$Genus[tax$ASV == "ASV38"] <- "Halomonadaceae"
tax$Genus[tax$ASV == "ASV383"] <- "Enterobacterales"
tax$Family[tax$ASV == "ASV383"] <- "Enterobacterales"
tax$Genus[tax$ASV == "ASV4"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV42"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV48"] <- "Atopobiaceae"
tax$Genus[tax$ASV == "ASV498"] <- "Comamonadaceae"
tax$Genus[tax$ASV == "ASV5"] <- "Bifidobacterium"
tax$Genus[tax$ASV == "ASV51"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV57"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV578"] <- "Sphingomonas"
tax$Genus[tax$ASV == "ASV6"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV60"] <- "Bifidobacterium"
tax$Genus[tax$ASV == "ASV627"] <- "Affinibrenneria"
tax$Genus[tax$ASV == "ASV63"] <- "Pediococcus"
tax$Order[tax$ASV == "ASV642"] <- "Hyphomicrobiales"
tax$Family[tax$ASV == "ASV642"] <- "Kichenibacteriaceae"
tax$Genus[tax$ASV == "ASV642"] <- "Lichenibacterium"
tax$Family[tax$ASV == "ASV663"] <- "Erythrobacteraceae"
tax$Genus[tax$ASV == "ASV663"] <- "Erythrobacter"
tax$Genus[tax$ASV == "ASV72"] <- "Enterobacteriaceae"
tax$Genus[tax$ASV == "ASV758"] <- "Nocardioidaceae"
tax$Genus[tax$ASV == "ASV78"] <- "Bifidobacteriaceae"
tax$Genus[tax$ASV == "ASV81"] <- "Agrobacterium"
tax$Genus[tax$ASV == "ASV829"] <- "Bifidobacteriaceae"
tax$Genus[tax$ASV == "ASV84"] <- "Rhizobiaceae"
tax$Genus[tax$ASV == "ASV86"] <- "Lactobacillaceae"
tax$Genus[tax$ASV == "ASV87"] <- "Sinorhizobium"
tax$Genus[tax$ASV == "ASV91"] <- "Enterobacteriaceae"
tax_table(ps.taxa) <- as.matrix(tax) # Replace with new taxonomy table
View(ps.taxa@tax_table)

# Rarefaction -------------------------------------------------------------
# Rarefy samples to 90% of the minimum sample depth
## Full Rare Curve
png("fullrare.png", units="in", width=4.5, height=4, res=300)
rarecurve(otu_table(ps.taxa), step=100, cex=0.5, label=TRUE, xlim=c(0, 100000)) 
dev.off()

## Partial Rare Curve
## X-axis max of 2000
png("partialrare.png", units="in", width=4.5, height=4, res=300)
rarecurve(otu_table(ps.taxa), step=100, cex=0.5, label=FALSE, xlim=c(0, 2000)) 
dev.off()

## Rarefied phyloseq object
ps.rarefied = rarefy_even_depth(ps.taxa, rngseed=1, sample.size=0.9*min(sample_sums(ps.taxa)), replace=F)

# Check Reads and Make Final Merged Table ---------------------------------
# Check number of reads kept after each filtering step
depth.ps1 <- data.frame(sample_sums(ps1)) # After removing chloroplast, mitochondria, unassigned phyla
# ps.clean.sum: After removing contaminant taxa
depth.taxa <- sample_sums(ps.taxa) # After presence-absence filtering and removing nonbacterial ASVs ID'd via BLAST
depth.rare <- sample_sums(ps.rarefied) # After rarefying
filter.reads <- data.frame(depth.ps1, ps.clean.sum, depth.taxa, depth.rare) 
write.csv(filter.reads, "./filter.reads.csv") # Download

# If the vectors are different lengths (b/c of samples being removed between filtering steps), you won't be able to make a table
# Solution: Configure the short vectors to have same length as the longest vector 
# (Add NA's to bottom of the short columns so vectors are the same length)
# Will need to fix the downloaded table so that the NAs are in the correct place
## Ex. length(depth.clean) <- length(depth.ps1) # 2nd shortest vector <- long vector
## Ex. length(depth.microbe) <- length(depth.ps1) # shortest vector <- long vector

# On "filter.reads.csv" sheet, check percentage of reads kept after each filtering step
# Want to keep >90% reads between filtering steps
## filter1 <- ps.clean.sum/depth.clean
## filter2 <- depth.taxa/ps.clean.sum
## filterdf <- data.frame(filter1, filter2)
## write.csv(filterdf, "./filterdf.csv")

# Make the final taxonomy and ASV tables into data frames and merge them
taxtab.rare <- data.frame(tax_table(ps.rarefied))
seqtab.rare <- data.frame(t(otu_table(ps.rarefied)))
merge.rare <- merge(seqtab.rare, taxtab.rare, by="row.names", all=FALSE)
write.csv(merge.rare, "./merge.rare1.csv")
 