# Load packages
library(phyloseq)
library(ggplot2)
library(svglite)

# Reorganizing axis labels by nest latitude and ordering names ------------
raredf <-ps.rarefied@sam_data
raredf$NestOrder = factor(raredf$Nest_name, levels=c("AUS","SA1", "SA2", "SA3", 
                                                     "SA4", "SA5", "SA6", "SA7", 
                                                     "SA8","PLV", "KV", "Har1", 
                                                     "Har2", "Ind.", "Individual")) 
                                           # changing order
raredf$WaspOrder = factor(raredf$Wasp.Type, levels=c("Honey Wasp", "Other Wasp", "Honey Bee"))
raredf$New_ID = factor(raredf$New_ID, levels=c("BrachA5", "BrachA4", "BrachA3", 
                                               "BrachA2", "BrachA1", "Brach16", 
                                               "Brach15", "Brach14", "Brach13", 
                                               "Brach12", "Brach11", "Brach23", 
                                               "Brach22", "Brach21", "Brach36", 
                                               "Brach35", "Brach34", "Brach33", 
                                               "Brach32", "Brach31", "Brach46", 
                                               "Brach45", "Brach44", "Brach43", 
                                               "Brach42", "Brach41", "Brach53", 
                                               "Brach52", "Brach51", "Brach63", 
                                               "Brach62", "Brach61", "Brach73", 
                                               "Brach72", "Brach71", "Brach83", 
                                               "Brach82", "Brach81", "BrachP16", 
                                               "BrachP15", "BrachP14", "BrachP13", 
                                               "BrachP12", "BrachP11", "BrachK16", 
                                               "BrachK15", "BrachK14", "BrachK13", 
                                               "BrachK12", "BrachK11", "BrachH15", 
                                               "BrachH14", "BrachH13", "BrachH12", 
                                               "BrachH11", "BrachH25", "BrachH24", 
                                               "BrachH23", "BrachH22", "BrachH21",
                                               "BrachIn5", "BrachIn4", "BrachIn3", 
                                               "BrachIn2", "BrachIn1",
                                               "Mel2","Mel1", 
                                               "Ich1",
                                               "Sphec4","Sphec2",
                                               "Crab7", "Crab5", "Crab2", "Crab1",
                                               "Thynn2", "Thynn1",
                                               "Mut2", "Mut1",
                                               "Vesp16", "Vesp11","Vesp17","Vesp6",
                                               "Vesp13","Vesp7","Vesp10","Vesp8",
                                               "Vesp5","Vesp2","Vesp1","Vesp15", 
                                               "Vesp12","Vesp4","Vesp3"))
sample_data(ps.rarefied) <- raredf # Replace with new metadata
View(ps.rarefied@sam_data)

# Top 25 relative abundance -----------------------------------------------
## Determine top 25
top25 <- names(sort(taxa_sums(ps.rarefied), decreasing=TRUE))[1:25]
ps.top25 <- transform_sample_counts(ps.rarefied, function(OTU) OTU/sum(OTU))
ps.top25 <- prune_taxa(top25, ps.top25)

## Relative abundance
library(microbiome); packageVersion("microbiome")
ps.top25.rel <- microbiome::transform(ps.top25, "compositional")

## Agglomerate bacterial taxa to the genus level to remove black lines within bars
genusGlommed = tax_glom(ps.top25.rel, "Genus")

## Color palette
mypalette8 <- c("Acinetobacter"="#6F7302", "Bifidobacterium"="#7570B3", 
                "Bombilactobacillus"="#CBC3E3", "Ca. Schmidhempelia"= "#4F7A8E", 
                "Ca. Stammerula" = "#D76865", "Convivina"="#D1DBBD", 
                "Enterobacteriaceae" = "#e77732", "Exiguobacterium" = "#F0D787", 
                "Fructilactobacillus"="#013F5B", "Gilliamella"="#B00068", 
                "Gluconobacter"= "#A65C32", "Lactobacillus"="#3CB0CC", 
                "Lactococcus"="#e5ba2e", "Neisseriaceae"= "#9CB4BF", 
                "Pediococcus" ="#A62F03", "Wolbachia" = "#C0E9F1")

## Plot
full_genus2 <- plot_bar(genusGlommed, "New_ID", fill="Genus") + coord_flip()+
  scale_fill_manual(values= mypalette8, na.value ="black")+
  facet_grid(NestOrder~., scales="free_y", space = "free_y", switch="y")+ 
  # Splitting samples into different panels by nests
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank())+
  theme(legend.key.size = unit(0.5,"line")) +
  theme(legend.text = element_text(face="italic", size = 25)) +
  theme(text = element_text(size=30)) + 
  theme(legend.position="bottom")+
  theme(axis.title.y = element_text(size=30))+
  ylab("Relative Abundance") +
  xlab("") +
  scale_y_continuous(expand = c(0, 0)) # Remove space between axis and bars
ggsave(file="full_genus9.svg", plot=full_genus2, width=15, height=20)


# Top 50 relative abundance -----------------------------------------------
## Determine Top 50
top50 <- names(sort(taxa_sums(ps.rarefied), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps.rarefied, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
ps.top50.rel <- microbiome::transform(ps.top50, "compositional") # Relative abundance

## Agglomerate bacterial taxa to the genus level to remove black lines within bars
genusGlommed50 = tax_glom(ps.top50.rel, "Genus")

## Plot
full_genus50 <- plot_bar(genusGlommed50, "New_ID", fill="Genus") + coord_flip()+
  facet_grid(NestOrder~., scales="free_y", space = "free_y", switch="y")+ 
  # Splitting samples into different panels by nests
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank())+
  theme(legend.key.size = unit(0.5,"line")) +
  theme(legend.text = element_text(face="italic", size = 25)) +
  theme(text = element_text(size=30)) + 
  theme(legend.position="bottom")+
  theme(axis.title.y = element_text(size=30))+
  ylab("Relative Abundance") +
  xlab("") +
  scale_y_continuous(expand = c(0, 0)) # Remove space between axis and bars
ggsave(file="full_genus50.svg", plot=full_genus50, width=25, height=20)

# Heatmaps -----------------------------------------------------------------
## Add ASVs to genera (outputs ASV#:Genus)
tax.rare <- as.data.frame(tax_table(ps.rarefied))
tax.rare$ASVGenus <- paste(row.names(tax.rare), tax.rare$Genus, sep=": ") # Creating a new column called ASVGenus
tax_table(ps.rarefied) <- as.matrix(tax.rare)
View(ps.rarefied@tax_table) # Replace taxonomy table

## Plot - using the top 50 relative abundance phyloseq object
top50heat <-plot_heatmap(ps.top50.rel, sample.label = "New_ID", 
                         taxa.label = "ASVGenus", 
                         taxa.order = "ASVGenus", 
                         sample.order = rev(c("BrachA1", "BrachA2", "BrachA3", 
                                              "BrachA4", "BrachA5", "Brach11", 
                                              "Brach12", "Brach13", "Brach14", 
                                              "Brach15", "Brach16", "Brach21", 
                                              "Brach22", "Brach23", "Brach31", 
                                              "Brach32", "Brach33", "Brach34", 
                                              "Brach35", "Brach36", "Brach41", 
                                              "Brach42", "Brach43", "Brach44", 
                                              "Brach45", "Brach46", "Brach51", 
                                              "Brach52", "Brach53", "Brach61", 
                                              "Brach62", "Brach63", "Brach71", 
                                              "Brach72", "Brach73", "Brach81", 
                                              "Brach82", "Brach83", "BrachP11", 
                                              "BrachP12", "BrachP13", "BrachP14", 
                                              "BrachP15", "BrachP16", "BrachK11", 
                                              "BrachK12", "BrachK13", "BrachK14", 
                                              "BrachK15", "BrachK16", "BrachH11", 
                                              "BrachH12", "BrachH13", "BrachH14", 
                                              "BrachH15", "BrachH21", "BrachH22", 
                                              "BrachH23", "BrachH24", "BrachH25",
                                              "BrachIn1", "BrachIn2", "BrachIn3", 
                                              "BrachIn4", "BrachIn5", "Vesp3",
                                              "Vesp4", "Vesp12", "Vesp15", "Vesp1", 
                                              "Vesp2", "Vesp5", "Vesp8", "Vesp10", 
                                              "Vesp7", "Vesp13", "Vesp6", "Vesp17", 
                                              "Vesp11", "Vesp16", "Mut1", "Mut2",
                                              "Thynn1", "Thynn2",
                                              "Crab1", "Crab2", "Crab5", "Crab7",
                                              "Sphec2","Sphec4", "Ich1",  
                                              "Mel1", "Mel2")), 
                                              # Setting the order of the samples
                         na.value="black", high = "#D5FFFF", low = "#00314D") + 
                         # Setting the gradient colors
  theme_bw(base_size=30) +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95))+
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  coord_flip() # Flipping the axes so that samples are on the y-axis
ggsave(file="top50heat3.svg", plot=top50heat, width=20, height=20)

## Top 100 relative abundance
top100 <- names(sort(taxa_sums(ps.rarefied), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps.rarefied, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
ps.top100.rel <- microbiome::transform(ps.top100, "compositional") # Relative abundance

## Plot
top100heat <-plot_heatmap(ps.top100.rel, sample.label = "New_ID", 
                          taxa.label = "ASVGenus", 
                          taxa.order = "ASVGenus", 
                          sample.order = rev(c("BrachA1", "BrachA2", "BrachA3", 
                                               "BrachA4", "BrachA5", "Brach11", 
                                               "Brach12", "Brach13", "Brach14", 
                                               "Brach15", "Brach16", "Brach21", 
                                               "Brach22", "Brach23", "Brach31", 
                                               "Brach32", "Brach33", "Brach34", 
                                               "Brach35", "Brach36", "Brach41", 
                                               "Brach42", "Brach43", "Brach44", 
                                               "Brach45", "Brach46", "Brach51", 
                                               "Brach52", "Brach53", "Brach61", 
                                               "Brach62", "Brach63", "Brach71", 
                                               "Brach72", "Brach73", "Brach81", 
                                               "Brach82", "Brach83", "BrachP11", 
                                               "BrachP12", "BrachP13", "BrachP14", 
                                               "BrachP15", "BrachP16", "BrachK11", 
                                               "BrachK12", "BrachK13", "BrachK14", 
                                               "BrachK15", "BrachK16", "BrachH11", 
                                               "BrachH12", "BrachH13", "BrachH14", 
                                               "BrachH15", "BrachH21", "BrachH22", 
                                               "BrachH23", "BrachH24", "BrachH25",
                                               "BrachIn1", "BrachIn2", "BrachIn3", 
                                               "BrachIn4", "BrachIn5", "Vesp3",
                                               "Vesp4", "Vesp12", "Vesp15", "Vesp1", 
                                               "Vesp2", "Vesp5", "Vesp8", "Vesp10", 
                                               "Vesp7", "Vesp13", "Vesp6", "Vesp17", 
                                               "Vesp11", "Vesp16", "Mut1", "Mut2",
                                               "Thynn1", "Thynn2",
                                               "Crab1", "Crab2", "Crab5", "Crab7",
                                               "Sphec2","Sphec4", "Ich1",  
                                               "Mel1", "Mel2")), 
                          # Setting the order of the samples
                          na.value="black", high = "#D5FFFF", low = "#00314D") + 
  # Setting the gradient colors
  theme_bw(base_size=30) +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95))+
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  coord_flip() # Flipping the axes so that samples are on the y-axis
ggsave(file="top100heat.svg", plot=top100heat, width=30, height=30)

# Bray-Curtis NMDS by insect type -------------------------------------------
## Performing the ordination
bray <- ordinate(ps.rarefied, method="NMDS", distance="bray")

## Color palette
mypalette7 <- c( "Honey Wasp" = "#006480", "Other Wasp" ="#D76865", "Honey Bee" = "#F0d787")

## Plot
braywasp <-plot_ordination(ps.rarefied, bray, color="WaspOrder", shape="WaspOrder") +
  scale_color_manual(values = mypalette7, name="") +
  scale_shape_manual(values=c(17, 3, 4), name="") + # Specify shape by species type
  theme_bw() +
  geom_point(size=3, stroke=1) + 
  theme(text = element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file="braywasp7.svg", plot=braywasp, width=8, height=5)

# Bray-Curtis NMDS by nest type -------------------------------------------
## Subset to honey wasps only
SamplesToRemove <- c("Mel1", "Mel2", "Vesp1", "Vesp2", "Crab1", "Thynn1", 
                     "Sphec2", "Mut1", "Vesp3", "Sphec3", "Ich1", "Vesp17", 
                     "Pomp1", "Vesp4", "Vesp5", "Vesp6", "Vesp7", "Crab2", 
                     "Vesp8", "Vesp10", "Vesp11", "Sphec4", "Vesp12", "Vesp13", 
                     "Crab3", "Vesp15", "Scol1", "Vesp16", "Thynn2", "Crab5", 
                     "Crab7", "Mut2", "Pomp2", "Sphec5", "Sphec6") 
honeywasponlyps <- subset_samples(ps.rarefied, !(New_ID %in% SamplesToRemove))
# OR instead of lines 225-231, could do: 
# honeywasponlyps <- subset_samples(ps.rarefied, Genus=="Brachygastra") 

## Performing the ordination
bray1 <- ordinate(honeywasponlyps, method="NMDS", distance="bray")

## Color palette
mypalette6 <- c("SA1" = "#013F5B", "SA2"= "#1E556F", "SA3"="#4F7A8E", 
                "SA4" = "#9CB4BF", "SA5" ="#3CB0CC", "SA6" = "#90D8E5", 
                "SA7" =  "#A3DEE9", "SA8" = "#C0E9F1","AUS" = "#A62F03", 
                "PLV" = "#A65C32", "KV" = "#F2A03D", "Har1" ="#6F7302", 
                "Har2"= "#D1DBBD", "Individual" = "#40231C")

## Plot
braynest1<-plot_ordination(honeywasponlyps, bray1, color="Nest_name", shape="Wasp.Type") +
  scale_color_manual(values = mypalette6, name="Nest Name") +
  scale_shape_manual(values=17, labels = NULL, name = NULL, guide = 'none')+
  theme_bw() +
  geom_point(size=3, stroke=1) +  
  theme(text = element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(file="braynest5.svg", plot=braynest1, width=8, height=5)


# PERMANOVA with respect to insect type -------------------------------------
# Includes all samples

# PERMANOVA test
library(vegan)
## Calculate bray curtis distance matrix
hw_bray <- phyloseq::distance(ps.rarefied, method = "bray")
## Make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.rarefied))
## Adonis test 
adonis2(hw_bray ~ Wasp.Type, data = sampledf) # adonis has been deprecated, use adonis2

# Homogeneity of dispersion test
dispersion <- betadisper(hw_bray, sampledf$Wasp.Type)
permutest(dispersion) 
plot(dispersion)

# Pairwise comparisons
library(pairwiseAdonis)
pairwise <- pairwise.adonis(hw_bray, sampledf$Wasp.Type) # Return Bonferonni adjusted p-values
pairwise
write.csv(pairwise, "./pairwise_wasptype.csv")

# PERMANOVA with respect to nest ------------------------------------------
# Includes only honey wasps (both the nests & individuals)

# PERMANOVA test
## Calculate bray curtis distance matrix
hw_bray <- phyloseq::distance(honeywasponlyps, method = "bray") # Using the honeywasponlyps
## Make a data frame from the sample_data
sampledf <- data.frame(sample_data(honeywasponlyps))
## Adonis test 
adonis2(hw_bray ~ Nest_name, data = sampledf) # adonis has been deprecated, use adonis2

# Homogeneity of dispersion test
dispersion <- betadisper(hw_bray, sampledf$Nest_name)
permutest(dispersion) 
plot(dispersion)

# Pairwise comparisons
library(pairwiseAdonis)
pairwise <- pairwise.adonis(hw_bray, sampledf$Nest_name) # Return Bonferonni adjusted p-values
pairwise
write.csv(pairwise, "./pairwise_honeywasponly.csv")

# PERMANOVA with respect to nest ------------------------------------------
# Includes only honey wasps from nests (no individuals)

## Honey Wasps IN NESTS Only
SamplesToRemove <- c("Mel1", "Mel2", "Vesp1", "Vesp2", "Crab1", "Thynn1", 
                     "Sphec2", "Mut1", "Vesp3", "Sphec3", "Ich1", "Vesp17", 
                     "Pomp1", "Vesp4", "Vesp5", "Vesp6", "Vesp7", "Crab2", 
                     "Vesp8", "Vesp10", "Vesp11", "Sphec4", "Vesp12", "Vesp13", 
                     "Crab3", "Vesp15", "Scol1", "Vesp16", "Thynn2", "Crab5", 
                     "Crab7", "Mut2", "Pomp2", "Sphec5", "Sphec6",
                     "BrachIn1", "BrachIn2", "BrachIn3", "BrachIn4", "BrachIn5")
nestwasponlyps <- subset_samples(ps.rarefied, !(New_ID %in% SamplesToRemove))

# PERMANOVA test
## Calculate bray curtis distance matrix
hw_bray <- phyloseq::distance(nestwasponlyps, method = "bray")
## Make a data frame from the sample_data
sampledf <- data.frame(sample_data(nestwasponlyps))
## Adonis test 
adonis2(hw_bray ~ Nest_name, data = sampledf) # adonis has been deprecated, use adonis2

# Homogeneity of dispersion test
dispersion <- betadisper(hw_bray, sampledf$Nest_name)
permutest(dispersion) 
plot(dispersion)

# Pairwise comparisons
library(pairwiseAdonis)
pairwise <- pairwise.adonis(hw_bray, sampledf$Nest_name) # Return Bonferonni adjusted p-values
pairwise
write.csv(pairwise, "./pairwise_nestwasponly.csv")
