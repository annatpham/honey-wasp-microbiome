qpcr <- read.csv("./qpcr.csv") # Load qPCR data

library(dplyr)
# Calculate the median 16S rRNA gene copy number for each insect type
medians <- qpcr %>% 
  group_by(type) %>%
  summarise(med_copies = median(avg))
# Add column to qpcr dataframe using the medians dataframe
qpcr$median <- medians$med_copies[match(qpcr$type, medians$type)] 

library(scales)
library(ggplot2)
library(svglite)

# Specify the order of the labels on the x-axis
qpcr$type <- factor(qpcr$type, levels = c("Brachygastra", "Polistes", "Solitary Wasps", "Native Bees"))

# Save plot as a png and svg file
# Remove "plot <-" on line 21 to save as a png
png("qpcr.png", units="in", width=6.875, height=7, res=600)
plot <- ggplot(qpcr, aes(x = type, y = avg, color = family, shape=type)) +
  scale_y_log10("Bacterial 16S rRNA gene copies",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + # Scientific notation
  geom_jitter(size=3, stroke = 1, width = 0.2) + # Spread out points
  scale_shape_manual(values=c(21, 19, 19, 19), name="type", guide="none") + # Shape of the points
  scale_color_manual(values=c("Andrenidae"="#6F7302", "Apidae"="#B00068", 
                              "Braconidae"="#CBC3E3", "Crabronidae"="#7570B3", 
                              "Halictidae"="#D76865", "Megachilidae"="#D1DBBD", 
                              "Mutillidae"="#e77732", "Pompilidae"="#F0D787", 
                              "Sphecidae"="#1E556F", "Thynnidae"="#A65C32", "Vespidae"="#3CB0CC")) +
  theme_bw(base_size=14) + # Font size
  theme(legend.title=element_blank(),
        legend.text=element_text(face="italic"),
        axis.title.x=element_blank(),
        panel.grid.major = element_blank()) +
  geom_errorbar(aes(ymax=median, ymin=median), color="#000000", linetype="dashed", size=1) + # Add dashed line at median 
  geom_rect(aes(ymin = 0, 
                ymax = 10^3, 
                xmin =Inf, 
                xmax = -Inf), fill = "gray", color = NA, alpha = 0.01) # Add rectangle for y-axis values below 10^3
dev.off()
ggsave(file="qpcr.svg", plot=plot, width=6.875, height=7)

# Check if data is normally distributed
shapiro.test(qpcr$avg)
# Result: W = 0.098994, p-value < 2.2e-16, data is not normally distributed

# Non-parametric Kruskal-Wallis test
kruskal.test(avg ~ type, data = qpcr)
# Result: Kruskal-Wallis chi-squared = 27.726, df = 3, p-value = 4.146e-06

library(rstatix)
# Pairwise comparisons between the insect types
pwc <- qpcr %>% 
  dunn_test(avg ~ type, p.adjust.method = "bonferroni") 
pwc
write.csv(pwc, "./pwc.csv") # Download pairwise comparison data
