# To create Google Maps API key:
# https://developers.google.com/maps/documentation/embed/get-api-key
register_google(key="XXX") # Change to API key

# Load packages and metadata
library(ggmap)
library(ggplot2)
library(ggsn) # To modify scale bar

collection <- read.csv("./collection.csv") # Load metadata

centermap <- get_googlemap(center = c(lon = mean(collection$lon), 28.6), 
                      zoom = 7, maptype = "terrain", scale = 2,
                      style = 'feature:road|element:all|visibility:off&style=feature:administrative.locality|element:labels|visibility:simplified')

# Specify the order of the labels in the legend
collection$insect <- factor(collection$insect, levels = c("Honey Wasp", "Other Wasp", "Honey Bee"))

# Create map and save as a png and svg file
library(svglite)

png("centermap.png", units="in", width=7, height=7, res=600)
plot <- ggmap(centermap) + # Remove "plot <-" to save as a png
  geom_point(data = collection, aes(x = lon, y = lat, fill = insect), 
             size = 5, stroke = 1, shape = 21) +
  scale_fill_manual(values = c( "Honey Wasp" = "#006480", "Other Wasp" ="#D76865", "Honey Bee" = "#F0d787")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.justification = c(0.95, 0.2), legend.position = c(0.95, 0.2)) +
  ggsn::scalebar(x.min = -95.8, x.max = -94.8, 
                 y.min = 25.9, y.max = 26.1, 
                 dist = 100, dist_unit = "km", transform = TRUE, 
                 model = "WGS84", height = 0.5, 
                 st.size = 4, st.dist = 0.5) 
ggsave(file="centermap.svg", plot=plot, width=7, height=7)
dev.off()
