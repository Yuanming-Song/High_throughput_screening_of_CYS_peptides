# Custom theme for text sizes
text_theme <- theme(
  plot.title = element_text(size = plot_title_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.title = element_text(size = legend_title_size),
  legend.text = element_text(size = legend_text_size)
)
pltthemetemp<-plttheme+text_theme

# Load required libraries
library(igraph)
library(ggplot2)
library(cowplot)

# Create peptide-per-node network
peptide_edges <- matrix(c(
  "A1", "A2",  # Disulfide bond
  "A1", "B1",  # A1's 5 unique partners
  "A1", "C1",
  "A1", "D1",
  "A1", "E1",
  "A1", "F1",
  "A2", "G2",  # A2's 3 unique partners
  "A2", "H2",
  "A2", "I2"
), ncol = 2, byrow = TRUE)

peptide_g <- graph_from_edgelist(peptide_edges, directed = FALSE)
V(peptide_g)$color <- ifelse(grepl("^A", V(peptide_g)$name), "red", "gray")
# Set edge colors using indices
E(peptide_g)$color <- "black"  # Default color for all edges
E(peptide_g)[1]$color <- "red"  # First edge (A1-A2) is red
E(peptide_g)$width <- 2  # Default width
E(peptide_g)[1]$width <- 3  # First edge (A1-A2) is thicker

# Create molecule-per-node network
molecule_edges <- matrix(c(
  "A", "B",  # All 8 edges from the combined A1-A2
  "A", "C",
  "A", "D",
  "A", "E",
  "A", "F",
  "A", "G",
  "A", "H",
  "A", "I"
), ncol = 2, byrow = TRUE)

molecule_g <- graph_from_edgelist(molecule_edges, directed = FALSE)
V(molecule_g)$color <- ifelse(V(molecule_g)$name == "A", "red", "gray")
E(molecule_g)$width <- 2

# Set up the plotting area
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))

# Plot peptide-per-node network
plot(peptide_g, 
     layout = layout_with_fr,
     vertex.label = V(peptide_g)$name,
     vertex.size = 30,  # Increased node size
     vertex.label.color = "black",
     vertex.label.cex = 0.8,  # Reduced label size for peptide network
     edge.width = E(peptide_g)$width,
     edge.curved = 0.2,  # Make edges slightly curved to reduce overlap
     main = "")

# Plot molecule-per-node network
plot(molecule_g, 
     layout = layout_with_fr,
     vertex.label = V(molecule_g)$name,
     vertex.size = 30,  # Increased node size
     vertex.label.color = "black",
     vertex.label.cex = 1.2,  # Keep original size for molecule network
     edge.width = E(molecule_g)$width,
     edge.curved = 0.2,  # Make edges slightly curved to reduce overlap
     main = "")

# Save the plot if save_plots is TRUE
if (exists("save_plots") && save_plots) {
  # Save as PNG
  png(paste0(pltsavedir, "/Edge_definition_comparison.png"),
      width = 6.5, height = 3, units = "in", res = 1100)
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
  
  # Plot peptide-per-node network
  plot(peptide_g, 
       layout = layout_with_fr,
       vertex.label = V(peptide_g)$name,
       vertex.size = 30,  # Increased node size
       vertex.label.color = "black",
       vertex.label.cex = 0.8,  # Reduced label size for peptide network
       edge.width = E(peptide_g)$width,
       edge.curved = 0.2,  # Make edges slightly curved to reduce overlap
       main = "")
  
  # Plot molecule-per-node network
  plot(molecule_g, 
       layout = layout_with_fr,
       vertex.label = V(molecule_g)$name,
       vertex.size = 30,  # Increased node size
       vertex.label.color = "black",
       vertex.label.cex = 1.2,  # Keep original size for molecule network
       edge.width = E(molecule_g)$width,
       edge.curved = 0.2,  # Make edges slightly curved to reduce overlap
       main = "")
  
  dev.off()
} 