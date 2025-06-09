# Load required libraries
library(network)
library(ergm)
library(ggplot2)
library(gridExtra)

# Create a random graph with 5 nodes
n_nodes <- 5
n_edges <- sample(5:10, 1)  # Random number of edges between 5 and 10
edges <- matrix(sample(1:n_nodes, 2*n_edges, replace=TRUE), ncol=2)  # Random edges
net <- network(edges, directed=FALSE, loops=FALSE, matrix.type="edgelist", num.vertices=n_nodes)

# Calculate statistics
stats <- summary(net ~ edges + nsp(1) + esp(0) + esp(1))

# Create the graph plot with minimal margins
par(mar = c(1, 1, 1, 8))  # Reduced margins, keeping right margin for text
plot(net, vertex.col = "lightblue", vertex.cex = 4)  # Increased node size

# Create a text grob for statistics
stats_text <- paste(
  "Network Statistics:",
  paste("Edges:", stats[1]),
  paste("NSP(1):", stats[2]),
  paste("ESP(0):", stats[3]),
  paste("ESP(1):", stats[4]),
  sep = "\n"
)

# Add statistics as a text box
text(x = par("usr")[2] - 1, y = mean(par("usr")[3:4]), 
     labels = stats_text, 
     adj = c(0, 0.5),
     cex = 1,  # This will give approximately size 12 text
     xpd = TRUE)  # Allow text to be drawn outside the plot area

# Save the plot
if (exists("save_plots") && save_plots) {
  png(filename = file.path(pltsavedir, "poster_network_example.png"),
      width = 6.5, height = 2, units = "in", res = 1100)
  
  # Recreate the plot
  par(mar = c(1, 1, 1, 8))  # Reduced margins, keeping right margin for text
  plot(net, vertex.col = "lightblue", vertex.cex = 6)  # Increased node size
  text(x = par("usr")[2] - 1, y = mean(par("usr")[3:4]), 
       labels = stats_text, 
       adj = c(0, 0.5),
       cex = 1,  # This will give approximately size 12 text
       xpd = TRUE)
  
  dev.off()
} 