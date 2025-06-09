# Main execution
# Set parameters
sequence <- "CPV"  # Change this to your sequence of interest
single_node <- 1 # Whether to use single_node version for dimer
# Load required libraries
step<-50
library(sna)
library(network)

# Function to plot network with nice layout
plot_peptide_network <- function(net, sequence, state, frame) {
  # Fruchterman-Reingold with stiffer springs
  layout_coords <- gplot.layout.fruchtermanreingold(
    net,
    layout.par = list(
      area = 10,      # smaller area = nodes pushed closer together (tighter springs)
      repulse = 10000,   # lower repulsion for tighter clusters
      niter = 9000,      # more iterations for convergence
      arrowhead.cex = 0
    )
  )
  # Plot with Fruchterman-Reingold layout
  gplot(net,
        main = paste(sequence, "-", state, "\nFrame:", frame),
        vertex.cex = 1,      # Larger nodes
        vertex.col = "lightblue",
        edge.col = "darkgrey",
        edge.width = 2,
        mode = "fruchtermanreingold",  # Use Fruchterman-Reingold layout
        displaylabels = FALSE,  # Show node labels
        displayisolates = TRUE,# Show isolated nodes
        coord = layout_coords,
        
        usearrows = FALSE)  
}

# Function to determine peptide type and build path
get_data_path <- function(sequence, single_node = FALSE) {
  n_chars <- nchar(sequence)
  if (n_chars == 2) {
    type <- "Dipeptide"
  } else if (n_chars == 3) {
    type <- "Tripeptide"
  } else if (n_chars == 4) {
    type <- "Tetrapeptide"
  } else {
    stop("Sequence must be 2, 3, or 4 amino acids long")
  }
  
  # Convert sequence to path format (e.g., "ACD" -> "A_C_D")
  path_seq <- paste(strsplit(sequence, "")[[1]], collapse = "_")
  
  # Add single_node suffix only for dimer if specified
  dimer_suffix <- if(single_node) "_single_node" else ""
  
  base_path <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist"
  
  return(list(
    monomer = file.path(base_path, type, "monomer", paste0(path_seq, ".rda")),
    dimer = file.path(base_path, type, "dimer", paste0(path_seq, dimer_suffix, ".rda"))
  ))
}

# Function to create network from edgelist
create_network_from_edgelist <- function(el, frame) {
  # Get number of unique nodes
  nodes <- unique(c(el[[frame]][,1], el[[frame]][,2]))
  n_nodes <- max(nodes)
  
  # Create network object
  net <- network.initialize(n_nodes, directed = FALSE)
  
  # Add edges
  for (i in 1:nrow(el[[frame]])) {
    net[el[[frame]][i,1], el[[frame]][i,2]] <- 1
  }
  
  return(net)
}

# Function to animate network evolution
animate_networks <- function(monomer_el, dimer_el, sequence, frames = NULL, delay = 0.5,step=1) {
  
    frames <- seq(1,length(monomer_el),step)  # Use all frames if not specified
  
  
  # Set up plot parameters
  par(mfrow = c(1,2), mar = c(1,1,4,1))
  
  for (frame in frames) {
    # Create and plot monomer network
    monomer_net <- create_network_from_edgelist(monomer_el, frame)
    plot_peptide_network(monomer_net, sequence, "monomer", frame)
    
    # Create and plot dimer network
    dimer_net <- create_network_from_edgelist(dimer_el, frame)
    plot_peptide_network(dimer_net, sequence, "dimer", frame)
    
    # Print current frame statistics
    cat("\rFrame", frame, "- Monomer edges:", network.edgecount(monomer_net),
        "Dimer edges:", network.edgecount(dimer_net), "     ")
    
    # Pause between frames
    Sys.sleep(delay)
  }
  cat("\n")  # New line after progress updates
}

frames_to_view <- 1:10  # Which frames to animate (NULL for all frames)
frame_delay <- 0.5  # Delay between frames in seconds

# Get data paths
paths <- get_data_path(sequence, single_node)
# cat("Analyzing sequence:", sequence, "\n")
# cat("Loading data from:\n")
# cat("Monomer:", paths$monomer, "\n")
# cat("Dimer:", paths$dimer, "\n\n")

# Load data
load(paths$monomer)
monomer_edgelist <- edgelist
load(paths$dimer)
dimer_edgelist <- edgelist

# Animate networks
#cat("Animating networks...\n")
animate_networks(monomer_edgelist, dimer_edgelist, sequence, frames_to_view, frame_delay,step=step)

# Reset plot parameters
par(mfrow = c(1,1)) 