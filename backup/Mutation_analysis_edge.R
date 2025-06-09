# Load the network stats data
load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/tripeptide_network_stats.rda")

# Combine monomer and dimer data
combined_df <- rbind(final_results[[1]], final_results[[2]])

# Define residue groups
acidic <- c("D", "E")
basic <- c("R", "K", "H")
polar <- c("S", "T", "N", "Q")
aromatic <- c("F", "W", "Y")
aliphatic <- c("A", "G", "I", "L", "M", "P", "V")

# Create ordered residue list
ordered_residues <- c(acidic, basic, polar, aliphatic, aromatic)

# Initialize matrix for mutation effects
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")
mutation_effects <- matrix(NA, nrow = length(residues), ncol = length(residues),
                          dimnames = list(from = residues, to = residues))

# Pre-process the data: create a lookup table for edge differences
# First, filter for only edge statistics
edge_data <- combined_df[combined_df$stat_name == "edges", ]

# Create a data frame with sequence, state, and mean value
edge_means <- aggregate(value ~ seqname + state, data = edge_data, FUN = mean)

# Reshape to wide format for easy difference calculation
edge_wide <- reshape(edge_means, 
                    idvar = "seqname", 
                    timevar = "state", 
                    direction = "wide")
edge_wide$edge_diff <- edge_wide$value.dimer - edge_wide$value.monomer

# Create a lookup table for quick access
edge_lookup <- setNames(edge_wide$edge_diff, edge_wide$seqname)

# Function to get edge difference (now just a lookup)
get_edge_difference <- function(seq) {
  return(edge_lookup[seq])
}

# Pre-compute all possible mutations
all_seqs <- unique(edge_wide$seqname)
mutation_pairs <- expand.grid(from = residues, to = residues)
mutation_pairs <- mutation_pairs[mutation_pairs$from != mutation_pairs$to, ]

# Calculate mutation effects using parallel processing
library(parallel)
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)

# Export necessary objects to the cluster
clusterExport(cl, c("residues", "all_seqs", "edge_lookup", "get_edge_difference"))

# Function to process a single mutation pair
process_mutation <- function(mutation) {
  from_res <- mutation[1]
  to_res <- mutation[2]
  deltas <- c()
  
  # Get all sequences containing from_res
  target_seqs <- all_seqs[grep(from_res, all_seqs)]
  
  for (seq in target_seqs) {
    seq_res <- strsplit(seq, "")[[1]]
    for (pos in 1:3) {
      if (seq_res[pos] == from_res) {
        mutated_seq <- seq_res
        mutated_seq[pos] <- to_res
        mutated_seq <- paste(mutated_seq, collapse = "")
        
        orig_diff <- get_edge_difference(seq)
        mut_diff <- get_edge_difference(mutated_seq)
        
        if (!is.na(orig_diff) && !is.na(mut_diff)) {
          deltas <- c(deltas, mut_diff - orig_diff)
        }
      }
    }
  }
  
  if (length(deltas) > 0) {
    return(mean(deltas))
  } else {
    return(NA)
  }
}

# Process mutations in parallel
mutation_results <- parApply(cl, mutation_pairs, 1, process_mutation)

# Stop the cluster
stopCluster(cl)

# Fill the mutation effects matrix
for (i in 1:nrow(mutation_pairs)) {
  mutation_effects[mutation_pairs$from[i], mutation_pairs$to[i]] <- mutation_results[i]
}

# Reorder the mutation effects matrix
mutation_effects <- mutation_effects[ordered_residues, ordered_residues]

# Convert to data frame (no filtering for triangle)
heatmap_df <- as.data.frame(as.table(mutation_effects))
colnames(heatmap_df) <- c("From", "To", "DeltaDeltaEdges")

# Add a column to indicate if it's on the diagonal or symmetric
heatmap_df$From_idx <- match(heatmap_df$From, ordered_residues)
heatmap_df$To_idx <- match(heatmap_df$To, ordered_residues)
heatmap_df$is_diagonal <- heatmap_df$From == heatmap_df$To
heatmap_df$is_symmetric <- heatmap_df$From_idx > heatmap_df$To_idx

# Plotting
mut_edge_plt <- ggplot(heatmap_df, aes(x = From, y = To, fill = DeltaDeltaEdges)) +
  geom_tile() +
  # Add diagonal line
  geom_abline(intercept = 20, slope = -1, color = "black", linetype = "dashed") +
  scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0, na.value = "white") +
  scale_x_discrete(limits = ordered_residues) +
  scale_y_discrete(limits = rev(ordered_residues)) +
  theme_minimal() +
  labs(title = "Average Change in ΔEdges by Mutation", fill = "ΔΔEdges") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Add text annotations for symmetric pairs
  geom_text(data = subset(heatmap_df, is_symmetric), 
            aes(label = "*"), 
            color = "black", 
            size = 3)

# Create interactive plot
print(ggplotly(mut_edge_plt))

# Save high-resolution plot
if (save_plots) {
  ggsave(paste0(APplotDir, "/Tripeptide_mutation_edges.png"), 
         plot = mut_edge_plt,
         dpi = 1100, 
         width = 24, 
         height = 16,
         units = "cm")
} 
