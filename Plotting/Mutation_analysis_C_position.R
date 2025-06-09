# Load the network stats data
if (loadnetstat) {
  base_dir <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide"
  if (single_node) {
    load(file.path(base_dir, "tripeptide_network_stats_single_node.rda"))
    # Combine monomer and dimer data
    combined_df_ddedge_tripep_mon <- c()
    state <- "monomer"
    for (i in 1:(length(final_results[[1]]))) {
      combined_df_ddedge_tripep_mon <- rbind(combined_df_ddedge_tripep_mon, final_results[[1]][[i]])
    }
    combined_df_ddedge_tripep_mon <- cbind(combined_df_ddedge_tripep_mon, state)
    
    combined_df_ddedge_tripep_dim <- c()
    state <- "dimer"
    for (i in 1:length(final_results[[2]])) {
      combined_df_ddedge_tripep_dim <- rbind(combined_df_ddedge_tripep_dim, final_results[[2]][[i]])
    }
    combined_df_ddedge_tripep_dim <- cbind(combined_df_ddedge_tripep_dim, state)
  } else {
    load(file.path(base_dir, "tripeptide_network_stats.rda"))
    combined_df_ddedge_tripep_mon <- final_results[[1]]
    combined_df_ddedge_tripep_dim <- final_results[[2]]
  }
  combined_df_ddedge_tripep <- rbind(combined_df_ddedge_tripep_dim, combined_df_ddedge_tripep_mon)
}

# Define residue groups
acidic <- c("D", "E")
basic <- c("R", "K", "H")
polar <- c("S", "T", "N", "Q")
aromatic <- c("F", "W", "Y")
aliphatic <- c("A", "G", "I", "L", "M", "P", "V")

# Create ordered residue list
ordered_residues <- c(acidic, basic, polar, aliphatic, aromatic)

# Initialize list for mutation effects by position
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Pre-process the data: create a lookup table for edge differences
# First, filter for only edge statistics
edge_data_tripep <- combined_df_ddedge_tripep[combined_df_ddedge_tripep$stat_name == "edges", ]

# Create a data frame with sequence, state, and mean value
edge_means_tripep <- aggregate(value ~ seqname + state, data = edge_data_tripep, FUN = mean)

# Reshape to wide format for easy difference calculation
edge_wide_tripep <- reshape(edge_means_tripep, 
                    idvar = "seqname", 
                    timevar = "state", 
                    direction = "wide")
edge_wide_tripep$edge_diff <- edge_wide_tripep$value.dimer - edge_wide_tripep$value.monomer

# Create a lookup table for quick access
edge_lookup_tripep <- setNames(edge_wide_tripep$edge_diff, edge_wide_tripep$seqname)

# Function to get edge difference (now just a lookup)
get_edge_difference <- function(seq) {
  return(edge_lookup_tripep[seq])
}

# Pre-compute all possible mutations
all_seqs_tripep <- unique(edge_wide_tripep$seqname)
mutation_pairs <- expand.grid(from = residues, to = residues)
mutation_pairs <- mutation_pairs[mutation_pairs$from != mutation_pairs$to, ]

# Calculate mutation effects using parallel processing
library(parallel)
library(ggplot2)
library(plotly)
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)

# Export necessary objects to the cluster
clusterExport(cl, c("residues", "all_seqs_tripep", "edge_lookup_tripep", "get_edge_difference"))

# Function to process a single mutation for sequences with C at specific position
process_mutation_c_pos <- function(mutation, c_position) {
  from_res <- mutation[1]
  to_res <- mutation[2]
  deltas <- c()
  
  # Get all sequences containing from_res and C at specified position
  target_seqs <- all_seqs_tripep[sapply(all_seqs_tripep, function(seq) {
    seq_res <- strsplit(seq, "")[[1]]
    seq_res[c_position] == "C"  # Check if C is at the specified position
  })]
  
  for (seq in target_seqs) {
    seq_res <- strsplit(seq, "")[[1]]
    for (pos in 1:3) {
      if (pos != c_position && seq_res[pos] == from_res) {  # Don't mutate C position
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
    return(c(mean(deltas), length(deltas)))
  } else {
    return(c(NA, 0))
  }
}

# Process mutations for each C position
results_by_c_pos <- list()
for (c_pos in 1:3) {
  # Process mutations in parallel
  mutation_results <- parApply(cl, mutation_pairs, 1, process_mutation_c_pos, c_position = c_pos)
  
  # Create a data frame for plotting
  plot_df <- data.frame(
    From = mutation_pairs$from,
    To = mutation_pairs$to,
    DeltaDeltaEdges = mutation_results[1,],
    Count = mutation_results[2,]
  )
  
  results_by_c_pos[[c_pos]] <- plot_df
}

# Stop the cluster
stopCluster(cl)

# Combine results
plot_df <- do.call(rbind, results_by_c_pos)
plot_df$C_Position <- rep(paste0("C at pos ", 1:3), each = nrow(mutation_pairs))

# Create the plot
mut_detail_plt <- ggplot(plot_df, aes(x = From, y = To, fill = DeltaDeltaEdges)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0, na.value = "white") +
  scale_x_discrete(limits = ordered_residues) +
  scale_y_discrete(limits = rev(ordered_residues)) +
  theme_minimal() +
  labs(title = "Average Change in ΔEdges by Mutation with C Position",
       x = "Original Residue",
       y = "Mutated To",
       fill = "ΔΔEdges") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create interactive plot
print(ggplotly(mut_detail_plt))

# Save high-resolution plot
if (exists("save_plots") && save_plots) {
  ggsave(paste0(APplotDir, "/Tripeptide_mutation_edges_C_position.png"), 
         plot = mut_detail_plt,
         dpi = 1100, 
         width = 32,
         height = 16,
         units = "cm")
} 