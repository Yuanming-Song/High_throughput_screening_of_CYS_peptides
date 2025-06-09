# Load the network stats data
load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/tripeptide_network_stats_single_node.rda")

# Combine monomer and dimer data
#combined_df <- rbind(final_results[[1]], final_results[[2]])
combined_df_mon<-c()
state<-"monomer"
for (i in 1:(length(final_results[[1]]))) {
  combined_df_mon<-rbind(combined_df_mon,final_results[[1]][[i]])
}
combined_df_mon<-cbind(combined_df_mon,state)

combined_df_dim<-c()
state<-"dimer"
for (i in 1:length(final_results[[2]])) {
  combined_df_dim<-rbind(combined_df_dim,final_results[[2]][[i]])
}
combined_df_dim<-cbind(combined_df_dim,state)
combined_df <- rbind(combined_df_dim, combined_df_mon)

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
mutation_effects_by_pos <- list()

# Pre-process the data: create a lookup table for edge differences
# First, filter for only edge statistics
edge_data <- combined_df[combined_df$stat_name == "edges" & combined_df$frame == 202, ]

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

# Pre-compute all possible mutations for each position
all_seqs <- unique(edge_wide$seqname)
mutation_pairs <- expand.grid(from = residues, to = residues, pos = 1:3)
mutation_pairs <- mutation_pairs[mutation_pairs$from != mutation_pairs$to, ]

# Calculate mutation effects using parallel processing
library(parallel)
library(ggplot2)
library(plotly)
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)

# Export necessary objects to the cluster
clusterExport(cl, c("residues", "all_seqs", "edge_lookup", "get_edge_difference"))

# Function to process a single mutation with position
process_mutation_pos <- function(mutation) {
  from_res <- mutation[[1]]
  to_res <- mutation[[2]]
  pos <- as.numeric(mutation[3])
  deltas <- c()
  
  # Get all sequences containing from_res at the specified position
  target_seqs <- all_seqs[sapply(all_seqs, function(seq) {
    seq_res <- strsplit(seq, "")[[1]]
    seq_res[pos] == from_res
  })]
  
  for (seq in target_seqs) {
    seq_res <- strsplit(seq, "")[[1]]
    mutated_seq <- seq_res
    mutated_seq[pos] <- as.character(to_res)
    mutated_seq <- paste(mutated_seq, collapse = "")
    
    orig_diff <- get_edge_difference(seq)
    mut_diff <- get_edge_difference(mutated_seq)
    
    if (!is.na(orig_diff) && !is.na(mut_diff)) {
      deltas <- c(deltas, mut_diff - orig_diff)
    }
  }
  
  if (length(deltas) > 0) {
    return(c(mean(deltas), length(deltas)))
  } else {
    return(c(NA, 0))
  }
}

# Process mutations in parallel
mutation_results <- parApply(cl, mutation_pairs, 1, process_mutation_pos)

# Stop the cluster
stopCluster(cl)

# Create a data frame for plotting
plot_df <- data.frame(
  From = mutation_pairs$from,
  To = mutation_pairs$to,
  Position = mutation_pairs$pos,
  DeltaDeltaEdges = mutation_results[1,],
  Count = mutation_results[2,]
)

# Create simplified position labels and ensure ordering
plot_df$pos_label <- paste0(plot_df$Position, plot_df$From)
plot_df$From <- factor(plot_df$From, levels = ordered_residues)
plot_df$To <- factor(plot_df$To, levels = ordered_residues)
plot_df$pos_label <- factor(plot_df$pos_label, 
                           levels = unlist(lapply(1:3, function(pos) paste0(pos, ordered_residues))))

# Create the plot
mut_detail_plt <- ggplot(plot_df, aes(x = pos_label, y = To, fill = DeltaDeltaEdges/300)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "grey", "yellow", "red"),
    values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
    limits = c(min(plot_df$DeltaDeltaEdges)/300, max(plot_df$DeltaDeltaEdges)/300),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  theme_minimal() +
  labs(title = "Average Change in ΔEdges by Mutation and Position",
       x = "Position and Original Residue",
       y = "Mutated To",
       fill = "ΔΔEdges") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(. ~ Position, scales = "free_x", space = "free_x")

# Create interactive plot
print(ggplotly(mut_detail_plt))

# Save high-resolution plot
if (exists("save_plots") && save_plots) {
  ggsave(paste0(APplotDir, "/Tripeptide_mutation_edges_detailed.png"), 
         plot = mut_detail_plt,
         dpi = 1100, 
         width = 32,  # Increased width to accommodate the additional detail
         height = 16,
         units = "cm")
} 