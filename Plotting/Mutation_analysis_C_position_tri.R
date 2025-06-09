# Custom theme for text sizes
text_theme <- theme(
  plot.title = element_text(size = plot_title_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.title = element_text(size = legend_title_size),
  legend.text = element_text(size = legend_text_size)
)
pltthemetemp<-plttheme+text_theme
if (loadnetstat) {
  #Load the network stats data
  load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/tripeptide_network_stats_single_node.rda")
  
  # Combine monomer and dimer data
  combined_df_ddedge_tripep_mon<-c()
  state<-"monomer"
  for (i in 1:(length(final_results[[1]]))) {
    combined_df_ddedge_tripep_mon<-rbind(combined_df_ddedge_tripep_mon,final_results[[1]][[i]])
  }
  combined_df_ddedge_tripep_mon<-cbind(combined_df_ddedge_tripep_mon,state)
  
  combined_df_ddedge_tripep_dim<-c()
  state<-"dimer"
  for (i in 1:length(final_results[[2]])) {
    combined_df_ddedge_tripep_dim<-rbind(combined_df_ddedge_tripep_dim,final_results[[2]][[i]])
  }
  combined_df_ddedge_tripep_dim<-cbind(combined_df_ddedge_tripep_dim,state)
  combined_df_ddedge_tripep <- rbind(combined_df_ddedge_tripep_dim, combined_df_ddedge_tripep_mon)
}
# Define residue groups and colors
acidic <- c("D", "E")
basic <- c("R", "K", "H")
polar <- c("S", "T", "N", "Q")
aromatic <- c("F", "W", "Y")
aliphatic <- c("A", "G", "I", "L", "M", "P", "V")

# Define colors for each group
group_colors <- c(
  setNames(rep("red", length(acidic)), acidic),
  setNames(rep("blue", length(basic)), basic),
  setNames(rep("orange", length(polar)), polar),
  setNames(rep("darkgreen", length(aliphatic)), aliphatic),
  setNames(rep("black", length(aromatic)), aromatic)
)

# Create ordered residue list
ordered_residues <- c(acidic, basic, polar, aliphatic, aromatic)

# Initialize list for mutation effects by position
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")
mutation_effects_by_pos <- list()

# Pre-process the data: create a lookup table for edge differences
# First, filter for only edge statistics
edge_data <- combined_df_ddedge_tripep[combined_df_ddedge_tripep$stat_name == "edges" & combined_df_ddedge_tripep$frame == max(combined_df_ddedge_tripep$frame), ]

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

# Function to get position of C
get_C_position <- function(seq) {
  seq_chars <- strsplit(seq, "")[[1]]
  return(which(seq_chars == "C"))
}

# Pre-compute all possible mutations
all_seqs <- unique(edge_wide$seqname)
mutation_pairs <- expand.grid(from = residues, to = residues, c_pos = 1:3)
mutation_pairs <- mutation_pairs[mutation_pairs$from != mutation_pairs$to, ]

# Calculate mutation effects using parallel processing
library(parallel)
library(ggplot2)
library(plotly)
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)

# Export necessary objects to the cluster
clusterExport(cl, c("residues", "all_seqs", "edge_lookup", "get_edge_difference", "get_C_position"))

# Function to process a single mutation based on C position
process_mutation_pos <- function(mutation) {
  from_res <- mutation[[1]]
  to_res <- mutation[[2]]
  target_c_pos <- as.numeric(mutation[3])
  deltas <- c()
  
  # Get all sequences where C is at the target position
  target_seqs <- all_seqs[sapply(all_seqs, function(seq) {
    c_pos <- get_C_position(seq)
    return(c_pos == target_c_pos && any(strsplit(seq, "")[[1]] == from_res))
  })]
  
  for (seq in target_seqs) {
    seq_res <- strsplit(seq, "")[[1]]
    mut_positions <- which(seq_res == from_res)
    
    for (pos in mut_positions) {
      mutated_seq <- seq_res
      mutated_seq[pos] <- as.character(to_res)
      mutated_seq <- paste(mutated_seq, collapse = "")
      
      orig_diff <- get_edge_difference(seq)
      mut_diff <- get_edge_difference(mutated_seq)
      
      if (!is.na(orig_diff) && !is.na(mut_diff)) {
        deltas <- c(deltas, mut_diff - orig_diff)
      }
    }
  }
  
  if (length(deltas) > 0) {
    return(c(mean(deltas), sd(deltas), min(deltas), max(deltas), median(deltas), length(deltas)))
  } else {
    return(c(NA, NA, NA, NA, NA, 0))
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
  C_Position = mutation_pairs$c_pos,
  DeltaDeltaEdges = mutation_results[1,],
  StdDevEdges = mutation_results[2,],
  MinEdges = mutation_results[3,],
  MaxEdges = mutation_results[4,],
  MedianEdges = mutation_results[5,],
  Count = mutation_results[6,]
)

# Create position labels and ensure ordering
plot_df$From <- factor(plot_df$From, levels = ordered_residues)
plot_df$To <- factor(plot_df$To, levels = ordered_residues)
plot_df$C_Position <- factor(plot_df$C_Position, 
                             levels = 1:3,
                             labels = paste0("C at position ", 1:3))

# Define common theme elements
colored_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, color = group_colors[ordered_residues]),
  axis.text.y = element_text(color = group_colors[ordered_residues]),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# Create the mean plot
mut_detail_plt <- ggplot(plot_df, aes(x = From, y = To, fill = DeltaDeltaEdges/300)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "grey", "yellow", "red"),
    values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
    limits = c(min(plot_df$DeltaDeltaEdges)/300, max(plot_df$DeltaDeltaEdges)/300),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Average Change in ΔEdges by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "ΔΔEdges") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create the standard deviation plot
mut_detail_std_plt <- ggplot(plot_df, aes(x = From, y = To, fill = StdDevEdges/300)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "yellow", "red"),
    values = scales::rescale(c(0, 1, 2), from = c(0, max(plot_df$StdDevEdges/300, na.rm=TRUE))),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Standard Deviation of ΔEdges by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Std Dev") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Get overall min and max for consistent scaling
overall_min <- min(c(plot_df$MinEdges, plot_df$MaxEdges, plot_df$MedianEdges), na.rm = TRUE)
overall_max <- max(c(plot_df$MinEdges, plot_df$MaxEdges, plot_df$MedianEdges), na.rm = TRUE)

# Create the minimum plot
mut_detail_min_plt <- ggplot(plot_df, aes(x = From, y = To, fill = MinEdges/300)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "grey", "yellow", "red"),
    values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
    limits = c(overall_min/300, overall_max/300),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Minimum Change in ΔEdges by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Min ΔΔEdges") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create the maximum plot
mut_detail_max_plt <- ggplot(plot_df, aes(x = From, y = To, fill = MaxEdges/300)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "grey", "yellow", "red"),
    values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
    limits = c(overall_min/300, overall_max/300),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Maximum Change in ΔEdges by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Max ΔΔEdges") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create the median plot
mut_detail_median_plt <- ggplot(plot_df, aes(x = From, y = To, fill = MedianEdges/300)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "grey", "yellow", "red"),
    values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
    limits = c(overall_min/300, overall_max/300),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Median Change in ΔEdges by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Median ΔΔEdges") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create interactive plots
print(ggplotly(mut_detail_plt))
print(ggplotly(mut_detail_std_plt))
print(ggplotly(mut_detail_min_plt))
print(ggplotly(mut_detail_max_plt))
print(ggplotly(mut_detail_median_plt))

# Save high-resolution plots
if (exists("save_plots") && save_plots) {
  ggsave(paste0(pltsavedir, "/Tripeptide_mutation_edges_C_absolute.png"), 
         plot = mut_detail_plt+theme(legend.position = "bottom",
                                     legend.direction = "horizontal"
                                     ),
         dpi = 1100, 
         width = 6.5,
         height = 3.5,
         units = "in")
  if (otherstatpltsave) {
    
  ggsave(paste0(pltsavedir, "/Tripeptide_mutation_edges_C_std_dev.png"), 
         plot = mut_detail_std_plt,
         dpi = 1100, 
         width = 32,
         height = 16,
         units = "cm")
  
  ggsave(paste0(pltsavedir, "/Tripeptide_mutation_edges_C_min.png"), 
         plot = mut_detail_min_plt,
         dpi = 1100, 
         width = 32,
         height = 16,
         units = "cm")
  
  ggsave(paste0(pltsavedir, "/Tripeptide_mutation_edges_C_max.png"), 
         plot = mut_detail_max_plt,
         dpi = 1100, 
         width = 32,
         height = 16,
         units = "cm")
  
  ggsave(paste0(pltsavedir, "/Tripeptide_mutation_edges_C_median.png"), 
         plot = mut_detail_median_plt,
         dpi = 1100, 
         width = 32,
         height = 16,
         units = "cm")
  } 
}