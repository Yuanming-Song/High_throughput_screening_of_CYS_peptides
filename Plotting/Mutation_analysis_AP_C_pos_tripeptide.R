# Custom theme for text sizes
text_theme <- theme(
  plot.title = element_text(size = plot_title_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.title = element_text(size = legend_title_size),
  legend.text = element_text(size = legend_text_size)
)
pltthemetemp<-plttheme+text_theme
if (readAP) {
  SASA_files <- list()
  for (pos in 1:3) {
    SASA_files[[pos]] <- read.table(
      paste0("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Tripeptide/SASA_score/SASA_result_with_common_mon_ini_C", pos, ".txt"),
      header = FALSE,
      col.names = c("Amino1", "Amino2", "Amino3", "ratio_mon_fin", "ratio_dim_fin")
    )
    SASA_files[[pos]]$label <- paste(SASA_files[[pos]]$Amino1, SASA_files[[pos]]$Amino2, SASA_files[[pos]]$Amino3, sep="_")
  }
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

# Initialize variables
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Pre-process the data: Calculate ΔAP for each sequence
for (pos in 1:3) {
  SASA_files[[pos]]$delta_AP <- SASA_files[[pos]]$ratio_dim_fin - SASA_files[[pos]]$ratio_mon_fin
}

# Create lookup tables for each C position
AP_lookup <- list()
for (pos in 1:3) {
  AP_lookup[[pos]] <- setNames(SASA_files[[pos]]$delta_AP, SASA_files[[pos]]$label)
}

# Function to get AP difference
get_AP_difference <- function(seq, c_pos) {
  return(AP_lookup[[c_pos]][seq])
}

# Pre-compute all possible mutations
mutation_pairs <- expand.grid(from = residues, to = residues, c_pos = 1:3)
mutation_pairs <- mutation_pairs[mutation_pairs$from != mutation_pairs$to, ]

# Calculate mutation effects using parallel processing
library(parallel)
library(ggplot2)
library(plotly)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export necessary objects to the cluster
clusterExport(cl, c("residues", "AP_lookup", "get_AP_difference"))

# Function to process a single mutation based on C position
process_mutation_pos <- function(mutation) {
  from_res <- mutation[[1]]
  to_res <- mutation[[2]]
  target_c_pos <- as.numeric(mutation[3])
  deltas <- c()
  
  # Get all sequences from the current C position file
  current_seqs <- names(AP_lookup[[target_c_pos]])
  target_seqs <- current_seqs[sapply(strsplit(current_seqs, "_"), function(x) any(x == from_res))]
  
  for (seq in target_seqs) {
    seq_res <- strsplit(seq, "_")[[1]]
    mut_positions <- which(seq_res == from_res)
    
    for (pos in mut_positions) {
      mutated_seq <- seq_res
      mutated_seq[pos] <- to_res
      mutated_seq <- paste(mutated_seq, collapse="_")
      
      orig_diff <- get_AP_difference(seq, target_c_pos)
      mut_diff <- get_AP_difference(mutated_seq, target_c_pos)
      
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
  DeltaDeltaAP = mutation_results[1,],
  StdDevAP = mutation_results[2,],
  MinAP = mutation_results[3,],
  MaxAP = mutation_results[4,],
  MedianAP = mutation_results[5,],
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
mut_detail_plt <- ggplot(plot_df, aes(x = From, y = To, fill = DeltaDeltaAP)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("blue", "grey", "yellow", "red"),
    values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
    limits = c(min(plot_df$DeltaDeltaAP), max(plot_df$DeltaDeltaAP)),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Average Change in ΔAP by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "ΔΔAP") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create the standard deviation plot
mut_detail_std_plt <- ggplot(plot_df, aes(x = From, y = To, fill = StdDevAP)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "yellow", "red"),
    values = scales::rescale(c(0, 1, 2), from = c(0, max(plot_df$StdDevAP, na.rm=TRUE))),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Standard Deviation of ΔAP by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Std Dev") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Get overall min and max for consistent scaling
overall_min <- min(c(plot_df$MinAP, plot_df$MaxAP, plot_df$MedianAP), na.rm = TRUE)
overall_max <- max(c(plot_df$MinAP, plot_df$MaxAP, plot_df$MedianAP), na.rm = TRUE)

# Create the minimum plot
mut_detail_min_plt <- ggplot(plot_df, aes(x = From, y = To, fill = MinAP)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 0,
    limits = c(overall_min, overall_max),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Minimum Change in ΔAP by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Min ΔΔAP") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create the maximum plot
mut_detail_max_plt <- ggplot(plot_df, aes(x = From, y = To, fill = MaxAP)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 0,
    limits = c(overall_min, overall_max),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Maximum Change in ΔAP by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Max ΔΔAP") +
  colored_theme +
  facet_grid(. ~ C_Position, scales = "free_x", space = "free_x")

# Create the median plot
mut_detail_median_plt <- ggplot(plot_df, aes(x = From, y = To, fill = MedianAP)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 0,
    limits = c(overall_min, overall_max),
    na.value = "black") +
  scale_y_discrete(limits = rev(ordered_residues)) +
  pltthemetemp +
  labs(title = "Median Change in ΔAP by Mutation (Grouped by C Position)",
       x = "Original Residue",
       y = "Mutated To",
       fill = "Median ΔΔAP") +
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
  ggsave(paste0(pltsavedir, "/Tripeptide_mutation_AP_C_absolute.png"), 
         plot = mut_detail_plt,
         dpi = 1100, 
         width = 6.5,
         height = 3.5,
         units = "in")
  if(0) {
    ggsave(paste0(pltsavedir, "/Tripeptide_mutation_AP_C_std_dev.png"), 
           plot = mut_detail_std_plt,
           dpi = 1100, 
           width = 32,
           height = 16,
           units = "cm")
    
    ggsave(paste0(pltsavedir, "/Tripeptide_mutation_AP_C_min.png"), 
           plot = mut_detail_min_plt,
           dpi = 1100, 
           width = 32,
           height = 16,
           units = "cm")
    
    ggsave(paste0(pltsavedir, "/Tripeptide_mutation_AP_C_max.png"), 
           plot = mut_detail_max_plt,
           dpi = 1100, 
           width = 32,
           height = 16,
           units = "cm")
    
    ggsave(paste0(pltsavedir, "/Tripeptide_mutation_AP_C_median.png"), 
           plot = mut_detail_median_plt,
           dpi = 1100, 
           width = 32,
           height = 16,
           units = "cm")
  } 
}