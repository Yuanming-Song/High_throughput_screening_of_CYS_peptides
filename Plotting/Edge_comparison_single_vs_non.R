# Custom theme for text sizes
source("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Plotting/WW_scale_base.R")
text_theme <- theme(
  plot.title = element_text(size = plot_title_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.title = element_text(size = legend_title_size),
  legend.text = element_text(size = legend_text_size)
)
pltthemetemp<-plttheme+text_theme

# Define amino acid groups
hydrophobic <- c("A",  "I", "L", "M", "V","P")
polar_uncharged <- c("S", "T", "N", "Q","G")
positive_charged <- c("R", "K", "H")
negative_charged <- c("D", "E")
aromatic <- c("F", "W", "Y")

# Function to classify sequences based on amino acids
classify_sequence <- function(a1, a2, a3) {
  # Convert inputs to character if they aren't already
  a1 <- as.character(a1)
  a2 <- as.character(a2)
  a3 <- as.character(a3)
  
  amino_acids <- c(a1, a2, a3)
  amino_acids <- amino_acids[amino_acids != "C"]  # Ignore C in classification
  
  hydrophobic_count <- sum(amino_acids %in% hydrophobic)
  polar_count <- sum(amino_acids %in% polar_uncharged)
  charged_count <- sum(amino_acids %in% c(positive_charged, negative_charged))
  aromatic_count <- sum(amino_acids %in% aromatic)
  
  # Return a single string value
  if (hydrophobic_count == 2) {
    return("2 Hydrophobic Side Chains")
  } else if (polar_count == 2) {
    return("2 Polar Uncharged Side Chains")
  } else if (charged_count == 2) {
    return("2 Charged Side Chains")
  } else if (aromatic_count == 2) {
    return("2 Aromatic")
  } else if (hydrophobic_count == 1) {
    if (aromatic_count == 1) {
      return("Hydrophobic + Aromatic")
    } else if (polar_count == 1) {
      return("Hydrophobic + Polar")
    } else if (charged_count == 1) {
      return("Hydrophobic + Charged")
    }
  } else if (aromatic_count == 1) {
    if (polar_count == 1) {
      return("Aromatic + Polar")
    } else if (charged_count == 1) {
      return("Aromatic + Charged")
    }
  } else if (polar_count == 1 && charged_count == 1) {
    return("Polar + Charged")
  }
  
  # Default return if no other conditions met
  return("Other")
}

# Load both versions of the data
base_dir <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide"

# Load single node version
load(file.path(base_dir, "tripeptide_network_stats_single_node.rda"))
combined_df_single_mon <- c()
state <- "monomer"
for (i in 1:(length(final_results[[1]]))) {
  combined_df_single_mon <- rbind(combined_df_single_mon, final_results[[1]][[i]])
}
combined_df_single_mon <- cbind(combined_df_single_mon, state)

combined_df_single_dim <- c()
state <- "dimer"
for (i in 1:length(final_results[[2]])) {
  combined_df_single_dim <- rbind(combined_df_single_dim, final_results[[2]][[i]])
}
combined_df_single_dim <- cbind(combined_df_single_dim, state)

combined_df_single <- rbind(combined_df_single_dim, combined_df_single_mon)

# Load non-single node version
load(file.path(base_dir, "tripeptide_network_stats.rda"))
combined_df_non_single <- rbind(final_results[[1]], final_results[[2]])

# Process the edge data for both versions
edge_data_single <- combined_df_single[combined_df_single$stat_name == "edges" & 
                                       combined_df_single$frame == max(combined_df_single$frame), ]

edge_data_non_single <- combined_df_non_single[combined_df_non_single$stat_name == "edges" & 
                                              combined_df_non_single$frame == max(combined_df_non_single$frame), ]

# Function to process data
process_data <- function(data, value_cols) {
  # Create a data frame with sequence, state, and mean value
  means <- aggregate(value ~ seqname + state, data = data, FUN = mean)
  
  # Reshape to wide format
  wide <- reshape(means, 
                  idvar = "seqname", 
                  timevar = "state", 
                  direction = "wide")
  colnames(wide) <- c("seqname", value_cols)
  
  # Add sequence information
  wide$amino1 <- substr(wide$seqname, 1, 1)
  wide$amino2 <- substr(wide$seqname, 2, 2)
  wide$amino3 <- substr(wide$seqname, 3, 3)
  wide$Category <- mapply(classify_sequence, 
                          wide$amino1, 
                          wide$amino2, 
                          wide$amino3)
  wide$label <- wide$seqname
  
  # Ensure Category is a proper factor
  wide$Category <- as.character(wide$Category)
  wide$Category[is.na(wide$Category)] <- "Other"
  wide$Category <- as.factor(wide$Category)
  
  return(wide)
}

# Process both versions
edge_wide_single <- process_data(edge_data_single, c("dimer_edges", "monomer_edges"))
edge_wide_non_single <- process_data(edge_data_non_single, c("dimer_edges", "monomer_edges"))

# Merge the two versions
comparison_data <- merge(edge_wide_single, edge_wide_non_single, 
                        by = c("seqname", "Category", "label"),
                        suffixes = c("_single", "_non_single"))

# Create the comparison plot
edge_comparison_plot <- ggplot(comparison_data, 
            aes(y = dimer_edges_single/300, 
                x = dimer_edges_non_single/150, 
                col = Category,
                text = label)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0.5, linetype = "dashed", color = "blue") +
  labs(x = "Edges per Peptide Dimer Node",
       y = "Edges per Peptide Node",
       col = "Side Chain Types") +
  pltthemetemp +
  theme(legend.position = "none") +
  guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(0, 6)) +
  scale_y_continuous(limits = c(0, 6))

# Print interactive plot
print(ggplotly(edge_comparison_plot))

# Save high-resolution plot
if (exists("save_plots") && save_plots) {
  ggsave(paste0(pltsavedir, "/Edge_comparison_single_vs_non.png"), 
         plot = edge_comparison_plot ,
         dpi = 1100, 
         width = 3.25,
         height = 3,
         units = "in")
} 