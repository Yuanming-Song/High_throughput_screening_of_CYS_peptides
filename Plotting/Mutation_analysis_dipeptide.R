# Load data if needed
if (loadnetstat) {
  #Load the network stats data for dipeptide
  load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Dipeptide/dipeptide_network_stats_single_node.rda")
  
  # Combine monomer and dimer data
  combined_df_net_di_mon <- c()
  state <- "monomer"
  for (i in 1:(length(final_results[[1]]))) {
    combined_df_net_di_mon <- rbind(combined_df_net_di_mon, final_results[[1]][[i]])
  }
  combined_df_net_di_mon <- cbind(combined_df_net_di_mon, state)
  
  combined_df_net_di_dim <- c()
  state <- "dimer"
  for (i in 1:length(final_results[[2]])) {
    combined_df_net_di_dim <- rbind(combined_df_net_di_dim, final_results[[2]][[i]])
  }
  combined_df_net_di_dim <- cbind(combined_df_net_di_dim, state)
  combined_df_net_di <- rbind(combined_df_net_di_dim, combined_df_net_di_mon)
}

if (readAP) {
  # Load single SASA file for all dipeptide sequences
  SASA_data_dip <- read.table(
    paste0("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Dipeptide/SASA_score_M21/SASA_result_with_common_mon_ini.txt"),
    header = FALSE,
    col.names = c("Amino1", "Amino2", "ratio_mon_fin", "ratio_dim_fin")
  )
  SASA_data_dip$label <- paste(SASA_data_dip$Amino1, SASA_data_dip$Amino2, sep="_")
  SASA_data_dip$delta_AP <- SASA_data_dip$ratio_dim_fin - SASA_data_dip$ratio_mon_fin
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

# Process Edge data
if ( plotedge) {
  edge_data <- combined_df_net_di[combined_df_net_di$stat_name == "edges" & combined_df_net_di$frame == max(combined_df_net_di$frame), ]
  edge_means <- aggregate(value ~ seqname + state, data = edge_data, FUN = mean)
  edge_wide <- reshape(edge_means, 
                      idvar = "seqname", 
                      timevar = "state", 
                      direction = "wide")
  edge_wide$edge_diff <- edge_wide$value.dimer - edge_wide$value.monomer
  
  # Create data frame for edge heatmap
  edge_matrix <- matrix(NA, nrow = 2, ncol = length(residues),
                       dimnames = list(Position = c("1", "2"),
                                     Residue = residues))
  
  # Fill edge matrix
  for (pos in 1:2) {
    for (res in residues) {
      # Find sequences where the residue is at the specified position
      if (pos == 1) {
        seq <- paste0(res, "C")
      } else {
        seq <- paste0("C", res)
      }
      if (seq %in% edge_wide$seqname) {
        edge_matrix[pos, res] <- edge_wide$edge_diff[edge_wide$seqname == seq]
      }
    }
  }
  
  # Create Edge heatmap
  library(ggplot2)
  library(reshape2)
  
  edge_df <- melt(edge_matrix)
  colnames(edge_df) <- c("Position", "Residue", "DeltaEdge")
  
  # Convert Position to factor with specific labels
  edge_df$Position <- factor(edge_df$Position, 
                           levels = c("1", "2"),
                           labels = c("C at position 2", "C at position 1"))
  edge_df$Residue<-factor(edge_df$Residue,ordered_residues)
  edge_plt <- ggplot(edge_df, aes(x = Residue, y = Position, fill = DeltaEdge/300)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("blue", "grey", "yellow", "red"),
      # values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
      limits = c(min(edge_df$DeltaEdge/300), max(edge_df$DeltaEdge/300))) +
    plttheme +
    labs(title = "Change in Edges by Position and Residue in Dipeptides",
         x = "Residue",
         y = "Position of C",
         fill = "ΔEdges") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = group_colors[ordered_residues]),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(ggplotly(edge_plt))
  
  if (exists("save_plots") && save_plots) {
    ggsave(paste0(pltsavedir, "/Dipeptide_edges_by_position.png"), 
           plot = edge_plt,
           dpi = 1100, 
           width = 6.5, 
           height = 3,
           units = "cm")
  }
}

# Process AP data
if (plotap) {
  # Create data frame for AP heatmap
  ap_matrix <- matrix(NA, nrow = 2, ncol = length(residues),
                     dimnames = list(Position = c("1", "2"),
                                   Residue = residues))
  
  # Fill AP matrix
  for (pos in 1:2) {
    for (res in residues) {
      if (pos == 1) {
        seq_label <- paste0(res, "_C")
      } else {
        seq_label <- paste0("C_", res)
      }
      if (seq_label %in% SASA_data_dip$label) {
        ap_matrix[pos, res] <- SASA_data_dip$delta_AP[SASA_data_dip$label == seq_label]
      }
    }
  }
  
  # Create AP heatmap
  library(ggplot2)
  library(reshape2)
  
  ap_df <- melt(ap_matrix)
  colnames(ap_df) <- c("Position", "Residue", "DeltaAP")
  
  # Convert Position to factor with specific labels
  ap_df$Position <- factor(ap_df$Position, 
                         levels = c("1", "2"),
                         labels = c("C at position 2", "C at position 1"))
  ap_df$Residue<-factor(ap_df$Residue,ordered_residues)
  ap_plt <- ggplot(ap_df, aes(x = Residue, y = Position, fill = DeltaAP)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("blue", "grey", "yellow", "red"),
     # values = scales::rescale(c(-2, 0, 1, 2), from = c(-2, 2)),
        limits = c(min(ap_df$DeltaAP), max(ap_df$DeltaAP)),
      na.value = "black") +
    plttheme +
    labs(title = "Change in AP by Position and Residue in Dipeptides",
         x = "Residue",
         y = "Position of C",
         fill = "ΔAP") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = group_colors[ordered_residues]),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(ggplotly(ap_plt))
  
  if (exists("save_plots") && save_plots) {
    ggsave(paste0(pltsavedir, "/Dipeptide_AP_by_position.png"), 
           plot = ap_plt,
           dpi = 1100, 
           width = 24, 
           height = 12,
           units = "cm")
  }
} 
