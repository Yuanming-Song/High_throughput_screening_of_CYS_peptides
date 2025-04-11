if (readAP) {
  SASA_files <- list()
  for (pos in 1:3) {
    SASA_files[[pos]] <- read.table(
      paste0("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Tripeptide/SASA_score/SASA_result_with_common_mon_ini_C", pos, ".txt"),
      header = FALSE,
      col.names = c("Amino1", "Amino2", "Amino3", "ratio_mon_fin", "ratio_dim_fin")
    )
    SASA_files[[pos]]$label <- paste0(SASA_files[[pos]]$Amino1, SASA_files[[pos]]$Amino2, SASA_files[[pos]]$Amino3)
    SASA_files[[pos]]$Category <- mapply(classify_sequence, SASA_files[[pos]]$Amino1, SASA_files[[pos]]$Amino2, SASA_files[[pos]]$Amino3)
    SASA_files[[pos]]$label <- paste(SASA_files[[pos]]$Amino1, SASA_files[[pos]]$Amino2, SASA_files[[pos]]$Amino3, sep="_")
  }
}


get_delta_AP <- function(df) {
  df$delta_AP <- df$ratio_dim_fin - df$ratio_mon_fin
  return(df)
}

# Compute delta_AP for each position
SASA_files <- lapply(SASA_files, get_delta_AP)

residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")
mutation_effects <- matrix(NA, nrow = length(residues), ncol = length(residues),
                           dimnames = list(from = residues, to = residues))

for (from_res in residues) {
  for (to_res in residues) {
    if (from_res == to_res) next
    deltas <- c()
    for (pos in 1:3) {
      df <- SASA_files[[pos]]
      target_seqs <- df %>% filter(Amino1 == from_res | Amino2 == from_res | Amino3 == from_res)
      for (i in seq_len(nrow(target_seqs))) {
        orig_row <- target_seqs[i, ]
        for (j in 1:3) {
          if (orig_row[[paste0("Amino", j)]] == from_res) {
            mutated <- orig_row
            mutated[[paste0("Amino", j)]] <- to_res
            mutated_label <- paste(mutated$Amino1, mutated$Amino2, mutated$Amino3, sep = "_")
            mutated_row <- df %>% filter(label == mutated_label)
            if (nrow(mutated_row) == 1) {
              delta_orig <- orig_row$delta_AP
              delta_mut <- mutated_row$delta_AP
              deltas <- c(deltas, delta_mut - delta_orig)
            }
          }
        }
      }
    }
    if (length(deltas) > 0) {
      mutation_effects[from_res, to_res] <- mean(deltas)
    }
  }
}

# Define residue groups
acidic <- c("D", "E")
basic <- c("R", "K", "H")
polar <- c("S", "T", "N", "Q")
aromatic <- c("F", "W", "Y")
aliphatic <- c("A", "G", "I", "L", "M", "P", "V")

# Create ordered residue list
ordered_residues <- c(acidic, basic, polar,  aliphatic,aromatic)

# Reorder the mutation effects matrix
mutation_effects <- mutation_effects[ordered_residues, ordered_residues]

# Convert to data frame and filter for lower triangle
heatmap_df <- as.data.frame(as.table(mutation_effects))
colnames(heatmap_df) <- c("From", "To", "DeltaDeltaAP")

# Filter for lower triangle (where To comes before From in the ordered list)
#heatmap_df$From_idx <- match(heatmap_df$From, ordered_residues)
#heatmap_df$To_idx <- match(heatmap_df$To, ordered_residues)
#heatmap_df <- heatmap_df[heatmap_df$To_idx > heatmap_df$From_idx, ]
# Add a column to indicate if it's on the diagonal or symmetric
heatmap_df$From_idx <- match(heatmap_df$From, ordered_residues)
heatmap_df$To_idx <- match(heatmap_df$To, ordered_residues)
heatmap_df$is_diagonal <- heatmap_df$From == heatmap_df$To
heatmap_df$is_symmetric <- heatmap_df$From_idx > heatmap_df$To_idx
# Plotting
mut_ap_plt<-ggplot(heatmap_df, aes(x = From, y = To, fill = DeltaDeltaAP)) +
  geom_tile() +
  # Add diagonal line
  geom_abline(intercept = 20, slope = -1, color = "black", linetype = "dashed") +
  scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0, na.value = "white") +
  scale_x_discrete(limits = ordered_residues) +
  scale_y_discrete(limits = rev(ordered_residues)) +
  theme_minimal() +
  labs(title = "Average Change in ΔAP by Mutation", fill = "ΔΔAP") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Add text annotations for symmetric pairs
  geom_text(data = subset(heatmap_df, is_symmetric), 
            aes(label = "*"), 
            color = "black", 
            size = 3)

# Create interactive plot
print(ggplotly(mut_ap_plt))

if (save_plots) {
  ggsave(paste0(APplotDir, "/Tripeptide_mutation_AP.png"), 
         plot = mut_ap_plt,
         dpi = 1100, 
         width = 24, 
         height = 16,
         units = "cm")
}