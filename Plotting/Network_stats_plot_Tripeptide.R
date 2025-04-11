load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/tripeptide_network_stats.rda")

# Assume final_results is a list with two elements:
#   final_results[[1]]: monomer_results 
#   final_results[[2]]: dimer_results
# Each has columns: frame, value, stat_name, seqname, ratio, state.
# Also assume seqlist is defined, e.g.:
seqlist <- c("C_P_V",
  "Q_R_C",
              "C_G_I", "T_A_C"#, "C_A_S"
)
# Remove underscores from seqlist for matching:
seqlist_clean <- gsub("_", "", seqlist)

# Combine monomer and dimer data
combined_df <- rbind(final_results[[1]], final_results[[2]])
# Create a cleaned sequence column (remove underscores) for matching
#combined_df <- combined_df %>% mutate(seqname_clean = gsub("_", "", seqname))

# Initialize list to store plots
plots_list <- list()

# Loop over each sequence in seqlist_clean
for (seq in seqlist_clean) {
  # Subset the data for the current sequence
  seq_df <- combined_df %>% filter(seqname == seq)
  # Get unique ratio values for monomer and dimer
  ratio_mono <- unique(seq_df$ratio[seq_df$state == "monomer"])
  ratio_dimer <- unique(seq_df$ratio[seq_df$state == "dimer"])
  
  plots_list[[seq]]<-list()
  # For each unique stat_name in this subset, generate a plot
  for (stat in unique(seq_df$stat_name)) {
    stat_df <- seq_df %>% filter(stat_name == stat)
    
    
    # Build title: "seq: ratio for monomer -> ratio for dimer, stat_name"
    title_str <- paste(seq, ": ", ratio_mono, " -> ", ratio_dimer, ", ", stat, sep = "")
    
    p <- ggplot(stat_df, aes(x = frame, y = value, color = state)) +
      geom_line() +
      ggtitle(title_str) +
      plttheme
    # Use a unique key combining seq and stat for storage
    key <- paste(seq, stat, sep = "_")
    plots_list[[seq]][[stat]] <- p
  }
}

# Print each plot (you can later combine them as needed)
for (seq in seqlist_clean) {
  for (stat in unique(seq_df$stat_name)) {
    print(ggplotly(plots_list[[seq]][[stat]] ))
  }
}

for (state in unique(combined_df$state))
for (stat in unique(combined_df$stat_name)) {
  stat_df <- combined_df[which(combined_df$stat_name==stat & combined_df$state== state),]

  stat_df$ratio[which(stat_df$ratio==100)]<-4
  # Build title: "seq: ratio for monomer -> ratio for dimer, stat_name"
  title_str <- paste(stat,state)
  
  #p <- 
    print(ggplot(stat_df, aes(x = frame, y = value, type=seqname, col=ratio)) +
    geom_line() +
    ggtitle(title_str) +
      scale_color_gradientn(colors = c("blue", "yellow", "red")) +
    plttheme)
  #plots_list_2[[state]][[stat]] <- p
}
statpltlist<-list()
# Remove legends and adjust axis labels for all plots
base_theme <- theme(legend.position = "none",
                    axis.title.x = element_blank()
                    #axis.title.y = element_blank()
)
for (state in unique(combined_df$state)) {
  
  statpltlist[[state]] <-list()
  for (stat in unique(combined_df$stat_name)) {
    stat_df <- combined_df[which(combined_df$stat_name==stat & combined_df$state== state & combined_df$frame== 202),]
    
    stat_df$ratio[which(stat_df$ratio==100)]<-4
    # Build title: "seq: ratio for monomer -> ratio for dimer, stat_name"
    title_str <- paste(stat,state)
    
    p <- ggplot(stat_df, aes(x = ratio, y = value, col=seqname)) +
      geom_point() +
      #ggtitle(title_str) +
      #scale_color_gradientn(colors = c("blue", "yellow", "red")) +
      labs(y=stat)+
      plttheme
    print(ggplotly(p+       labs(x="APs")))
    statpltlist[[state]][[stat]] <- p+base_theme
  }
  
  # Combine the plots in a 2x2 grid
  combined <- plot_grid(plotlist = statpltlist[[state]], ncol = 1, align = "v")
  
  # Add common y-axis title
  combined_final <- ggdraw() +
    draw_plot(combined, x = 0.02, y = 0, width = 0.98, height = 1) +
    draw_label("APs", x = 0.5, y = 0.02, angle = 0, size = 12)
  
  
  
  # Print the combined plot
  print(combined_final)
  
  # Save plot if save_plots is TRUE
  if (save_plots) {
    ggsave(paste0(APplotDir, "Tripeptide_sufficient_stat",state, "_combined.png"), 
           plot = combined_final,
           dpi = 1100,
           width = 24,
           height = 16,
           units = "cm")
  }
}

# Combine plots by statistics
stat_plots <- list()
all_stats <- unique(combined_df$stat_name)

# Process stats in groups of 3
for (i in seq(1, length(all_stats), by = 3)) {
  # Get current group of stats
  current_stats <- all_stats[i:min(i+2, length(all_stats))]
  
  # Create a list to store combined plots for this group
  group_plots <- list()
  
  for (stat in current_stats) {
    # Create a list to store all sequence plots for this stat
    seq_plots <- list()
    
    for (seq in seqlist_clean) {
      # Get the plot for this sequence and stat
      p <- plots_list[[seq]][[stat]]
      
      # Remove title and adjust theme
      p <- p + 
        theme(legend.position = "none",
              axis.title.y = element_blank(),
              axis.title.x = element_blank()) 
      
      seq_plots[[seq]] <- p
    }
    
    # Combine all sequence plots for this stat
    combined_stat <- plot_grid(plotlist = seq_plots, ncol = 1, align = "v")
    
    # Add stat label
    combined_stat_final <- ggdraw() +
      draw_plot(combined_stat, x = 0.02, y = 0, width = 0.98, height = 1) +
      draw_label(stat, x = 0.02, y = 0.5, angle = 90, size = 12)
    
    group_plots[[stat]] <- combined_stat_final
  }
  
  # Combine the group of stats horizontally
  combined_group <- plot_grid(plotlist = group_plots, ncol = length(current_stats), align = "h")
  
  # Print the combined plot
  print(combined_group)
  
  # Save plot if save_plots is TRUE
  if (save_plots) {
    group_name <- paste(current_stats, collapse = "_")
    ggsave(paste0(APplotDir, "/Tripeptide_stats_", group_name, "_combined.png"), 
           plot = combined_group,
           dpi = 1100,
           width = 24 * length(current_stats),  # Adjust width based on number of stats
           height = 16,
           units = "cm")
  }
}
