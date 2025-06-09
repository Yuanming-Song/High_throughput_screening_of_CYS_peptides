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

if (0) {
  # Create a cleaned sequence column (remove underscores) for matching
  #combined_df_ddedge_tripep <- combined_df_ddedge_tripep %>% mutate(seqname_clean = gsub("_", "", seqname))
  # Assume final_results is a list with two elements:
  #   final_results[[1]]: monomer_results 
  #   final_results[[2]]: dimer_results
  # Each has columns: frame, value, stat_name, seqname, ratio, state.
  # Also assume seqlist is defined, e.g.:

  for (state in unique(combined_df_ddedge_tripep$state))
    for (stat in unique(combined_df_ddedge_tripep$stat_name)) {
      stat_df <- combined_df_ddedge_tripep[which(combined_df_ddedge_tripep$stat_name==stat & combined_df_ddedge_tripep$state== state),]
      
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
}
statpltlist<-list()
# Remove legends and adjust axis labels for all plots
base_theme <- theme(legend.position = "none",
                    axis.title.x = element_blank()
                    #axis.title.y = element_blank()
)

# Initialize list to store all plots
all_plots <- list()

for (state in unique(combined_df_ddedge_tripep$state)) {
  
  statpltlist[[state]] <-list()
  for (stat in unique(combined_df_ddedge_tripep$stat_name)) {
    stat_df <- combined_df_ddedge_tripep[which(combined_df_ddedge_tripep$stat_name==stat & combined_df_ddedge_tripep$state== state & combined_df_ddedge_tripep$frame== 202),]
    
    stat_df$ratio[which(stat_df$ratio==100)]<-4
    # Build title: "seq: ratio for monomer -> ratio for dimer, stat_name"
    title_str <- paste(stat,state)
    
    p <- ggplot(stat_df, aes(x = ratio, y = value, col=seqname)) +
      geom_point() +
      #ggtitle(title_str) +
      #scale_color_gradientn(colors = c("blue", "yellow", "red")) +
      labs(y=stat)+
      plttheme
    # Store the interactive plot
    all_plots[[paste(state, stat, "interactive", sep="_")]] <- ggplotly(p + labs(x="APs"))
    statpltlist[[state]][[stat]] <- p+base_theme
  }
  
  # Combine the plots in a 2x2 grid
  combined <- plot_grid(plotlist = statpltlist[[state]], ncol = 1, align = "v")
  
  # Add common y-axis title
  combined_final <- ggdraw() +
    draw_plot(combined, x = 0.02, y = 0, width = 0.98, height = 1) +
    draw_label("APs", x = 0.5, y = 0.02, angle = 0, size = 12)
  
  # Store the combined plot
  all_plots[[paste(state, "combined", sep="_")]] <- combined_final
  
  # Save plot if save_plots is TRUE
  if (save_plots) {
    ggsave(paste0(APplotDir, "Tripeptide_sufficient_stat_",state, "_combined.png"), 
           plot = combined_final,
           dpi = 1100,
           width = 6.5/2,
           height = 8,
           units = "in")
  }
}

# Print all stored plots after the loop is complete
for (plot_name in names(all_plots)) {
  print(all_plots[[plot_name]])
}

# Only proceed with individual sequence plotting if plot_individual_seqs is TRUE
if (exists("plot_individual_seqs") && plot_individual_seqs) {
  seqlist <- c("C_P_V",
               #"Q_R_C",
               "C_G_I", "T_A_C", "C_A_S"
  )
  # Remove underscores from seqlist for matching:
  seqlist_clean <- gsub("_", "", seqlist)

  # First calculate y-axis ranges for each statistic
  stat_ranges <- list()
  for (stat in unique(combined_df_ddedge_tripep$stat_name)) {
    stat_data <- combined_df_ddedge_tripep %>% 
      filter(seqname %in% seqlist_clean, stat_name == stat)
    stat_ranges[[stat]] <- c(min(stat_data$value), max(stat_data$value))
  }

  # Initialize list to store plots
  plots_list <- list()

  # Loop over each sequence in seqlist_clean
  for (seq in seqlist_clean) {
    # Subset the data for the current sequence
    seq_df <- combined_df_ddedge_tripep %>% filter(seqname == seq)
    # Get unique ratio values for monomer and dimer
    ratio_mono <- unique(seq_df$ratio[seq_df$state == "monomer"])
    ratio_dimer <- unique(seq_df$ratio[seq_df$state == "dimer"])
    
    plots_list[[seq]]<-list()
    # For each unique stat_name in this subset, generate a plot
    for (stat in unique(seq_df$stat_name)) {
      stat_df <- seq_df %>% filter(stat_name == stat)
      
      # Build title: "seq: ratio for monomer -> ratio for dimer, stat_name"
      title_str <- paste0(seq, ":\n", ratio_mono, " -> \n", ratio_dimer)
      
      p <- ggplot(stat_df, aes(x = frame, y = value, color = state)) +
        geom_line() +
        labs(y=stat)+
        base_theme+
        plttheme+
        theme(legend.position = "none")+
        scale_x_continuous(breaks = pretty_breaks(n = 3)) +
        scale_y_continuous(limits = stat_ranges[[stat]], breaks = pretty_breaks(n = 4))
      
      # Use a unique key combining seq and stat for storage
      key <- paste(seq, stat, sep = "_")
      if (stat ==unique(seq_df$stat_name)[1]) {
        p<-p+ggtitle(title_str) 
      } 
      if (seq !=seqlist_clean[1]) {
        p<-p+
          theme(axis.title.y = element_blank())
      }
      plots_list[[seq]][[stat]] <- p
    }
  }

  combseq<-list()
  for (seq in seqlist_clean) {
    # Combine the plots in a 2x2 grid
    combseq[[seq]] <- plot_grid(plotlist = plots_list[[seq]], ncol = 1, align = "v")
  }
  seqcombined <- plot_grid(plotlist = combseq, ncol = 4, align = "vh")
  comb_seq_legend<-get_legend(plots_list[[1]][[1]]+
                                theme(legend.position = "right",
                                      legend.direction = "horizontal"))


  # Add common y-axis title
  combined_final_seq <- ggdraw() +
    draw_plot(seqcombined, x = 0.02, y = 0.02, width = 0.98, height = 0.98) +
    draw_label("Time (ns)", x = 0.5, y = 0.02, angle = 0, size = 12)
  combined_final_seq <- plot_grid(plotlist = list(combined_final_seq,comb_seq_legend), nrow = 2,rel_heights=c(0.95,0.05))
  
  # Store the final sequence plot
  all_plots[["sequence_combined"]] <- combined_final_seq
  
  # Save the plot
  ggsave(paste0(APplotDir, "Tripeptide_sufficient_stat_exp_seq.png"), 
         plot = combined_final_seq,
         dpi = 1100,
         width = 6.5,
         height = 8,
         units = "in")
}
if (0) {
# Combine plots by statistics
stat_plots <- list()
all_stats <- unique(combined_df_ddedge_tripep$stat_name)

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
}