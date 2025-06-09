# Code Block 2: Load data and create plots
load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/csize_SIRAH/SIRAH_validation_boxsize_csize_t.rda")


avogadro <- 6.02214076e23  # Avogadro's number
sum_size_t$box<-as.numeric(sum_size_t$box)
# Get unique box sizes and calculate concentrations
unique_boxes <- sort(unique(sum_size_t$box))
concentrations <- round(300 * 1000 / (avogadro * 1e-24 * (unique_boxes^3 )))

# Create mapping between box sizes and concentration labels
box_to_conc <- setNames(concentrations, factor(unique_boxes))

# Add concentration column to the data frame
sum_size_t$concentration <- round(300 * 1000 / (avogadro * 1e-24 * (sum_size_t$box^3 )))

# Define colors for each concentration using a more diverse palette
conc_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
    "#A65628", "#F781BF", "#999999", "#66C2A5",
    "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"),
  factor(concentrations)
)

# Define colors for each box size using a more diverse palette
box_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
    "#A65628", "#F781BF", "#999999", "#66C2A5",
    "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"),
  factor(unique_boxes)
)

# Check for missing box sizes in each sequence
for (s in unique(sum_size_t$seq)) {
  seq_boxes <- unique(subset(sum_size_t, seq == s)$box)
  missing_boxes <- setdiff(unique_boxes, seq_boxes)
  if (length(missing_boxes) > 0) {
    cat(sprintf("Sequence %s is missing box sizes: %s\n", 
                s, 
                paste(missing_boxes, collapse = ", ")))
  } else {
    goodseq<-s
  }
}

plot_list1 <- list()
plot_list2 <- list()

for (s in unique_seqs) {
  df_seq <- subset(sum_size_t, seq == s)
  if (concentration) {
    p1 <- ggplot(df_seq, aes(x = V1, y = V2, 
                             color =factor(concentration)))+
      scale_color_manual(values = conc_colors)
  }else {
    p1<-ggplot(df_seq, aes(x = V1, y = V2, 
                           color =factor(box)))+
      scale_color_manual(values = box_colors)
  }
  p1<-p1+
    geom_line() +
    labs(title = paste0("n-", gsub("_", "", s), "-c"),
         x = "",
         y = "",
         color = if(concentration) "Concentration (mM)" else "Box Size (nm)") +
    plttheme +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")   # Remove legend from individual plots
  
  
  # Add log scale for y-axis in p2
  ymin <- min(log10(df_seq$V3[which(df_seq$V3 != 0)]))
  ymax <- max(log10(df_seq$V3))
  
  if (concentration) {
    p2 <- ggplot(df_seq, aes(x = V1, y = V3, 
                             color = factor(concentration))) +
      scale_color_manual(values = conc_colors)
  } else {
    p2 <- ggplot(df_seq, aes(x = V1, y = V3, 
                             color = factor(box))) +
      scale_color_manual(values = box_colors)
  }
  p2 <- p2 +
    geom_line() +
    labs(title = paste0("n-", gsub("_", "", s), "-c"),
         x = "",
         y = "",
         color = if(concentration) "Concentration (mM)" else "Box Size (nm)") +
    plttheme +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +  # Remove legend from individual plots
    scale_y_continuous(
      trans = "log10",
      limits = c(10^ymin, 10^ymax),
      breaks = c(1, 10, 50, 100, 300),
      labels = c("1", "10", "50", "100", "300")
    )
  
  plot_list1[[s]] <- p1
  plot_list2[[s]] <- p2
}

# Create a separate plot just for the legend with 3x3 layout
legend_plot <- plot_list1[[goodseq]] +
  plttheme +
  theme(legend.direction = "vertical")+
  guides(color = guide_legend(nrow = 3, ncol = 3)) +
  labs(col = if(concentration) "Concentration (mM)" else "Box Size (nm)")

# Extract legend using ggplotGrob
legend_grob <- ggplotGrob(legend_plot)
legend <- legend_grob$grobs[which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")][[1]]



# Combine plots from plot_list1 with legend in upper right
combined_plot1 <- plot_grid(
  plot_list1[[1]], legend, plot_list1[[2]], 
  plot_list1[[3]], 
  ncol = 2, 
  align = "v",
  axis = "l"  # Add axis parameter for vertical alignment
)

# Add common labels
combined_plot1 <- ggdraw() +
  draw_plot(combined_plot1, x = 0.02, y = 0.02, width = 0.98, height = 1) +
  draw_label("Largest Cluster Size", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Time (ns)", x = 0.5, y = 0.02, size = 12)

# Combine plots from plot_list2 with legend in upper right
combined_plot2 <- plot_grid(
  plot_list2[[1]], legend, plot_list2[[2]], 
  plot_list2[[3]], 
  ncol = 2, 
  align = "v",
  axis = "l"  # Add axis parameter for vertical alignment
)

# Add common labels
combined_plot2 <- ggdraw() +
  draw_plot(combined_plot2, x = 0.02, y = 0, width = 0.98, height = 1) +
  draw_label("Number of Clusters", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Time (ns)", x = 0.5, y = 0.02, size = 12)

# Explicitly print the plots
print(combined_plot1)
print(combined_plot2)

# Save plots if save_plots is TRUE
if (save_plots) {
  ggsave(paste0(APplotDir,"SIRAH_plt/SIRAH_largest_cluster_size_", 
                if(concentration) "concentration" else "boxsize", ".png"), 
         plot = combined_plot1,
         dpi = 1100, 
         width = 24, 
         height = 16,
         units = "cm")
  
  ggsave(paste0(APplotDir,"SIRAH_plt/SIRAH_number_of_clusters_", 
                if(concentration) "concentration" else "boxsize", ".png"), 
         plot = combined_plot2,
         dpi = 1100, 
         width = 24, 
         height = 16,
         units = "cm")
}
