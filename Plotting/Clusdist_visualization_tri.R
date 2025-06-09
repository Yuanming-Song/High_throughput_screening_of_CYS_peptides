# File paths and corresponding data names
files <- list(
  dimer = "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Csizedist/dimer_cdist_tripeptide_single_node.rda",
  monomer = "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Csizedist/monomer_cdist_tripeptide.rda"
)

pltxscale <- {scale_x_continuous(
  limits = log10(c(1, 300)),  # Set x-axis range in log base 10
  breaks = log10(c(1, 2, 10, 100, 300)),  # Set major breaks
  labels = c("1", "2","10", "100", "300")  # Label in numeric scale
)} 

# Remove axis labels from individual plots
pltlabs <- labs(x = NULL, y = NULL, color = NULL)
if (loaddata) {
# Loop through files and generate plots
for (state in names(files)) {
  # Load the file
  load(files[[state]])
  # Select the appropriate data based on the state
  data <- if (state == "dimer") dimer_sizehis else monomer_sizehis
  
  # Prepare data for plotting
  plot_data <- data.frame()
  x_values <- data[, 1]  # Use the first column for x-axis values
  
  for (i in 2:ncol(data)) {
    y_values <- data[, i]  # Select the current column for y-axis
    temp_data <- data.frame(
      x = x_values,
      y = y_values,
      Sequence = colnames(data)[i]  # Add the column name as Sequence
    )
    plot_data <- rbind(plot_data, temp_data)
  }
  plot_data$state <- ifelse(state == "dimer", "Dimer", "Monomer")
  assign(paste0(ifelse((state == "dimer"), "dimer_sizehis", "monomer_sizehis"), "_plot_data"), plot_data)
}
}
if (!exists("seqlist")) {
  seqlist <- c("C_P_V", "C_G_I", "T_A_C", "C_A_S")
}
plotlist <- list()
for (seq in seqlist) {
  residues <- unlist(strsplit(seq, "_"))
  pos <- which(residues == "C")
  
  # Read SASA data
  SASA_common_ini <- read.table(
    paste0("/Users/song/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Tripeptide/SASA_score/SASA_result_with_common_mon_ini_C", pos, ".txt"), 
    header = FALSE, 
    col.names = c("Amino1", "Amino2", "Amino3", "ratio_mon_fin", "ratio_dim_fin")
  )
  
  # Create label and get AP values
  SASA_common_ini$label <- paste(SASA_common_ini$Amino1, SASA_common_ini$Amino2, SASA_common_ini$Amino3, sep="_")
  APmon <- round(SASA_common_ini[which(SASA_common_ini$label==seq), 4], digits = 2)
  APdim <- round(SASA_common_ini[which(SASA_common_ini$label==seq), 5], digits = 2)
  
  # Combine plot data
  plotdata_comb <- c()
  for (state in names(files)) {
    data <- if (state == "dimer") dimer_sizehis_plot_data else monomer_sizehis_plot_data
    plotdata_comb <- rbind(plotdata_comb, data[which(data$Sequence==seq),])
  }
  
  pltyscale <- {
    scale_y_continuous(
      limits = c(min(log10(plotdata_comb$y[which(plotdata_comb$y!=0)])), 0),
      breaks = log10(c(.001, .1, 1)),
      labels = c(".001", ".1", "1")
    ) 
  }
  
  # Create and display plot with no axis titles and no legend
  clusdist_plot <- ggplot(plotdata_comb, aes(x = log10(x), y = log10(y), color = state)) +
    geom_point() +
    plttheme +
    pltxscale +
    pltyscale +
    pltlabs +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(title = paste(seq, "APs:", paste(APmon, APdim, sep="->")))
  
  print(ggplotly(clusdist_plot))
  plotlist[[seq]] <- clusdist_plot
} 

# Create combined plot
seqcombined <- plot_grid(plotlist = plotlist, ncol = 2, nrow=2, align = "hv")

# Extract legend from one plot (temporarily add legend back to get it)
legend_plot <- plotlist[[1]] + 
  theme(legend.position = "right", 
        legend.direction = "horizontal",
        legend.title = element_blank())
comb_seq_legend <- get_legend(legend_plot)

# Create the final combined plot with shared axis labels
# First create plot with y-axis label
y_axis_plot <- ggdraw() +
  draw_plot(seqcombined,x = 0.02,y=0.02,height = 0.98) +
  draw_label("Frequency", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Cluster Size", x = 0.5, y = 0.02, size = 12)

# Then add x-axis label and legend
combined_final_seq <- plot_grid(
  y_axis_plot,
  comb_seq_legend,
  nrow = 2,
  rel_heights = c(0.90, 0.05)
) 
print(combined_final_seq)

# Save the combined plot
if (save_plots) {
  
  
  ggsave(paste0(pltsavedir, "Tripeptide_clusdist_exp_seq.png"), 
         plot = combined_final_seq,
         dpi = 1100,
         width = 6.5,
         height = 3,
         units = "in")
  
  # Save individual plots if requested
  if (0) {
    if (!exists("pltsavedir")) {
      pltsavedir <- "~/Documents/Research/manuscript/ML_MD_pep/ML_MD Peptide/MARTINI_plt/"
    }
    for (seq in seqlist) {
      ggsave(
        paste0(pltsavedir, "/Tripeptide_clusdist_", seq, ".png"),
        plot = plotlist[[seq]],
        dpi = 1100,
        width = 6.5,
        height = 3,
        units = "in"
      )
    }
  }
}