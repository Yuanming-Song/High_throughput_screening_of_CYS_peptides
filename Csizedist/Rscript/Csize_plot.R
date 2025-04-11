seqlist<-c(#"L_V_C",
  "C_P_V","C_G_I","T_A_C","C_A_S"
  
)
# File paths and corresponding data names
files <- list(
  dimer = "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Csizedist/dimer_cdist_tripeptide.rda",
  monomer = "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Csizedist/monomer_cdist_tripeptide.rda"
)

# Define rebinning bins
monomer_rebins <- c(round(exp(seq(0,5,1))),300)
dimer_rebins <- unique(round(monomer_rebins/2))
dimer_rebins<-c(1,dimer_rebins[which(dimer_rebins>1)])
# Define colors for bins (with transparency)
bin_colors <- c(rgb(1, 0, 0, 0.2), rgb(0, 1, 0, 0.2), rgb(0, 0, 1, 0.2), 
                rgb(1, 1, 0, 0.2), rgb(1, 0, 1, 0.2), rgb(0, 1, 1, 0.2), 
                rgb(1, 0.5, 0, 0.2), rgb(0.5, 0, 0.5, 0.2), rgb(0, 0.5, 0, 0.2),
                rgb(0.5, 0, 0, 0.2), rgb(0, 0, 0.5, 0.2), rgb(0.5, 0.5, 0, 0.2))

# Define colors for points (monomer and dimer)
point_colors <- c("Monomer" = "blue", "Dimer" = "red")

pltxscale<- {scale_x_continuous(
  limits = log10(c(1, 300)),  # Set x-axis range in log base 10
  breaks = log10(c(1, 2, 10, 100, 300)),  # Set major breaks
  labels = c("1", "2","10", "100", "300")  # Label in numeric scale
)} 
pltlabs<- labs(x = "Cluster Size", y = "Frequency", color = "")

# Function to rebin data
rebin_data <- function(data, bins, is_dimer = FALSE) {
# Adjust bins for dimer data (multiply by 2)

  
  # Create empty result vector
  result <- numeric(length(bins) - 1)
  
  # For each bin
  for (i in 1:(length(bins) - 1)) {
    # Find indices of values within this bin
    idx <- which(data[, 1] >= bins[i] & data[, 1] < bins[i + 1])
    # Sum the frequencies in this bin
    result[i] <- sum(data[idx, 2])
  }
  
  # Normalize
  #result <- result / sum(result)
    if (is_dimer) {
    bins <- c(1,bins[-1] * 2)
  }
  # Calculate bin centers on the log scale
  bin_centers <- 10^((log10(bins[-length(bins)]) + log10(bins[-1])) / 2)


  
    return(data.frame(x = bin_centers, y = result))
}

# Load all SASA files at the beginning
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
if (readsizehis) {
  # Loop through files and generate plots
  for (state in names(files)) {
    # Load the file
    load(files[[state]])
    # Select the appropriate data based on the state
    data <- if (state == "dimer") dimer_sizehis else monomer_sizehis
    
    # Prepare data for plotting
    plot_data <- data.frame()
    x_values <- data[, 1]  # Use the first column for x-axis values
    
    for (i in 2:ncol(data)
       ) {
      y_values <- data[, i]  # Select the current column for y-axis
      temp_data <- data.frame(
        x = x_values,
        y = y_values,
        Sequence = colnames(data)[i],  # Add the column name as Sequence
        type = "original",
        state = ifelse(state == "dimer", "Dimer", "Monomer")  # Keep only necessary columns
      )
      plot_data <- rbind(plot_data, temp_data)
      
      # Add rebinned data
      seq_data <- data.frame(x = x_values, y = y_values)
      # Use appropriate bins based on state
      current_bins <- if (state == "dimer") dimer_rebins else monomer_rebins
      rebinned <- rebin_data(seq_data, current_bins, is_dimer = (state == "dimer"))
      rebinned$Sequence <- colnames(data)[i]
      rebinned$type <- "rebinned"
      rebinned$state <- ifelse(state == "dimer", "Dimer", "Monomer")  # Keep only necessary columns
      plot_data <- rbind(plot_data, rebinned)
    }
        
    plot_data$state <- ifelse(state == "dimer", "Dimer", "Monomer")
    assign(paste0(ifelse((state == "dimer"), "dimer_sizehis", "monomer_sizehis"), "_plot_data"), plot_data)
  }
}


for (seq in seqlist) {
  residues <- unlist(strsplit(seq, "_"))
  pos <- which(residues == "C")
  
  # Get AP values from pre-loaded SASA data
  APmon <- round(SASA_files[[pos]]$ratio_mon_fin[which(SASA_files[[pos]]$label==seq)], digits = 2)
  APdim <- round(SASA_files[[pos]]$ratio_dim_fin[which(SASA_files[[pos]]$label==seq)], digits = 2)
  
  # Create separate plots for monomer and dimer
  for (state in c("monomer", "dimer")) {
    if (state == "dimer") {
    data <-  dimer_sizehis_plot_data 
    plotdata <- data[which(data$Sequence == seq), ]
    plotdata[ plotdata$type=="original" & plotdata$x <= 150, ]$x<-plotdata[ plotdata$type=="original" & plotdata$x <= 150, ]$x*2
    
    } else {
      data <- monomer_sizehis_plot_data
      plotdata <- data[which(data$Sequence == seq), ]
      
    }
    # Determine ymin and ymax in log scale
    ymin <- min(log10(plotdata$y[which(plotdata$y != 0)]))
    ymax <- 0
    
    # Use appropriate bins for the state
    current_bins <- if (state == "dimer") dimer_rebins else monomer_rebins
    bin_starts <- current_bins[-length(current_bins)]
    bin_ends <- current_bins[-1]
    num_bins <- length(bin_starts)
    extended_bin_colors <- rep(bin_colors, length.out = num_bins)
    bin_rects <- data.frame(
      xmin = if (state == "dimer") log10(c(1,bin_starts[-1]*2)) else log10(bin_starts),
      xmax = log10(bin_ends*ifelse(state == "dimer",2,1)),
      ymin = ymin,
      ymax = ymax,
      fill_color = extended_bin_colors
    )
    
    pltyscale <- {
      scale_y_continuous(
        limits = c(ymin, ymax),
        breaks = log10(c(.001, .1, 1)),
        labels = c(".001", ".1", "1")
      )
    }
    
    plot_name <- paste0("plot_", seq, "_", state, "_separate")
    assign(plot_name, ggplot() +
             # Add background rectangles for the corresponding bins
             geom_rect(data = bin_rects, aes(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax), fill = bin_rects$fill_color, alpha = 0.2) +
             # Add points colored by type (original/rebinned)
             geom_point(data = plotdata, aes(x = log10(x), y = log10(y), color = type)) +
             plttheme +
             pltxscale +
             pltyscale +
             pltlabs +
             scale_color_manual(values = c("original" = "black", "rebinned" = "red")) +
             theme(legend.position = "right", legend.key = element_blank()) +
             labs(title = paste(seq, state, "APs:", paste(APmon, APdim, sep = "->")))
    )
  }

  # Create combined plot with only rebinned data
  plotdata_comb <- rbind(
    monomer_sizehis_plot_data[which(monomer_sizehis_plot_data$Sequence == seq & monomer_sizehis_plot_data$type == "rebinned"), ],
    dimer_sizehis_plot_data[which(dimer_sizehis_plot_data$Sequence == seq & dimer_sizehis_plot_data$type == "rebinned"), ]
  )
  
  # Determine ymin, ymax, and midpoint in log scale for combined plot
  ymin_comb <- min(log10(plotdata_comb$y[which(plotdata_comb$y != 0)]))
  ymax_comb <- 0
  midpoint_comb <- (ymin_comb + ymax_comb) / 2
  
  # Use monomer_rebins for monomer background rectangles
  bin_starts_monomer <- monomer_rebins[-length(monomer_rebins)]
  bin_ends_monomer <- monomer_rebins[-1]
  num_bins_monomer <- length(bin_starts_monomer)
  extended_bin_colors_monomer <- rep(bin_colors, length.out = num_bins_monomer)
  bin_rects_monomer <- data.frame(
    xmin = log10(bin_starts_monomer),
    xmax = log10(bin_ends_monomer),
    ymin = midpoint_comb,  # Top half of the plot
    ymax = ymax_comb,
    fill_color = extended_bin_colors_monomer
  )
  
  # Use dimer_rebins for dimer background rectangles
  bin_starts_dimer <- dimer_rebins[-length(dimer_rebins)]
  bin_ends_dimer <- dimer_rebins[-1]
  num_bins_dimer <- length(bin_starts_dimer)
  extended_bin_colors_dimer <- rep(bin_colors, length.out = num_bins_dimer)
  bin_rects_dimer <- data.frame(
    xmin = log10(c(1,bin_starts_dimer[-1]*2)),
    xmax = log10(bin_ends_dimer*2),
    ymin = ymin_comb,  # Bottom half of the plot
    ymax = midpoint_comb,
    fill_color = extended_bin_colors_dimer
  )
  
  pltyscale_comb <- {
    scale_y_continuous(
      limits = c(ymin_comb, ymax_comb),
      breaks = log10(c(.001, .1, 1)),
      labels = c(".001", ".1", "1")
    )
  }
  
  plot_name <- paste0("plot_", seq, "_combined_rebinned")
  assign(plot_name, ggplot() +
           # Add background rectangles for monomer bins
           geom_rect(data = bin_rects_monomer, aes(xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax), fill = bin_rects_monomer$fill_color, alpha = 0.2) +
           # Add background rectangles for dimer bins
           geom_rect(data = bin_rects_dimer, aes(xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax), fill = bin_rects_dimer$fill_color, alpha = 0.2) +
           # Add points colored by state (Monomer/Dimer)
           geom_point(data = plotdata_comb, aes(x = log10(x), y = log10(y), color = state)) +
           plttheme +
           pltxscale +
           pltyscale_comb +
           pltlabs +
           scale_color_manual(values = point_colors) +
           theme(legend.position = "right", legend.key = element_blank()) +
           labs(title = paste(seq, "Combined (Rebinned)", "APs:", paste(APmon, APdim, sep = "->")))
  )
  
  # Create original combined plot with both original data
  plotdata_orig <- rbind(
    monomer_sizehis_plot_data[which(monomer_sizehis_plot_data$Sequence == seq & monomer_sizehis_plot_data$type == "original"), ],
    dimer_sizehis_plot_data[which(dimer_sizehis_plot_data$Sequence == seq & dimer_sizehis_plot_data$type == "original"), ]
  )
  plotdata_orig[ plotdata_orig$type=="original" & plotdata_orig$x <= 150 & plotdata_orig$state =="Dimer", ]$x <- plotdata_orig[ plotdata_orig$type=="original" & plotdata_orig$x <= 150 & plotdata_orig$state =="Dimer", ]$x*2
  
  # Determine ymin, ymax, and midpoint in log scale for combined plot
  ymin_orig <- min(log10(plotdata_orig$y[which(plotdata_orig$y != 0)]))
  ymax_orig <- 0
  midpoint_orig <- (ymin_orig + ymax_orig) / 2
  
  # Use monomer_rebins for monomer background rectangles
  bin_starts_monomer <- monomer_rebins[-length(monomer_rebins)]
  bin_ends_monomer <- monomer_rebins[-1]
  num_bins_monomer <- length(bin_starts_monomer)
  extended_bin_colors_monomer <- rep(bin_colors, length.out = num_bins_monomer)
  bin_rects_monomer <- data.frame(
    xmin = log10(bin_starts_monomer),
    xmax = log10(bin_ends_monomer),
    ymin = midpoint_orig,  # Top half of the plot
    ymax = ymax_orig,
    fill_color = extended_bin_colors_monomer
  )
  
  # Use dimer_rebins for dimer background rectangles
  bin_starts_dimer <- dimer_rebins[-length(dimer_rebins)]
  bin_ends_dimer <- dimer_rebins[-1]
  num_bins_dimer <- length(bin_starts_dimer)
  extended_bin_colors_dimer <- rep(bin_colors, length.out = num_bins_dimer)
  bin_rects_dimer <- data.frame(
    xmin = log10(c(1,bin_starts_dimer[-1]*2)),
    xmax = log10(bin_ends_dimer*2),
    ymin = ymin_orig,  # Bottom half of the plot
    ymax = midpoint_orig,
    fill_color = extended_bin_colors_dimer
  )
  
  pltyscale_orig <- {
    scale_y_continuous(
      limits = c(ymin_orig, ymax_orig),
      breaks = log10(c(.001, .1, 1)),
      labels = c(".001", ".1", "1")
    )
  }
  
  plot_name <- paste0("plot_", seq, "_combined_original")
  assign(plot_name, ggplot() +
           # Add background rectangles for monomer bins
           geom_rect(data = bin_rects_monomer, aes(xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax), fill = bin_rects_monomer$fill_color, alpha = 0.2) +
           # Add background rectangles for dimer bins
           geom_rect(data = bin_rects_dimer, aes(xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax), fill = bin_rects_dimer$fill_color, alpha = 0.2) +
           # Add points colored by state (Monomer/Dimer)
           geom_point(data = plotdata_orig, aes(x = log10(x), y = log10(y), color = state)) +
           plttheme +
           pltxscale +
           pltyscale_orig +
           pltlabs +
           scale_color_manual(values = point_colors) +
           theme(legend.position = "right", legend.key = element_blank()) +
           labs(title = paste(seq, "Combined (Original)", "APs:", paste(APmon, APdim, sep = "->")))
  )
  
  # Print all plots
  #print(get(paste0("plot_", seq, "_monomer_separate")))
  #print(get(paste0("plot_", seq, "_dimer_separate")))
  #print(get(paste0("plot_", seq, "_combined_original")))
  #print(get(paste0("plot_", seq, "_combined_rebinned")))
}

# New section for combined plots (after line 313)
# Set save option
save_plots <- TRUE
plot_dir <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Tripeptide/plots/"

# Remove legends and adjust axis labels for all plots
base_theme <- theme(legend.position = "none",
                   
                   axis.title.y = element_blank())

# Create four separate combined plots (one for each type)
plot_types <- list(
  list(name = "monomer", suffix = "_monomer_separate"),
  list(name = "dimer", suffix = "_dimer_separate"),
  list(name = "combined_original", suffix = "_combined_original"),
  list(name = "combined_rebinned", suffix = "_combined_rebinned")
)

for (plot_type in plot_types) {
  # Get all plots of this type for all sequences
  plots <- lapply(seqlist, function(seq) {
    plot <- get(paste0("plot_", seq, plot_type$suffix))
    # Apply base theme
    plot + base_theme
  })
  
  # Remove x-axis labels for first two plots
  for (i in 1:length(plots)) {
    if (i <= 2) {
      plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    }
  }
  
  # Get legend from the first plot
  legend <- get_legend(plots[[1]] + theme(legend.position = "right", legend.direction = "vertical"))
  
  # Combine the plots in a 2x2 grid
  combined <- plot_grid(plotlist = plots, ncol = 2, align = "v")
  
  # Add common y-axis title
  combined_final <- ggdraw() +
    draw_plot(combined, x = 0.02, y = 0, width = 0.98, height = 1) +
    draw_label("Relative Frequency", x = 0.02, y = 0.5, angle = 90, size = 12)
  
  
  # Combine with legend
  combined_final <- plot_grid(combined_final, legend, ncol = 2, rel_widths = c(0.87, 0.13))
  
  # Print the combined plot
  print(combined_final)
  
  # Save plot if save_plots is TRUE
  if (save_plots) {
    ggsave(paste0(APplotDir, "Tripeptide_csize",plot_type$name, "_combined.png"), 
           plot = combined_final,
           dpi = 1100,
           width = 24,
           height = 16,
           units = "cm")
  }
}
