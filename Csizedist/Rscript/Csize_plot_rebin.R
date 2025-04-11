# File paths and corresponding data names
files <- list(
  dimer = "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Csizedist/dimer_cdist_tripeptide.rda",
  monomer = "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Csizedist/monomer_cdist_tripeptide.rda"
)

# Define rebinning bins
rebins <- c(seq(1,10,1), seq(30,100,20), 150, 300)

pltxscale <- {scale_x_continuous(
  limits = log10(c(1, 300)),  # Set x-axis range in log base 10
  breaks = log10(c(1, 2, 10, 100, 300)),  # Set major breaks
  labels = c("1", "2","10", "100", "300")  # Label in numeric scale
)} 

pltlabs <- labs(x = "Cluster Size", y = "Frequency", color = "")

# Function to rebin data
rebin_data <- function(data, bins) {
  # Create empty result vector
  result <- numeric(length(bins))
  
  # For each bin
  for (i in 1:(length(bins)-1)) {
    # Find indices of values within this bin
    idx <- which(data[,1] >= bins[i] & data[,1] < bins[i+1])
    # Sum the frequencies in this bin
    result[i] <- sum(data[idx,2])
  }
  
  # For the last bin
  idx <- which(data[,1] >= bins[length(bins)-1] & data[,1] <= bins[length(bins)])
  result[length(bins)-1] <- sum(data[idx,2])
  
  # Normalize
  result <- result/sum(result)
  
  return(data.frame(x = bins[-length(bins)], y = result[-length(result)]))
}

# Load the data
for (state in names(files)) {
  load(files[[state]])
  data <- if (state == "dimer") dimer_sizehis else monomer_sizehis
  
  # Prepare original data for plotting
  plot_data <- data.frame()
  x_values <- data[, 1]
  
  for (i in 2:ncol(data)) {
    y_values <- data[, i]
    temp_data <- data.frame(
      x = x_values,
      y = y_values,
      Sequence = colnames(data)[i],
      type = "original"
    )
    plot_data <- rbind(plot_data, temp_data)
  }
  
  # Rebin the data
  for (i in 2:ncol(data)) {
    seq_data <- data.frame(x = data[,1], y = data[,i])
    rebinned <- rebin_data(seq_data, rebins)
    rebinned$Sequence <- colnames(data)[i]
    rebinned$type <- "rebinned"
    plot_data <- rbind(plot_data, rebinned)
  }
  
  plot_data$x <- if (state == "dimer") plot_data$x*2 else plot_data$x
  plot_data$state <- ifelse(state == "dimer", "Dimer", "Monomer")
  assign(paste0(ifelse((state == "dimer"), "dimer_sizehis", "monomer_sizehis"), "_plot_data"), plot_data)
}

# Select one sequence to plot
seq <- "C_P_V"  # You can change this to any sequence you want

# Create separate plots for monomer and dimer
for (state in c("monomer", "dimer")) {
  data <- if (state == "dimer") dimer_sizehis_plot_data else monomer_sizehis_plot_data
  plotdata <- data[which(data$Sequence == seq),]
  
  pltyscale <- {
    scale_y_continuous(
      limits = c(min(log10(plotdata$y[which(plotdata$y!=0)])), 0),
      breaks = log10(c(.001, .1, 1)),
      labels = c(".001", ".1", "1")
    ) 
  }
  
  p <- ggplot(plotdata, aes(x = log10(x), y = log10(y), color = type)) +
    geom_point() +
    plttheme +
    pltxscale +
    pltyscale +
    pltlabs +
    theme(legend.position = "right") +
    labs(title = paste(seq, state))
  
  print(ggplotly(p))
} 