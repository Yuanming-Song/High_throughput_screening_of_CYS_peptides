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
if (loadnetstat) {
  #Load the network stats data
  base_dir <- "~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide"
  if (single_node) {
    load(file.path(base_dir, "tripeptide_network_stats_single_node.rda"))
    # Combine monomer and dimer data
    combined_df_ddedge_tripep_mon<-c()
    state<-"monomer"
    for (i in 1:(length(final_results[[1]]))) {
      combined_df_ddedge_tripep_mon<-rbind(combined_df_ddedge_tripep_mon,final_results[[1]][[i]])
    }
    combined_df_ddedge_tripep_mon<-cbind(combined_df_ddedge_tripep_mon,state)
    
    combined_df_ddedge_tripep_dim<-c()
    state<-"dimer"
    for (i in 1:length(final_results[[2]])) {
      combined_df_ddedge_tripep_dim<-rbind(combined_df_ddedge_tripep_dim,final_results[[2]][[i]])
    }
    combined_df_ddedge_tripep_dim<-cbind(combined_df_ddedge_tripep_dim,state)
    

  } else {
    load(file.path(base_dir, "tripeptide_network_stats.rda"))
    combined_df_ddedge_tripep_mon<-final_results[[1]]
    combined_df_ddedge_tripep_dim<-final_results[[2]]
  }
  combined_df_ddedge_tripep <- rbind(combined_df_ddedge_tripep_dim, combined_df_ddedge_tripep_mon)
 
}

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

# Function to get position of C
get_C_position <- function(seq) {
  seq_chars <- strsplit(seq, "")[[1]]
  return(which(seq_chars == "C"))
}

# Process the edge data
edge_data <- combined_df_ddedge_tripep[combined_df_ddedge_tripep$stat_name == "edges" & 
                                         combined_df_ddedge_tripep$frame == max(combined_df_ddedge_tripep$frame), ]

# Process the AP score data (using ratio column)
ap_data <- combined_df_ddedge_tripep[combined_df_ddedge_tripep$frame == max(combined_df_ddedge_tripep$frame), ]
ap_data$value <- ap_data$ratio  # Use ratio as the value for aggregation

# Function to process data (either edges or AP scores)
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
  wide$C_position <- sapply(wide$seqname, get_C_position)
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
  
  # Calculate delta values
  wide$delta <- wide[[value_cols[1]]] - wide[[value_cols[2]]]
  
  return(wide)
}

# Process both edge and AP data
edge_wide <- process_data(edge_data, c("dimer_edges", "monomer_edges"))
ap_wide <- process_data(ap_data, c("dimer_ap", "monomer_ap"))

# Function to create correlation plot
create_correlation_plot <- function(data, x_col, y_col, x_label, y_label, title, ratio=1) {
  # Handle ratio input
  if (length(ratio) == 1) {
    x_ratio <- ratio
    y_ratio <- ratio
  } else {
    x_ratio <- ratio[1]
    y_ratio <- ratio[2]
  }
  
  # Calculate data ranges
  x_values <- data[[x_col]]/x_ratio
  y_values <- data[[y_col]]/y_ratio
  
  # Calculate 5th and 95th percentiles
  x_range <- quantile(x_values, probs = c(0.05, 0.95), na.rm = TRUE)
  y_range <- quantile(y_values, probs = c(0.05, 0.95), na.rm = TRUE)
  
  # Check if data is within default range (0.5, 6)
  if (x_range[1] >= 0.5 && x_range[2] <= 6 && 
      y_range[1] >= 0.5 && y_range[2] <= 6) {
    min_val <- 0.5
    max_val <- 6
  } else {
    # Determine common range for both axes
    min_val <- min(x_range[1], y_range[1])
    max_val <- max(x_range[2], y_range[2])
    
    # Round to nice numbers
    min_val <- floor(min_val * 10) / 10
    max_val <- ceiling(max_val * 10) / 10
  }
  
  plot <- ggplot(data, 
                 aes(x = get(x_col)/x_ratio, 
                     y = get(y_col)/y_ratio, 
                     col = Category,
                     text = label)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(x = x_label, 
         y = y_label,
         col = "Side Chain Types") +
    pltthemetemp +
    scale_x_continuous(limits = c(min_val, max_val)) +
    scale_y_continuous(limits = c(min_val, max_val)) +
    coord_fixed(ratio = 1) +
    scale_color_discrete() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
    facet_wrap(~C_position, labeller = labeller(C_position = function(x) paste("C at position", x)))
  
  return(plot)
}

# Create edge correlation plot
edge_correlation_plot <- create_correlation_plot(
  edge_wide, 
  "monomer_edges", "dimer_edges",
  "Monomer Edges", "Dimer Edges",
  "Edge Correlation in Tripeptides",
  if (single_node) 300 else c(300, 150)
)

# Create AP correlation plot
ap_correlation_plot <- create_correlation_plot(
  ap_wide, 
  "monomer_ap", "dimer_ap",
  "Monomer AP Score", "Dimer AP Score",
  "AP Score Correlation in Tripeptides"
)

# Create interactive plots with simpler legends
edge_correlation_interactive <- edge_correlation_plot +
  theme(legend.position = "right") +
  guides(color = guide_legend(nrow = NULL))

ap_correlation_interactive <- ap_correlation_plot +
  theme(legend.position = "right") +
  guides(color = guide_legend(nrow = NULL))

# Print interactive plots
print(ggplotly(edge_correlation_interactive))
print(ggplotly(ap_correlation_interactive))

# Save high-resolution plots
if (exists("save_plots") && save_plots) {
  ggsave(paste0(pltsavedir, "/Edge_correlation_by_C_position_tripep.png"), 
         plot = edge_correlation_plot+theme(
           legend.key.spacing.y  = unit(-10,"pt")
         ),
         dpi = 1100, 
         width = 6.5,
         height = 3,
         units = "in")
  
  ggsave(paste0(pltsavedir, "/AP_correlation_by_C_position_tripep.png"), 
         plot = ap_correlation_plot+theme(
           legend.key.spacing.y  = unit(-10,"pt")
         ),
         dpi = 1100, 
         width = 6.5,
         height = 3.5,
         units = "in")
  
  
  # Get top 10 candidates
  top_10_edge <- edge_wide %>% 
    arrange(desc(delta)) %>% 
    select(seqname, dimer_edges, monomer_edges, delta, Category) %>% 
    head(10)
  
  # Calculate WW scale values for top 10 edges
  top_10_edge$ww_absolute <- sapply(top_10_edge$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$total_hydrophobicity)
  })
  
  top_10_edge$ww_normalized <- sapply(top_10_edge$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$normalized_hydrophobicity)
  })
  
  # Get top by category for edges
  top_by_category_edge <- edge_wide %>%
    group_by(Category) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(seqname, dimer_edges, monomer_edges, delta, Category)
  
  # Calculate WW scale values for top by category edges
  top_by_category_edge$ww_absolute <- sapply(top_by_category_edge$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$total_hydrophobicity)
  })
  
  top_by_category_edge$ww_normalized <- sapply(top_by_category_edge$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$normalized_hydrophobicity)
  })
  
  # Get top 10 candidates for AP
  top_10_ap <- ap_wide %>% 
    arrange(desc(delta)) %>% 
    select(seqname, dimer_ap, monomer_ap, delta, Category) %>% 
    head(10)
  
  # Calculate WW scale values for top 10 AP
  top_10_ap$ww_absolute <- sapply(top_10_ap$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$total_hydrophobicity)
  })
  
  top_10_ap$ww_normalized <- sapply(top_10_ap$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$normalized_hydrophobicity)
  })
  
  # Get top by category for AP
  top_by_category_ap <- ap_wide %>%
    group_by(Category) %>%
    arrange(desc(delta), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(seqname, dimer_ap, monomer_ap, delta, Category)
  
  # Calculate WW scale values for top by category AP
  top_by_category_ap$ww_absolute <- sapply(top_by_category_ap$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$total_hydrophobicity)
  })
  
  top_by_category_ap$ww_normalized <- sapply(top_by_category_ap$seqname, function(seq) {
    result <- calculate_ww_hydrophobicity(seq)
    return(result$normalized_hydrophobicity)
  })
  
  # Create LaTeX tables for edges
  edge_latex_table <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
    "        \\hline",
    "        Sequence & Dimer Edges & Monomer Edges & dEdge & logP & logP$^\\prime$ & Category \\\\",
    "        \\hline"
  )
  
  # Add rows for edge data
  for(i in 1:10) {
    row <- sprintf("        %s & %.2f & %.2f & %.2f & %.3f & %.3f & %s \\\\",
                   top_10_edge$seqname[i], 
                   top_10_edge$dimer_edges[i]/300, 
                   top_10_edge$monomer_edges[i]/300,
                   top_10_edge$delta[i]/300,
                   top_10_edge$ww_absolute[i],
                   top_10_edge$ww_normalized[i],
                   as.character(top_10_edge$Category[i]))
    edge_latex_table <- c(edge_latex_table, row)
  }
  
  edge_latex_table <- c(edge_latex_table,
                        "        \\hline",
                        "    \\end{tabular}",
                        "    \\caption{Top 10 tripeptide sequences ranked by the highest difference in network edges (dEdge) between dimeric and monomeric forms. logP$^\\prime$ is the normalized Wimley-White hydrophobicity scale calculated as (logP - logP$_{min}$)/(logP$_{max}$ - logP$_{min}$), where logP$_{min}$ and logP$_{max}$ are the theoretical minimum and maximum values possible for a tripeptide.}",
                        "    \\label{tab:tri_top10_dedge}",
                        "\\end{table}")
  
  # Create LaTeX table for top_by_category_edge
  edge_by_category_latex_table <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
    "        \\hline",
    "        Sequence & Dimer Edges & Monomer Edges & dEdge & logP & logP$^\\prime$ & Category \\\\",
    "        \\hline"
  )
  
  for(i in 1:nrow(top_by_category_edge)) {
    row <- sprintf("        %s & %.2f & %.2f & %.2f & %.3f & %.3f & %s \\\\",
                   top_by_category_edge$seqname[i], 
                   top_by_category_edge$dimer_edges[i]/300, 
                   top_by_category_edge$monomer_edges[i]/300,
                   top_by_category_edge$delta[i]/300,
                   top_by_category_edge$ww_absolute[i],
                   top_by_category_edge$ww_normalized[i],
                   as.character(top_by_category_edge$Category[i]))
    edge_by_category_latex_table <- c(edge_by_category_latex_table, row)
  }
  
  edge_by_category_latex_table <- c(edge_by_category_latex_table,
                                    "        \\hline",
                                    "    \\end{tabular}",
                                    "    \\caption{Top tripeptide sequence with the highest dEdge from each sidechain chemical category. logP$^\\prime$ is the normalized Wimley-White hydrophobicity scale calculated as (logP - logP$_{min}$)/(logP$_{max}$ - logP$_{min}$), where logP$_{min}$ and logP$_{max}$ are the theoretical minimum and maximum values possible for a tripeptide.}",
                                    "    \\label{tab:tri_top_bycat_dedge}",
                                    "\\end{table}")
  
  # Create LaTeX tables for AP
  ap_latex_table <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
    "        \\hline",
    "        Sequence & Dimer AP & Monomer AP & dAP & logP & logP$^\\prime$ & Category \\\\",
    "        \\hline"
  )
  
  # Add rows for AP data
  for(i in 1:10) {
    row <- sprintf("        %s & %.3f & %.3f & %.3f & %.3f & %.3f & %s \\\\",
                   top_10_ap$seqname[i], 
                   top_10_ap$dimer_ap[i], 
                   top_10_ap$monomer_ap[i],
                   top_10_ap$delta[i],
                   top_10_ap$ww_absolute[i],
                   top_10_ap$ww_normalized[i],
                   as.character(top_10_ap$Category[i]))
    ap_latex_table <- c(ap_latex_table, row)
  }
  
  ap_latex_table <- c(ap_latex_table,
                      "        \\hline",
                      "    \\end{tabular}",
                      "    \\caption{Top 10 tripeptide sequences ranked by the highest difference in Aggregation Propensity (dAP) between dimeric and monomeric states. logP$^\\prime$ is the normalized Wimley-White hydrophobicity scale calculated as (logP - logP$_{min}$)/(logP$_{max}$ - logP$_{min}$), where logP$_{min}$ and logP$_{max}$ are the theoretical minimum and maximum values possible for a tripeptide.}",
                      "    \\label{tab:tri_top10_dap}",
                      "\\end{table}")
  
  # Create LaTeX table for top_by_category_ap
  ap_by_category_latex_table <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
    "        \\hline",
    "        Sequence & Dimer AP & Monomer AP & dAP & logP & logP$^\\prime$ & Category \\\\",
    "        \\hline"
  )
  
  for(i in 1:nrow(top_by_category_ap)) {
    row <- sprintf("        %s & %.3f & %.3f & %.3f & %.3f & %.3f & %s \\\\",
                   top_by_category_ap$seqname[i],
                   top_by_category_ap$dimer_ap[i],
                   top_by_category_ap$monomer_ap[i],
                   top_by_category_ap$delta[i],
                   top_by_category_ap$ww_absolute[i],
                   top_by_category_ap$ww_normalized[i],
                   as.character(top_by_category_ap$Category[i]))
    ap_by_category_latex_table <- c(ap_by_category_latex_table, row)
  }
  
  ap_by_category_latex_table <- c(ap_by_category_latex_table,
                                  "        \\hline",
                                  "    \\end{tabular}",
                                  "    \\caption{Top tripeptide sequence with the highest dAP from each sidechain chemical category. logP$^\\prime$ is the normalized Wimley-White hydrophobicity scale calculated as (logP - logP$_{min}$)/(logP$_{max}$ - logP$_{min}$), where logP$_{min}$ and logP$_{max}$ are the theoretical minimum and maximum values possible for a tripeptide.}",
                                  "    \\label{tab:tri_top_bycat_dap}",
                                  "\\end{table}")
  
  writeLines(edge_latex_table,
             paste0(pltsavedir, "/top_10_dedge_tripep.tex"))
  writeLines(edge_by_category_latex_table,
             paste0(pltsavedir, "/top_by_category_dedge_tripep.tex"))
  writeLines(ap_latex_table,
             paste0(pltsavedir, "/top_10_dap_tripep.tex"))
  writeLines(ap_by_category_latex_table,
             paste0(pltsavedir, "/top_by_category_dap_tripep.tex"))
} 
