# Custom theme for text sizes
#if (APonly) loadnetstat <- FALSE
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
hydrophobic <- c("A", "I", "L", "M", "V", "P")
polar_uncharged <- c("S", "T", "N", "Q", "G")
positive_charged <- c("R", "K", "H")
negative_charged <- c("D", "E")
aromatic <- c("F", "W", "Y")

# Function to get unique amino acid types in sequence
get_aa_types <- function(amino_acids) {
  types <- c()
  if(any(amino_acids %in% hydrophobic)) types <- c(types, "Hydrophobic")
  if(any(amino_acids %in% polar_uncharged)) types <- c(types, "Polar")
  if(any(amino_acids %in% c(positive_charged, negative_charged))) types <- c(types, "Charged")
  if(any(amino_acids %in% aromatic)) types <- c(types, "Aromatic")
  return(paste(sort(types), collapse = " + "))
}


# Function to get position of C
get_C_position <- function(seq) {
  seq_chars <- strsplit(seq, "")[[1]]
  return(which(seq_chars == "C"))
}
# Function to classify sequences based on amino acid types present
classify_sequence_tetra <- function(a1, a2, a3, a4) {
  # Convert inputs to character if they aren't already
  amino_acids <- c(as.character(a1), as.character(a2), 
                   as.character(a3), as.character(a4))
  amino_acids <- amino_acids[amino_acids != "C"]  # Ignore C in classification
  return(get_aa_types(amino_acids))
}
if (loadnetstat && !APonly) {
  # Initialize empty dataframe for combined results
  combined_df_netstat_tetra <- data.frame()
  
  # List of available C positions (currently only 1, add others when available)
  c_positions <- c(1,2,3,4)  # Add 2,3,4 when those files are available
  
  # Load and combine data for each C position
  for(pos in c_positions) {
    rda_file <- paste0("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/tetrapeptide_network_stats_C", pos, "_single_node.rda")
    if(file.exists(rda_file)) {
      load(rda_file)
      pos_data <- rbind(final_results[[1]], final_results[[2]])
      combined_df_netstat_tetra <- rbind(combined_df_netstat_tetra, pos_data)
    }
  }
}

if (APonly && loadnetstat) {
  c_positions <- 1:4
  combined_ap_data_tetra <- data.frame()
  for (pos in c_positions) {
    file_path <- paste0("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/AP_analysis/Tetrapeptide/SASA_score/SASA_result_with_common_mon_ini_tetra_C", pos, ".txt")
    if (file.exists(file_path)) {
      tryCatch({
        temp_data <- read.table(file_path, header = FALSE)
        temp_data$C_position <- pos
        combined_ap_data_tetra <- rbind(combined_ap_data_tetra, temp_data)
      }, error = function(e) {
        message(sprintf("Skipping position %d: %s", pos, e$message))
      })
    }
  }
  colnames(combined_ap_data_tetra)[1:7] <- c("amino1", "amino2", "amino3", "amino4", "monomer_ap", "dimer_ap","C_position")
  combined_ap_data_tetra$seqname <- paste0(combined_ap_data_tetra$amino1,
                                     combined_ap_data_tetra$amino2,
                                     combined_ap_data_tetra$amino3,
                                     combined_ap_data_tetra$amino4)
  combined_ap_data_tetra$Category <- mapply(classify_sequence_tetra, 
                                      combined_ap_data_tetra$amino1, 
                                      combined_ap_data_tetra$amino2, 
                                      combined_ap_data_tetra$amino3,
                                      combined_ap_data_tetra$amino4)
  combined_ap_data_tetra$label <- combined_ap_data_tetra$seqname
  combined_ap_data_tetra$Category <- as.factor(combined_ap_data_tetra$Category)
  combined_ap_data_tetra$Category <- as.factor(combined_ap_data_tetra$Category)
  combined_ap_data_tetra$delta<-combined_ap_data_tetra$dimer_ap-combined_ap_data_tetra$monomer_ap
}

if (!APonly) {
  # Process the edge data
  edge_data <- combined_df_netstat_tetra[combined_df_netstat_tetra$stat_name == "edges" & 
                                       combined_df_netstat_tetra$frame == max(combined_df_netstat_tetra$frame), ]
}

if (!APonly) {
  # Process the AP score data (using ratio column)
  ap_data <- combined_df_netstat_tetra[combined_df_netstat_tetra$frame == max(combined_df_netstat_tetra$frame), ]
  ap_data$value <- ap_data$ratio  # Use ratio as the value for aggregation
}

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
  wide$amino4 <- substr(wide$seqname, 4, 4)
  wide$Category <- mapply(classify_sequence_tetra, 
                         wide$amino1, 
                         wide$amino2, 
                         wide$amino3,
                         wide$amino4)
  wide$label <- wide$seqname
  
  # Ensure Category is a proper factor
  wide$Category <- as.character(wide$Category)
  wide$Category[is.na(wide$Category)] <- "Other"
  wide$Category <- as.factor(wide$Category)
  
  return(wide)
}

if (!APonly) {
  # Process edge data
  edge_wide <- process_data(edge_data, c("dimer_edges", "monomer_edges"))
}

if (APonly) {
  ap_wide <- combined_ap_data_tetra

} else {
  ap_wide <- process_data(ap_data, c("dimer_ap", "monomer_ap"))
}

# Function to create correlation plot
create_correlation_plot <- function(data, x_col, y_col, x_label, y_label, title,ratio=1) {
  plot <- ggplot(data, 
         aes(x = get(x_col)/ratio, 
             y = get(y_col)/ratio, 
             col = Category,
             text = label)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(x = x_label, 
         y = y_label,
         col = "Side Chain Types"
         ) +
    pltthemetemp +
    scale_x_continuous(limits = c(0.5, 6)) +
    scale_y_continuous(limits = c(0.5, 6)) +
    coord_fixed(ratio = 1) +
    scale_color_discrete() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
    facet_wrap(~C_position, labeller = labeller(C_position = function(x) paste("C at position", x)))
  
  return(plot)
}

if (!APonly) {
  # Create edge correlation plot
  edge_correlation_plot <- create_correlation_plot(
    edge_wide, 
    "monomer_edges", "dimer_edges",
    "Monomer Edges", "Dimer Edges",
    "Edge Correlation in Tetrapeptides",300
  )
}

# Create AP correlation plot
ap_correlation_plot <- create_correlation_plot(
  ap_wide, 
  "monomer_ap", "dimer_ap",
  "Monomer AP Score", "Dimer AP Score",
  "AP Score Correlation in Tetrapeptides"
)

if (!APonly) {
  # Create interactive plots with simpler legends
  edge_correlation_interactive <- edge_correlation_plot +
    theme(legend.position = "right") +
    guides(color = guide_legend(nrow = NULL))
}

ap_correlation_interactive <- ap_correlation_plot +
  theme(legend.position = "right") +
  guides(color = guide_legend(nrow = NULL))

if (!APonly) {
  # Print interactive edge plot
  print(ggplotly(edge_correlation_interactive))
}

print(ggplotly(ap_correlation_interactive))

# Save high-resolution plots
if (exists("save_plots") && save_plots) {
  if (!APonly) {
    ggsave(paste0(pltsavedir, "/Edge_correlation_by_C_position_tetra.png"), 
           plot = edge_correlation_plot+theme(
             legend.key.spacing.y  = unit(-10,"pt")
           ),
           dpi = 1100, 
           width = 6.5,
           height = 3,
           units = "in")
  }
  
  ggsave(paste0(pltsavedir, "/AP_correlation_by_C_position_tetra.png"), 
         plot = ap_correlation_plot+theme(
           legend.key.spacing.y  = unit(-10,"pt")
         ),
         dpi = 1100, 
         width = 6.5,
         height = 5.5,
         units = "in")
} 

if (!APonly) {
  edge_wide$delta<-edge_wide$dimer_edges-edge_wide$monomer_edges
  # Get top 10 candidates for edges
  top_10_edge <- edge_wide %>% 
    arrange(desc(delta)) %>% 
    select(seqname, dimer_edges, monomer_edges, delta, Category) %>% 
    head(20)
  
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
}
ap_wide$delta<-ap_wide$dimer_ap-ap_wide$monomer_ap
# Get top 10 candidates for AP (always)
top_10_ap <- ap_wide %>%
  arrange(desc(delta)) %>%
  select(seqname, dimer_ap, monomer_ap, delta, Category) %>%
  head(10)

# Calculate WW scale values for top 10
top_10_ap$ww_absolute <- sapply(top_10_ap$seqname, function(seq) {
  result <- calculate_ww_hydrophobicity(seq)
  return(result$total_hydrophobicity)
})

top_10_ap$ww_normalized <- sapply(top_10_ap$seqname, function(seq) {
  result <- calculate_ww_hydrophobicity(seq)
  return(result$normalized_hydrophobicity)
})

# Get top candidate by Category for AP (always)
top_by_category_ap <- ap_wide %>%
  group_by(Category) %>%
  arrange(desc(delta), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(seqname, dimer_ap, monomer_ap, delta, Category)

# Calculate WW scale values for top by category
top_by_category_ap$ww_absolute <- sapply(top_by_category_ap$seqname, function(seq) {
  result <- calculate_ww_hydrophobicity(seq)
  return(result$total_hydrophobicity)
})

top_by_category_ap$ww_normalized <- sapply(top_by_category_ap$seqname, function(seq) {
  result <- calculate_ww_hydrophobicity(seq)
  return(result$normalized_hydrophobicity)
})

if (!APonly) {
  # Create LaTeX tables for edges
  edge_latex_table <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
    "        \\hline",
    "        Sequence & Dimer Edges & Monomer Edges & dEdge & WW Abs & WW Norm & Category \\\\",
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
  
  # Add closing lines for edge table
  edge_latex_table <- c(edge_latex_table,
                       "        \\hline",
                       "    \\end{tabular}",
                       "    \\caption{Top 10 tetrapeptide sequences ranked by the highest difference in network edges (dEdge) between dimeric and monomeric forms.}",
                       "    \\label{tab:tetra_top10_dedge}",
                       "\\end{table}")
  
  # Create LaTeX table for top_by_category_edge
  edge_by_category_latex_table <- c(
    "\\begin{table}[H]",
    "    \\centering",
    "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
    "        \\hline",
    "        Sequence & Dimer Edges & Monomer Edges & dEdge & WW Abs & WW Norm & Category \\\\",
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
                       "    \\caption{Top tetrapeptide sequence with the highest dEdge from each sidechain chemical category.}",
                       "    \\label{tab:tetra_top_bycat_dedge}",
                       "\\end{table}")
}

# Create LaTeX tables for AP (always)
ap_latex_table <- c(
  "\\begin{table}[H]",
  "    \\centering",
  "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
  "        \\hline",
  "        Sequence & Dimer AP & Monomer AP & dAP & WW Abs & WW Norm & Category \\\\",
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

# Add closing lines for AP table
ap_latex_table <- c(ap_latex_table,
                   "        \\hline",
                   "    \\end{tabular}",
                   "    \\caption{Top 10 tetrapeptide sequences ranked by the highest difference in Aggregation Propensity (dAP) between dimeric and monomeric states.}",
                   "    \\label{tab:tetra_top10_dap}",
                   "\\end{table}")

# Create LaTeX table for top_by_category_ap
ap_by_category_latex_table <- c(
  "\\begin{table}[H]",
  "    \\centering",
  "    \\begin{tabular}{|c|c|c|c|c|c|c|}",
  "        \\hline",
  "        Sequence & Dimer AP & Monomer AP & dAP & WW Abs & WW Norm & Category \\\\",
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
                 "    \\caption{Top tetrapeptide sequence with the highest dAP from each sidechain chemical category.}",
                 "    \\label{tab:tetra_top_bycat_dap}",
                 "\\end{table}")

# Save the LaTeX tables if save_plots is true
if (exists("save_plots") && save_plots) {
  if (!APonly) {
    writeLines(edge_latex_table,
              paste0(pltsavedir, "/top_10_dedge_tetra.tex"))
    writeLines(edge_by_category_latex_table,
               paste0(pltsavedir, "/top_by_category_dedge_tetra.tex"))
  }
  writeLines(ap_latex_table,
            paste0(pltsavedir, "/top_10_dap_tetra.tex"))
  writeLines(ap_by_category_latex_table,
             paste0(pltsavedir, "/top_by_category_dap_tetra.tex"))
}

