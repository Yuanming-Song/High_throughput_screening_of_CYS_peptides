text_theme <- theme(
  plot.title = element_text(size = plot_title_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.title = element_text(size = legend_title_size),
  legend.text = element_text(size = legend_text_size)
)
pltthemetemp<-plttheme+text_theme
# Set base directory and file path
base_dir <- "/Users/song/Documents/Research/MRSEC/ML-MD Peptide /Experimental data"
xlsx_file <- file.path(base_dir, "tripeps-data-01.xlsx")

# Read all sheets
sheets <- excel_sheets(xlsx_file)
data_list <- list()

special_sheets <- c(1, 3, 5)

for (i in seq_along(sheets)) {
  if (i %in% special_sheets) {
    temp <- suppressMessages(read_excel(xlsx_file, sheet = sheets[i], skip = 1, col_names = FALSE))
    if (i == 1) {
      first_col_name <- temp[2,1,drop=TRUE]
      other_col_names <- temp[2,-1]
      colnames(temp) <- c(first_col_name, other_col_names)
      temp <- temp[-c(1,2), ]
    } else {
      colnames(temp) <- c(temp[1,1], temp[2,-1])
      temp <- temp[-c(1,2), ]
    }
    data_list[[sheets[i]]] <- temp
  } else {
    data_list[[sheets[i]]] <- suppressMessages(read_excel(xlsx_file, sheet = sheets[i]))
  }
}

# Extract sheet 5 (data) and sheet 6 (info)
data_sheet5 <- data_list[[sheets[5]]]
info_sheet6 <- data_list[[sheets[6]]]
colnames(info_sheet6) <- c("colname", "seqname", "something", "pH")

# Filter for specific sequences and pH values
target_pairs <- list(
  list(seq = "CAS", pH = "9"),
  list(seq = "CGI", pH = "9"),
  list(seq = "CPV", pH = "3"),
  list(seq = "TAC", pH = "3")
)

# Filter info sheet for target sequences and pH values
filtered_info <- info_sheet6[info_sheet6$seqname %in% sapply(target_pairs, function(x) x$seq) & 
                            info_sheet6$pH %in% sapply(target_pairs, function(x) x$pH), ]

# Determine global y-axis range from only the used columns in sheet6
cols_used <- unique(filtered_info$colname)
all_y_values <- unlist(data_sheet5[, cols_used, drop = FALSE])
global_y_range <- range(as.numeric(all_y_values), na.rm = TRUE)

# Create a list to store the exact plots we want
plot_list <- list()
dEdge_values <- numeric(length(target_pairs))

# Create plots for each target pair
for (i in seq_along(target_pairs)) {
  pair <- target_pairs[[i]]
  seqname <- pair$seq
  pH_val <- pair$pH
  
  # Get the data for this specific sequence and pH
  sub_info <- filtered_info[filtered_info$seqname == seqname & filtered_info$pH == pH_val, ]
  
  if (nrow(sub_info) > 0) {
    cols_to_plot <- sub_info$colname
    legend_labels <- sub_info[[5]] # Use the 5th column for legends
    
    plot_data <- data_sheet5[, c("Hours", cols_to_plot), drop = FALSE]
    colnames(plot_data)[1] <- "X"
    colnames(plot_data)[-1] <- legend_labels
    
    plot_data_long <- melt(plot_data, id.vars = "X")
    
    plot_data_long$variable <- factor(plot_data_long$variable, levels = legend_labels, labels = legend_labels)

    x_range <- range(as.numeric(plot_data_long$X), na.rm = TRUE)
    
    # Get dEdge value for this sequence
    seq_data <- combined_df_ddedge_tripep[combined_df_ddedge_tripep$seq == seqname & 
                                         combined_df_ddedge_tripep$stat_name == "edges" & 
                                         combined_df_ddedge_tripep$frame == max(combined_df_ddedge_tripep$frame), ]
    means <- aggregate(value ~ state, data = seq_data, FUN = mean)
    dEdge <- (means$value[means$state == "dimer"] - means$value[means$state == "monomer"]) / 300
    dEdge_values[i] <- dEdge
    
    p <- ggplot(plot_data_long, aes(x = as.numeric(X), y = as.numeric(value), color = variable)) +
      geom_line() +
      labs(title = paste0("Seq", i, ", pH=", pH_val, ", dEdge=", sprintf("%.2f", dEdge)),
           x = NULL, y = NULL, color = expression("["~H[2]~O[2]~"] (mM)")) +
      coord_cartesian(xlim = x_range, ylim = global_y_range) +
      scale_y_continuous(labels = scales::scientific) +
      pltthemetemp +
      theme(axis.title = element_blank(),
            legend.position = "none")
    
    plot_list[[i]] <- p
  }
}

# Sort plots by dEdge values
plot_order <- order(dEdge_values)
plot_list <- plot_list[plot_order]

# Update sequence numbers in titles to reflect the sorted order
for (i in seq_along(plot_list)) {
  old_title <- plot_list[[i]]$labels$title
  # Extract pH and dEdge values from old title
  parts <- strsplit(old_title, ", ")[[1]]
  pH_part <- parts[2]
  dEdge_part <- parts[3]
  # Create new title with updated sequence number
  plot_list[[i]]$labels$title <- paste0("Seq", i, ", ", pH_part, ", ", dEdge_part)
}

# Create the plot grid with exactly our four plots
combined_body <- plot_grid(plotlist = plot_list, ncol = 2)

combined_with_labels <- ggdraw() +
  draw_plot(combined_body, x = 0.05, y = 0.01, width = 0.95, height = 0.98) +
  draw_label("RNU", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Time (Hours)", x = 0.5, y = 0.02, size = 12)

## Extract shared legend once with correct theme for extraction
first_plot <- plot_list[[1]]
legend_sheet56 <- get_legend(first_plot + 
                               theme(legend.position = "right", 
                                     legend.direction = "horizontal",
                                     legend.title.align = 0.5,
                                     legend.spacing.x = unit(0.2, 'cm'),
                                     legend.title = element_text(margin = margin(r = 10))))

final_plot <- plot_grid(legend_sheet56, combined_with_labels, ncol = 1, rel_heights = c(0.05, 0.95))

# Print final plot
print(final_plot)

# Optional saving
if (exists("save_plots") && save_plots) {
  ggsave(filename = file.path(pltsavedir, "Experimental_data_tripeptide_poster.png"), 
         plot = final_plot, 
         width = 6.5, 
         height = 3,
         dpi = 1100)
}
