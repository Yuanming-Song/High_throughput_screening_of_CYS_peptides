# Load required libraries
library(readxl)
library(ggplot2)
library(cowplot)
library(reshape2)
if (reloaddata) {
# Set base directory and file path
base_dir <- "/Users/song/Documents/Research/MRSEC/ML-MD Peptide /Experimental data"
xlsx_file <- file.path(base_dir, "tripeps-data-01.xlsx")

# Read all sheets
sheets <- excel_sheets(xlsx_file)
data_list <- list()

special_sheets <- c(1, 3, 5)

for (i in seq_along(sheets)) {
  if (i %in% special_sheets) {
    temp <- read_excel(xlsx_file, sheet = sheets[i], skip = 1, col_names = FALSE)
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
    data_list[[sheets[i]]] <- read_excel(xlsx_file, sheet = sheets[i])
  }
}
}
# Extract sheet 5 (data) and sheet 6 (info)
data_sheet5 <- data_list[[sheets[5]]]
info_sheet6 <- data_list[[sheets[6]]]
colnames(info_sheet6) <- c("colname", "seqname", "something", "pH")

# Determine global y-axis range from only the used columns in sheet6
cols_used <- unique(info_sheet6$colname)
all_y_values <- unlist(data_sheet5[, cols_used, drop = FALSE])
global_y_range <- range(as.numeric(all_y_values), na.rm = TRUE)

# Create plots
plot_list <- list()
combined_plots <- list()

# Split info by pH
pH_groups <- split(info_sheet6, info_sheet6$pH)

for (pH_val in names(pH_groups)) {
  seq_groups <- split(pH_groups[[pH_val]], pH_groups[[pH_val]]$seqname)
  
  plot_list[[pH_val]] <- list()
  
  for (seqname in names(seq_groups)) {
    sub_info <- seq_groups[[seqname]]
    
    cols_to_plot <- sub_info$colname
    legend_labels <- sub_info[[5]] # Use the 5th column for legends
    
    plot_data <- data_sheet5[, c("Hours", cols_to_plot), drop = FALSE]
    colnames(plot_data)[1] <- "X"
    colnames(plot_data)[-1] <- legend_labels
    
    plot_data_long <- melt(plot_data, id.vars = "X")
    
    plot_data_long$variable <- factor(plot_data_long$variable, levels = legend_labels, labels = legend_labels)

    x_range <- range(as.numeric(plot_data_long$X), na.rm = TRUE)
    
    p <- ggplot(plot_data_long, aes(x = as.numeric(X), y = as.numeric(value), color = variable)) +
      geom_line() +
      labs(title = paste0(seqname, ", pH=", pH_val),
           x = NULL, y = NULL, color = expression("["~H[2]~O[2]~"] (mM)")) +
      coord_cartesian(xlim = x_range, ylim = global_y_range) +
      scale_y_continuous(labels = scales::scientific) +
      plttheme +
      theme(axis.title = element_blank(),
            legend.position = "none")
    
    plot_list[[pH_val]][[seqname]] <- p
  }
}

## Correctly order the plots into a matrix: rows = sequences, columns = pH
all_seqs <- sort(unique(info_sheet6$seqname))
all_pHs <- sort(unique(info_sheet6$pH))

ordered_plots <- list()
for (seq in all_seqs) {
  for (pH in all_pHs) {
    ordered_plots[[paste(seq, pH, sep = "_")]] <- plot_list[[as.character(pH)]][[seq]]
  }
}

combined_body <- plot_grid(plotlist = ordered_plots, ncol = length(all_pHs))

combined_with_labels <- ggdraw() +
  draw_plot(combined_body, x = 0.05, y = 0.01, width = 0.95, height = 0.99) +
  draw_label("RNU", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Time (Hours)", x = 0.5, y = 0.01, size = 12)

## Extract shared legend once with correct theme for extraction
first_pH <- names(plot_list)[1]
first_seq <- names(plot_list[[first_pH]])[1]
legend_sheet56 <- get_legend(plot_list[[first_pH]][[first_seq]] + 
                               theme(legend.position = "right", 
                                     legend.direction = "horizontal",
                                     legend.title.align = 0.5,
                                     legend.spacing.x = unit(0.2, 'cm'),
                                     legend.title = element_text(margin = margin(r = 10))))

final_combined_sheet56 <- plot_grid(combined_with_labels, legend_sheet56, ncol = 1, rel_heights = c(0.95, 0.05))

# Example: view one combined plot
# print(combined_plots[["7"]])

## ---- Now handle Sheets 1-2 (TAC) ----
data_sheet1 <- data_list[[sheets[1]]]
data_sheet1[,1]<-rbind(data_sheet1[-1,1],0)
data_sheet1<-data_sheet1[-nrow(data_sheet1),]
colnames(data_sheet1)[1]<-"Hours"

info_sheet2 <- data_list[[sheets[2]]]
colnames(info_sheet2) <- c("colname", "buffer_conc", "pH", "peptide_conc")

# Find unique combinations
info_sheet2 <- unique(info_sheet2[, c("buffer_conc", "peptide_conc", "colname", "pH")])

# Build plots
plots_TAC <- list()
combo_TAC <- info_sheet2[, c("buffer_conc", "peptide_conc")]
combo_TAC <- unique(combo_TAC)
combo_TAC <- combo_TAC[order(combo_TAC$peptide_conc, combo_TAC$buffer_conc), ]

tac_all_values <- unlist(data_sheet1[, unique(info_sheet2$colname)])
tac_y_range <- range(as.numeric(tac_all_values), na.rm = TRUE)

for (i in seq_len(nrow(combo_TAC))) {
  buf <- combo_TAC$buffer_conc[i]
  pep <- combo_TAC$peptide_conc[i]
  subset_info <- subset(info_sheet2, buffer_conc == buf & peptide_conc == pep)
  cols_to_plot <- subset_info$colname
  pH_labels <- subset_info$pH
  
  plot_data <- data_sheet1[, c("Hours", cols_to_plot)]
  colnames(plot_data)[1] <- "X"
  colnames(plot_data)[-1] <- pH_labels
  
  plot_data_long <- melt(plot_data, id.vars = "X")
  
  p <- ggplot(plot_data_long, aes(x = as.numeric(X), y = as.numeric(value), color = factor(variable))) +
    geom_line() +
    labs(x = NULL, y = NULL, color = "pH", title = paste0("[P]=", pep, ",[PB]=", buf)) +
    plttheme +
    coord_cartesian(ylim = tac_y_range) +
          scale_y_continuous(labels = scales::scientific) +
    theme(legend.position = "none")
  
  plots_TAC[[i]] <- p
}

plt <- ggplot(plot_data_long, aes(x = as.numeric(X), y = as.numeric(value), color = factor(variable))) +
  geom_line() +
  labs(x = NULL, y = NULL, color = "pH", title = paste0("[P]=", pep, ",[PB]=", buf)) + plttheme + 
  theme(legend.position = "right", legend.direction = "horizontal")

# Combine TAC plots
legend_TAC <- get_legend(plt)
combined_TAC <- plot_grid(plotlist = plots_TAC, ncol = 3)

combined_TAC_final <- ggdraw() +
  draw_plot(combined_TAC, x = 0.05, y = 0.01, width = 0.95, height = 0.99) +
  draw_label("RNU", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Time (Hours)", x = 0.5, y = 0.01, size = 12)

final_TAC <- plot_grid(combined_TAC_final, legend_TAC, ncol = 1, rel_heights = c(0.95, 0.05))

## ---- Now handle Sheets 3-4 (CAS) ----
data_sheet3 <- data_list[[sheets[3]]]
info_sheet4 <- data_list[[sheets[4]]]
colnames(info_sheet4) <- c("colname", "buffer_conc", "pH", "peptide_conc")

info_sheet4 <- unique(info_sheet4[, c("buffer_conc", "peptide_conc", "colname", "pH")])

plots_CAS <- list()
combo_CAS <- info_sheet4[, c("buffer_conc", "peptide_conc")]
combo_CAS <- unique(combo_CAS)
combo_CAS <- combo_CAS[order(combo_CAS$peptide_conc, combo_CAS$buffer_conc), ]

cas_all_values <- unlist(data_sheet3[, unique(info_sheet4$colname)])
cas_y_range <- range(as.numeric(cas_all_values), na.rm = TRUE)

for (i in seq_len(nrow(combo_CAS))) {
  buf <- combo_CAS$buffer_conc[i]
  pep <- combo_CAS$peptide_conc[i]
  subset_info <- subset(info_sheet4, buffer_conc == buf & peptide_conc == pep)
  cols_to_plot <- subset_info$colname
  pH_labels <- subset_info$pH
  
  plot_data <- data_sheet3[, c("Hours", cols_to_plot)]
  colnames(plot_data)[1] <- "X"
  colnames(plot_data)[-1] <- pH_labels
  
  plot_data_long <- melt(plot_data, id.vars = "X")
  
  p <- ggplot(plot_data_long, aes(x = as.numeric(X), y = as.numeric(value), color = factor(variable))) +
    geom_line() +
    labs(x = NULL, y = NULL, color = "pH", title = paste0("[P]=", pep, ",[PB]=", buf)) +
    plttheme +
    coord_cartesian(ylim = cas_y_range) +
          scale_y_continuous(labels = scales::scientific) +
    theme(legend.position = "none")
  
  plots_CAS[[i]] <- p
}

plt <- ggplot(plot_data_long, aes(x = as.numeric(X), y = as.numeric(value), color = factor(variable))) +
  geom_line() +
  labs(x = NULL, y = NULL, color = "pH") + 
  plttheme + 
  theme(legend.position = "right", legend.direction = "horizontal")

# Combine CAS plots
legend_CAS <- get_legend(plt)
combined_CAS <- plot_grid(plotlist = plots_CAS, ncol = 3)

combined_CAS_final <- ggdraw() +
  draw_plot(combined_CAS, x = 0.05, y = 0.01, width = 0.95, height = 0.99) +
  draw_label("RNU", x = 0.02, y = 0.5, angle = 90, size = 12) +
  draw_label("Time (Hours)", x = 0.5, y = 0.01, size = 12)

final_CAS <- plot_grid(combined_CAS_final, legend_CAS, ncol = 1, rel_heights = c(0.95, 0.05))

# Print final plots
print(final_TAC)
print(final_CAS)
print(final_combined_sheet56)

# Optional saving
if (exists("saveplt") && saveplt) {
  ggsave(filename = file.path(base_path_plt, "Experimental_data_TAC_combined_plot.png"), plot = final_TAC, width = 6.5, height = 4)
  ggsave(filename = file.path(base_path_plt , "Experimental_data_CAS_combined_plot.png"), plot = final_CAS, width = 6.5, height = 4)
  ggsave(filename = file.path(base_path_plt, "Experimental_data_Tripeptide_combined_plot.png"), plot = final_combined_sheet56, width = 6.5, height = 5.6)
}
