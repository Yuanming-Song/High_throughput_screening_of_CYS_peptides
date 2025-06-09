# dEdge_dAP_distribution_plots.R
# This script assumes you have three data frames in your environment:
#   - combined_d_dipep
#   - combined_d_tripep
#   - combined_d_tetrapeptide
# Each should have columns for dAP and dEdge (named as such, or as 'delta' for each property)
text_theme <- theme(
  plot.title = element_text(size = plot_title_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.title = element_text(size = legend_title_size),
  legend.text = element_text(size = legend_text_size)
)
pltthemetemp<-plttheme+text_theme
# Helper to calculate dEdge and dAP for the final frame
calculate_dEdge_dAP <- function(df) {
  # Filter for the final frame
  final_frame <- max(df$frame)
  df_final <- df[df$frame == final_frame, ]
  
  # For edge data, compute dEdge
  edge_data <- df_final[df_final$stat_name == "edges", ]
  edge_wide <- reshape(aggregate(value ~ seqname + state, data = edge_data, FUN = mean),
                       idvar = "seqname", timevar = "state", direction = "wide")
  names(edge_wide) <- c("seqname", "dimer_edges", "monomer_edges")
  edge_wide$dEdge <- edge_wide$dimer_edges - edge_wide$monomer_edges
  edge_wide$dEdge <- edge_wide$dEdge/300
  # For AP data, compute dAP
  ap_data <- df_final
  ap_data$value <- ap_data$ratio
  ap_wide <- reshape(aggregate(value ~ seqname + state, data = ap_data, FUN = mean),
                     idvar = "seqname", timevar = "state", direction = "wide")
  names(ap_wide) <- c("seqname", "dimer_ap", "monomer_ap")
  ap_wide$dAP <- ap_wide$dimer_ap - ap_wide$monomer_ap
  
  # Merge edge and AP data
  result <- merge(edge_wide, ap_wide, by = "seqname")
  return(result)
}

# ---- Data loading and processing ----
# If canonical data frames are not defined, load and process from .rda files

# Dipeptide loading and processing
if (!exists("combined_df_ddedge_dipep")) {
  load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Dipeptide/dipeptide_network_stats_single_node.rda")
  combined_df_ddedge_dipep_mon <- c()
  state <- "monomer"
  for (i in 1:(length(final_results[[1]]))) {
    combined_df_ddedge_dipep_mon <- rbind(combined_df_ddedge_dipep_mon, final_results[[1]][[i]])
  }
  combined_df_ddedge_dipep_mon <- cbind(combined_df_ddedge_dipep_mon, state)
  
  combined_df_ddedge_dipep_dim <- c()
  state <- "dimer"
  for (i in 1:length(final_results[[2]])) {
    combined_df_ddedge_dipep_dim <- rbind(combined_df_ddedge_dipep_dim, final_results[[2]][[i]])
  }
  combined_df_ddedge_dipep_dim <- cbind(combined_df_ddedge_dipep_dim, state)
  combined_df_ddedge_dipep <- rbind(combined_df_ddedge_dipep_dim, combined_df_ddedge_dipep_mon)
}

# Tripeptide loading and processing
if (!exists("combined_df_ddedge_tripep")) {
  load("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tripeptide/tripeptide_network_stats_single_node.rda")
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
  combined_df_ddedge_tripep <- rbind(combined_df_ddedge_tripep_dim, combined_df_ddedge_tripep_mon)
}

# Tetrapeptide loading and processing
if (!exists("combined_df_netstat_tetra")) {
  combined_df_netstat_tetra <- data.frame()
  c_positions <- c(1, 2, 3, 4)
  for (pos in c_positions) {
    rda_file <- paste0("~/Documents/Research/HPC/dfs2/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/tetrapeptide_network_stats_C", pos, "_single_node.rda")
    if (file.exists(rda_file)) {
      load(rda_file)
      pos_data <- rbind(final_results[[1]], final_results[[2]])
      combined_df_netstat_tetra <- rbind(combined_df_netstat_tetra, pos_data)
    }
  }
}

# ---- Data processing for each peptide type ----
# Calculate dEdge and dAP for each peptide type
di_df <- calculate_dEdge_dAP(combined_df_ddedge_dipep)
tri_df <- calculate_dEdge_dAP(combined_df_ddedge_tripep)
tetra_df <- calculate_dEdge_dAP(combined_df_netstat_tetra)

# Add peptide length labels
di_df$PeptideLength <- "Dipeptide"
tri_df$PeptideLength <- "Tripeptide"
tetra_df$PeptideLength <- "Tetrapeptide"

# Stack all data
all_df <- bind_rows(di_df, tri_df, tetra_df)

# Convert PeptideLength to factor with specific order
all_df$PeptideLength <- factor(all_df$PeptideLength, 
                              levels = c("Dipeptide", "Tripeptide", "Tetrapeptide"))

# 1. Distribution of dAP (normalized)
g1 <- ggplot(all_df, aes(x = dAP, fill = PeptideLength)) +
  geom_density(alpha = 0.5, stat = "density") +
  labs(title = "Normalized dAP Distribution by Peptide Length", x = "dAP", y = "Density") +
  pltthemetemp

# 2. Distribution of dEdge (normalized)
g2 <- ggplot(all_df, aes(x = dEdge, fill = PeptideLength)) +
  geom_density(alpha = 0.5, stat = "density") +
  labs(title = "Normalized dEdge Distribution by Peptide Length", x = "dEdge", y = "Density") +
  pltthemetemp

# 3. Scatter plots: dEdge vs dAP for each peptide length (not colored)
# Option 1: Facet by peptide length
g3 <- ggplot(all_df, aes(x = dAP, y = dEdge)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~PeptideLength, ncol = 3) +
  labs(title = "dEdge vs dAP by Peptide Length", x = "dAP", y = "dEdge") +
  pltthemetemp

# 4. Edge dimer vs Edge monomer scatter plot
g4 <- ggplot(all_df, aes(x = monomer_edges/300, y = dimer_edges/300)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~PeptideLength, ncol = 3) +
  labs(title = "", 
       x = "Monomer Edges", 
       y = "Dimer Edges") +
  coord_fixed(ratio = 1) +
  pltthemetemp

# 5. Edge vs AP plot, separated by Monomer/Dimer
# Create separate data frames for monomer and dimer
monomer_data <- all_df %>%
  select(monomer_edges, monomer_ap, PeptideLength) %>%
  rename(Edges = monomer_edges, AP = monomer_ap) %>%
  mutate(State = "Monomer") %>%
  filter(AP != 100)  # Remove points where AP is 100

dimer_data <- all_df %>%
  select(dimer_edges, dimer_ap, PeptideLength) %>%
  rename(Edges = dimer_edges, AP = dimer_ap) %>%
  mutate(State = "Dimer") %>%
  filter(AP != 100)  # Remove points where AP is 100

# Combine the data
edge_ap_data <- bind_rows(monomer_data, dimer_data)
# Ensure Monomer appears on the left
edge_ap_data$State <- factor(edge_ap_data$State, levels = c("Monomer", "Dimer"))

g5 <- ggplot(edge_ap_data, aes(x = AP, y = Edges/300)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~State, ncol = 2) +
  labs(title = "", 
       x = "AP Score", 
       y = "Edges") +
  pltthemetemp

# Create data frame for network statistics vs AP
network_stats_data <- all_df %>%
  select(seqname, monomer_ap, dimer_ap, PeptideLength) %>%
  pivot_longer(
    cols = c(monomer_ap, dimer_ap),
    names_to = "State",
    values_to = "AP"
  ) %>%
  filter(AP != 100) %>%  # Remove points where AP is 100
  mutate(State = ifelse(State == "monomer_ap", "Monomer", "Dimer"))

# Add network statistics
network_stats <- combined_df_ddedge_tripep %>%
  filter(frame == max(frame)) %>%
  group_by(seqname, stat_name) %>%
  summarise(value = mean(value), .groups = 'drop') %>%
  pivot_wider(names_from = stat_name, values_from = value)

# Merge with AP data
network_stats_data <- merge(network_stats_data, network_stats, by = "seqname")

# Create plots for each network statistic
g6_nsp1 <- ggplot(network_stats_data, aes(x = AP, y = nsp1, color = State)) +
  geom_point(alpha = 0.6) +
  labs(x = "", y = "NSP(1)") +
  pltthemetemp +
  theme(legend.position = "none")

g6_esp0 <- ggplot(network_stats_data, aes(x = AP, y = esp0, color = State)) +
  geom_point(alpha = 0.6) +
  labs(x = "", y = "ESP(0)") +
  pltthemetemp +
  theme(legend.position = "none")

g6_esp1 <- ggplot(network_stats_data, aes(x = AP, y = esp1, color = State)) +
  geom_point(alpha = 0.6) +
  labs(x = "", y = "ESP(1)") +
  pltthemetemp +
  theme(legend.position = "none")

# Extract legend from esp0 plot
legend_plot <- ggplot(network_stats_data, aes(x = AP, y = esp0, color = State)) +
  geom_point(alpha = 0.6) +
  pltthemetemp +
  labs(col="") +
  theme(legend.position = "right")
legend <- get_legend(legend_plot)

# Combine the plots
g6 <- plot_grid(legend, 
                plot_grid(g6_nsp1, g6_esp0, g6_esp1, ncol = 3),
                ncol = 1, rel_heights = c(0.1, 1))

# Add common x-axis label
g6 <- g6 + 
  draw_label("AP", x = 0.5, y = 0.02, size = axis_title_size)

# Print the combined plot
print(g6)
# Show plots
print(ggplotly(g1))
print(ggplotly(g2))
print(ggplotly(g3))
print(ggplotly(g4))
print(ggplotly(g5))

if (exists("save_plots") && save_plots) {
  # Optionally, save plots
  # ggsave(file.path(pltsavedir, "dAP_distribution_by_length.png"), g1, width=6.5, height=4, dpi=1100)
  # ggsave(file.path(pltsavedir, "dEdge_distribution_by_length.png"), g2, width=6.5, height=4, dpi=1100)
  # ggsave(file.path(pltsavedir, "dEdge_vs_dAP_scatter_by_length.png"), g3, width=6.5, height=4, dpi=1100)
  ggsave(file.path(pltsavedir, "Edge_correlation_by_length.png"), g4, width=6.5, height=3, dpi=1100)
  ggsave(file.path(pltsavedir, "Edge_vs_AP_by_state.png"), g5, width=6.5, height=3, dpi=1100)
  ggsave(file.path(pltsavedir, "Network_stats_vs_AP.png"), g6, width=6.5, height=2, dpi=1100)

}


