# Tetrapeptide multi-stat lookup utility

library(reshape2)

# Build a lookup-ready data frame for tetrapeptides
build_tetra_lookup_df <- function(df) {
  # Only use last frame for each sequence/state/stat
  df <- df[df$frame == ave(df$frame, df$seqname, df$state, df$stat_name, FUN = max), ]
  # Helper to get wide format for a given value column
  get_wide <- function(value_col) {
    dcast(df, seqname ~ state + stat_name, value.var = value_col)
  }
  # Get wide tables for value and ratio
  wide_value <- get_wide("value")
  wide_ratio <- get_wide("ratio")
  # Merge on seqname
  wide <- merge(wide_value, wide_ratio, by = "seqname", suffixes = c("_value", "_ratio"))
  # Calculate deltas and APs (example for 'edges' stat)
  if (all(c("dimer_edges", "monomer_edges") %in% colnames(wide))) {
    wide$dEdge <- wide$dimer_edges - wide$monomer_edges
  }
  if (all(c("dimer_edges_ratio", "monomer_edges_ratio") %in% colnames(wide))) {
    wide$dAP <- wide$dimer_edges_ratio - wide$monomer_edges_ratio
    wide$monomer_ap <- wide$monomer_edges_ratio
    wide$dimer_ap <- wide$dimer_edges_ratio
  }
  # Add more stats as needed
  wide
}

# General multi-stat lookup function
lookup_multi_stat <- function(df, ranges, rank_by = NULL, decreasing = TRUE) {
  for (stat in names(ranges)) {
    min_val <- ranges[[stat]][1]
    max_val <- ranges[[stat]][2]
    df <- df[df[[stat]] >= min_val & df[[stat]] <= max_val, ]
  }
  if (!is.null(rank_by)) {
    df <- df[order(df[[rank_by]], decreasing = decreasing), ]
  }
  print(df$seqname)
  return(df)
}

# Example usage:
 lookup_df <- build_tetra_lookup_df(combined_df_netstat_tetra)
lookup_multi_stat(lookup_df, list(dEdge = c(2, 3), monomer_edges = c(100, 200), dAP = c(0.5, 1.0), dimer_ap = c(1.2, 2.0)), rank_by = "dAP") 