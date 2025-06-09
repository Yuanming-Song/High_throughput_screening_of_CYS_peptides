#Treat each monomer as single node?
single_node<-FALSE
# define maindir
main_dir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Dipeptide/"
missing_log<-"missing_edge_Dipeptide"
# Define residues excluding "C"
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Generate all possible dipeptide sequences with exactly one "C"
all_dipeptides <- c()
for (pos in 1:2) {
  if (pos == 1) {
    for (r2 in residues) {
      all_dipeptides <- c(all_dipeptides, paste("C", r2, sep = "_"))
    }
  } else if (pos == 2) {
    for (r1 in residues) {
      all_dipeptides <- c(all_dipeptides, paste(r1, "C", sep = "_"))
    }
  }
}

# Verify total count is 2*19
expected_count <- 2 * 19
if (length(all_dipeptides) != expected_count) {
  stop("Error: The number of generated sequences does not match 2*19.")
}

# Define the states to process (e.g., dimer and monomer)
states <- c("dimer", "monomer"
	    )

# Loop through each state
for (state in states) {
  tempedgelist <- list()             # Temporary list for edgelists for this state
  log_file <- file.path(main_dir, paste0(missing_log, state, ".log"))
  log_con <- file(log_file, open = "wt")
  
  # Loop through each dipeptide sequence
  for (seq in all_dipeptides) {
    # Construct the full simulation directory for this sequence
    # (Assuming that under main_dir/state/ the file for a sequence is named "SEQ.rda")
    rda_file <- file.path(main_dir, state, paste0(seq, 
        if (single_node && state == "dimer") "_single_node" else "",
        ".rda"))
    
    # If the .rda file does not exist, log the missing file and continue
    if (!file.exists(rda_file)) {
      writeLines(paste(state, seq, "missing"), log_con)
      next
    }
    
    # Load the .rda file; assume it loads an object named "edgelist"
    load(rda_file)
    
    # Store the edgelist in the temporary list, keyed by the sequence
    tempedgelist[[seq]] <- edgelist
  }
  close(log_con)
  
  # Save the cumulative edgelist for this state as "dipeptide_edgelist_<state>.rda" in main_dir
  out_rda <- file.path(main_dir, paste0("dipeptide_edgelist_", state, 
      if (single_node && state == "dimer") "_single_node" else "",
      ".rda"))
  edgelist <- tempedgelist
  save(edgelist, file = out_rda)
} 