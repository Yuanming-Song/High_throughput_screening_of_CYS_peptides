# define maindir
main_dir <- "/dfs9/tw/yuanmis1/mrsec/ML-MD-Peptide/Edgelist/Tetrapeptide/"
outdir<-file.path(main_dir,"SubData/")
missing_log<-"missing_edge_Tetrapeptide"
outpre<-"Tetrapeptide_edgelist_"

# Treat each monomer as single node?
single_node <- FALSE

# Define residues excluding "C"
residues <- c("A", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "E", "D")

# Read command-line arguments: pos (position of "C"), inpres (first non-C residue), and optionally single_node
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript |this script| <pos> <inpres> [single_node]")
}
pos <- as.numeric(args[1])
inpres <- args[2]
if (length(args) >= 3) {
  single_node <- as.logical(args[3])
}

# Generate sequences for this batch
all_peptides <- c()

if (pos == 1) {
  # "C" in position 1 and r2 = inpres: C, inpres, r3, r4
  for (r3 in residues) {
    for (r4 in residues) {
      all_peptides <- c(all_peptides, paste("C", inpres, r3, r4, sep = "_"))
    }
  }
} else if (pos == 2) {
  # "C" in position 2 and r1 = inpres: inpres, C, r3, r4
  for (r3 in residues) {
    for (r4 in residues) {
      all_peptides <- c(all_peptides, paste(inpres, "C", r3, r4, sep = "_"))
    }
  }
} else if (pos == 3) {
  # "C" in position 3 and r1 = inpres: inpres, r2, C, r4
  for (r2 in residues) {
    for (r4 in residues) {
      all_peptides <- c(all_peptides, paste(inpres, r2, "C", r4, sep = "_"))
    }
  }
} else if (pos == 4) {
  # "C" in position 4 and r1 = inpres: inpres, r2, r3, C
  for (r2 in residues) {
    for (r3 in residues) {
      all_peptides <- c(all_peptides, paste(inpres, r2, r3, "C", sep = "_"))
    }
  }
}

# Define the states to process (e.g., dimer and monomer)
states <- c("dimer", "monomer")

# Loop through each state
for (state in states) {
  tempedgelist <- list()             # Temporary list for edgelists for this state
  # Only add single_node suffix for dimer state
  log_file <- file.path(outdir, paste0(missing_log, "_", state, "_C", pos, "_", inpres, 
                                        if(state == "dimer" && single_node) "_single_node" else "", ".log"))
  log_con <- file(log_file, open = "wt")
  
  # Loop through each Tetrapeptide sequence
  for (seq in all_peptides) {
    # Construct the full simulation directory for this sequence
    # Only add single_node suffix for dimer state
    rda_file <- file.path(main_dir, state, paste0(seq, if(state == "dimer" && single_node) "_single_node" else "", ".rda"))
    
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
  
  # Save the cumulative edgelist for this batch
  # Only add single_node suffix for dimer state
  out_rda <- file.path(outdir, paste0(outpre, state, "_C", pos, "_", inpres, 
                                       if(state == "dimer" && single_node) "_single_node" else "", ".rda"))
  edgelist <- tempedgelist
  save(edgelist, file = out_rda)
}